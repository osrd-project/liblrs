//! Tools to make it easier to build an LRS
//! It also avoids the need to manipulate flatbuffer data

use std::collections::HashMap;
use std::path::Path;

use flatbuffers::{ForwardsUOffset, Vector, WIPOffset};
use geo::Coord;

use crate::curves::{Curve, SphericalLineStringCurve};

use crate::lrs_generated::{self, *};
use crate::osm_helpers::sort_edges;

/// A key-value `HashMap` to add metadata to the objects.
pub type Properties = std::collections::HashMap<String, String>;

#[macro_export]
/// Build a properties map:
/// `properties!("source" => "openstreetmap", "licence" => "ODbL")`.
macro_rules! properties {
    ($($k:expr => $v:expr),* $(,)?) => {{
        core::convert::From::from([$(($k.to_owned(), $v.to_owned()),)*])
    }};
}

/// The linear position of an [`Anchor`] doesn’t always match the measured distance.
/// For example if a road was transformed into a bypass, resulting in a longer road,
/// but measurements are kept the same.
/// The start of the [`Curve`] might also be different from the `0` of the LRM.
pub struct AnchorOnLrm {
    /// Index of the considered [`Anchor`]. Use the value returned by [`Builder::add_anchor`].
    pub anchor_index: usize,
    /// The distance from the start of the LRM.
    /// It can be different from the measured distance.
    pub distance_along_lrm: f64,
}

#[derive(Copy, Clone)]
/// A [`Traversal`] is composed by many [`Segment`]s.
pub struct SegmentOfTraversal {
    /// Index of the considered [`Segment`]. Use the value returned by [`Builder::add_segment`]
    pub segment_index: usize,
    /// When integrating the [`Segment`] in the [`Traversal`], if we should consider the coordinates in reverse order.
    pub reversed: bool,
}

impl From<SegmentOfTraversal> for lrs_generated::SegmentOfTraversal {
    fn from(val: SegmentOfTraversal) -> Self {
        lrs_generated::SegmentOfTraversal::new(
            val.segment_index as u64,
            match val.reversed {
                true => Direction::Decreasing,
                false => Direction::Increasing,
            },
        )
    }
}

/// Helper structure to help building an LRS file.
/// It holds all the temporary structures and is called to append more data.
#[derive(Default)]
pub struct Builder<'fbb> {
    fbb: flatbuffers::FlatBufferBuilder<'fbb>,

    // Temporary geometry of [`Segment`]s, we use them to project [`Anchor`]s.
    segment_geom: Vec<Vec<Coord>>,
    // Temporary geometry of [`Traversal`] [`Curve`]s, we use them to project [`Anchor`]s and compute length.
    traversal_curve: Vec<SphericalLineStringCurve>,
    // Temporary [`Anchor`]s because we need to project them on the [`Traversal`] of each LRM they belong to.
    anchors_coordinates: Vec<Coord>,

    // Structures that allows to find indices from their id
    traversal_map: HashMap<String, usize>,

    // Final objects that will be in the binary file.
    nodes: Vec<WIPOffset<Node<'fbb>>>,
    segments: Vec<WIPOffset<Segment<'fbb>>>,
    traversals: Vec<WIPOffset<Traversal<'fbb>>>,
    anchors: Vec<WIPOffset<Anchor<'fbb>>>,
    lrms: Vec<WIPOffset<LinearReferencingMethod<'fbb>>>,
}

impl<'fbb> Builder<'fbb> {
    /// Private helper function to transform an `HashMap` into a flatbuffer vector of `Property`
    fn build_properties(
        &mut self,
        properties: Properties,
    ) -> Option<WIPOffset<Vector<'fbb, ForwardsUOffset<Property<'fbb>>>>> {
        let properties_vec: Vec<_> = properties
            .iter()
            .map(|(k, v)| {
                let key = Some(self.fbb.create_string(k));
                let value = Some(self.fbb.create_string(v));
                Property::create(&mut self.fbb, &PropertyArgs { key, value })
            })
            .collect();
        Some(self.fbb.create_vector(&properties_vec))
    }

    /// Create a new [`Builder`] (default `size = 1024`).
    pub fn new() -> Self {
        Self {
            fbb: flatbuffers::FlatBufferBuilder::with_capacity(1024),
            ..Default::default()
        }
    }

    /// Add a new [`Node`].
    pub fn add_node(&mut self, id: &str, properties: Properties) -> usize {
        let properties = self.build_properties(properties);
        let id = Some(self.fbb.create_string(id));
        let node = Node::create(&mut self.fbb, &NodeArgs { id, properties });
        self.nodes.push(node);
        self.nodes.len() - 1
    }

    /// Add a new [`Anchor`].
    pub fn add_anchor(
        &mut self,
        id: &str,
        name: &str,
        coord: Coord,
        properties: Properties,
    ) -> usize {
        let properties = self.build_properties(properties);
        let geometry = Point::new(coord.x, coord.y);
        let anchor_arg = AnchorArgs {
            id: Some(self.fbb.create_string(id)),
            name: Some(self.fbb.create_string(name)),
            geometry: Some(&geometry),
            properties,
            ..Default::default()
        };

        self.anchors
            .push(Anchor::create(&mut self.fbb, &anchor_arg));
        self.anchors_coordinates.push(coord);

        self.anchors.len() - 1
    }

    /// Add a new [`Segment`].
    pub fn add_segment(
        &mut self,
        id: &str,
        geometry: &[Coord],
        start_node_index: usize,
        end_node_index: usize,
    ) -> usize {
        let points_iter = geometry.iter().map(|c| Point::new(c.x, c.y));
        let points = self.fbb.create_vector_from_iter(points_iter);
        let args = SegmentArgs {
            id: Some(self.fbb.create_string(id)),
            properties: None,
            geometry: Some(points),
            start_node_index: start_node_index as u64,
            end_node_index: end_node_index as u64,
        };
        self.segments.push(Segment::create(&mut self.fbb, &args));
        self.segment_geom.push(geometry.to_vec());
        self.segments.len() - 1
    }

    /// Add a new [`Traversal`], created from the [`Segment`]s provided through `Builder::add_segment_to_traversal`.
    /// The existing [`Segment`]s are consumed and will not be accessible anymore.
    pub fn add_traversal(&mut self, traversal_id: &str, segments: &[SegmentOfTraversal]) -> usize {
        let mut coords = vec![];
        for segment in segments {
            if segment.reversed {
                for &coord in self.segment_geom[segment.segment_index].iter().rev() {
                    coords.push(coord)
                }
            } else {
                for &coord in self.segment_geom[segment.segment_index].iter() {
                    coords.push(coord)
                }
            }
        }
        let curve = SphericalLineStringCurve::new(geo::LineString::new(coords), 100.);
        self.traversal_curve.push(curve);

        let segments = self.fbb.create_vector_from_iter(
            segments
                .iter()
                .map(|s| Into::<lrs_generated::SegmentOfTraversal>::into(*s)),
        );

        let args = TraversalArgs {
            id: Some(self.fbb.create_string(traversal_id)),
            segments: Some(segments),
            properties: None,
        };
        self.traversals
            .push(Traversal::create(&mut self.fbb, &args));
        let idx = self.traversals.len() - 1;
        self.traversal_map.insert(traversal_id.to_owned(), idx);
        idx
    }

    /// Create a linear referencing method where the distance is provided.
    /// If the distance is not known, use `Builder::add_lrm` to compute the distance.
    /// The [`Anchor`]s will be projected on the [`Curve`].
    pub fn add_lrm_with_distances(
        &mut self,
        id: &str,
        traversal_index: usize,
        anchors: &[AnchorOnLrm],
        properties: Properties,
    ) {
        let id = Some(self.fbb.create_string(id));
        let properties = self.build_properties(properties);
        let anchor_indices = anchors.iter().map(|a| a.anchor_index as u64);
        let distances = anchors.iter().map(|a| a.distance_along_lrm);
        let args = LinearReferencingMethodArgs {
            id,
            properties,
            traversal_index: traversal_index as u32,
            anchor_indices: Some(self.fbb.create_vector_from_iter(anchor_indices)),
            distances: Some(self.fbb.create_vector_from_iter(distances)),
            projected_anchors: Some(self.project_anchors(anchors, traversal_index)),
            ..Default::default()
        };
        self.lrms
            .push(LinearReferencingMethod::create(&mut self.fbb, &args));
    }

    /// Private helper that projects [`Anchor`]s onto a [`Curve`].
    fn project_anchors(
        &mut self,
        anchors: &[AnchorOnLrm],
        traversal: usize,
    ) -> WIPOffset<Vector<'fbb, ForwardsUOffset<lrs_generated::ProjectedAnchor<'fbb>>>> {
        let curve = &self.traversal_curve[traversal];
        let projected_anchors: Vec<_> = anchors
            .iter()
            .map(|anchor| self.anchors_coordinates[anchor.anchor_index].into())
            .map(|coord| curve.project(coord).expect("could not projets anchor"))
            .map(|projection| (projection.projected_coords, projection.distance_along_curve))
            .map(|(coords, distance)| (lrs_generated::Point::new(coords.x(), coords.y()), distance))
            .map(|(geometry, distance_along_curve)| {
                let args = ProjectedAnchorArgs {
                    geometry: Some(&geometry),
                    distance_along_curve,
                };
                ProjectedAnchor::create(&mut self.fbb, &args)
            })
            .collect();
        self.fbb.create_vector(&projected_anchors)
    }

    /// Save the flatbuffer to the given file.
    pub fn save<P: AsRef<Path>>(&mut self, out_file: &P, properties: Properties) {
        std::fs::write(out_file, self.build_data(properties)).unwrap();
    }

    /// Return the binary data.
    pub fn build_data(&mut self, properties: Properties) -> &[u8] {
        let properties = self.build_properties(properties);
        let lrs_args = LrsArgs {
            properties,
            nodes: Some(self.fbb.create_vector(&self.nodes)),
            segments: Some(self.fbb.create_vector(&self.segments)),
            traversals: Some(self.fbb.create_vector(&self.traversals)),
            anchors: Some(self.fbb.create_vector(&self.anchors)),
            linear_referencing_methods: Some(self.fbb.create_vector(&self.lrms)),
            geometry_type: GeometryType::Geographic,
        };

        let lrs = Lrs::create(&mut self.fbb, &lrs_args);
        self.fbb.finish(lrs, None);
        self.fbb.finished_data()
    }

    /// Return the mapping between a traversal id and its index in the builder.
    pub fn get_traversal_indexes(&self) -> HashMap<String, usize> {
        self.traversal_map.clone()
    }

    /// Read the topology from an OpenStreetMap source.
    /// It will read incoming [`Node`]s and [`Segment`]s to create the [`Traversal`]s.
    pub fn read_from_osm(
        &mut self,
        input_file: &str,
        lrm_tag: &str,
        required: Vec<(String, String)>,
        to_reject: Vec<(String, String)>,
    ) {
        let mut reader = osm4routing::Reader::new().merge_ways().read_tag(lrm_tag);

        for (key, value) in required.iter() {
            reader = reader.require(key, value)
        }

        for (key, value) in to_reject.iter() {
            reader = reader.reject(key, value)
        }

        let (nodes, edges) = reader
            .read(input_file)
            .expect("could not read the osm file");

        let mut edges_map = HashMap::<_, _>::new();
        let mut traversals = HashMap::<String, Vec<_>>::new();
        let nodes_index: HashMap<_, _> = nodes
            .iter()
            .map(|n| (n.id, self.add_node(&n.id.0.to_string(), properties!())))
            .collect();

        let mut edge_index = 0;
        for edge in edges.iter() {
            if let Some(srv_ref) = edge.tags.get(lrm_tag) {
                edges_map.insert(edge.id.clone(), edge_index);
                traversals
                    .entry(srv_ref.clone())
                    .or_default()
                    .push(edge.clone());

                let start_node_index = nodes_index[&edge.source];
                let end_node_index = nodes_index[&edge.target];
                self.add_segment(&edge.id, &edge.geometry, start_node_index, end_node_index);
                edge_index += 1;
            }
        }

        // Sort the traversals
        for (srv_ref, edges) in traversals.into_iter() {
            let segments: Vec<_> = sort_edges(edges, &srv_ref)
                .into_iter()
                .map(|(edge, reversed)| SegmentOfTraversal {
                    segment_index: edges_map[&edge.id],
                    reversed,
                })
                .collect();
            self.add_traversal(&srv_ref, &segments);
        }
    }
}
