//! Tools to make it easier to build an LRS
//! It also avoids the need to manipulate flatbuffer data

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use flatbuffers::{ForwardsUOffset, Vector, WIPOffset};
use geo::Coord;

use crate::curves::{Curve, CurveError, CurveProjection, SphericalLineStringCurve};

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
#[derive(Clone, Copy, Debug)]
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

enum AnchorPosition {
    Geographical(Coord),
    Curve(f64),
}

struct TempSegment {
    id: String,
    geometry: Vec<Coord>,
    start_node_index: u64,
    end_node_index: u64,
}

struct TempTraversal {
    id: String,
    curve: SphericalLineStringCurve,
    segments: Vec<SegmentOfTraversal>,
}

impl TempTraversal {
    fn reverse(&mut self) {
        self.curve.reverse();
        self.segments.reverse();
        for segment_of_traversal in &mut self.segments {
            segment_of_traversal.reversed = !segment_of_traversal.reversed;
        }
    }
}

/// Helper structure to help building an LRS file.
/// It holds all the temporary structures and is called to append more data.
#[derive(Default)]
pub struct Builder<'fbb> {
    fbb: flatbuffers::FlatBufferBuilder<'fbb>,

    // Temporary geometry of [`Segment`]s, we use them to project [`Anchor`]s.
    temp_segments: Vec<TempSegment>,
    // Temporary geometry of [`Traversal`]' with [`Curve`] and list of [`Segment`], we use them to project [`Anchor`]s and compute length.
    temp_traversal: Vec<TempTraversal>,
    // Temporary [`Anchor`]s because we need to project them on the [`Traversal`] of each LRM they belong to.
    temp_anchors: Vec<AnchorPosition>,
    // Position of every node
    nodes_coords: Vec<Coord>,
    // Every node of a given traversal
    nodes_of_traversal: Vec<Vec<usize>>,

    // Final objects that will be in the binary file.
    nodes: Vec<WIPOffset<Node<'fbb>>>,
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

    /// Gives the indexes of all the nodes of a traversal
    pub fn get_nodes_of_traversal(&self, traversal_idx: usize) -> &[usize] {
        &self.nodes_of_traversal[traversal_idx]
    }

    /// Add a new [`Node`].
    pub fn add_node(&mut self, id: &str, coord: Coord, properties: Properties) -> usize {
        let point = Point::new(coord.x, coord.y);
        let args = NodeArgs {
            id: Some(self.fbb.create_string(id)),
            geometry: Some(&point),
            properties: self.build_properties(properties),
        };

        self.nodes_coords.push(coord);
        self.nodes.push(Node::create(&mut self.fbb, &args));
        self.nodes.len() - 1
    }

    /// Add a new [`Anchor`] based on the coordinates.
    pub fn add_anchor(
        &mut self,
        id: &str,
        name: Option<&str>,
        coord: Coord,
        properties: Properties,
    ) -> usize {
        let properties = self.build_properties(properties);
        let geometry = Point::new(coord.x, coord.y);
        let anchor_arg = AnchorArgs {
            id: Some(self.fbb.create_string(id)),
            name: name.map(|n| self.fbb.create_string(n)),
            geometry: Some(&geometry),
            properties,
            ..Default::default()
        };

        self.anchors
            .push(Anchor::create(&mut self.fbb, &anchor_arg));
        self.temp_anchors.push(AnchorPosition::Geographical(coord));

        self.temp_anchors.len() - 1
    }

    /// A new [`Anchor`] based on its position along the curve.
    pub fn add_projected_anchor(
        &mut self,
        id: &str,
        name: Option<&str>,
        position_on_curve: f64,
        properties: Properties,
    ) -> usize {
        let properties = self.build_properties(properties);
        let anchor_arg = AnchorArgs {
            id: Some(self.fbb.create_string(id)),
            name: name.map(|n| self.fbb.create_string(n)),
            properties,
            ..Default::default()
        };

        self.anchors
            .push(Anchor::create(&mut self.fbb, &anchor_arg));
        self.temp_anchors
            .push(AnchorPosition::Curve(position_on_curve));

        self.temp_anchors.len() - 1
    }

    /// Add a new [`Segment`].
    pub fn add_segment(
        &mut self,
        id: &str,
        geometry: &[Coord],
        start_node_index: usize,
        end_node_index: usize,
    ) -> usize {
        self.temp_segments.push(TempSegment {
            id: id.to_owned(),
            geometry: geometry.to_vec(),
            start_node_index: start_node_index as u64,
            end_node_index: end_node_index as u64,
        });
        self.temp_segments.len() - 1
    }

    /// Add a new [`Traversal`], created from the [`Segment`]s provided through `Builder::add_segment`.
    /// The existing [`Segment`]s are consumed and will not be accessible anymore.
    pub fn add_traversal(&mut self, traversal_id: &str, segments: &[SegmentOfTraversal]) -> usize {
        let mut coords = vec![];
        let mut nodes_of_traversal = vec![];
        for segment in segments {
            let start_node = self.temp_segments[segment.segment_index].start_node_index as usize;
            let end_node = self.temp_segments[segment.segment_index].end_node_index as usize;
            if nodes_of_traversal.is_empty() {
                nodes_of_traversal.push(end_node);
            }
            if segment.reversed {
                nodes_of_traversal.push(start_node);
                for &coord in self.temp_segments[segment.segment_index]
                    .geometry
                    .iter()
                    .rev()
                {
                    coords.push(coord);
                }
            } else {
                if nodes_of_traversal.is_empty() {
                    nodes_of_traversal.push(start_node);
                }
                nodes_of_traversal.push(end_node);
                for &coord in self.temp_segments[segment.segment_index].geometry.iter() {
                    coords.push(coord)
                }
            }
        }

        self.temp_traversal.push(TempTraversal {
            id: traversal_id.to_owned(),
            curve: SphericalLineStringCurve::new(geo::LineString::new(coords), 100.),
            segments: segments.to_vec(),
        });
        self.nodes_of_traversal.push(nodes_of_traversal);

        self.temp_traversal.len() - 1
    }

    /// Create a linear referencing method where the distance is provided.
    /// The [`Anchor`]s will be projected on the [`Curve`].
    pub fn add_lrm(
        &mut self,
        id: &str,
        traversal_index: usize,
        anchors: &[AnchorOnLrm],
        properties: Properties,
    ) {
        let id = Some(self.fbb.create_string(id));
        let properties = self.build_properties(properties);
        let mut anchors = anchors.to_vec();
        anchors.sort_by_key(|anchor| (anchor.distance_along_lrm * 10e6) as i64);
        let anchor_indices = anchors.iter().map(|a| a.anchor_index as u64);
        let distances = anchors.iter().map(|a| a.distance_along_lrm);

        let args = LinearReferencingMethodArgs {
            id,
            properties,
            traversal_index: traversal_index as u32,
            anchor_indices: Some(self.fbb.create_vector_from_iter(anchor_indices)),
            distances: Some(self.fbb.create_vector_from_iter(distances)),
            projected_anchors: Some(self.project_anchors(&anchors, traversal_index)),
            ..Default::default()
        };
        self.lrms
            .push(LinearReferencingMethod::create(&mut self.fbb, &args));
    }

    /// Private helper that projects [`Anchor`]s onto a [`Curve`].
    fn project_anchors(
        &mut self,
        anchors: &[AnchorOnLrm],
        traversal_idx: usize,
    ) -> WIPOffset<Vector<'fbb, ForwardsUOffset<lrs_generated::ProjectedAnchor<'fbb>>>> {
        let curve = &self.temp_traversal[traversal_idx].curve;

        let projected_anchors: Vec<_> = anchors
            .iter()
            .map(|anchor| match self.temp_anchors[anchor.anchor_index] {
                AnchorPosition::Curve(distance_along_curve) => (None, distance_along_curve),
                AnchorPosition::Geographical(coord) => {
                    let projected = curve
                        .project(coord.into())
                        .expect("could not projets anchor");
                    let geometry = lrs_generated::Point::new(
                        projected.projected_coords.x(),
                        projected.projected_coords.y(),
                    );

                    (Some(geometry), projected.distance_along_curve)
                }
            })
            .map(|(geom, distance_along_curve)| {
                ProjectedAnchor::create(
                    &mut self.fbb,
                    &ProjectedAnchorArgs {
                        geometry: geom.as_ref(),
                        distance_along_curve,
                    },
                )
            })
            .collect();
        self.fbb.create_vector(&projected_anchors)
    }

    /// Save the flatbuffer to the given file.
    pub fn save<P: AsRef<Path>>(&mut self, out_file: &P, properties: Properties) {
        std::fs::write(out_file, self.build_data(properties)).unwrap();
    }

    /// Private function that builds the segments data for serialization
    fn build_segments(&mut self) -> Vec<WIPOffset<Segment<'fbb>>> {
        let segments: Vec<_> = self
            .temp_segments
            .iter()
            .map(|segment| {
                let points_iter = segment.geometry.iter().map(|c| Point::new(c.x, c.y));
                (
                    self.fbb.create_string(&segment.id),
                    self.fbb.create_vector_from_iter(points_iter),
                    segment.start_node_index,
                    segment.end_node_index,
                )
            })
            .collect();
        segments
            .into_iter()
            .map(|(id, points, start_node_index, end_node_index)| {
                Segment::create(
                    &mut self.fbb,
                    &SegmentArgs {
                        id: Some(id),
                        properties: None,
                        geometry: Some(points),
                        start_node_index,
                        end_node_index,
                    },
                )
            })
            .collect()
    }

    /// Private function that builds the traversal data for serialization
    fn build_traversals(&mut self) -> Vec<WIPOffset<Traversal<'fbb>>> {
        self.temp_traversal
            .iter()
            .map(|traversal| {
                let segments_of_traversal = self.fbb.create_vector_from_iter(
                    traversal
                        .segments
                        .iter()
                        .map(|s| Into::<lrs_generated::SegmentOfTraversal>::into(*s)),
                );

                let args = TraversalArgs {
                    id: Some(self.fbb.create_string(&traversal.id)),
                    segments: Some(segments_of_traversal),
                    properties: None,
                };
                Traversal::create(&mut self.fbb, &args)
            })
            .collect()
    }

    /// Return the binary data.
    pub fn build_data(&mut self, properties: Properties) -> &[u8] {
        let segments = self.build_segments();
        let traversals = self.build_traversals();

        let lrs_args = LrsArgs {
            properties: self.build_properties(properties),
            nodes: Some(self.fbb.create_vector(&self.nodes)),
            segments: Some(self.fbb.create_vector(&segments)),
            traversals: Some(self.fbb.create_vector(&traversals)),
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
        self.temp_traversal
            .iter()
            .enumerate()
            .map(|(idx, traversal)| (traversal.id.to_owned(), idx))
            .collect()
    }

    /// Read the topology from an OpenStreetMap source.
    /// It will read incoming [`Node`]s and [`Segment`]s to create the [`Traversal`]s.
    pub fn read_from_osm(
        &mut self,
        input_file: &PathBuf,
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
            .map(|n| {
                (
                    n.id,
                    self.add_node(&n.id.0.to_string(), n.coord, properties!()),
                )
            })
            .collect();

        for edge in edges.iter() {
            if let Some(srv_ref) = edge.tags.get(lrm_tag) {
                traversals
                    .entry(srv_ref.clone())
                    .or_default()
                    .push(edge.clone());

                let start_node_idx = nodes_index[&edge.source];
                let end_node_idx = nodes_index[&edge.target];
                let idx = self.add_segment(&edge.id, &edge.geometry, start_node_idx, end_node_idx);
                edges_map.insert(edge.id.clone(), idx);
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

    /// Gives the euclidian distance between two traversals
    /// While working on spherical coordinates, this usually doesn’t make much sense,
    /// this is good enough to sort curves by distance
    pub fn euclidean_distance(&self, lrm_index_a: usize, lrm_index_b: usize) -> f64 {
        let a = &self.temp_traversal[lrm_index_a].curve.geom;
        let b = &self.temp_traversal[lrm_index_b].curve.geom;
        geo::EuclideanDistance::euclidean_distance(a, b)
    }

    /// Returns the position along the curve of the traversal
    /// The value will be between 0.0 and 1.0, both included
    pub fn project(
        &self,
        lrm_index: usize,
        point: geo::Point,
    ) -> Result<CurveProjection, CurveError> {
        self.temp_traversal[lrm_index].curve.project(point)
    }

    /// Reverses the direction of the traversal
    pub fn reverse(&mut self, lrm_index: usize) {
        self.temp_traversal[lrm_index].reverse();
        self.nodes_of_traversal[lrm_index].reverse();
    }

    /// Returns the coordinates of a node
    pub fn get_node_coord(&self, node_index: usize) -> Coord {
        self.nodes_coords[node_index]
    }
}
