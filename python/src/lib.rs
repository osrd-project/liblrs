//! High level extensions meant for an easy usage
//! Those functions are exposed in wasm-bindings

use std::path::PathBuf;

use liblrs::lrs_ext::*;
use liblrs::{builder::Properties, lrs::LrsBase};
use pyo3::{exceptions::PyTypeError, prelude::*};

/// Holds the whole Linear Referencing System.
#[pyclass]
pub struct Lrs {
    lrs: ExtLrs,
}

#[pymodule]
fn liblrs_python(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Lrs>()?;
    m.add_class::<LrmScaleMeasure>()?;
    m.add_class::<Anchor>()?;
    m.add_class::<Point>()?;
    m.add_class::<AnchorOnLrm>()?;
    m.add_class::<SegmentOfTraversal>()?;
    m.add_class::<Builder>()?;
    Ok(())
}

#[derive(Clone, Copy)]
/// A geographical [`Point`], it can be either a projected or spherical coordinates.
#[pyclass]
pub struct Point {
    /// Position on x-axis or `longitude`.
    #[pyo3(get, set)]
    pub x: f64,
    /// Position on y-axis or `latitude`.
    #[pyo3(get, set)]
    pub y: f64,
}

#[pymethods]
impl Point {
    #[new]
    /// Build a new geographical point.
    ///
    /// When using spherical coordinates, longitude is x and latitude y
    fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
}

impl From<geo_types::Point> for Point {
    fn from(value: geo_types::Point) -> Self {
        Self {
            x: value.x(),
            y: value.y(),
        }
    }
}

impl From<geo_types::Coord> for Point {
    fn from(value: geo_types::Coord) -> Self {
        Self {
            x: value.x,
            y: value.y,
        }
    }
}

impl From<Point> for geo_types::Coord {
    fn from(val: Point) -> Self {
        geo_types::Coord { x: val.x, y: val.y }
    }
}

#[pyclass]
/// Represent a position on an [`LrmScale`] relative as an `offset` to an [`Anchor`].
pub struct LrmScaleMeasure {
    #[pyo3(get, set)]
    /// `name` of the reference [`Anchor`].
    anchor_name: String,
    #[pyo3(get, set)]
    /// `offset` to the reference [`Anchor`].
    scale_offset: f64,
}

impl From<&liblrs::lrm_scale::LrmScaleMeasure> for LrmScaleMeasure {
    fn from(value: &liblrs::lrm_scale::LrmScaleMeasure) -> Self {
        Self {
            anchor_name: value.anchor_name.clone(),
            scale_offset: value.scale_offset,
        }
    }
}

impl From<&LrmScaleMeasure> for liblrs::lrm_scale::LrmScaleMeasure {
    fn from(val: &LrmScaleMeasure) -> Self {
        liblrs::lrm_scale::LrmScaleMeasure {
            anchor_name: val.anchor_name.clone(),
            scale_offset: val.scale_offset,
        }
    }
}

#[pymethods]
impl LrmScaleMeasure {
    #[new]
    /// Build a new [`LrmMeasure`] from an [`Anchor`] `name` and the `offset` on the [`LrmScale`].
    fn new(anchor_name: String, scale_offset: f64) -> Self {
        Self {
            anchor_name,
            scale_offset,
        }
    }
}

#[derive(Clone, Copy)]
#[pyclass]
/// A traversal is composed by segments
pub struct SegmentOfTraversal {
    /// Index of the considered segment. Use the value returned by [`Builder::add_segment`]
    #[pyo3(get, set)]
    pub segment_index: usize,
    /// When integrating the segment in the traversal, should we consider the coordinates in the reverse order
    #[pyo3(get, set)]
    pub reversed: bool,
}

#[pymethods]
impl SegmentOfTraversal {
    #[new]
    fn new(segment_index: usize, reversed: bool) -> Self {
        Self {
            segment_index,
            reversed,
        }
    }
}

impl From<SegmentOfTraversal> for liblrs::builder::SegmentOfTraversal {
    fn from(val: SegmentOfTraversal) -> Self {
        liblrs::builder::SegmentOfTraversal {
            segment_index: val.segment_index,
            reversed: val.reversed,
        }
    }
}

#[derive(Clone, Copy)]
#[pyclass]
/// The linear position of an anchor doesn’t always match the measured distance
/// For example if a road was transformed into a bypass, resulting in a longer road,
/// but measurements are kept the same
/// The start of the curve might also be different from the `0` of the LRM
pub struct AnchorOnLrm {
    /// Index of the considered anchor. Use the value returned by [`Builder::add_anchor`]
    #[pyo3(get, set)]
    pub anchor_index: usize,
    /// The distance from the start of the LRM.
    /// It can be different from the measured distance
    #[pyo3(get, set)]
    pub distance_along_lrm: f64,
}

#[pymethods]
impl AnchorOnLrm {
    #[new]
    fn new(anchor_index: usize, distance_along_lrm: f64) -> Self {
        Self {
            anchor_index,
            distance_along_lrm,
        }
    }
}

impl From<AnchorOnLrm> for liblrs::builder::AnchorOnLrm {
    fn from(val: AnchorOnLrm) -> Self {
        liblrs::builder::AnchorOnLrm {
            anchor_index: val.anchor_index,
            distance_along_lrm: val.distance_along_lrm,
        }
    }
}

#[pyclass]
/// An `Anchor` is a reference point for a given [`Curve`].
/// It can be a milestone, a bridge…
pub struct Anchor {
    /// `name` of the [`Anchor`].
    #[pyo3(get, set)]
    pub name: String,
    /// Projected position on the [`Curve`] (the reference point isn’t always on the curve).
    #[pyo3(get, set)]
    pub position: Option<Point>,
    /// Position on the [`Curve`].
    #[pyo3(get, set)]
    pub curve_position: f64,
    /// Position on the scale.
    #[pyo3(get, set)]
    pub scale_position: f64,
}

impl From<&liblrs::lrm_scale::Anchor> for Anchor {
    fn from(value: &liblrs::lrm_scale::Anchor) -> Self {
        Self {
            name: value.clone().id.unwrap_or_else(|| "-".to_owned()),
            position: value.point.map(|p| p.into()),
            curve_position: value.curve_position,
            scale_position: value.scale_position,
        }
    }
}

#[pymethods]
impl Lrs {
    /// Load the data.
    #[new]
    pub fn load(data: &[u8]) -> PyResult<Lrs> {
        ExtLrs::load(data)
            .map(|lrs| Self { lrs })
            .map_err(|e| PyTypeError::new_err(e.to_string()))
    }

    /// How many LRMs compose the LRS.
    pub fn lrm_len(&self) -> usize {
        self.lrs.lrm_len()
    }

    /// Return the geometry of the LRM.
    pub fn get_lrm_geom(&self, index: usize) -> PyResult<Vec<Point>> {
        self.lrs
            .get_lrm_geom(index)
            .map(|coords| coords.into_iter().map(|coord| coord.into()).collect())
            .map_err(|e| PyTypeError::new_err(e.to_string()))
    }

    /// `id` of the [`LrmScale`].
    pub fn get_lrm_scale_id(&self, index: usize) -> String {
        self.lrs.get_lrm_scale_id(index)
    }

    /// All the [`Anchor`]s of a LRM.
    pub fn get_anchors(&self, lrm_index: usize) -> Vec<Anchor> {
        self.lrs
            .get_anchors(lrm_index)
            .iter()
            .map(Anchor::from)
            .collect()
    }

    /// Get the position given a [`LrmScaleMeasure`].
    pub fn resolve(&self, lrm_index: usize, measure: &LrmScaleMeasure) -> PyResult<Point> {
        self.lrs
            .resolve(lrm_index, &measure.into())
            .map(Point::from)
            .map_err(|e| PyTypeError::new_err(e.to_string()))
    }

    /// Get the positon along the curve given a [`LrmScaleMeasure`]
    /// The value will be between 0.0 and 1.0, both included
    pub fn locate_point(&self, lrm_index: usize, measure: &LrmScaleMeasure) -> PyResult<f64> {
        self.lrs.lrs.lrms[lrm_index]
            .scale
            .locate_point(&measure.into())
            .map_err(|e| PyTypeError::new_err(e.to_string()))
    }

    /// Given two [`LrmScaleMeasure`]s, return a range of [`Point`] that represent a line string.
    pub fn resolve_range(
        &self,
        lrm_index: usize,
        from: &LrmScaleMeasure,
        to: &LrmScaleMeasure,
    ) -> PyResult<Vec<Point>> {
        self.lrs
            .resolve_range(lrm_index, &from.into(), &to.into())
            .map(|coords| coords.into_iter().map(|coord| coord.into()).collect())
            .map_err(|e| PyTypeError::new_err(e.to_string()))
    }

    /// Given a ID returns the corresponding lrs index (or None if not found)
    pub fn find_lrm(&self, lrm_id: &str) -> Option<usize> {
        self.lrs.lrs.get_lrm(lrm_id).map(|handle| handle.0)
    }
}

#[pyclass]
struct Builder {
    inner: liblrs::builder::Builder<'static>,
}

#[pymethods]
impl Builder {
    #[new]
    /// Instantiate a new builder
    fn new() -> Self {
        Self {
            inner: liblrs::builder::Builder::new(),
        }
    }

    /// Add a new topological node (e.g. a railway switch)
    pub fn add_node(&mut self, id: &str, coord: Point, properties: Properties) -> usize {
        self.inner.add_node(id, coord.into(), properties)
    }

    /// Add a new anchor by its cooordinates
    #[pyo3(signature = (id, coord, properties, name=None))]
    pub fn add_anchor(
        &mut self,
        id: &str,
        coord: Point,
        properties: Properties,
        name: Option<&str>,
    ) -> usize {
        self.inner.add_anchor(id, name, coord.into(), properties)
    }

    /// Add a new anchor by its position along the curve
    #[pyo3(signature = (id, position_on_curve, properties, name=None))]
    pub fn add_projected_anchor(
        &mut self,
        id: &str,
        position_on_curve: f64,
        properties: Properties,
        name: Option<&str>,
    ) -> usize {
        self.inner
            .add_projected_anchor(id, name, position_on_curve, properties)
    }

    /// Add a new segment
    ///
    /// The geometry represents the curve
    /// start_node_index and end_node_index are the topological extremeties returned by `add_node`
    pub fn add_segment(
        &mut self,
        id: &str,
        geometry: Vec<Point>,
        start_node_index: usize,
        end_node_index: usize,
    ) -> usize {
        let geometry: Vec<_> = geometry.into_iter().map(|point| point.into()).collect();
        self.inner
            .add_segment(id, &geometry, start_node_index, end_node_index)
    }

    /// Add a traversal
    ///
    /// segments represent the curve of the traversal
    pub fn add_traversal(&mut self, traversal_id: &str, segments: Vec<SegmentOfTraversal>) {
        let segments: Vec<_> = segments.into_iter().map(|segment| segment.into()).collect();
        self.inner.add_traversal(traversal_id, &segments);
    }

    /// Add a linear referencing model
    ///
    /// It is composed by the traversal identified by traversa_index (that represents the curve)
    /// and the anchors (that represent the milestones)
    pub fn add_lrm(
        &mut self,
        id: &str,
        traversal_index: usize,
        anchors: Vec<AnchorOnLrm>,
        properties: Properties,
    ) {
        let anchors: Vec<_> = anchors.into_iter().map(|anchor| anchor.into()).collect();
        self.inner
            .add_lrm(id, traversal_index, &anchors, properties)
    }

    /// List all the traversals by their id and index
    pub fn get_traversal_indexes(&mut self) -> std::collections::HashMap<String, usize> {
        self.inner.get_traversal_indexes()
    }

    /// Read the topology from an OpenStreetMap source
    ///
    /// It reads the nodes, segments and traversals.
    pub fn read_from_osm(
        &mut self,
        input_osm_file: PathBuf,
        lrm_tag: String,
        required: Vec<(String, String)>,
        to_reject: Vec<(String, String)>,
    ) {
        self.inner
            .read_from_osm(&input_osm_file, &lrm_tag, required, to_reject)
    }

    /// Save the lrs to a file
    pub fn save(&mut self, out_file: PathBuf, properties: Properties) {
        self.inner.save(&out_file, properties)
    }

    /// Compute the euclidean distance between two lrms
    pub fn euclidean_distance(&self, lrm_index_a: usize, lrm_index_b: usize) -> f64 {
        self.inner.euclidean_distance(lrm_index_a, lrm_index_b)
    }

    /// List all the node indices of a traversal
    pub fn get_nodes_of_traversal(&self, lrm_index: usize) -> Vec<usize> {
        self.inner.get_nodes_of_traversal(lrm_index).to_vec()
    }

    /// Get the coordinates of a node identified by its index
    pub fn get_node_coord(&self, node_index: usize) -> Point {
        self.inner.get_node_coord(node_index).into()
    }

    /// Project a point on a the curve of an lrm
    ///
    /// Return a value between 0 and 1, both included
    /// Return None if the curve of the traversal is not defined
    pub fn project(&self, lrm_index: usize, point: Point) -> Option<f64> {
        self.inner
            .project(lrm_index, geo_types::point! {x: point.x, y: point.y})
            .map(|p| p.distance_along_curve)
            .ok()
    }

    /// Reverse the orientation of the lrm
    ///
    /// If it is composed by the segments (a, b)-(b, c) it will be (c, b)-(b, a)
    pub fn reverse(&mut self, lrm_index: usize) {
        self.inner.reverse(lrm_index)
    }
}
