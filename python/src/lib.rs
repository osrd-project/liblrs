//! High level extensions meant for an easy usage
//! Those functions are exposed in wasm-bindings

use liblrs::builder::Properties;
use liblrs::lrs_ext::*;
use pyo3::{exceptions::PyTypeError, prelude::*};

/// Holds the whole Linear Referencing System.
#[pyclass]
pub struct Lrs {
    lrs: ExtLrs,
}

#[pymodule]
fn liblrs_python<'py>(_py: Python, m: &Bound<'py, PyModule>) -> PyResult<()> {
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

impl Into<geo_types::Coord> for Point {
    fn into(self) -> geo_types::Coord {
        geo_types::Coord {
            x: self.x,
            y: self.y,
        }
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

impl Into<liblrs::lrm_scale::LrmScaleMeasure> for &LrmScaleMeasure {
    fn into(self) -> liblrs::lrm_scale::LrmScaleMeasure {
        liblrs::lrm_scale::LrmScaleMeasure {
            anchor_name: self.anchor_name.clone(),
            scale_offset: self.scale_offset,
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
    pub segment_index: usize,
    /// When integrating the segment in the traversal, should we consider the coordinates in the reverse order
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

impl Into<liblrs::builder::SegmentOfTraversal> for SegmentOfTraversal {
    fn into(self) -> liblrs::builder::SegmentOfTraversal {
        liblrs::builder::SegmentOfTraversal {
            segment_index: self.segment_index,
            reversed: self.reversed,
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
    pub anchor_index: usize,
    /// The distance from the start of the LRM.
    /// It can be different from the measured distance
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

impl Into<liblrs::builder::AnchorOnLrm> for AnchorOnLrm {
    fn into(self) -> liblrs::builder::AnchorOnLrm {
        liblrs::builder::AnchorOnLrm {
            anchor_index: self.anchor_index,
            distance_along_lrm: self.distance_along_lrm,
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
    pub position: Point,
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
            position: value.point.into(),
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

    /// Given two [`LrmScaleMeasure`]s, return a range of [`LineString`].
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
}

#[pyclass]
struct Builder {
    inner: liblrs::builder::Builder<'static>,
}

#[pymethods]
impl Builder {
    #[new]
    fn new() -> Self {
        Self {
            inner: liblrs::builder::Builder::new(),
        }
    }

    pub fn add_node(&mut self, id: &str, properties: Properties) -> usize {
        self.inner.add_node(id, properties)
    }

    pub fn add_anchor(
        &mut self,
        id: &str,
        name: &str,
        coord: Point,
        properties: Properties,
    ) -> usize {
        self.inner.add_anchor(id, name, coord.into(), properties)
    }

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

    pub fn add_traversal(&mut self, traversal_id: &str, segments: Vec<SegmentOfTraversal>) {
        let segments: Vec<_> = segments.into_iter().map(|segment| segment.into()).collect();
        self.inner.add_traversal(traversal_id, &segments);
    }

    pub fn add_lrm_with_distances(
        &mut self,
        id: &str,
        traversal_index: usize,
        anchors: Vec<AnchorOnLrm>,
        properties: Properties,
    ) {
        let anchors: Vec<_> = anchors.into_iter().map(|anchor| anchor.into()).collect();
        self.inner
            .add_lrm_with_distances(id, traversal_index, &anchors, properties)
    }

    pub fn get_traversal_indexes(&mut self) -> std::collections::HashMap<String, usize> {
        self.inner.get_traversal_indexes()
    }

    pub fn read_from_osm(
        &mut self,
        input_osm_file: String,
        lrm_tag: String,
        required: Vec<(String, String)>,
        to_reject: Vec<(String, String)>,
    ) {
        self.inner
            .read_from_osm(&input_osm_file, &lrm_tag, required, to_reject)
    }

    pub fn save(&mut self, out_file: String, properties: Properties) {
        self.inner.save(&out_file, properties)
    }
}
