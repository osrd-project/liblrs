//! High level extensions meant for an easy usage
//! Those functions are exposed in wasm-bindings

use liblrs::lrs_ext::*;
use pyo3::{exceptions::PyTypeError, prelude::*};

/// Hold the whole Linear Referencing System.
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

impl From<PositionnedAnchor> for Anchor {
    fn from(value: PositionnedAnchor) -> Self {
        Self {
            name: value.name,
            position: value.position.into(),
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
    pub fn get_anchors(&self, lrm_index: usize) -> PyResult<Vec<Anchor>> {
        self.lrs
            .get_anchors(lrm_index)
            .map(|anchor| anchor.into_iter().map(Anchor::from).collect())
            .map_err(|e| PyTypeError::new_err(e.to_string()))
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
