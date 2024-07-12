//! High level extensions meant for an easy usage
//! Those functions are exposed in wasm-bindings

use liblrs::lrs_ext::*;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
/// Struct exposed to js
pub struct Lrs {
    lrs: ExtLrs,
}

#[derive(Clone, Copy)]
#[wasm_bindgen]
/// A geographical [`Point`], it can be either a projected or spherical coordinates.
pub struct Point {
    /// Position on x-axis or `longitude`.
    pub x: f64,
    /// Position on y-axis or `latitude`.
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

#[wasm_bindgen(getter_with_clone)]
/// Represent a position on an [`LrmScale`] relative as an `offset` to an [`Anchor`].
pub struct LrmScaleMeasure {
    /// `name` of the reference [`Anchor`].
    anchor_name: String,
    /// `offset` to the reference [`Anchor`].
    scale_offset: f64,
}

#[wasm_bindgen]
impl LrmScaleMeasure {
    /// Build a new [`LrmMeasure`] from an [`Anchor`] `name` and the `offset` on the [`LrmScale`].
    #[wasm_bindgen(constructor)]
    pub fn new(anchor_name: &str, scale_offset: f64) -> Self {
        Self {
            anchor_name: anchor_name.to_owned(),
            scale_offset,
        }
    }
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
        Self {
            anchor_name: val.anchor_name.clone(),
            scale_offset: val.scale_offset,
        }
    }
}

#[wasm_bindgen]
/// An `Anchor` is a reference point for a given [`Curve`].
pub struct Anchor {
    #[wasm_bindgen(getter_with_clone)]
    /// `name` of the [`Anchor`].
    pub name: String,
    /// Projected position on the [`Curve`] (the reference point isn’t always on the curve).
    pub position: Point,
    /// Position on the [`Curve`].
    pub curve_position: f64,
    /// Position on the scale.
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

#[wasm_bindgen]
impl Lrs {
    /// Load the data.
    pub fn load(data: &[u8]) -> Result<Lrs, String> {
        ExtLrs::load(data).map(|lrs| Self { lrs })
    }

    /// How many LRMs compose the LRS.
    pub fn lrm_len(&self) -> usize {
        self.lrs.lrm_len()
    }

    /// Return the geometry of the LRM.
    pub fn get_lrm_geom(&self, index: usize) -> Result<Vec<Point>, String> {
        self.lrs
            .get_lrm_geom(index)
            .map(|coords| coords.into_iter().map(|coord| coord.into()).collect())
            .map_err(|e| e.to_string())
    }

    /// `id` of the [`LrmScale`].
    pub fn get_lrm_scale_id(&self, index: usize) -> String {
        self.lrs.get_lrm_scale_id(index)
    }

    /// All the [`Anchor`]s of a LRM.
    pub fn get_anchors(&self, lrm_index: usize) -> Result<Vec<Anchor>, String> {
        self.lrs
            .get_anchors(lrm_index)
            .map(|anchor| anchor.into_iter().map(Anchor::from).collect())
            .map_err(|e| e.to_string())
    }

    /// Get the position given a [`LrmScaleMeasure`].
    pub fn resolve(&self, lrm_index: usize, measure: &LrmScaleMeasure) -> Result<Point, String> {
        self.lrs
            .resolve(lrm_index, &measure.into())
            .map(Point::from)
            .map_err(|e| e.to_string())
    }

    /// Given two [`LrmScaleMeasure`]s, return a range of [`LineString`].
    pub fn resolve_range(
        &self,
        lrm_index: usize,
        from: &LrmScaleMeasure,
        to: &LrmScaleMeasure,
    ) -> Result<Vec<Point>, String> {
        self.lrs
            .resolve_range(lrm_index, &from.into(), &to.into())
            .map(|coords| coords.into_iter().map(|coord| coord.into()).collect())
    }
}

#[wasm_bindgen]
/// Display stacktrace in case of a panic.
pub fn set_panic_hook() {
    // When the `console_error_panic_hook` feature is enabled, we can call the
    // `set_panic_hook` function at least once during initialization, and then
    // we will get better error messages if our code ever panics.
    //
    // For more details see
    // https://github.com/rustwasm/console_error_panic_hook#readme
    #[cfg(feature = "console_error_panic_hook")]
    console_error_panic_hook::set_once();
}
