//! High level extensions meant for an easy usage
//! Those functions are exposed in wasm-bindings

use wasm_bindgen::prelude::*;

use crate::curves::{Curve, SphericalLineStringCurve};
use crate::lrm_scale::LrmScaleMeasure;
use crate::lrs::LrsBase;
use crate::lrs::{self, TraversalHandle, TraversalPosition};

type Lrs = lrs::Lrs<SphericalLineStringCurve>;

#[wasm_bindgen]
/// Struct exposed to js
pub struct ExtLrs {
    lrs: Lrs,
}

#[derive(Clone, Copy)]
#[wasm_bindgen]
/// A point
pub struct Point {
    /// x
    pub x: f64,
    /// y
    pub y: f64,
}

impl From<geo::Point> for Point {
    fn from(value: geo::Point) -> Self {
        Self {
            x: value.x(),
            y: value.y(),
        }
    }
}

impl From<geo::Coord> for Point {
    fn from(value: geo::Coord) -> Self {
        Self {
            x: value.x,
            y: value.y,
        }
    }
}

#[wasm_bindgen]
/// An anchor
pub struct Anchor {
    /// Name
    #[wasm_bindgen(getter_with_clone)]
    pub name: String,
    /// Project position on the curve
    pub position: Point,
    /// bla
    pub curve_position: f64,
    /// foo
    pub scale_position: f64,
}

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);
}

#[wasm_bindgen]
impl ExtLrs {
    /// Load the data
    pub fn load(data: &[u8]) -> Result<ExtLrs, String> {
        Lrs::from_bytes(data)
            .map(|lrs| Self { lrs })
            .map_err(|err| err.to_string())
    }

    /// How many lrms are there
    pub fn lrm_len(&self) -> usize {
        self.lrs.lrm_len()
    }

    /// Returns the geometry of the lrm
    pub fn get_lrm_geom(&self, index: usize) -> Result<Vec<Point>, String> {
        self.lrs
            .get_linestring(TraversalHandle(index))
            .map_err(|err| err.to_string())
            .map(|linestring| {
                linestring
                    .0
                    .iter()
                    .map(|coord| Point {
                        x: coord.x,
                        y: coord.y,
                    })
                    .collect()
            })
    }

    /// `id` of the [`LrmScale`]
    pub fn get_lrm_scale_id(&self, index: usize) -> String {
        self.lrs.lrms[index].scale.id.clone()
    }

    /// All the [`Anchor`]s of a LRM
    pub fn get_anchors(&self, lrm_index: usize) -> Result<Vec<Anchor>, String> {
        let lrm = &self.lrs.lrms[lrm_index];
        let curve = &self.lrs.traversals[lrm.reference_traversal.0].curve;
        lrm.scale
            .anchors
            .iter()
            .map(|a| make_anchor(curve, a))
            .collect()
    }

    /// Get the position given a [`LrmScaleMeasure`]
    pub fn resolve(&self, lrm_index: usize, measure: &LrmScaleMeasure) -> Result<Point, String> {
        let curve_position = self.lrs.lrms[lrm_index]
            .scale
            .locate_point(measure)
            .map_err(|e| e.to_string())?;

        let traversal_position = TraversalPosition {
            distance_from_start: curve_position,
            traversal: TraversalHandle(lrm_index),
        };
        let position = self
            .lrs
            .locate_traversal(traversal_position)
            .map_err(|e| e.to_string())?;

        Ok(position.into())
    }

    /// Given two [`LrmScaleMeasure`]s, returns a range of [`LineString`]
    pub fn resolve_range(
        &self,
        lrm_index: usize,
        from: &LrmScaleMeasure,
        to: &LrmScaleMeasure,
    ) -> Result<Vec<Point>, String> {
        let scale = &self.lrs.lrms[lrm_index].scale;
        let curve = &self.lrs.traversals[lrm_index].curve;
        let from = scale.locate_point(from).map_err(|e| e.to_string())?;
        let to = scale.locate_point(to).map_err(|e| e.to_string())?;

        log(&format!("from {}, to {}", from, to));
        match curve.sublinestring(from / curve.length(), to / curve.length()) {
            Some(linestring) => Ok(linestring.into_iter().map(|p| p.into()).collect()),
            None => Err("Could not find sublinestring".to_string()),
        }
    }
}

fn make_anchor(
    curve: &SphericalLineStringCurve,
    anchor: &crate::lrm_scale::Anchor,
) -> Result<Anchor, String> {
    let position: Point = curve
        .resolve(crate::curves::CurveProjection {
            distance_along_curve: anchor.curve_position,
            offset: 0.,
        })
        .map_err(|err| err.to_string())?
        .into();

    let name = anchor.id.clone().unwrap_or("-".to_string());
    Ok(Anchor {
        name,
        position,
        curve_position: anchor.curve_position,
        scale_position: anchor.scale_position,
    })
}

#[wasm_bindgen]
/// Display stacktrace in case of a panic
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
