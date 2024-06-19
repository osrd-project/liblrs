//! High level extensions meant for an easy usage
//! Those functions are exposed in wasm-bindings

use geo::{Coord, Point};

use crate::curves::{Curve, CurveError, SphericalLineStringCurve};
use crate::lrm_scale::Anchor;
use crate::lrm_scale::LrmScaleMeasure;
use crate::lrs::{self, TraversalHandle, TraversalPosition};
use crate::lrs::{LrsBase, LrsError};

type Lrs = lrs::Lrs<SphericalLineStringCurve>;

/// Struct exposed to js
pub struct ExtLrs {
    lrs: Lrs,
}

/// And [`Anchor`] with its coordinates
pub struct PositionnedAnchor {
    /// Name
    pub name: String,
    /// Projected position on the curve
    pub position: Coord,
    /// Position on the curve
    pub curve_position: f64,
    /// Position on the scale
    pub scale_position: f64,
}

impl PositionnedAnchor {
    fn new(anchor: &Anchor, position: Coord) -> PositionnedAnchor {
        PositionnedAnchor {
            name: anchor.id.clone().unwrap_or("-".to_owned()),
            curve_position: anchor.curve_position,
            scale_position: anchor.curve_position,
            position,
        }
    }
}

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
    pub fn get_lrm_geom(&self, index: usize) -> Result<Vec<geo::Coord>, String> {
        self.lrs
            .get_linestring(TraversalHandle(index))
            .map_err(|err| err.to_string())
            .map(|linestring| linestring.0)
    }

    /// `id` of the [`LrmScale`]
    pub fn get_lrm_scale_id(&self, index: usize) -> String {
        self.lrs.lrms[index].scale.id.clone()
    }

    /// All the [`Anchor`]s of a LRM
    pub fn get_anchors(&self, lrm_index: usize) -> Result<Vec<PositionnedAnchor>, CurveError> {
        let lrm = &self.lrs.lrms[lrm_index];
        let curve = &self.lrs.traversals[lrm.reference_traversal.0].curve;
        lrm.scale
            .anchors
            .iter()
            .map(|a| make_anchor(curve, a))
            .collect()
    }

    /// Get the position given a [`LrmScaleMeasure`]
    pub fn resolve(&self, lrm_index: usize, measure: &LrmScaleMeasure) -> Result<Point, LrsError> {
        let curve_position = self.lrs.lrms[lrm_index].scale.locate_point(measure)?;

        let traversal_position = TraversalPosition {
            distance_from_start: curve_position,
            traversal: TraversalHandle(lrm_index),
        };
        self.lrs.locate_traversal(traversal_position)
    }

    /// Given two [`LrmScaleMeasure`]s, returns a range of [`LineString`]
    pub fn resolve_range(
        &self,
        lrm_index: usize,
        from: &LrmScaleMeasure,
        to: &LrmScaleMeasure,
    ) -> Result<Vec<Coord>, String> {
        let scale = &self.lrs.lrms[lrm_index].scale;
        let curve = &self.lrs.traversals[lrm_index].curve;
        let from = scale.locate_point(from).map_err(|e| e.to_string())?;
        let to = scale.locate_point(to).map_err(|e| e.to_string())?;

        match curve.sublinestring(from / curve.length(), to / curve.length()) {
            Some(linestring) => Ok(linestring.0),
            None => Err("Could not find sublinestring".to_string()),
        }
    }
}

fn make_anchor(
    curve: &SphericalLineStringCurve,
    anchor: &crate::lrm_scale::Anchor,
) -> Result<PositionnedAnchor, CurveError> {
    let position = curve
        .resolve(crate::curves::CurveProjection {
            distance_along_curve: anchor.curve_position,
            offset: 0.,
        })
        .map(|point| point.into())?;
    Ok(PositionnedAnchor::new(anchor, position))
}
