//! High level extensions meant for an easy usage
//! Those functions are exposed in wasm-bindings

use geo::{Coord, Point};

use crate::curves::{Curve, SphericalLineStringCurve};
use crate::lrm_scale::Anchor;
use crate::lrm_scale::LrmScaleMeasure;
use crate::lrs::{self, TraversalPosition};
use crate::lrs::{LrsBase, LrsError};

type Lrs = lrs::Lrs<SphericalLineStringCurve>;

/// Struct exposed to js.
pub struct ExtLrs {
    /// The linear referencing system
    pub lrs: Lrs,
}

impl ExtLrs {
    /// Load the data.
    pub fn load(data: &[u8]) -> Result<ExtLrs, String> {
        Lrs::from_bytes(data)
            .map(|lrs| Self { lrs })
            .map_err(|err| err.to_string())
    }

    /// How many LRMs compose the LRS.
    pub fn lrm_len(&self) -> usize {
        self.lrs.lrm_len()
    }

    /// Return the geometry of the LRM.
    pub fn get_lrm_geom(&self, index: usize) -> Result<Vec<geo::Coord>, String> {
        let lrm = self.lrs.lrms.get(index).ok_or("Invalid index")?;
        self.lrs
            .get_linestring(lrm.reference_traversal)
            .map_err(|err| err.to_string())
            .map(|linestring| linestring.0)
    }

    /// `id` of the [`LrmScale`].
    pub fn get_lrm_scale_id(&self, index: usize) -> String {
        self.lrs.lrms[index].scale.id.clone()
    }

    /// All the [`Anchor`]s of a LRM.
    pub fn get_anchors(&self, lrm_index: usize) -> Vec<Anchor> {
        self.lrs.lrms[lrm_index].scale.anchors.to_vec()
    }

    /// Get the position given a [`LrmScaleMeasure`].
    pub fn resolve(&self, lrm_index: usize, measure: &LrmScaleMeasure) -> Result<Point, LrsError> {
        let lrm = &self.lrs.lrms[lrm_index];
        let curve_position = lrm.scale.locate_point(measure)?.clamp(0., 1.0);

        let traversal_position = TraversalPosition {
            curve_position,
            traversal: lrm.reference_traversal,
        };
        self.lrs.locate_traversal(traversal_position)
    }

    /// Given two [`LrmScaleMeasure`]s, return a range of [`LineString`].
    pub fn resolve_range(
        &self,
        lrm_index: usize,
        from: &LrmScaleMeasure,
        to: &LrmScaleMeasure,
    ) -> Result<Vec<Coord>, String> {
        let lrm = &self.lrs.lrms[lrm_index];
        let scale = &lrm.scale;
        let curve = &self.lrs.traversals[lrm.reference_traversal.0].curve;
        let from = scale
            .locate_point(from)
            .map_err(|e| e.to_string())?
            .clamp(0., 1.);
        let to = scale
            .locate_point(to)
            .map_err(|e| e.to_string())?
            .clamp(0., 1.);

        match curve.sublinestring(from, to) {
            Some(linestring) => Ok(linestring.0),
            None => Err("Could not find sublinestring".to_string()),
        }
    }
}
