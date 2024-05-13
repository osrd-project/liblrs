//! A Linear Reference System ([`Lrs`]) is the combination of multiple [`LrmScale`] and [`Traversal`]s.
//!
//! For instance a highway could have a [`Lrs`], with an [`LrmScale`] for every direction.
//!
//! A traversal is a chain of `Segment` that builds a [`Curve`]. A segment could be the road between two intersections.

use std::cmp::Ordering;

use geo::orient::Direction;
use thiserror::Error;

use crate::curves::{Curve, CurveError, CurveProjection};
use crate::lrm_scale::{
    Anchor, CurvePosition, LrmScale, LrmScaleError, LrmScaleMeasure, ScalePosition,
};
use crate::lrs_generated;
use geo::{coord, point, LineString, Point};

/// Used as handle to identify a [`LrmScale`] within a specific [`Lrs`].
#[derive(Copy, Clone, Debug, Eq, Hash, PartialEq)]
pub struct LrmHandle(pub usize);
/// Used as handle to identify a [`Traversal`] within a specific [`Lrs`].
#[derive(Copy, Clone, Debug, Eq, Hash, PartialEq)]
pub struct TraversalHandle(pub usize);

/// Represents an Linear Reference Method (LRM).
/// It is the combination of one (or more) [`Traversal`]s for one [`LrmScale`].
pub struct Lrm {
    /// The scale of this [`Lrm`].
    pub scale: LrmScale,
    /// The [`Traversal`] that is the reference of this [`Lrm`].
    pub reference_traversal: TraversalHandle,
    /// All the [`Traversal`]s where this [`Lrm`] applies.
    pub traversals: Vec<TraversalHandle>,
}

/// A [`Traversal`] is path in the network that ends [`Curve`].
/// That [`Traversal`]s can be used for many different [`Lrm`]s.
pub struct Traversal<CurveImpl: Curve> {
    /// Identifies this [`Traversal`].
    pub id: String,
    /// The geometrical [`Curve`] of this [`Traversal`].
    pub curve: CurveImpl,
    /// All the [`Lrm`]s that use this [`Traversal`].
    pub lrms: Vec<LrmHandle>,
}

/// The Linear Reference System. It must be specified for a given implementation
/// of [Curve], such as [crate::curves::LineStringCurve].
pub struct Lrs<CurveImpl: Curve> {
    /// All the [`Lrm`] of this Lrs
    pub lrms: Vec<Lrm>,
    /// All the [`Traversal`] of this Lrs
    pub traversals: Vec<Traversal<CurveImpl>>,
}

/// The result of a projection onto an [`LrmScale`].
pub struct LrmProjection {
    /// Contains `measure` ([`LrmScaleMeasure`]) and `lrm` ([`LrmHandle`]).
    pub measure: LrmMeasure,
    /// How far from the [`Lrm`] is the [`Point`] that has been projected.
    pub orthogonal_offset: f64,
}

/// Identifies a [`ScalePosition`] on an [`LrmScale`] by the distance from the start of the scale.
#[derive(Clone, Copy, Debug)]
pub struct LrmPosition {
    /// The distance from that of the scale.
    pub distance_from_start: ScalePosition,
    /// Identifies the [`LrmScale`].
    pub lrm: LrmHandle,
}

/// Identifies a position on an [LrmScale] by distance from an [`Anchor`] of the scale.
#[derive(Clone, Debug)]
pub struct LrmMeasure {
    /// Contains `anchor_name` and `scale_offset` ([`ScalePosition`]).
    pub measure: LrmScaleMeasure,
    /// Identifies the [`LrmScale`].
    pub lrm: LrmHandle,
}

/// The result of a projection an a [`Traversal`].
#[derive(Clone, Copy, Debug)]
pub struct TraversalProjection {
    /// Distance from the start of the [`Curve`] to the [`Traversal`].
    pub distance_from_start: CurvePosition,
    /// How far from the [`Traversal`] is the [`Point`] that has been projected.
    pub orthogonal_offset: f64,
    /// Identifies the [`Traversal`].
    pub traversal: TraversalHandle,
}

/// Identifies a position on a [`Traversal`].
#[derive(Clone, Copy, Debug)]
pub struct TraversalPosition {
    /// Distance from the start of the [`Curve`] to the [`Traversal`].
    pub distance_from_start: CurvePosition,
    /// Identifies the [`Traversal`].
    pub traversal: TraversalHandle,
}

impl From<TraversalPosition> for CurveProjection {
    fn from(val: TraversalPosition) -> Self {
        CurveProjection {
            distance_along_curve: val.distance_from_start,
            offset: 0.,
        }
    }
}

/// Describes an interval (= range) on a [`Traversal`].
/// The borders are [`CurvePosition`]s.
/// It can be used to identify a speed limit zone for instance.
pub struct TraversalRange {
    /// Identifies the [`Traversal`].
    pub traversal: TraversalHandle,
    /// Begin of the range.
    pub begin: CurvePosition,
    /// End of the range.
    pub end: CurvePosition,
    /// [`Direction`] of the range.
    pub direction: Direction,
}

/// Describes an interval (= range) on a [`LrmScale`].
/// The borders are [`LrmScaleMeasure`]s.
/// It can be used to identify a speed limit zone for instance.
pub struct LrmRange {
    /// Identifies the [`Lrm`].
    pub lrm: LrmHandle,
    /// Begin of the range.
    pub begin: LrmScaleMeasure,
    /// End of the range.
    pub end: LrmScaleMeasure,
    /// [`Direction`] of the range.
    pub direction: Direction,
}

impl<CurveImpl: Curve> Lrs<CurveImpl> {
    /// Number of lrms
    pub fn lrm_len(&self) -> usize {
        self.lrms.len()
    }

    /// Loads an [`Lrs`] from an byte array
    pub fn from_bytes(buf: &[u8]) -> Result<Self, LrsError> {
        let lrs = lrs_generated::root_as_lrs(buf).map_err(LrsError::InvalidArchive)?;

        let network = lrs
            .networks()
            .ok_or(LrsError::IncompleteArchive("network".to_string()))?
            .get(0);
        let geometry_view = lrs
            .views()
            .ok_or(LrsError::IncompleteArchive("geometry_view".to_string()))?
            .get(0);
        let segments_geometry = geometry_view
            .networks()
            .ok_or(LrsError::IncompleteArchive(
                "geometry view’s network".to_string(),
            ))?
            .get(0)
            .segments();

        let mut result = Self {
            lrms: vec![],
            traversals: vec![],
        };

        let source_anchors = lrs
            .anchors()
            .ok_or(LrsError::IncompleteArchive("anchors".to_owned()))?;
        // Read the traversals and build the curves
        for traversal in network.traversals() {
            let mut coords = vec![];
            for idx in 0..traversal.segments().len() {
                let segment_idx = traversal.segments().get(idx) as usize;
                let direction = traversal.directions().get(idx);
                let mut geom: Vec<_> = segments_geometry
                    .get(segment_idx)
                    .points()
                    .iter()
                    .map(|p| (coord! {x: p.x(),y: p.y()}))
                    .collect();
                if direction == lrs_generated::Direction::Decreasing {
                    geom.reverse();
                }
                coords.append(&mut geom);
            }
            let line_string = geo::LineString::new(coords);
            result.traversals.push(Traversal {
                id: traversal.id().to_owned(),
                curve: CurveImpl::new(line_string, 1000.),
                lrms: vec![],
            });
        }

        // Read the lrm scales
        for (lrm_idx, raw_lrm) in lrs.linear_referencing_methods().unwrap().iter().enumerate() {
            let reference_traversal_idx = raw_lrm
                .traversal_index()
                .ok_or(LrsError::IncompleteArchive(format!(
                    "reference traversal for LRM {lrm_idx}"
                )))?
                .traversal_index() as usize;
            let curve = &result
                .traversals
                .get(reference_traversal_idx)
                .ok_or(LrsError::IncompleteArchive(format!(
                    "traversal {reference_traversal_idx} from lrm {lrm_idx}"
                )))?
                .curve;

            let anchors: Vec<_> = raw_lrm
                .anchor_indices()
                .iter()
                .enumerate()
                .map(|(idx, anchor_idx)| {
                    let anchor = source_anchors.get(anchor_idx as usize);
                    let scale_position = raw_lrm.distances().get(idx);
                    let p = geometry_view
                        .anchors()
                        .expect("No anchors")
                        .get(anchor_idx as usize)
                        .geom()
                        .map(|p| point! {x: p.x(), y: p.y()})
                        .expect("Anchor without geometry");
                    let curve_position = curve
                        .project(p)
                        .expect("Could not project anchor to curve")
                        .distance_along_curve;

                    match anchor.name() {
                        Some(name) => Anchor::new(name, scale_position, curve_position),
                        None => Anchor::new_unnamed(scale_position, curve_position),
                    }
                })
                .collect();

            let mut lrm = Lrm {
                scale: LrmScale {
                    id: raw_lrm.id().to_owned(),
                    anchors,
                },
                reference_traversal: TraversalHandle(reference_traversal_idx),
                traversals: vec![],
            };

            for traversal_ref in raw_lrm
                .used_on()
                .ok_or(LrsError::IncompleteArchive(format!(
                    "used_on for lrm {lrm_idx}"
                )))?
            {
                let traversal_idx = traversal_ref.traversal_index() as usize;
                result.traversals[traversal_idx]
                    .lrms
                    .push(LrmHandle(lrm_idx));
                lrm.traversals.push(TraversalHandle(traversal_idx));
            }

            result.lrms.push(lrm);
        }

        Ok(result)
    }

    /// Loads an [`Lrs`] from the file system
    pub fn new<P: AsRef<std::path::Path>>(filename: P) -> Result<Self, LrsError> {
        use std::io::Read;
        let mut f = std::fs::File::open(filename).map_err(|_| LrsError::OpenFileError)?;
        let mut buf = Vec::new();
        f.read_to_end(&mut buf)
            .map_err(|_| LrsError::ReadFileError)?;

        Self::from_bytes(&buf)
    }
}

/// Errors when manipulating [`Lrs`].
#[derive(Error, Debug, PartialEq)]
pub enum LrsError {
    /// The [`LrmHandle`] is not valid. It might have been built manually or the structure mutated.
    #[error("invalid handle")]
    InvalidHandle,
    /// An error occured while manipulating a [`Curve`] of the [`Lrs`].
    #[error("curve error")]
    CurveError(#[from] CurveError),
    /// An error occured while manipulating a [`LrmScale`] of the [`Lrs`].
    #[error("curve error")]
    LrmScaleError(#[from] LrmScaleError),
    /// Could not open the LRS file.
    #[error("open file error")]
    OpenFileError,
    /// Could not read the LRS file.
    #[error("read file error")]
    ReadFileError,
    /// Could not parse the LRS file.
    #[error("invalid flatbuffer content {0}")]
    InvalidArchive(#[from] flatbuffers::InvalidFlatbuffer),
    /// The archive does not have all the required data
    #[error("the archive does not have all the required data: {0} is missing")]
    IncompleteArchive(String),
}

/// The basic functions to manipulate the [`Lrs`].
pub trait LrsBase {
    /// Returns the [`LrmHandle`] (if it exists) of the [`Lrm`] identified by its `lrm_id`.
    fn get_lrm(&self, lrm_id: &str) -> Option<LrmHandle>;
    /// Returns the [`TraversalHandle`] (if it exists) of the [`Traversal`] identified by its `traversal_id`.
    fn get_traversal(&self, traversal_id: &str) -> Option<TraversalHandle>;

    /// Returns the [`Curve`] as a [`LineString`]
    /// If the implementation uses an other format (e.g. splines),
    /// it will be segmentized as a [`LineString`] and might not be as acurate as the underlying representation
    fn get_linestring(&self, traversal: TraversalHandle) -> Result<LineString, LrsError>;

    /// Projects a [`Point`] on all applicable [`Traversal`]s to a given [`Lrm`].
    /// The [`Point`] must be in the bounding box of the [`Curve`] of the [`Traversal`].
    /// The result is sorted by `orthogonal_offset`: the nearest [`Lrm`] to the [`Point`] is the first item.
    fn lookup(&self, point: Point, lrm: LrmHandle) -> Vec<LrmProjection>;
    /// Projects a [`Point`] on all [`Lrm`] where the [`Point`] is in the bounding box.
    /// The result is sorted by `orthogonal_offset`: the nearest [`Lrm`] to the [`Point`] is the first item.
    fn lookup_lrms(&self, point: Point) -> Vec<LrmProjection>;
    /// Projects a [`Point`] on all [`Traversal`]s where the [`Point`] is in the bounding box.
    /// The result is sorted by `orthogonal_offset`: the nearest [`Lrm`] to the [`Point`] is the first item.
    fn lookup_traversals(&self, point: Point) -> Vec<TraversalProjection>;

    /// Given a [`TraversalPosition`], returns it geographical position ([`Point`]).
    fn locate_traversal(&self, position: TraversalPosition) -> Result<Point, LrsError>;

    /// And [`Lrm`] can be used on many [`Traversal`]s.
    /// For example, for both directions of a highway.
    /// This method returns all applicable [`TraversalHandle`]s to that [`Traversal`].
    fn get_lrm_applicable_traversals(&self, lrm: LrmHandle) -> &[TraversalHandle];
    /// An [`Lrm`] has a reference [`Traversal`]
    /// For example, for the centerline of a highway, or a specific track.
    /// This methods returns the [`TraversalHandle`].
    fn get_lrm_reference_traversal(&self, lrm: LrmHandle) -> TraversalHandle;

    /// A [`Traversal`] can be use for multiple [`Lrm`]s.
    /// For example, a highway could have milestones referenced in `miles` AND `kilometers`.
    fn get_traversal_lrms(&self, traversal: TraversalHandle) -> &[LrmHandle];

    /// Projects a [`TraversalPosition`] on a [`Traversal`] onto an other [`Traversal`],
    /// e.g. when placing a point on both sides of the highway.
    fn traversal_project(
        &self,
        position: TraversalPosition,
        onto: TraversalHandle,
    ) -> Result<TraversalProjection, LrsError>;

    /// Projects a [`TraversalRange`] on a [`Traversal`] onto an other [`Traversal`],
    /// e.g. when placing a stretch where wild animals cross.
    fn traversal_project_range(
        &self,
        range: TraversalRange,
        onto: TraversalHandle,
    ) -> Result<TraversalRange, LrsError>;

    /// Given the [`TraversalPosition`] on a [`Traversal`], projects that [`TraversalPosition`] onto an [`Lrm`].
    fn lrm_project(
        &self,
        position: TraversalPosition,
        onto: LrmHandle,
    ) -> Result<LrmProjection, LrsError>;

    /// Projects a [`TraversalRange`] on a [`Traversal`] onto an [LrmScale],
    /// e.g. when placing a stretch where wild animals cross.
    fn lrm_project_range(
        &self,
        range: TraversalRange,
        onto: LrmHandle,
    ) -> Result<LrmRange, LrsError>;

    /// Given a [`LrmPosition`], returns its [`LrmMeasure`].
    /// It will find the nearest [`Anchor`] that gives a positive `offset`.
    fn lrm_get_measure(&self, position: LrmPosition) -> Result<LrmMeasure, LrsError>;
    /// Given an [`LrmMeasure`], returns its [`LrmPosition`].
    fn lrm_get_position(&self, measure: LrmMeasure) -> Result<LrmPosition, LrsError>;

    // TODO
    // fn traversal_get_segment(position: TraversalPosition) -> SegmentPosition;
    // fn traversal_range_get_segments(range: TraversalRange) -> Vec<SegmentRange>;
}

impl<CurveImpl: Curve> LrsBase for Lrs<CurveImpl> {
    fn get_lrm(&self, lrm_id: &str) -> Option<LrmHandle> {
        self.lrms
            .iter()
            .position(|lrm| lrm.scale.id == lrm_id)
            .map(LrmHandle)
    }

    fn get_traversal(&self, traversal_id: &str) -> Option<TraversalHandle> {
        self.traversals
            .iter()
            .position(|traversal| traversal.id == traversal_id)
            .map(TraversalHandle)
    }

    fn lookup(&self, point: Point, lrm_handle: LrmHandle) -> Vec<LrmProjection> {
        let lrm = &self.lrms[lrm_handle.0];
        let mut result: Vec<_> = self
            .get_lrm_applicable_traversals(lrm_handle)
            .iter()
            .flat_map(|t| self.traversals[t.0].curve.project(point))
            .flat_map(|projection| {
                let measure = lrm.scale.locate_anchor(projection.distance_along_curve)?;
                Ok::<LrmProjection, LrsError>(LrmProjection {
                    measure: LrmMeasure {
                        lrm: lrm_handle,
                        measure,
                    },
                    orthogonal_offset: projection.offset,
                })
            })
            .collect();
        result.sort_by(|a, b| {
            a.orthogonal_offset
                .partial_cmp(&b.orthogonal_offset)
                .unwrap_or(Ordering::Equal)
        });
        result
    }

    fn lookup_lrms(&self, point: Point) -> Vec<LrmProjection> {
        let mut result: Vec<_> = self
            .lrms
            .iter()
            .enumerate()
            .flat_map(|(lrm_idx, _lrm)| self.lookup(point, LrmHandle(lrm_idx)))
            .collect();
        result.sort_by(|a, b| {
            a.orthogonal_offset
                .partial_cmp(&b.orthogonal_offset)
                .unwrap_or(Ordering::Equal)
        });
        result
    }

    fn lookup_traversals(&self, point: Point) -> Vec<TraversalProjection> {
        let mut result: Vec<_> = self
            .traversals
            .iter()
            .enumerate()
            .flat_map(|(idx, traversal)| traversal.curve.project(point).map(|proj| (idx, proj)))
            .map(|(idx, proj)| TraversalProjection {
                traversal: TraversalHandle(idx),
                orthogonal_offset: proj.offset,
                distance_from_start: proj.distance_along_curve,
            })
            .collect();
        result.sort_by(|a, b| {
            a.orthogonal_offset
                .partial_cmp(&b.orthogonal_offset)
                .unwrap_or(Ordering::Equal)
        });
        result
    }

    fn locate_traversal(&self, position: TraversalPosition) -> Result<Point, LrsError> {
        Ok(self
            .get_curve(position.traversal)?
            .resolve(position.into())?)
    }

    fn get_lrm_applicable_traversals(&self, lrm: LrmHandle) -> &[TraversalHandle] {
        &self.lrms[lrm.0].traversals
    }

    fn get_lrm_reference_traversal(&self, lrm: LrmHandle) -> TraversalHandle {
        self.lrms[lrm.0].reference_traversal
    }

    fn get_traversal_lrms(&self, traversal: TraversalHandle) -> &[LrmHandle] {
        &self.traversals[traversal.0].lrms
    }

    fn traversal_project(
        &self,
        position: TraversalPosition,
        onto: TraversalHandle,
    ) -> Result<TraversalProjection, LrsError> {
        let segment = self.orthogonal_segment(position.traversal, position.distance_from_start)?;
        let onto_curve = self.get_curve(onto)?;

        let point = onto_curve
            .intersect_segment(segment)
            .ok_or(CurveError::NotOnTheCurve)?;
        let projected = onto_curve.project(point)?;

        Ok(TraversalProjection {
            distance_from_start: projected.distance_along_curve,
            orthogonal_offset: projected.offset,
            traversal: onto,
        })
    }

    fn traversal_project_range(
        &self,
        range: TraversalRange,
        onto: TraversalHandle,
    ) -> Result<TraversalRange, LrsError> {
        let begin_pos = TraversalPosition {
            traversal: range.traversal,
            distance_from_start: range.begin,
        };
        let end_pos = TraversalPosition {
            traversal: range.traversal,
            distance_from_start: range.end,
        };
        let begin_projection = self.traversal_project(begin_pos, onto)?;
        let end_position = self.traversal_project(end_pos, onto)?;

        Ok(TraversalRange {
            traversal: onto,
            begin: begin_projection.distance_from_start,
            end: end_position.distance_from_start,
            direction: range.direction,
        })
    }

    fn lrm_project(
        &self,
        position: TraversalPosition,
        onto: LrmHandle,
    ) -> Result<LrmProjection, LrsError> {
        let lrm = self.lrms.get(onto.0).ok_or(LrsError::InvalidHandle)?;
        let measure = lrm.scale.locate_anchor(position.distance_from_start)?;
        Ok(LrmProjection {
            measure: LrmMeasure { lrm: onto, measure },
            orthogonal_offset: 0.,
        })
    }

    fn lrm_project_range(
        &self,
        range: TraversalRange,
        onto: LrmHandle,
    ) -> Result<LrmRange, LrsError> {
        let begin_pos = TraversalPosition {
            traversal: range.traversal,
            distance_from_start: range.begin,
        };

        let end_pos = TraversalPosition {
            traversal: range.traversal,
            distance_from_start: range.end,
        };

        let begin_projection = self.lrm_project(begin_pos, onto)?;
        let end_projection = self.lrm_project(end_pos, onto)?;

        Ok(LrmRange {
            lrm: onto,
            begin: begin_projection.measure.measure,
            end: end_projection.measure.measure,
            direction: range.direction,
        })
    }

    fn lrm_get_measure(&self, position: LrmPosition) -> Result<LrmMeasure, LrsError> {
        let scale = self.get_lrm_by_handle(position.lrm)?;
        Ok(LrmMeasure {
            measure: scale.get_measure(position.distance_from_start)?,
            lrm: position.lrm,
        })
    }

    fn lrm_get_position(&self, measure: LrmMeasure) -> Result<LrmPosition, LrsError> {
        let scale = self.get_lrm_by_handle(measure.lrm)?;
        Ok(LrmPosition {
            distance_from_start: scale.get_position(measure.measure)?,
            lrm: measure.lrm,
        })
    }

    fn get_linestring(&self, traversal: TraversalHandle) -> Result<LineString, LrsError> {
        self.get_curve(traversal).map(|c| c.as_linestring())
    }
}

impl<CurveImpl: Curve> Lrs<CurveImpl> {
    fn get_curve(&self, handle: TraversalHandle) -> Result<&CurveImpl, LrsError> {
        self.traversals
            .get(handle.0)
            .map(|traversal| &traversal.curve)
            .ok_or(LrsError::InvalidHandle)
    }

    fn get_lrm_by_handle(&self, handle: LrmHandle) -> Result<&LrmScale, LrsError> {
        self.lrms
            .get(handle.0)
            .map(|lrm| &lrm.scale)
            .ok_or(LrsError::InvalidHandle)
    }

    fn orthogonal_segment(
        &self,
        handle: TraversalHandle,
        from_start: CurvePosition,
    ) -> Result<geo::Line, LrsError> {
        let from_curve = self.get_curve(handle)?;
        let normal = from_curve.get_normal(from_start)?;

        let position = from_curve.resolve(CurveProjection {
            distance_along_curve: from_start,
            offset: 0.,
        })?;
        let start = geo::coord! {
            x: position.x() + normal.0 * from_curve.max_extent(),
            y: position.y() + normal.1 * from_curve.max_extent()
        };
        let end = geo::coord! {
            x: position.x() - normal.0 * from_curve.max_extent(),
            y: position.y() - normal.1 * from_curve.max_extent()
        };

        Ok(geo::Line::new(start, end))
    }
}

#[cfg(test)]
mod tests {
    use geo::line_string;

    use crate::curves::PlanarLineStringCurve;

    use super::*;

    fn lrs() -> Lrs<PlanarLineStringCurve> {
        let traversal = Traversal {
            curve: PlanarLineStringCurve::new(line_string![(x: 0., y:0.), (x: 200., y:0.)], 1.),
            id: "curve".to_owned(),
            lrms: vec![LrmHandle(0), LrmHandle(1)],
        };

        let traversal2 = Traversal {
            curve: PlanarLineStringCurve::new(line_string![(x: 0., y:-1.), (x: 200., y:-1.)], 1.),
            id: "curve".to_owned(),
            lrms: vec![LrmHandle(1)],
        };

        let lrm = Lrm {
            reference_traversal: TraversalHandle(0),
            scale: crate::lrm_scale::tests::scale(),
            traversals: vec![TraversalHandle(0)],
        };

        let mut lrm2 = Lrm {
            reference_traversal: TraversalHandle(0),
            scale: crate::lrm_scale::tests::scale(),
            traversals: vec![TraversalHandle(0), TraversalHandle(1)],
        };
        "id2".clone_into(&mut lrm2.scale.id);

        Lrs {
            lrms: vec![lrm, lrm2],
            traversals: vec![traversal, traversal2],
        }
    }

    #[test]
    fn get_lrm() {
        assert!(lrs().get_lrm("id").is_some());
        assert!(lrs().get_lrm("ideology").is_none());
    }

    #[test]
    fn get_traversal() {
        assert!(lrs().get_traversal("curve").is_some());
        assert!(lrs().get_traversal("Achtung, die Kurve!").is_none());
    }

    #[test]
    fn lookup_single_lrm() {
        let lrs = lrs();
        let result = lrs.lookup(point! {x: 50., y:0.5}, lrs.get_lrm("id").unwrap());
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].orthogonal_offset, 0.5);
        assert_eq!(result[0].measure.measure.scale_offset, 5.);
    }

    #[test]
    fn lookup_multiple_lrm() {
        let lrs = lrs();
        let result = lrs.lookup(point! {x: 50., y:0.5}, lrs.get_lrm("id2").unwrap());
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].orthogonal_offset, 0.5);
        assert_eq!(result[0].measure.measure.scale_offset, 5.);
        assert_eq!(result[0].measure.measure.anchor_name, "a");
        assert_eq!(result[1].orthogonal_offset, 1.5);
        assert_eq!(result[1].measure.measure.scale_offset, 5.);
    }

    #[test]
    fn lookup_lrms() {
        let result = lrs().lookup_lrms(point! {x: 50., y:0.5});
        assert_eq!(result.len(), 3);
        assert_eq!(result[0].orthogonal_offset, 0.5);
        assert_eq!(result[0].measure.measure.scale_offset, 5.);
        assert_eq!(result[1].orthogonal_offset, 0.5);
        assert_eq!(result[1].measure.measure.scale_offset, 5.);
        assert_eq!(result[2].orthogonal_offset, 1.5);
        assert_eq!(result[2].measure.measure.scale_offset, 5.);
    }

    #[test]
    fn lookup_traversals() {
        let result = lrs().lookup_traversals(point! {x: 50., y:0.5});
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn locate_traversal() {
        let result = lrs()
            .locate_traversal(TraversalPosition {
                distance_from_start: 10.,
                traversal: TraversalHandle(0),
            })
            .unwrap();
        assert_eq!(result, point! {x: 10., y: 0.});

        let result = lrs()
            .locate_traversal(TraversalPosition {
                distance_from_start: 10.,
                traversal: TraversalHandle(1),
            })
            .unwrap();
        assert_eq!(result, point! {x: 10., y: -1.});
    }

    #[test]
    fn get_lrm_applicable_traversals() {
        let lrs = lrs();
        let result = lrs.get_lrm_applicable_traversals(LrmHandle(0));
        assert_eq!(result.len(), 1);

        let result = lrs.get_lrm_applicable_traversals(LrmHandle(1));
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn get_lrm_reference_traversal() {
        let result = lrs().get_lrm_reference_traversal(LrmHandle(0));
        assert_eq!(TraversalHandle(0), result);
        let result = lrs().get_lrm_reference_traversal(LrmHandle(1));
        assert_eq!(TraversalHandle(0), result);
    }

    #[test]
    fn get_traversal_lrms() {
        let lrs = lrs();
        let result = lrs.get_traversal_lrms(TraversalHandle(0));
        assert_eq!(result, &[LrmHandle(0), LrmHandle(1)]);

        let result = lrs.get_traversal_lrms(TraversalHandle(1));
        assert_eq!(result, &[LrmHandle(1)]);
    }

    #[test]
    fn traversal_project() {
        let position = TraversalPosition {
            distance_from_start: 10.,
            traversal: TraversalHandle(0),
        };
        let result = lrs()
            .traversal_project(position, TraversalHandle(0))
            .unwrap();
        assert_eq!(result.distance_from_start, 10.);
    }

    #[test]
    fn traversal_project_range() {
        let range = TraversalRange {
            begin: 10.,
            end: 20.,
            direction: Direction::Default,
            traversal: TraversalHandle(0),
        };

        let result = lrs()
            .traversal_project_range(range, TraversalHandle(0))
            .unwrap();
        assert_eq!(result.begin, 10.);
        assert_eq!(result.end, 20.);
    }

    #[test]
    fn lrm_project() {
        let mut position = TraversalPosition {
            distance_from_start: 50.,
            traversal: TraversalHandle(0),
        };
        let result = lrs().lrm_project(position, LrmHandle(0)).unwrap();
        assert_eq!(result.orthogonal_offset, 0.);
        assert_eq!(result.measure.measure.scale_offset, 5.);
        assert_eq!(result.measure.measure.anchor_name, "a");

        position.distance_from_start = 130.;
        let result = lrs().lrm_project(position, LrmHandle(0)).unwrap();
        assert_eq!(result.orthogonal_offset, 0.);
        assert_eq!(result.measure.measure.scale_offset, 3.);
        assert_eq!(result.measure.measure.anchor_name, "b");
    }

    #[test]
    fn lrm_project_range() {
        let range = TraversalRange {
            begin: 50.,
            end: 130.,
            direction: Direction::Default,
            traversal: TraversalHandle(0),
        };
        let result = lrs().lrm_project_range(range, LrmHandle(0)).unwrap();
        assert_eq!(result.begin.scale_offset, 5.);
        assert_eq!(result.begin.anchor_name, "a");
        assert_eq!(result.end.scale_offset, 3.);
        assert_eq!(result.end.anchor_name, "b");
    }

    #[test]
    fn lrm_get_measure() {
        let position = LrmPosition {
            distance_from_start: 5.,
            lrm: LrmHandle(0),
        };
        let result = lrs().lrm_get_measure(position).unwrap();
        assert_eq!(result.measure.anchor_name, "a");
        assert_eq!(result.measure.scale_offset, 5.);

        let position = LrmPosition {
            distance_from_start: 25.,
            lrm: LrmHandle(0),
        };
        let result = lrs().lrm_get_measure(position).unwrap();
        assert_eq!(result.measure.anchor_name, "b");
        assert_eq!(result.measure.scale_offset, 15.);
    }

    #[test]
    fn lrm_get_position() {
        let measure = LrmMeasure {
            lrm: LrmHandle(0),
            measure: LrmScaleMeasure {
                anchor_name: "a".to_owned(),
                scale_offset: 5.,
            },
        };
        let result = lrs().lrm_get_position(measure).unwrap();

        assert_eq!(result.distance_from_start, 5.);
    }
}
