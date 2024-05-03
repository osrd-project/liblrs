//! This module defines a curve and primitive functions on them
//! Most common manipulations are projecting a point on the curve
//! and otherway round find the coordinates of a point along the curve
//! For now the implementation is based on a [`LineString`],
//! but other implementations could be considered such as splines.

use geo::kernels::RobustKernel;
use geo::prelude::*;
use geo::{coord, Line, LineString, Point, Rect};
use thiserror::Error;

/// A [`Curve`] is the fundamental building block for an LRM.
/// It provides basic primitives to locate/project points on it.
/// A [`Curve`] can be part of a larger [`Curve`] (e.g. for optimisation purposes and to have better bounding boxes).
/// The [`Curve`] can be implemented.
pub trait Curve {
    /// Builds a new [`Curve`] from a [`LineString`].
    /// `max_extent` is the maximum distance that is considered to be “on the curve”.
    /// `max_extent` plays a role in the bounding box.
    fn new(geom: LineString, max_extent: f64) -> Self;

    /// The length of the [`Curve`].
    fn length(&self) -> f64;

    /// How far from the [`Curve`] could be considered to be still on the [`Curve`].
    fn max_extent(&self) -> f64;

    /// Is the geometry valid. Depending on the representation.
    /// It must have at least two [`Coord`]s.
    /// If there are exactly two [`Coord`]s, they must be different.
    fn is_valid(&self) -> bool;

    /// Projects the [`Point`] to the closest position on the [`Curve`].
    /// Will fail if the [`Curve`] is invalid (e.g. no [`Point`] on it)
    /// or if the [`Point`] is too far away.
    /// If the [`Curve`] is a piece of a larger [`Curve`] (`start_offset > 0`)
    /// then the `distance_along_curve` if from the whole [`Curve`], not just the current piece.
    fn project(&self, point: Point) -> Result<CurveProjection, CurveError>;

    /// Returns the geographical position of a [`Point`] on the [`Curve`].
    /// Will return an error if the [`CurveProjection`] is not on this [`Curve`].
    fn resolve(&self, projection: CurveProjection) -> Result<Point, CurveError>;

    /// Bounding box of the [`Curve`] with a buffer of `max_extent`.
    fn bbox(&self) -> Rect;

    /// Computes the normal at a given offset on the [`Curve`].
    /// Will return an error if the [`Curve`] is invalid or the `offset` is outside of the [`Curve`].
    /// Points to the positive side (left).
    fn get_normal(&self, offset: f64) -> Result<(f64, f64), CurveError>;

    /// Returns the [`Point`] where the [`Curve`] and the segment ([`Line`]) intersect.
    /// If the segment intersects the [`Curve`] multiple times, an intersection is chosen randomly.
    /// TODO: implement return of all the points of intersection.
    /// When the segment is collinear with the [`Curve`] it is ignored.
    fn intersect_segment(&self, segment: Line) -> Option<Point>;
}

/// Errors when manipulating the [`Curve`]s.
#[derive(Error, Debug, PartialEq)]
pub enum CurveError {
    /// The condition of validity might differ depending on the [`Curve`] implementation.
    #[error("the curve geometry is not valid")]
    InvalidGeometry,
    /// At least one coordinate is non a finite number (`NaN`, `infinite`).
    #[error("the coordinates are not finite")]
    NotFiniteCoordinates,
    /// The considered [`Point`] is not on the [`Curve`].
    #[error("the point is not on the curve")]
    NotOnTheCurve,
}

/// Implementation based on [`LineString`]:
/// the [`Curve`] is a string of continous [`Line`]s.
/// Each [`Line`] made up of 2 [`Coord`]s.
/// This implementation doesn't take in account the ellipsoidal model of the earth.
/// The coordinates are reprensented by `f64`.
/// That means a precison of about 1_000_000th of a mm for a [`Curve`] that spans around the Earth.
pub struct PlanarLineStringCurve {
    /// When a [`Curve`] might be a piece of a longer [`Curve`]
    /// then the `start_offset` allows to know how far along the longer [`Curve`] we are.
    pub start_offset: f64,

    /// The max distance that is considered of being part of the [`Curve`].
    /// It is used to compute the bounding box.
    pub max_extent: f64,

    /// The coordinates are considered to be planar.
    /// All distance and length calculations are expressed in the same units as coordinates.
    pub geom: LineString,

    length: f64,
}

impl PlanarLineStringCurve {
    /// Splits the [`LineString`] into smaller [`Curve`]s of at most `max_len` length.
    /// If the initial geometry is invalid, it returns an empty vector.
    pub fn new_fragmented(geom: LineString, max_len: f64, max_extent: f64) -> Vec<Self> {
        let n = (geom.euclidean_length() / max_len).ceil() as usize;
        geom.line_segmentize(n)
            .map(|multi| {
                multi
                    .0
                    .into_iter()
                    .map(|geom| Self::new(geom, max_extent))
                    .collect()
            })
            .unwrap_or_default()
    }
}

impl Curve for PlanarLineStringCurve {
    fn new(geom: LineString, max_extent: f64) -> Self {
        let length = geom.euclidean_length();
        Self {
            start_offset: 0.,
            max_extent,
            geom,
            length,
        }
    }

    fn length(&self) -> f64 {
        self.length
    }

    fn max_extent(&self) -> f64 {
        self.max_extent
    }

    fn is_valid(&self) -> bool {
        self.geom.coords_count() >= 2 && (self.geom.coords_count() > 2 || !self.geom.is_closed())
    }

    fn project(&self, point: Point) -> Result<CurveProjection, CurveError> {
        if !self.is_valid() {
            return Err(CurveError::InvalidGeometry);
        }

        match self.geom.line_locate_point(&point) {
            Some(location) => {
                let distance_along_curve = location * self.length + self.start_offset;

                let begin = self.geom.coords().next().unwrap();
                let end = self.geom.coords().next_back().unwrap();

                let sign = match RobustKernel::orient2d(point.into(), *end, *begin) {
                    Orientation::Clockwise => 1.,
                    _ => -1.,
                };
                let offset = point.euclidean_distance(&self.geom) * sign;

                Ok(CurveProjection {
                    distance_along_curve,
                    offset,
                })
            }
            None => Err(CurveError::NotFiniteCoordinates),
        }
    }

    fn resolve(&self, projection: CurveProjection) -> Result<Point, CurveError> {
        let fraction = (projection.distance_along_curve - self.start_offset) / self.length;
        if !(0. ..=1.).contains(&fraction) || fraction.is_nan() {
            Err(CurveError::NotOnTheCurve)
        } else {
            Ok(self.geom.line_interpolate_point(fraction).unwrap())
        }
    }

    fn bbox(&self) -> Rect {
        let bounding_rect = self.geom.bounding_rect().unwrap();
        Rect::new(
            coord! {
                x: bounding_rect.min().x - self.max_extent,
                y: bounding_rect.min().y - self.max_extent,
            },
            coord! {
                x: bounding_rect.max().x + self.max_extent,
                y: bounding_rect.max().y + self.max_extent,
            },
        )
    }

    fn intersect_segment(&self, segment: Line) -> Option<Point> {
        self.geom
            .lines()
            .flat_map(|curve_line| {
                match geo::line_intersection::line_intersection(segment, curve_line) {
                    Some(LineIntersection::SinglePoint {
                        intersection,
                        is_proper: _,
                    }) => Some(intersection.into()),
                    Some(LineIntersection::Collinear { intersection: _ }) => None,
                    None => None,
                }
            })
            .next()
    }

    fn get_normal(&self, offset: f64) -> Result<(f64, f64), CurveError> {
        // We find the Point where the normal is computed
        let point = self.resolve(CurveProjection {
            distance_along_curve: offset,
            offset: 0.,
        })?;

        let line = self
            .geom
            .lines_iter()
            .find(|line| line.contains(&point))
            .ok_or(CurveError::NotFiniteCoordinates)?;

        // translate to (0, 0) and normalize by the length of the curve to get unit vector of tangent
        let tangent = (
            (line.end.x - line.start.x) / self.length,
            (line.end.y - line.start.y) / self.length,
        );

        // 90° clockwise rotation
        let normal = (-tangent.1, tangent.0);

        Ok(normal)
    }
}

/// Implementation based on [`LineString`]:
/// the [`Curve`] is a string of continous [`Line`]s.
/// Each [`Line`] made up of 2 [`Coord`]s.
/// A spherical coordinates model takes 3 coordinates: longitude (`x`), latitude (`y`)
/// and radius as const ([`geo::MEAN_EARTH_RADIUS`]).
/// [`GeodesicLength`] for a Paris to New-York is about 5837.283 km, and
/// [`HaversineLength`] for a Paris to New-York is about 5852.970 km.
/// We've chosen here to stick with Haversine Formula.
/// The computations are made along the Great Circle.
/// When required, some methods use [geo::algorithm::densify_haversine::DensifyHaversine]
/// to get an approximation of the great circle by generating more [`Line`]s on the [`LineString`].
/// The coordinates are reprensented by `f64`.
/// That means a precison of about 1_000_000th of a mm for a [`Curve`] that spans around the Earth.
pub struct SphericalLineStringCurve {
    /// When a [`Curve`] might be a piece of a longer [`Curve`]
    /// then the `start_offset` allows to know how far along the longer [`Curve`] we are.
    pub start_offset: f64,

    /// The max distance that is considered of being part of the [`Curve`].
    /// It is used to compute the bounding box.
    pub max_extent: f64,

    /// The coordinates are considered to be spherical.
    /// All distance and length calculations are expressed in the same units as coordinates.
    pub geom: LineString,

    /// In meters.
    length: f64,

    /// In meters. Represents the minimum length by which the curve can be densified.
    pub densify_by: f64,
}

impl SphericalLineStringCurve {
    /// Splits the [`LineString`] into smaller [`Curve`]s of at most `max_len` length.
    /// If the initial geometry is invalid, it returns an empty vector.
    pub fn new_fragmented(geom: LineString, max_len: f64, max_extent: f64) -> Vec<Self> {
        let n = (geom.haversine_length() / max_len).ceil() as usize;

        geom.line_segmentize_haversine(n)
            .map(|multi| {
                multi
                    .0
                    .into_iter()
                    .map(|geom| Self::new(geom, max_extent))
                    .collect()
            })
            .unwrap_or_default()
    }
}

impl Curve for SphericalLineStringCurve {
    fn new(geom: LineString, max_extent: f64) -> Self {
        let length = geom.haversine_length();
        Self {
            start_offset: 0.,
            max_extent,
            geom,
            length,
            densify_by: 100., // arbitrary, maximum length of a curve will be 100m, otherwise it will be densified
        }
    }

    fn length(&self) -> f64 {
        self.length
    }

    fn max_extent(&self) -> f64 {
        self.max_extent
    }

    fn is_valid(&self) -> bool {
        if !(self.geom.coords_count() >= 2
            && (self.geom.coords_count() > 2 || !self.geom.is_closed()))
        {
            return false;
        }
        for coord in self.geom.coords() {
            if coord.x < -180.0 || coord.x > 180.0 || coord.y < -90.0 || coord.y > 90.0 {
                return false;
            }
        }
        true
    }

    fn project(&self, point: Point) -> Result<CurveProjection, CurveError> {
        if !self.is_valid() {
            return Err(CurveError::InvalidGeometry);
        }

        match self.geom.haversine_closest_point(&point) {
            geo::Closest::SinglePoint(closest_point) => {
                let distance_along_curve = closest_point
                    .haversine_distance(&self.geom.points().next().unwrap())
                    + self.start_offset;

                let begin = self.geom.coords().next().unwrap();
                let end = self.geom.coords().next_back().unwrap();

                let sign = match RobustKernel::orient2d(point.into(), *end, *begin) {
                    Orientation::Clockwise => 1.,
                    _ => -1.,
                };
                let offset = closest_point.haversine_distance(&point) * sign;

                Ok(CurveProjection {
                    distance_along_curve,
                    offset,
                })
            }
            geo::Closest::Intersection(_) => Err(CurveError::InvalidGeometry),
            geo::Closest::Indeterminate => Err(CurveError::NotFiniteCoordinates),
        }
    }

    fn resolve(&self, projection: CurveProjection) -> Result<Point, CurveError> {
        let fraction = (projection.distance_along_curve - self.start_offset) / self.length;
        if !(0. ..=1.).contains(&fraction) || fraction.is_nan() {
            return Err(CurveError::NotOnTheCurve);
        }

        let total_length = self.length;
        let fractional_length = total_length * fraction;
        let mut accumulated_length = 0.;

        // go through each segment and look for the one that frame the distance that we seek
        for segment in self.geom.lines() {
            let segment_length = segment.haversine_length();
            if accumulated_length + segment_length >= fractional_length {
                let segment_fraction = (fractional_length - accumulated_length) / segment_length;
                let segment_start = Point::from(segment.start);
                let segment_end = Point::from(segment.end);

                // get the Point at a haversine distance between two Points
                // of a certain fraction of length on this distance
                return Ok(segment_start.haversine_intermediate(&segment_end, segment_fraction));
            }
            accumulated_length += segment_length;
        }

        Err(CurveError::NotOnTheCurve)
    }

    fn bbox(&self) -> Rect {
        let bounding_rect = self.geom.bounding_rect().unwrap();
        Rect::new(
            coord! {
                x: bounding_rect.min().x - self.max_extent,
                y: bounding_rect.min().y - self.max_extent,
            },
            coord! {
                x: bounding_rect.max().x + self.max_extent,
                y: bounding_rect.max().y + self.max_extent,
            },
        )
    }

    // Important:
    // - the SphericalLineStringCurve is densified for long curves
    // to get the intersection(s) closer to the real closest path.
    fn intersect_segment(&self, segment: Line) -> Option<Point> {
        self.geom
            .densify_haversine(self.densify_by)
            .lines()
            .flat_map(|curve_line| {
                match geo::line_intersection::line_intersection(segment, curve_line) {
                    Some(LineIntersection::SinglePoint {
                        intersection,
                        is_proper: _,
                    }) => Some(intersection.into()),
                    Some(LineIntersection::Collinear { intersection: _ }) => None,
                    None => None,
                }
            })
            .next()
    }

    // Important :
    // - the output normal vector is normalized considering haversine formula,
    // and thus can be only used where it has been computed.
    // - the SphericalLineStringCurve is densified for long curves
    // to get the intersection(s) closer to the real closest path.
    fn get_normal(&self, offset: f64) -> Result<(f64, f64), CurveError> {
        // go through each segment and look for the one that frame the distance that we seek
        let mut accumulated_length = 0.;
        for segment in self.geom.densify_haversine(self.densify_by).lines() {
            let segment_length = segment.haversine_length();
            if accumulated_length + segment_length >= offset {
                // get Points from the segment that frame the point
                let start = geo::Point(segment.start);
                let end = geo::Point(segment.end);

                // get bearing from start Point and end Point of the segment, and add 90° clockwise rotation to it
                let normal_vector_bearing = start.haversine_bearing(end) + 90.;

                // get end Point from the end of the segment for the normal vector bearing value and 1m of (haversine) length
                let end_normal = end.haversine_destination(normal_vector_bearing, 1.);

                return Ok((end_normal.x() - end.x(), end_normal.y() - end.y()));
            }
            accumulated_length += segment_length;
        }
        Err(CurveError::NotFiniteCoordinates)
    }
}

/// Represents a [`Point`] in space projected on the [`Curve`].
#[derive(Clone, Copy)]
pub struct CurveProjection {
    /// How far from the [`Curve`] start is located the [`Point`].
    /// If the [`Curve`] is part of a larger [`Curve`], `start_offset` is strictly positive
    /// and the `start_offset` will be considered.
    pub distance_along_curve: f64,

    /// How far is the [`Point`] from the [`Curve`]:
    /// - [`geo::algorithm::euclidean_distance`] for [`PlanarLineStringCurve`]
    /// - [`geo::algorithm::haversine_distance`] for [`SphericalLineStringCurve`]
    ///
    /// The distance is `positive` if the [`Point`] is located on the left of the [`Curve`]
    /// and `negative` if the [`Point`] is on the right.
    pub offset: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    pub use geo::line_string;
    use geo::point;

    const PARIS_LON: f64 = 2.352565660016694;
    const PARIS_LAT: f64 = 48.85643268390663;
    const NEW_YORK_LON: f64 = -74.00599134051316;
    const NEW_YORK_LAT: f64 = 40.71274961837565;
    const LILLE_LON: f64 = 3.066667;
    const LILLE_LAT: f64 = 50.633333;
    const PERPIGNAN_LON: f64 = 2.8948332;
    const PERPIGNAN_LAT: f64 = 42.6886591;
    const BREST_LON: f64 = -4.486076;
    const BREST_LAT: f64 = 48.390394;
    const NANCY_LON: f64 = 6.184417;
    const NANCY_LAT: f64 = 48.692054;
    const REYKJAVIK_LON: f64 = -21.827774;
    const REYKJAVIK_LAT: f64 = 64.128288;

    #[test]
    fn planar_fragmented() {
        let framentation_max_length = 1.;
        let c = PlanarLineStringCurve::new_fragmented(
            line_string![(x: 0., y: 0.), (x: 2., y: 0.)],
            framentation_max_length,
            1.,
        );
        assert_eq!(2, c.len());
        assert_eq!(framentation_max_length, c[0].length());
    }

    #[test]
    fn planar_length() {
        let c = PlanarLineStringCurve::new(line_string![(x: 0., y: 0.), (x: 2., y: 0.)], 1.);
        assert_eq!(2., c.length());
    }

    #[test]
    fn planar_is_invalid() {
        // Valid curve
        let curve = PlanarLineStringCurve::new(line_string![(x: 0., y: 0.), (x: 1., y: 1.)], 1.);
        assert!(curve.is_valid());

        // Invalid curve: too few coordinates
        let curve = PlanarLineStringCurve::new(line_string![(x: 0., y: 0.)], 1.);
        assert!(!curve.is_valid());

        // Invalid curve: closed LineString with only 2 coordinates
        let curve = PlanarLineStringCurve::new(line_string![(x: 0., y: 0.), (x: 0., y: 0.)], 1.);
        assert!(!curve.is_valid());
    }

    #[test]
    fn planar_projection() {
        let mut c = PlanarLineStringCurve::new(line_string![(x: 0., y: 0.), (x: 2., y: 0.)], 1.);

        let projected = c.project(point! {x: 1., y: 1.}).unwrap();
        assert_eq!(1., projected.distance_along_curve);
        assert_eq!(1., projected.offset);

        let projected = c.project(point! {x: 1., y: -1.}).unwrap();
        assert_eq!(1., projected.distance_along_curve);
        assert_eq!(-1., projected.offset);

        c.start_offset = 1.;
        let projected = c.project(point! {x: 1., y: -1.}).unwrap();
        assert_eq!(2., projected.distance_along_curve);
    }

    #[test]
    fn planar_resolve() {
        let mut c = PlanarLineStringCurve::new(line_string![(x: 0., y: 0.), (x: 2., y: 0.)], 1.);

        let projection = CurveProjection {
            distance_along_curve: 1.,
            offset: 0.,
        };
        let p = c.resolve(projection).unwrap();
        assert_eq!(p.x(), 1.);
        assert_eq!(p.y(), 0.);

        c.start_offset = 1.;
        let p = c.resolve(projection).unwrap();
        assert_eq!(p.x(), 0.);
    }

    #[test]
    fn planar_bbox() {
        let c = PlanarLineStringCurve::new(line_string![(x: 0., y: 0.), (x: 2., y: 0.)], 1.);
        let bbox = c.bbox();

        assert_eq!(bbox.min(), coord! {x: -1., y: -1.});
        assert_eq!(bbox.max(), coord! {x: 3., y: 1.});
    }

    #[test]
    fn planar_intersect_segment() {
        // Right angle
        let c = PlanarLineStringCurve::new(line_string![(x: 0., y: 0.), (x: 2., y:0.)], 1.);
        let segment = Line::new(coord! {x: 1., y: 1.}, coord! {x: 1., y: -1.});
        assert_eq!(c.intersect_segment(segment), Some(point! {x: 1., y: 0.}));

        // No intersection
        let segment = Line::new(coord! {x: 10., y: 10.}, coord! {x:20., y: 10.});
        assert!(c.intersect_segment(segment).is_none());

        // Collinear
        let segment = Line::new(coord! {x: 0., y: 0.,}, coord! {x: 1., y:0.});
        assert!(c.intersect_segment(segment).is_none());

        // Multiple intersection
        let c = PlanarLineStringCurve::new(
            line_string![(x: 0., y: 0.), (x: 1., y: 2.), (x: 2., y: 0.)],
            1.,
        );
        let segment = Line::new(coord! {x: 0., y: 1.}, coord! {x: 2., y: 1.});
        assert!(c.intersect_segment(segment).is_some());
    }

    #[test]
    fn planar_normal() {
        let c = PlanarLineStringCurve::new(line_string![(x: 0., y: 0.), (x: 2., y: 0.)], 1.);
        let normal_c = c.get_normal(1.).unwrap();
        assert_relative_eq!(normal_c.0, 0.);
        assert_relative_eq!(normal_c.1, 1.);
    }

    #[test]
    fn spherical_fragmented() {
        let framentation_max_length = 1.;
        let paris_to_new_york = SphericalLineStringCurve::new_fragmented(
            line_string![(x: PARIS_LON, y: PARIS_LAT), (x: NEW_YORK_LON, y: NEW_YORK_LAT)],
            framentation_max_length,
            1.,
        );

        assert_eq!(5837284, paris_to_new_york.len());
        assert_relative_eq!(
            framentation_max_length,
            paris_to_new_york[0].length(),
            epsilon = 1e-7
        );
        // 1e-7 means we lose 0.1 micrometer per segment
    }

    #[test]
    fn spherical_length() {
        let paris_to_new_york = SphericalLineStringCurve::new(
            line_string![(x: PARIS_LON, y: PARIS_LAT), (x: NEW_YORK_LON, y: NEW_YORK_LAT)],
            1.,
        );
        assert_relative_eq!(5837283.441678336, paris_to_new_york.length()); // 5837283.441678336 using [`HaversineLength`] and 5852969.839293494 using [`GeodesicLength`]

        let lille_to_perpignan = SphericalLineStringCurve::new(
            line_string![(x: LILLE_LON, y: LILLE_LAT), (x: PERPIGNAN_LON, y: PERPIGNAN_LAT)],
            1.,
        );
        assert_relative_eq!(883505.2931188548, lille_to_perpignan.length()); // 883505.2931188548 using [`HaversineLength`] and 883260.051153502 using [`GeodesicLength`]

        let brest_to_nancy = SphericalLineStringCurve::new(
            line_string![(x: BREST_LON, y: BREST_LAT), (x: NANCY_LON, y: NANCY_LAT)],
            1.,
        );
        assert_relative_eq!(785636.8730262491, brest_to_nancy.length()); // 785636.8730262491 using [`HaversineLength`] and 787994.4363866252 using [`GeodesicLength`]
    }

    #[test]
    fn spherical_is_valid() {
        // Valid curve
        let curve = SphericalLineStringCurve::new(line_string![(x: 0., y: 0.), (x: 1., y: 1.)], 1.);
        assert!(curve.is_valid());

        // Invalid curve: too few coordinates
        let curve = SphericalLineStringCurve::new(line_string![(x: 0., y: 0.)], 1.);
        assert!(!curve.is_valid());

        // Invalid curve: closed LineString with only 2 coordinates
        let curve = SphericalLineStringCurve::new(line_string![(x: 0., y: 0.), (x: 0., y: 0.)], 1.);
        assert!(!curve.is_valid());

        // Invalid curve: longitude > 180.
        let curve =
            SphericalLineStringCurve::new(line_string![(x: 180.1, y: 0.), (x: 0., y: 0.)], 1.);
        assert!(!curve.is_valid());

        // Invalid curve: longitude < -180.
        let curve =
            SphericalLineStringCurve::new(line_string![(x: -180.1, y: 0.), (x: 0., y: 0.)], 1.);
        assert!(!curve.is_valid());

        // Invalid curve: latitude > 90.
        let curve =
            SphericalLineStringCurve::new(line_string![(x: 0., y: 90.1), (x: 0., y: 0.)], 1.);
        assert!(!curve.is_valid());

        // Invalid curve: latitude > 180.
        let curve =
            SphericalLineStringCurve::new(line_string![(x: 0., y: -90.1), (x: 0., y: 0.)], 1.);
        assert!(!curve.is_valid());
    }

    #[test]
    fn spherical_projection() {
        let mut paris_to_new_york = SphericalLineStringCurve::new(
            line_string![(x: PARIS_LON, y: PARIS_LAT), (x: NEW_YORK_LON, y: NEW_YORK_LAT)],
            1.,
        );

        // Point is located on the right (north) of the curve
        let projected = paris_to_new_york
            .project(point! {x: -6.705403880820967, y: 51.42135181702875})
            .unwrap();
        assert_eq!(701924.3809693493, projected.distance_along_curve);
        assert_eq!(-67157.93913531031, projected.offset);

        // Point is located on the left (south) of the curve
        let projected = paris_to_new_york
            .project(point! {x: -12.250890759346419, y: 45.857650969554356})
            .unwrap();
        assert_eq!(963365.3768036617, projected.distance_along_curve);
        assert_eq!(625592.3211438804, projected.offset);

        // Same point, but with an offset from the curve
        paris_to_new_york.start_offset = 1000000.;
        let projected = paris_to_new_york
            .project(point! {x: -12.250890759346419, y: 45.857650969554356})
            .unwrap();
        assert_eq!(1963365.3768036617, projected.distance_along_curve);
        assert_eq!(625592.3211438804, projected.offset);

        // ################################################################################
        let mut new_york_to_paris = SphericalLineStringCurve::new(
            line_string![(x: NEW_YORK_LON, y: NEW_YORK_LAT), (x: PARIS_LON, y: PARIS_LAT)],
            1.,
        );

        // Point is located on the left (north) of the curve
        let projected = new_york_to_paris
            .project(point! {x: -6.705403880820967, y: 51.42135181702875})
            .unwrap();
        assert_eq!(5135359.060708988, projected.distance_along_curve);
        assert_eq!(67157.93913531031, projected.offset);

        // Point is located on the right (south)  of the curve
        let projected = new_york_to_paris
            .project(point! {x: -12.250890759346419, y: 45.857650969554356})
            .unwrap();
        assert_eq!(4873918.064874676, projected.distance_along_curve);
        assert_eq!(-625592.3211438811, projected.offset); // Note: result is weird -> distance should remain the same than the other way curve, difference is 0.7mm

        // Same point, but with an offset from the curve
        new_york_to_paris.start_offset = 1000000.;
        let projected = new_york_to_paris
            .project(point! {x: -12.250890759346419, y: 45.857650969554356})
            .unwrap();
        assert_eq!(5873918.064874676, projected.distance_along_curve);
        assert_eq!(-625592.3211438811, projected.offset); // Note: same, difference is 0.7mm
    }

    #[test]
    fn spherical_resolve() {
        let mut paris_to_new_york = SphericalLineStringCurve::new(
            line_string![(x: PARIS_LON, y: PARIS_LAT), (x: NEW_YORK_LON, y: NEW_YORK_LAT)],
            1.,
        );

        let mut projection = CurveProjection {
            distance_along_curve: 1000000.,
            offset: 0.,
        };
        let paris_to_new_york_projection = paris_to_new_york.resolve(projection).unwrap();
        assert_eq!(paris_to_new_york_projection.x(), -11.113419713640527);
        assert_eq!(paris_to_new_york_projection.y(), 51.44320774918762);

        paris_to_new_york.start_offset = 300000.;
        let paris_to_new_york_projection = paris_to_new_york.resolve(projection).unwrap();
        assert_eq!(paris_to_new_york_projection.x(), -6.924269520648392);
        assert_eq!(paris_to_new_york_projection.y(), 50.83295845015173);

        projection.distance_along_curve = 8000000.;
        assert!(paris_to_new_york.resolve(projection).is_err());

        // ################################################################################
        let lille_to_perpignan = SphericalLineStringCurve::new(
            line_string![(x: LILLE_LON, y: LILLE_LAT), (x: PERPIGNAN_LON, y: PERPIGNAN_LAT)],
            1.,
        );

        let projection = CurveProjection {
            distance_along_curve: 500000.,
            offset: 0.,
        };
        let lille_to_perpignan_p = lille_to_perpignan.resolve(projection).unwrap();
        assert_eq!(lille_to_perpignan_p.x(), 2.963285977058639);
        assert_eq!(lille_to_perpignan_p.y(), 46.13725407237963);

        // ################################################################################
        let brest_to_nancy = SphericalLineStringCurve::new(
            line_string![(x: BREST_LON, y: BREST_LAT), (x: NANCY_LON, y: NANCY_LAT)],
            1.,
        );

        let projection = CurveProjection {
            distance_along_curve: 500000.,
            offset: 0.,
        };
        let brest_to_nancy_p = brest_to_nancy.resolve(projection).unwrap();
        assert_eq!(brest_to_nancy_p.x(), 2.292338879356053);
        assert_eq!(brest_to_nancy_p.y(), 48.69669359753767);
    }

    #[test]
    fn spherical_bbox() {
        let paris_to_new_york = SphericalLineStringCurve::new(
            line_string![(x: PARIS_LON, y: PARIS_LAT), (x: NEW_YORK_LON, y: NEW_YORK_LAT)],
            1.,
        );
        let bbox = paris_to_new_york.bbox();

        assert_eq!(
            bbox.min(),
            coord! {x: -75.00599134051316, y: 39.71274961837565}
        );
        assert_eq!(
            bbox.max(),
            coord! {x: 3.352565660016694, y: 49.85643268390663}
        );
    }

    #[test]
    fn spherical_intersect_segment() {
        // Note: following tests have been computed with a maximum length of curve of 100m, otherwise the curve is densified.

        // Intersection
        let paris_to_new_york = SphericalLineStringCurve::new(
            line_string![(x: PARIS_LON, y: PARIS_LAT), (x: NEW_YORK_LON, y: NEW_YORK_LAT)],
            1.,
        );
        let segment = Line::new(
            coord! {x: -36.76627263796084, y: 69.72980545457074},
            coord! {x: -53.52127629098692, y: 15.34337895024332},
        );
        assert_eq!(
            paris_to_new_york.intersect_segment(segment),
            Some(point! {x: -42.500669938830555, y: 51.11605974559634})
        );

        // No intersection
        let segment = Line::new(
            coord! {x: -88.45243862592235, y: 20.758717928501483},
            coord! {x:19.035989490700018, y: 41.32134615429521},
        );
        assert!(paris_to_new_york.intersect_segment(segment).is_none());

        // TODO: Collinear
        // Notes:
        // - because of the haversine densification, the geometry is slightly different and includes more points
        // than before, thus creating intersection(s) point(s).
        // - is very rare in reality

        // Multiple intersection
        let paris_to_reykjavik_to_new_york = SphericalLineStringCurve::new(
            line_string![(x: PARIS_LON, y: PARIS_LAT), (x: REYKJAVIK_LON, y: REYKJAVIK_LAT), (x: NEW_YORK_LON, y: NEW_YORK_LAT)],
            1.,
        );

        let segment = Line::new(
            coord! {x: -70.77775907909825, y: 47.835409180411006},
            coord! {x: 9.293636086504506, y: 54.83039737996501},
        );
        assert!(paris_to_reykjavik_to_new_york
            .intersect_segment(segment)
            .is_some());
    }

    #[test]
    fn spherical_normal() {
        // Earth radius is equal to 6371008.8m, considering geo::MEAN_EARTH_RADIUS.
        let earth_circumference = 6371008.8 * std::f64::consts::PI * 2.; // = 40030228.88407185 m
        let normalized_translation_on_earth = 360. / earth_circumference;

        let longitudinal_curve =
            SphericalLineStringCurve::new(line_string![(x: 0., y: 0.), (x: 1., y: 0.)], 1.);
        let longitudinal_normal = longitudinal_curve.get_normal(0.).unwrap();
        assert_relative_eq!(longitudinal_normal.0, 0., epsilon = 1e-7);
        assert_relative_eq!(
            longitudinal_normal.1,
            normalized_translation_on_earth,
            epsilon = 1e-4
        );

        let latitudinal_curve =
            SphericalLineStringCurve::new(line_string![(x: 0., y: 0.), (x: 0., y: 1.)], 1.);
        let latitudinal_normal = latitudinal_curve.get_normal(0.).unwrap();
        assert_relative_eq!(
            latitudinal_normal.0,
            normalized_translation_on_earth,
            epsilon = 1e-8
        );
        assert_relative_eq!(latitudinal_normal.1, 0., epsilon = 1e-7);
    }
}
