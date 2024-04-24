//! A LRM (linear reference model) is an abstract representation
//! where the geometry and real distances are not considered.

use thiserror::Error;

/// Measurement along the `Curve`. Typically in meters.
pub type CurvePosition = f64;

/// Measurement along the `LrmScale`. Often in meters, but it could be anything.
pub type ScalePosition = f64;

/// Errors when manipulating a `LrmScale`.
#[derive(Error, Debug, PartialEq)]
pub enum LrmScaleError {
    /// Returned when building a `LrmScale` from a builder and less than 2 [NamedAnchor] objects were provided.
    #[error("a scale needs at least two named anchor")]
    NoEnoughNamedAnchor,
    /// All the [NamedAnchor] objects must be unique within a same `LrmScale`.
    #[error("duplicated anchor: {0}")]
    DuplicatedAnchorName(String),
    /// Could not find the position on the `Curve` as the [Anchor] is not known.
    #[error("anchor is unknown in the LrmScale")]
    UnknownAnchorName,
    /// Could not find an [Anchor] that matches a given offset.
    #[error("no anchor found")]
    NoAnchorFound,
}

/// An `Anchor` is a reference point that is well known from which the location is computed.
#[derive(PartialEq, Debug)]
pub struct Anchor {
    /// Some `Anchor` objects might not be named,
    /// e.g. the first anchor of the LRM.
    pub id: Option<String>,

    /// Distance from the start of the scale in the scale space, can be negative.
    pub scale_position: ScalePosition,

    /// Real distance from the start of the `Curve`.
    /// The `Curve` might not start at the same 0 (e.g. the `Curve` is longer than the scale),
    /// or the `Curve` might not progress at the same rate (e.g. the `Curve` is a schematic representation that distorts distances).
    pub curve_position: CurvePosition,
}

impl Anchor {
    /// Builds a named `Anchor`.
    pub fn new(name: &str, scale_position: ScalePosition, curve_position: CurvePosition) -> Self {
        Self {
            id: Some(name.to_owned()),
            scale_position,
            curve_position,
        }
    }

    /// Builds an unnamed `Anchor`.
    pub fn new_unnamed(scale_position: ScalePosition, curve_position: CurvePosition) -> Self {
        Self {
            id: None,
            scale_position,
            curve_position,
        }
    }

    fn as_named(&self) -> Option<NamedAnchor> {
        self.id.as_ref().map(|id| NamedAnchor {
            id: id.to_owned(),
            scale_position: self.scale_position,
            curve_position: self.curve_position,
        })
    }
}

// Private struct to be used when we only deal with Anchor that has name.
struct NamedAnchor {
    id: String,
    scale_position: ScalePosition,
    curve_position: CurvePosition,
}

/// An helper to build a scale by adding consecutive [Anchor] objects with relative distances.
/// When having all the `Anchor` objects and their distances in both scale and real position,
/// it is simpler to directly build the [LrmScale] from an `Vec<Anchor>`.
/// Using the builder will however ensure that the scale is valid.
pub struct ScaleBuilder {
    anchors: Vec<Anchor>,
}

impl ScaleBuilder {
    /// Creates a new scale with an initial [Anchor].
    pub fn new(anchor: Anchor) -> Self {
        Self {
            anchors: vec![anchor],
        }
    }

    /// Builds a named [Anchor] and adds it to the `ScaleBuilder`.
    /// Distances are relative to previous `Anchor`.
    pub fn add_named(self, id: &str, scale_dist: ScalePosition, curve_dist: CurvePosition) -> Self {
        self.add(Some(id.to_owned()), scale_dist, curve_dist)
    }

    /// Builds an unnamed [Anchor] and adds it to the `ScaleBuilder`.
    /// Distances are relative to previous `Anchor`.
    pub fn add_unnamed(self, scale_dist: ScalePosition, curve_dist: CurvePosition) -> Self {
        self.add(None, scale_dist, curve_dist)
    }

    /// Builds an [Anchor] and adds it to the `ScaleBuilder`.
    /// Distances are relative to previous `Anchor`.
    pub fn add(
        mut self,
        id: Option<String>,
        scale_dist: ScalePosition,
        curve_dist: CurvePosition,
    ) -> Self {
        let last_anchor = self
            .anchors
            .last()
            .expect("The builder should have at least one anchor");

        self.anchors.push(Anchor {
            id,
            scale_position: last_anchor.scale_position + scale_dist,
            curve_position: last_anchor.curve_position + curve_dist,
        });
        self
    }

    /// Requires at least one named [Anchor].
    /// Will fail if none is present and if the `Anchor` names are duplicated.
    /// This will consume the `ScaleBuilder` that can not be used after.
    pub fn build(self, id: &str) -> Result<LrmScale, LrmScaleError> {
        let mut names = std::collections::HashSet::new();
        for anchor in self.anchors.iter() {
            if let Some(name) = &anchor.id {
                if !names.insert(name) {
                    return Err(LrmScaleError::DuplicatedAnchorName(name.to_string()));
                }
            }
        }

        if names.is_empty() {
            Err(LrmScaleError::NoEnoughNamedAnchor)
        } else {
            Ok(LrmScale {
                id: id.to_owned(),
                anchors: self.anchors,
            })
        }
    }
}

/// A measure defines a location on the [LrmScale].
/// It is given as an [Anchor] name and an `offset` on that scale.
/// It is often represented as `12+100` to say `“100 scale units after the Anchor 12`”.
#[derive(Clone, Debug)]
pub struct LrmScaleMeasure {
    /// `Name` of the [Anchor]. While it is often named after a kilometer position,
    /// it can be anything (a letter, a landmark).
    pub anchor_name: String,
    /// The `offset` from the anchor in the scale units.
    /// there is no guarantee that its value matches actual distance on the `Curve` and is defined in scale units.
    pub scale_offset: ScalePosition,
}

impl LrmScaleMeasure {
    /// Builds a new `LrmMeasure` from an [Anchor] `name` and the `offset` on the [LrmScale].
    pub fn new(anchor_name: &str, scale_offset: ScalePosition) -> Self {
        Self {
            anchor_name: anchor_name.to_owned(),
            scale_offset,
        }
    }
}

/// Represents an `LrmScale` and allows to map [Measure] to a position along a `Curve`.
#[derive(PartialEq, Debug)]
pub struct LrmScale {
    /// Unique identifier.
    pub id: String,
    /// The [Anchor] objects are reference points on the scale from which relative distances are used.
    pub anchors: Vec<Anchor>,
}

impl LrmScale {
    /// Locates a point along a `Curve` given an [Anchor] and an `offset`,
    /// which might be negative.
    pub fn locate_point(&self, measure: &LrmScaleMeasure) -> Result<CurvePosition, LrmScaleError> {
        let named_anchor = self
            .iter_named()
            .find(|anchor| anchor.id == measure.anchor_name)
            .ok_or(LrmScaleError::UnknownAnchorName)?;
        let nearest_anchor = self
            .next_anchor(&named_anchor.id)
            .or_else(|| self.previous_anchor(&named_anchor.id))
            .ok_or(LrmScaleError::NoAnchorFound)?;

        let scale_interval = named_anchor.scale_position - nearest_anchor.scale_position;
        let curve_interval = named_anchor.curve_position - nearest_anchor.curve_position;
        Ok(named_anchor.curve_position + curve_interval * measure.scale_offset / scale_interval)
    }

    /// Returns a measure given a distance along the `Curve`.
    /// The corresponding [Anchor] is the named `Anchor` that gives the smallest positive `offset`.
    /// If such an `Anchor` does not exists, the first named `Anchor` is used.
    pub fn locate_anchor(
        &self,
        curve_position: CurvePosition,
    ) -> Result<LrmScaleMeasure, LrmScaleError> {
        // First, we find the nearest named Anchor to the Curve.
        let named_anchor = self
            .nearest_named(curve_position)
            .ok_or(LrmScaleError::NoAnchorFound)?;

        // Then we search the nearest Anchor that will be the reference
        // to convert from Curve units to scale units.
        let nearest_anchor = if named_anchor.curve_position < curve_position {
            self.next_anchor(&named_anchor.id)
                .or(self.previous_anchor(&named_anchor.id))
        } else {
            self.previous_anchor(&named_anchor.id)
                .or(self.next_anchor(&named_anchor.id))
        }
        .ok_or(LrmScaleError::NoAnchorFound)?;

        let ratio = (nearest_anchor.scale_position - named_anchor.scale_position)
            / (nearest_anchor.curve_position - named_anchor.curve_position);

        Ok(LrmScaleMeasure {
            anchor_name: named_anchor.id,
            scale_offset: (curve_position - named_anchor.curve_position) * ratio,
        })
    }

    /// Returns a measure given a distance along the `LrmScale`.
    /// The corresponding [Anchor] is the named `Anchor` that gives the smallest positive `offset`.
    /// If such an `Anchor` does not exists, the first named `Anchor` is used.
    pub fn get_measure(
        &self,
        scale_position: ScalePosition,
    ) -> Result<LrmScaleMeasure, LrmScaleError> {
        let named_anchor = self
            .scale_nearest_named(scale_position)
            .ok_or(LrmScaleError::NoAnchorFound)?;

        Ok(LrmScaleMeasure {
            anchor_name: named_anchor.id,
            scale_offset: scale_position - named_anchor.scale_position,
        })
    }

    /// Locates a point along the scale given an [Anchor] and an `offset`,
    /// which might be negative.
    pub fn get_position(&self, measure: LrmScaleMeasure) -> Result<ScalePosition, LrmScaleError> {
        let named_anchor = self
            .iter_named()
            .find(|anchor| anchor.id == measure.anchor_name)
            .ok_or(LrmScaleError::UnknownAnchorName)?;

        Ok(named_anchor.scale_position + measure.scale_offset)
    }

    fn nearest_named(&self, curve_position: CurvePosition) -> Option<NamedAnchor> {
        // Tries to find the Anchor whose curve_position is the biggest possible, yet smaller than Curve position
        // Otherwise take the first named
        // Anchor names   ----A----B----
        // Curve positions    2    3
        // With Curve position = 2.1, we want A
        // With Curve position = 2.9, we want A
        //                       3.5, we want B
        //                       1.5, we want A
        self.iter_named()
            .rev()
            .find(|anchor| anchor.curve_position <= curve_position)
            .or_else(|| self.iter_named().next())
    }

    fn scale_nearest_named(&self, scale_position: ScalePosition) -> Option<NamedAnchor> {
        // Like nearest_named, but our position is along the scale
        self.iter_named()
            .rev()
            .find(|anchor| anchor.scale_position <= scale_position)
            .or_else(|| self.iter_named().next())
    }

    // Finds the closest Anchor before the Anchor having the name `name`
    fn previous_anchor(&self, name: &str) -> Option<&Anchor> {
        self.anchors
            .iter()
            .rev()
            .skip_while(|anchor| anchor.id.as_deref() != Some(name))
            .nth(1)
    }

    // Finds the closest Anchor after the Anchor having the name `name`
    fn next_anchor(&self, name: &str) -> Option<&Anchor> {
        self.anchors
            .iter()
            .skip_while(|anchor| anchor.id.as_deref() != Some(name))
            .nth(1)
    }

    // Iterates only on named Anchor objects
    fn iter_named(&self) -> impl DoubleEndedIterator<Item = NamedAnchor> + '_ {
        self.anchors.iter().filter_map(|anchor| anchor.as_named())
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;
    pub fn scale() -> LrmScale {
        ScaleBuilder::new(Anchor::new("a", 0., 0.))
            .add_named("b", 10., 100.)
            .build("id")
            .unwrap()
    }

    #[test]
    fn builder() {
        // Everything as planed
        assert_eq!(scale().anchors[0].curve_position, 0.);
        assert_eq!(scale().anchors[1].curve_position, 100.);

        // Missing named Anchor
        let b = ScaleBuilder::new(Anchor::new_unnamed(0., 0.));
        let scale = b.build("id");
        assert_eq!(scale, Err(LrmScaleError::NoEnoughNamedAnchor));

        // Duplicated names
        let scale = ScaleBuilder::new(Anchor::new("a", 0., 0.))
            .add_named("a", 100., 100.)
            .build("id");
        assert_eq!(
            scale,
            Err(LrmScaleError::DuplicatedAnchorName("a".to_string()))
        );
    }

    #[test]
    fn locate_point() {
        // Everything a usual
        assert_eq!(
            scale().locate_point(&LrmScaleMeasure::new("a", 5.)),
            Ok(50.)
        );
        assert_eq!(
            scale().locate_point(&LrmScaleMeasure::new("b", 5.)),
            Ok(150.)
        );

        // Negative offsets
        assert_eq!(
            scale().locate_point(&LrmScaleMeasure::new("a", -5.)),
            Ok(-50.)
        );

        // Unknown Anchor
        assert_eq!(
            scale().locate_point(&LrmScaleMeasure::new("c", 5.)),
            Err(LrmScaleError::UnknownAnchorName)
        );
    }

    #[test]
    fn nearest_named() {
        let scale = ScaleBuilder::new(Anchor::new("a", 0., 2.))
            .add_named("b", 10., 1.)
            .build("id")
            .unwrap();

        assert_eq!(scale.nearest_named(2.1).unwrap().id, "a");
        assert_eq!(scale.nearest_named(2.9).unwrap().id, "a");
        assert_eq!(scale.nearest_named(1.5).unwrap().id, "a");
        assert_eq!(scale.nearest_named(3.5).unwrap().id, "b");
    }

    #[test]
    fn locate_anchor() {
        let measure = scale().locate_anchor(40.).unwrap();
        assert_eq!(measure.anchor_name, "a");
        assert_eq!(measure.scale_offset, 4.);

        let measure = scale().locate_anchor(150.).unwrap();
        assert_eq!(measure.anchor_name, "b");
        assert_eq!(measure.scale_offset, 5.);

        let measure = scale().locate_anchor(-10.).unwrap();
        assert_eq!(measure.anchor_name, "a");
        assert_eq!(measure.scale_offset, -1.);
    }

    #[test]
    fn locate_anchor_with_unnamed() {
        // ----Unnamed(100)----A(200)----B(300)----Unnamed(400)---
        let scale = ScaleBuilder::new(Anchor::new_unnamed(0., 100.))
            .add_named("a", 1., 100.)
            .add_named("b", 1., 100.)
            .add_unnamed(1., 100.)
            .build("id")
            .unwrap();

        // Unnamed----position----Named
        let measure = scale.locate_anchor(150.).unwrap();
        assert_eq!(measure.anchor_name, "a");
        assert_eq!(measure.scale_offset, -0.5);

        // position----Unnamed----Named
        let measure = scale.locate_anchor(50.).unwrap();
        assert_eq!(measure.anchor_name, "a");
        assert_eq!(measure.scale_offset, -1.5);

        // Unnamed----Named----position----Unnamed
        let measure = scale.locate_anchor(350.).unwrap();
        assert_eq!(measure.anchor_name, "b");
        assert_eq!(measure.scale_offset, 0.5);

        // Unnamed----Named----Unnamed----position
        let measure = scale.locate_anchor(500.).unwrap();
        assert_eq!(measure.anchor_name, "b");
        assert_eq!(measure.scale_offset, 2.);
    }

    #[test]
    fn get_measure() {
        // a(scale 0)----measure(scale 5)----b(scale 10)
        let measure = scale().get_measure(5.).unwrap();
        assert_eq!(measure.anchor_name, "a");
        assert_eq!(measure.scale_offset, 5.);

        // a(scale 0)----b(scale 10)----measure(scale 25)
        let measure = scale().get_measure(25.).unwrap();
        assert_eq!(measure.anchor_name, "b");
        assert_eq!(measure.scale_offset, 15.);
    }

    #[test]
    fn get_position() {
        // a(scale 0)----position(scale a+5)----b(scale 10)
        let position = scale()
            .get_position(LrmScaleMeasure {
                anchor_name: "a".to_string(),
                scale_offset: 5.,
            })
            .unwrap();
        assert_eq!(position, 5.);

        // a(scale 0)----b(scale 10)----position(scale b+15)
        let position = scale()
            .get_position(LrmScaleMeasure {
                anchor_name: "b".to_string(),
                scale_offset: 15.,
            })
            .unwrap();
        assert_eq!(position, 25.);
    }
}
