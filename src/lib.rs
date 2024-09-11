#[allow(unused_imports)]
#[allow(clippy::all)]
#[allow(dead_code)]
#[rustfmt::skip]
mod lrs_generated;

#[deny(missing_docs)]
mod osm_helpers;

#[deny(missing_docs)]
pub mod curves;
#[deny(missing_docs)]
pub mod lrm_scale;
#[deny(missing_docs)]
pub mod lrs;
#[deny(missing_docs)]
pub mod lrs_ext;

#[deny(missing_docs)]
pub mod builder;

#[test]
fn read_and_write_lrs() {
    use builder::*;
    use curves::SphericalLineStringCurve;
    use geo::Coord;

    let mut builder = Builder::new();
    let anchor_index = builder.add_anchor(
        "Ancre",
        Some("12"),
        Coord { x: 0., y: 0. },
        properties!("some key" => "some value"),
    );
    let start_node = builder.add_node("a", Coord { x: 0., y: 0. }, properties!());
    let end_node = builder.add_node("b", Coord { x: 1., y: 1. }, properties!());
    let segment_geometry = &[Coord { x: 0., y: 0. }, Coord { x: 1., y: 1. }];
    let segment = SegmentOfTraversal {
        segment_index: builder.add_segment("segment", segment_geometry, start_node, end_node),
        reversed: false,
    };
    let traversal = builder.add_traversal("traversal", &[segment]);
    let anchor_on_lrm = AnchorOnLrm {
        anchor_index,
        distance_along_lrm: 12.0,
    };
    builder.add_lrm("lrm", traversal, &[anchor_on_lrm], properties!());

    let buffer = builder.build_data(properties!("source" => "example"));
    let lrs = lrs::Lrs::<SphericalLineStringCurve>::from_bytes(buffer).unwrap();
    let anchor = &lrs.lrms[0].scale.anchors[0];
    assert_eq!(anchor.id.as_ref().unwrap(), "12");
}
