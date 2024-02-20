extern crate flatbuffers;

<<<<<<< HEAD
#[allow(unused_imports)]
#[allow(clippy::all)]
#[rustfmt::skip]
mod lrs_generated;

=======
mod lrs_generated;

use flatbuffers::WIPOffset;
>>>>>>> 82e23bd (Initial use of the flatbuffer in ruste)
pub use lrs_generated::*;

#[test]
fn read_and_write_lrs() {
    let mut fbb = flatbuffers::FlatBufferBuilder::with_capacity(1024);
    let property = PropertyArgs {
        key: Some(fbb.create_string("some key")),
        value: Some(fbb.create_string("some value")),
    };
    let properties = &[Property::create(&mut fbb, &property)];
    let anchor_arg = AnchorArgs {
        id: Some(fbb.create_string("Ancre")),
        name: Some(fbb.create_string("12")),
        ..Default::default()
    };

    let anchor = Anchor::create(&mut fbb, &anchor_arg);
    let lrs_args = LrsArgs {
        properties: Some(fbb.create_vector(properties)),
        anchors: Some(fbb.create_vector(&[anchor])),
<<<<<<< HEAD
        ..Default::default()
=======
        networks: Some(fbb.create_vector::<WIPOffset<Network>>(&[])),
        linear_referencing_methods: Some(
            fbb.create_vector::<WIPOffset<LinearReferencingMethod>>(&[]),
        ),
        views: Some(fbb.create_vector::<WIPOffset<GeometryView>>(&[])),
>>>>>>> 82e23bd (Initial use of the flatbuffer in ruste)
    };
    let lrs = Lrs::create(&mut fbb, &lrs_args);

    fbb.finish(lrs, None);
    let buffer = fbb.finished_data();

    let read_rls = root_as_lrs(buffer).unwrap();
<<<<<<< HEAD
    assert_eq!(read_rls.anchors().unwrap().get(0).id(), "Ancre");
    assert!(read_rls.networks().is_none());
=======
    assert_eq!(read_rls.anchors().get(0).id(), "Ancre");
    assert_eq!(read_rls.networks().len(), 0);
>>>>>>> 82e23bd (Initial use of the flatbuffer in ruste)
}
