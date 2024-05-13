extern crate flatbuffers;

#[allow(unused_imports)]
#[allow(clippy::all)]
#[rustfmt::skip]
mod lrs_generated;

#[deny(missing_docs)]
pub mod curves;
#[deny(missing_docs)]
pub mod lrm_scale;
#[deny(missing_docs)]
pub mod lrs;
#[deny(missing_docs)]
pub mod lrs_ext;
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
        ..Default::default()
    };
    let lrs = Lrs::create(&mut fbb, &lrs_args);

    fbb.finish(lrs, None);
    let buffer = fbb.finished_data();

    let read_rls = root_as_lrs(buffer).unwrap();
    assert_eq!(read_rls.anchors().unwrap().get(0).id(), "Ancre");
    assert!(read_rls.networks().is_none());
}
