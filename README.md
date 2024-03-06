# Library for Linear Reference System

TODO: explain what this library is for

## Using liblrs

### Data serialization with FlatBuffers

The data that defines an LRS is serialized using the [FlatBuffers format](https://flatbuffers.dev/).

The schema is described in [schema/lrs.fbs](schema/lrs.fbs). The library is written in rust and the [generated file](src/lrs_generated.rs) is commited. This means there is no need to have the `flatc` executable to build and run this project.

If your contribution changes the schema, you will need to generate the file with flatc. The version must be the release 23.5.26. Do not use a version built from master.

`flatc -o src --rust schema/lrs.fbs`

## Norms

### Comment convention

See [How to write documentation in Rust](https://doc.rust-lang.org/rustdoc/how-to-write-documentation.html) to keep the code clean and clear (also [this](https://github.com/rust-lang/rfcs/blob/master/text/1574-more-api-documentation-conventions.md#appendix-a-full-conventions-text) for other examples).