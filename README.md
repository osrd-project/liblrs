# Library for Linear Reference System

TODO: explain what this library is for

## Using liblrs

### Data serialization with FlatBuffers

The data that defines an LRS is serialized using the [FlatBuffers format](https://flatbuffers.dev/).

The schema is described in [schema/lrs.fbs](schema/lrs.fbs). The library is written in rust and the [generated file](src/lrs_generated.rs) is commited. This means there is no need to have the `flatc` executable to build and run this project.

If your contribution changes the schema, you will need to generate the file with

`flatc -o src --rust schema/lrs.fbs`