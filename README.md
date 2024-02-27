# Library for Linear Reference System

TODO: explain what this library is for

## Using liblrs

### Data serialization with FlatBuffers

The data that defines an LRS is serialized using the [FlatBuffers format](https://flatbuffers.dev/).

The schema is described in [schema/lrs.fbs](schema/lrs.fbs). The library is written in rust and the [generated file](src/lrs_generated.rs) is commited. This means there is no need to have the `flatc` executable to build and run this project.

If your contribution changes the schema, you will need to generate the file with flatc. The version must be the release 23.5.26. Do not use a version built from master.

`flatc -o src --rust schema/lrs.fbs`

### Extracting geometry from OpenStreetMap

We provide a binary that extracts geometry data from OpenStreetMap and saves it in the FlatBuffer format.

For now we only handle Railway data. A tag that describes the LRM must be provided. In France, use the [ref:FR:SNCF_Reseau](https://wiki.openstreetmap.org/wiki/FR:Key:ref:FR:SNCF_Reseau).

Run the the binary:

`cargo run --release --bin geometry_from_osm -- -i france.rail.osm.pbf  -o osm.lrs.bin2 --lrm-tag=ref:FR:SNCF_Reseau`