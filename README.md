# Library for Linear Reference System

[Linear Reference](https://en.wikipedia.org/wiki/Linear_referencing) allows to identify a position along a curve (a road, a canal, a railroad…) relative to fixed reference points. 

We call:
- **anchor** a reference point: it can be milestones, a landmark, an intersection…
- **scale** the list of _anchors_ and the theoretical distance between them,
- **curve** the physical description of the object to reference,
- **linear referencing method** (LRM) the combination of a _scale_ and a _curve_,
- **linear referencing system** (LRS) the complete set of data.


While the logic is quite simple, multiple small subtleties make the LRS difficult to use.
Distances between milestone change (construction of a bypass around a town), the origin of curve is displaced (the railwail station moved), there is only one scale for a river and its parallel running canal…

This library aims to handle many edge cases and makes little assumptions about the data:
- anchors are not always numbers,
- anchors don’t need to be on the curve,
- distance between anchors are not fixed,
- distance between anchor might not match the measured distance,
- works on spherical and projected coordinates,
- a single scale can be used for many curves.

## Bindings and HTML demonstrator

The core library is written in rust. We expose javascript binding through [WebAssembly](https://webassembly.org/). Those binding can be built in the `wasm` directory.

```
cd wasm
npm install
npm build
```

### Demonstrator

A simple HTML demonstrator allows to test the data and the functions:

```
cd wasm
npm run serve
```

And open your browser at `http://localhost:8080`

You can customize the map background if you provide your own [maplibre](https://maplibre.org/) style:

```
MAPLIBRE_STYLE="https://your_tile_provider/style.json?key=42" npm run serve
```

## Using liblrs

### Data serialization with FlatBuffers

The data that defines an LRS is serialized using the [FlatBuffers format](https://flatbuffers.dev/).

The schema is described in [schema/lrs.fbs](schema/lrs.fbs). The library is written in rust and the [generated file](src/lrs_generated.rs) is commited. This means there is no need to have the `flatc` executable to build and run this project.

If your contribution changes the schema, you will need to generate the file with flatc. The version must be the release 23.5.26. Do not use a version built from master.

`flatc -o src --rust schema/lrs.fbs`

## Norms

### Comment convention

See [How to write documentation in Rust](https://doc.rust-lang.org/rustdoc/how-to-write-documentation.html) to keep the code clean and clear (also [this](https://github.com/rust-lang/rfcs/blob/master/text/1574-more-api-documentation-conventions.md#appendix-a-full-conventions-text) for other examples).

### Extracting geometry from OpenStreetMap

We provide a binary that extracts geometry data from OpenStreetMap and saves it in the FlatBuffer format.

For now we only handle Railway data. A tag that describes the LRM must be provided. In France, use the [ref:FR:SNCF_Reseau](https://wiki.openstreetmap.org/wiki/FR:Key:ref:FR:SNCF_Reseau).

Run the the binary:

`cargo run --release --bin geometry_from_osm -- -i france.rail.osm.pbf  -o osm.lrs.bin2 --lrm-tag=ref:FR:SNCF_Reseau`
