use std::path::PathBuf;

use clap::Parser;

use liblrs::{builder::Builder, properties};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
/// Arguments given by the command line interface.
struct Args {
    /// OpenStreetMap file to parse.
    #[arg(short, long)]
    input_osm_file: PathBuf,

    /// Output file where the [`Lrs`] will be written.
    #[arg(short, long)]
    output_lrs: PathBuf,

    /// OpenStreetMap tag identifying the LRM. The french railway network uses `ref:FR:SNCF_Reseau`.
    #[arg(short, long)]
    lrm_tag: String,
}

/// Example: to generate an LRS from an OpenStreetMap dump
///
/// `$ cargo run --release --bin geometry_from_osm -- -i france.osm.pbf  -o osm_83000.lrs.bin --lrm-tag=ref:fr:SNCF_Reseau`
fn main() {
    let cli_args = Args::parse();

    let required = properties!("railway" => "rail");
    let to_reject = properties!(
        "service"=> "siding",
        "service"=> "spur",
        "building"=> "*",
        "area"=> "yes",
        "gauge"=> "600",
        "roller_coaster"=> "*",
        "construction"=> "*"
    );

    let mut builder = Builder::new();
    builder.read_from_osm(
        &cli_args.input_osm_file,
        &cli_args.lrm_tag,
        required,
        to_reject,
    );

    builder.save(
        &cli_args.output_lrs,
        properties!("source" => "OpenStreetMap", "licence" => "OdBL"),
    );
}
