use std::collections::HashMap;

use clap::Parser;
use liblrs::*;

fn read_osm(input_file: &str, lrm_tag: &str) -> (Vec<osm4routing::Node>, Vec<osm4routing::Edge>) {
    osm4routing::Reader::new()
        .merge_ways()
        .require("railway", "rail")
        .reject("service", "yard")
        .reject("service", "siding")
        .reject("service", "spur")
        .reject("building", "*")
        .reject("area", "yes")
        .reject("gauge", "600")
        .reject("roller_coaster", "*")
        .reject("construction", "*")
        .read_tag(lrm_tag)
        .read(input_file)
        .expect("could not read the osm file")
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// OpenStreetMap file to parse
    #[arg(short, long)]
    input_osm_file: String,

    /// Output file where the LRS will be written
    #[arg(short, long)]
    output_lrs: String,

    /// OpenStreetMap tag identifying the lrm. For the french railway network use ref:FR:SNCF_Reseau
    #[arg(short, long)]
    lrm_tag: String,
}

fn main() {
    let cli_args = Args::parse();
    let (nodes, edges) = read_osm(&cli_args.input_osm_file, &cli_args.lrm_tag);
    println!(
        "In OpenStreetMap, we have {} nodes and {} edges",
        nodes.len(),
        edges.len()
    );

    let mut nodes_map = HashMap::<&osm4routing::NodeId, Vec<_>>::new();
    let mut traversal_indices = HashMap::<String, Vec<_>>::new();
    let mut traversal_directions = HashMap::<String, Vec<_>>::new();

    let mut fbb = flatbuffers::FlatBufferBuilder::with_capacity(1024);

    for (edge_index, edge) in edges.iter().enumerate() {
        let edge_index = edge_index as u64;
        let mut args = ConnectionArgs {
            endpoint: Some(Endpoint::Begin),
            segment_index: edge_index,
            ..Default::default()
        };
        let connection = Connection::create(&mut fbb, &args);
        nodes_map.entry(&edge.source).or_default().push(connection);

        args.endpoint = Some(Endpoint::End);
        let connection = Connection::create(&mut fbb, &args);
        nodes_map.entry(&edge.target).or_default().push(connection);

        if let Some(srv_ref) = edge.tags.get(&cli_args.lrm_tag) {
            traversal_indices
                .entry(srv_ref.clone())
                .or_default()
                .push(edge_index);

            // TODO compute the actual direction
            traversal_directions
                .entry(srv_ref.clone())
                .or_default()
                .push(Direction::Increasing);
        }
    }

    let traversals: Vec<_> = traversal_indices
        .into_iter()
        .map(|(lrm_ref, segments)| {
            let directions = fbb.create_vector(&traversal_directions.remove(&lrm_ref).unwrap());
            let args = TraversalArgs {
                id: Some(fbb.create_string(&lrm_ref)),
                directions: Some(directions),
                segments: Some(fbb.create_vector(&segments)),
                properties: None,
            };
            Traversal::create(&mut fbb, &args)
        })
        .collect();

    println!("In the LRS we have {} traversals", traversals.len());

    let nodes: Vec<_> = nodes
        .iter()
        .map(|n| {
            let args = NodeArgs {
                id: Some(fbb.create_string(&n.id.0.to_string())),
                properties: None,
                connections: Some(fbb.create_vector(&nodes_map.remove(&n.id).unwrap_or_default())),
            };
            Node::create(&mut fbb, &args)
        })
        .collect();

    let segments: Vec<_> = edges
        .iter()
        .map(|e| {
            let args = SegmentArgs {
                id: Some(fbb.create_string(&e.id)),
                properties: None,
            };
            Segment::create(&mut fbb, &args)
        })
        .collect();

    let network_args = NetworkArgs {
        id: Some(fbb.create_string("OpenStreetMap")),
        nodes: Some(fbb.create_vector(&nodes)),
        segments: Some(fbb.create_vector(&segments)),
        traversals: Some(fbb.create_vector_from_iter(traversals.into_iter())),
    };
    let network = Network::create(&mut fbb, &network_args);

    let segments_geometry: Vec<_> = edges
        .iter()
        .map(|e| {
            let points_iter = e.geometry.iter().map(|c| Point::new(c.lon, c.lat, 0.0));
            let points = Some(fbb.create_vector_from_iter(points_iter));
            SegmentGeometry::create(&mut fbb, &SegmentGeometryArgs { points })
        })
        .collect();

    let network_geometry_args = NetworkGeometryArgs {
        segments: Some(fbb.create_vector(&segments_geometry)),
    };

    let network_geometry = NetworkGeometry::create(&mut fbb, &network_geometry_args);

    let geometry_view_args = GeometryViewArgs {
        geometry_type: Some(GeometryType::Geographic),
        networks: Some(fbb.create_vector(&[network_geometry])),
        ..Default::default()
    };

    let view = GeometryView::create(&mut fbb, &geometry_view_args);

    let key = Some(fbb.create_string("source"));
    let value = Some(fbb.create_string("OpenStreetMap"));
    let source = Property::create(&mut fbb, &PropertyArgs { key, value });

    let lrs_args = LrsArgs {
        properties: Some(fbb.create_vector(&[source])),
        networks: Some(fbb.create_vector(&[network])),
        views: Some(fbb.create_vector(&[view])),
        ..Default::default()
    };
    let lrs = Lrs::create(&mut fbb, &lrs_args);

    fbb.finish(lrs, None);
    let buffer = fbb.finished_data();

    std::fs::write(&cli_args.output_lrs, buffer).unwrap();
}
