use std::collections::HashMap;

use clap::Parser;
use liblrs::*;
use osm4routing::Edge;

fn read_osm(input_file: &str, lrm_tag: &str) -> (Vec<osm4routing::Node>, Vec<osm4routing::Edge>) {
    osm4routing::Reader::new()
        .merge_ways()
        .require("railway", "rail")
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
/// Arguments given by the command line interface.
struct Args {
    /// OpenStreetMap file to parse.
    #[arg(short, long)]
    input_osm_file: String,

    /// Output file where the [`Lrs`] will be written.
    #[arg(short, long)]
    output_lrs: String,

    /// OpenStreetMap tag identifying the LRM. The french railway network uses `ref:FR:SNCF_Reseau`.
    #[arg(short, long)]
    lrm_tag: String,
}

// When sorting the edges, each candidate is tested to see if they match and if they need to be reversed.
#[derive(PartialEq, Eq, Debug)]
enum Candidate {
    Source,
    Target,
    Impossible,
}

// If we have a string of edges that ends with `end_edge`,
// the candidate `Edge(s, t)` can be joined at the `end_edge(source, target)`.
// Returns Forward if the `Edge` can be append in the same direction `(target == s)`,
// Backward if it must be reversed `(target == t)` or
// And NotCandidate if it can’t be appended.
// `consider_source` means that the source (or target if false) of the `last_edge` is considered.
fn can_be_appended(candidate: &Edge, last_edge: &Edge, consider_source: bool) -> Candidate {
    let last_node = if consider_source {
        last_edge.source
    } else {
        last_edge.target
    };
    if candidate.source == last_node {
        Candidate::Source
    } else if candidate.target == last_node {
        Candidate::Target
    } else {
        Candidate::Impossible
    }
}

fn insert(
    mut to_insert: Vec<Edge>,
    position: usize,
    reversed: bool,
    at_start: bool,
    mut sorted: Vec<(Edge, bool)>,
) -> (Vec<Edge>, Vec<(Edge, bool)>) {
    let to_insert_value = (to_insert.remove(position), reversed);
    if at_start {
        sorted.insert(0, to_insert_value);
    } else {
        sorted.push(to_insert_value)
    }
    sort_iteration(to_insert, sorted)
}

fn sort_iteration(
    to_insert: Vec<Edge>,
    sorted: Vec<(Edge, bool)>,
) -> (Vec<Edge>, Vec<(Edge, bool)>) {
    if sorted.is_empty() {
        insert(to_insert, 0, false, false, vec![])
    } else {
        for i in 0..to_insert.len() {
            let (begin_edge, begin_direction) = sorted.first().unwrap();
            let (end_edge, end_direction) = sorted.last().unwrap();

            let at_begin = can_be_appended(&to_insert[i], begin_edge, !begin_direction);
            let at_end = can_be_appended(&to_insert[i], end_edge, *end_direction);

            match (at_begin, at_end) {
                (Candidate::Target, _) => return insert(to_insert, i, false, true, sorted),
                (Candidate::Source, _) => return insert(to_insert, i, true, true, sorted),
                (_, Candidate::Source) => return insert(to_insert, i, false, false, sorted),
                (_, Candidate::Target) => return insert(to_insert, i, true, false, sorted),
                (Candidate::Impossible, Candidate::Impossible) => continue,
            }
        }
        (to_insert, sorted)
    }
}

fn sort_edges(edges: Vec<Edge>, traversal_ref: &str) -> Vec<(Edge, bool)> {
    let (to_insert, sorted) = sort_iteration(edges, vec![]);

    // Print some stats about edges that could not be matched
    if !to_insert.is_empty() {
        let ignored = to_insert.len();
        let total = ignored + sorted.len();
        let first = if sorted[0].1 {
            sorted[0].0.target
        } else {
            sorted[0].0.source
        };

        let last_edge = sorted.last().unwrap();
        let last = if last_edge.1 {
            last_edge.0.source
        } else {
            last_edge.0.target
        };
        println!("[WARN] on traversal {traversal_ref}, ignoring {ignored} edges out of {total}");
        println!(
            "       Sorted traversal from osm node: {} to: {}",
            first.0, last.0
        );
    }

    sorted
}

/// Example: to generate an LRS from an OpenStreetMap dump
///
/// `$ cargo run --release --bin geometry_from_osm -- -i france.osm.pbf  -o osm_83000.lrs.bin --lrm-tag=ref:fr:SNCF_Reseau`
fn main() {
    let cli_args = Args::parse();
    let (nodes, edges) = read_osm(&cli_args.input_osm_file, &cli_args.lrm_tag);
    let (nodes_len, edges_len) = (nodes.len(), edges.len());
    println!("In OpenStreetMap, we have {nodes_len} nodes and {edges_len} edges.");

    let mut edges_map = HashMap::<_, _>::new();

    let mut fbb = flatbuffers::FlatBufferBuilder::with_capacity(1024);

    let mut traversals = HashMap::<String, Vec<_>>::new();

    for (edge_index, edge) in edges.iter().enumerate() {
        if let Some(srv_ref) = edge.tags.get(&cli_args.lrm_tag) {
            edges_map.insert(edge.id.clone(), edge_index as u64);
            traversals
                .entry(srv_ref.clone())
                .or_default()
                .push(edge.clone());
        }
    }

    // Sort the traversals
    let mut fb_traversals = Vec::new();
    for (srv_ref, edges) in traversals.into_iter() {
        let segments =
            fbb.create_vector_from_iter(sort_edges(edges.clone(), &srv_ref).into_iter().map(
                |(edge, reversed)| {
                    SegmentOfTraversal::new(
                        edges_map[&edge.id],
                        match reversed {
                            true => Direction::Decreasing,
                            false => Direction::Increasing,
                        },
                    )
                },
            ));
        let args = TraversalArgs {
            id: Some(fbb.create_string(&srv_ref)),
            segments: Some(segments),
            properties: None,
        };
        fb_traversals.push(Traversal::create(&mut fbb, &args));
    }

    println!("In the LRS, we have {} traversals.", fb_traversals.len());

    let nodes_index = HashMap::<osm4routing::NodeId, usize>::from_iter(
        nodes
            .iter()
            .enumerate()
            .map(|(index, node)| (node.id, index)),
    );
    let nodes: Vec<_> = nodes
        .iter()
        .map(|n| {
            let args = NodeArgs {
                id: Some(fbb.create_string(&n.id.0.to_string())),
                properties: None,
            };
            Node::create(&mut fbb, &args)
        })
        .collect();

    let segments: Vec<_> = edges
        .iter()
        .map(|e| {
            let points_iter = e.geometry.iter().map(|c| Point::new(c.lon, c.lat));
            let points = Some(fbb.create_vector_from_iter(points_iter));
            let args = SegmentArgs {
                id: Some(fbb.create_string(&e.id)),
                properties: None,
                geometry: points,
                start_node_index: nodes_index[&e.source] as u64,
                end_node_index: nodes_index[&e.target] as u64,
            };
            Segment::create(&mut fbb, &args)
        })
        .collect();

    let key = Some(fbb.create_string("source"));
    let value = Some(fbb.create_string("OpenStreetMap"));
    let source = Property::create(&mut fbb, &PropertyArgs { key, value });

    let lrs_args = LrsArgs {
        properties: Some(fbb.create_vector(&[source])),
        nodes: Some(fbb.create_vector(&nodes)),
        segments: Some(fbb.create_vector(&segments)),
        traversals: Some(fbb.create_vector_from_iter(fb_traversals.into_iter())),
        ..Default::default()
    };

    let lrs = Lrs::create(&mut fbb, &lrs_args);

    fbb.finish(lrs, None);
    let buffer = fbb.finished_data();

    std::fs::write(&cli_args.output_lrs, buffer).unwrap();
}

#[cfg(test)]
pub mod tests {
    use osm4routing::{Edge, NodeId};

    use crate::{can_be_appended, sort_edges, Candidate};

    fn edge(source: i64, target: i64) -> Edge {
        Edge {
            source: NodeId(source),
            target: NodeId(target),
            ..Default::default()
        }
    }

    #[test]
    fn test_is_candidate() {
        // last_edge is constant
        assert_eq!(
            can_be_appended(&edge(0, 1), &edge(1, 2), true),
            Candidate::Target
        );

        assert_eq!(
            can_be_appended(&edge(0, 1), &edge(1, 2), false),
            Candidate::Impossible
        );

        assert_eq!(
            can_be_appended(&edge(1, 0), &edge(1, 2), true),
            Candidate::Source
        );

        assert_eq!(
            can_be_appended(&edge(1, 0), &edge(1, 2), false),
            Candidate::Impossible
        );

        // last_edge in opposite directions
        assert_eq!(
            can_be_appended(&edge(0, 1), &edge(2, 1), true),
            Candidate::Impossible
        );

        assert_eq!(
            can_be_appended(&edge(0, 1), &edge(2, 1), false),
            Candidate::Target
        );

        assert_eq!(
            can_be_appended(&edge(1, 0), &edge(2, 1), true),
            Candidate::Impossible
        );

        assert_eq!(
            can_be_appended(&edge(1, 0), &edge(2, 1), false),
            Candidate::Source
        );
    }

    #[test]
    fn sort_edges_simple() {
        let e = edge(0, 1);

        let sorted = sort_edges(vec![e.clone()], "");
        assert_eq!(sorted[0].0, e);
        assert!(!sorted[0].1);
    }

    #[test]
    fn sort_edges_two_edges_in_order() {
        let e1 = edge(0, 1);
        let e2 = edge(1, 2);

        let sorted = sort_edges(vec![e1.clone(), e2.clone()], "");
        assert_eq!(sorted[0].0, e1);
        assert_eq!(sorted[1].0, e2);
        assert!(!sorted[0].1);
        assert!(!sorted[1].1);
    }

    #[test]
    fn sort_edges_two_edges_first_reversed() {
        let e1 = edge(1, 0);
        let e2 = edge(1, 2);

        let sorted = sort_edges(vec![e1.clone(), e2.clone()], "");
        assert_eq!(sorted[0].0, e1);
        assert_eq!(sorted[1].0, e2);
        assert!(sorted[0].1);
        assert!(!sorted[1].1);
    }
}
