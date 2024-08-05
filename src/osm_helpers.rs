//! Helper functions to manipulate OpenStreetMap data
//! and extract an ordered topology

use osm4routing::Edge;

/// When sorting the edges, each candidate is tested to see if they match and if they need to be reversed.
#[derive(PartialEq, Eq, Debug)]
enum Candidate {
    Source,
    Target,
    Impossible,
}

/// If we have a string of edges that ends with `end_edge`,
/// the candidate `Edge(s, t)` can be joined at the `end_edge(source, target)`.
/// Returns Forward if the `Edge` can be append in the same direction `(target == s)`,
/// Backward if it must be reversed `(target == t)` or
/// And NotCandidate if it can’t be appended.
/// `consider_source` means that the source (or target if false) of the `last_edge` is considered.
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

/// Sort edges from OpenStreetMap to build continous traversals.
///
/// The traversals are identified by a tag that is used on many ways.
/// We try to build the longest continous chain of ways, but the is no guarantee to succeed.
/// The ways might not share nodes or they might represent a tree.
pub fn sort_edges(edges: Vec<Edge>, traversal_ref: &str) -> Vec<(Edge, bool)> {
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
        println!("[WARN] on traversal {traversal_ref}, ignoring {ignored} edges out of {total}. Sorted from {} to {}", first.0, last.0);
    }

    sorted
}

#[cfg(test)]
pub mod tests {
    use osm4routing::{Edge, NodeId};

    use super::*;

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
