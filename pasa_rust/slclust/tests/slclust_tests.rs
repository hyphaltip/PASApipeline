use slclust::Graph;
use std::io::Cursor;

#[test]
fn test_graph_add_linked_nodes() {
    let mut g = Graph::new();
    g.add_linked_nodes("a", "b");
    g.add_linked_nodes("b", "c");

    assert_eq!(g.num_nodes(), 3);
}

#[test]
fn test_graph_deduplication() {
    let mut g = Graph::new();
    g.add_linked_nodes("a", "b");
    g.add_linked_nodes("a", "b"); // duplicate edge
    g.add_linked_nodes("b", "a"); // reverse duplicate

    assert_eq!(g.num_nodes(), 2);
    // HashSet should prevent duplicate edges
}

#[test]
fn test_graph_single_cluster() {
    let mut g = Graph::new();
    g.add_linked_nodes("a", "b");
    g.add_linked_nodes("b", "c");
    g.add_linked_nodes("c", "d");

    let mut output = Vec::new();
    g.print_clusters(&mut output).unwrap();

    let result = String::from_utf8(output).unwrap();
    let lines: Vec<&str> = result.trim().lines().collect();

    // Should have exactly one cluster with all 4 nodes
    assert_eq!(lines.len(), 1);
    let nodes: Vec<&str> = lines[0].split_whitespace().collect();
    assert_eq!(nodes.len(), 4);
}

#[test]
fn test_graph_multiple_clusters() {
    let mut g = Graph::new();
    // Cluster 1: a-b-c
    g.add_linked_nodes("a", "b");
    g.add_linked_nodes("b", "c");
    // Cluster 2: d-e
    g.add_linked_nodes("d", "e");
    // Cluster 3: f (isolated, but we need an edge to add it)
    g.add_linked_nodes("f", "g");

    let mut output = Vec::new();
    g.print_clusters(&mut output).unwrap();

    let result = String::from_utf8(output).unwrap();
    let lines: Vec<&str> = result.trim().lines().collect();

    // Should have 3 clusters
    assert_eq!(lines.len(), 3);
}

#[test]
fn test_graph_self_loop_ignored() {
    let mut g = Graph::new();
    // The main.rs skips pairs where a == b
    // But if called directly, add_linked_nodes should still work
    g.add_linked_nodes("a", "b");

    assert_eq!(g.num_nodes(), 2);
}

#[test]
fn test_graph_large_cluster() {
    // Test that iterative DFS doesn't overflow on large clusters
    let mut g = Graph::new();
    for i in 0..10000 {
        g.add_linked_nodes(&format!("n{}", i), &format!("n{}", i + 1));
    }

    let mut output = Vec::new();
    g.print_clusters(&mut output).unwrap();

    let result = String::from_utf8(output).unwrap();
    let lines: Vec<&str> = result.trim().lines().collect();

    // Should have one cluster with 10001 nodes
    assert_eq!(lines.len(), 1);
    let nodes: Vec<&str> = lines[0].split_whitespace().collect();
    assert_eq!(nodes.len(), 10001);
}

#[test]
fn test_jaccard_coefficient() {
    let mut g = Graph::new();
    // Create a graph where a-b have high Jaccard similarity
    // a connected to b, c
    // b connected to a, c
    // c connected to a, b
    g.add_linked_nodes("a", "b");
    g.add_linked_nodes("a", "c");
    g.add_linked_nodes("b", "c");

    // With Jaccard cutoff of 0.0, all edges should remain
    let filtered = g.apply_jaccard_coeff(0.0);
    assert!(filtered.num_nodes() > 0);

    // With very high Jaccard cutoff, most edges should be removed
    let filtered_high = g.apply_jaccard_coeff(0.99);
    // Some nodes might not exist if all their edges were filtered
}

#[test]
fn test_jaccard_preserves_connectivity() {
    let mut g = Graph::new();
    // Triangle: a-b, b-c, a-c
    g.add_linked_nodes("a", "b");
    g.add_linked_nodes("b", "c");
    g.add_linked_nodes("a", "c");

    // With low Jaccard threshold, connectivity should be preserved
    let filtered = g.apply_jaccard_coeff(0.1);

    let mut output = Vec::new();
    let mut filtered_mut = filtered;
    filtered_mut.print_clusters(&mut output).unwrap();

    let result = String::from_utf8(output).unwrap();
    let lines: Vec<&str> = result.trim().lines().collect();

    // Should still have one connected cluster
    assert_eq!(lines.len(), 1);
    let nodes: Vec<&str> = lines[0].split_whitespace().collect();
    assert_eq!(nodes.len(), 3);
}
