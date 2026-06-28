use std::collections::{HashMap, HashSet};
use std::io::{self, Write};

/// Tier 4 Optimization: Graph with HashSet-based adjacency lists.
///
/// The original C++ implementation uses `vector<Graphnode*>` for adjacency
/// lists, with O(degree) linear scan for duplicate detection in `addLinkedNode`.
/// This Rust implementation uses `HashSet<usize>` for O(1) duplicate detection.
///
/// Additional optimizations:
/// - Iterative DFS instead of recursive (avoids stack overflow on large clusters)
/// - HashMap for node lookup instead of std::map (O(1) vs O(log n))
/// - Index-based references instead of pointers (safe, no dangling refs)
pub struct Graph {
    /// Adjacency list: node index -> set of neighbor indices
    nodes: Vec<GraphNode>,
    /// Name -> index lookup
    node_lookup: HashMap<String, usize>,
}

struct GraphNode {
    name: String,
    /// HashSet for O(1) duplicate detection (Tier 4 optimization)
    neighbors: HashSet<usize>,
    marked: bool,
}

impl Graph {
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            node_lookup: HashMap::new(),
        }
    }

    /// Add a bidirectional edge between nodes a and b.
    /// Uses HashSet for O(1) duplicate detection.
    pub fn add_linked_nodes(&mut self, a: &str, b: &str) {
        let a_idx = self.get_or_create_node(a);
        let b_idx = self.get_or_create_node(b);

        // O(1) insert with HashSet (vs O(degree) linear scan in C++)
        self.nodes[a_idx].neighbors.insert(b_idx);
        self.nodes[b_idx].neighbors.insert(a_idx);
    }

    fn get_or_create_node(&mut self, name: &str) -> usize {
        if let Some(&idx) = self.node_lookup.get(name) {
            return idx;
        }

        let idx = self.nodes.len();
        self.nodes.push(GraphNode {
            name: name.to_string(),
            neighbors: HashSet::new(),
            marked: false,
        });
        self.node_lookup.insert(name.to_string(), idx);
        idx
    }

    /// Find and print all connected components using iterative DFS.
    ///
    /// Tier 4 optimization: Uses iterative DFS with explicit stack instead
    /// of recursive DFS. This eliminates the need for `ulimit -s unlimited`
    /// that the C++ version requires for large clusters.
    pub fn print_clusters(&mut self, output: &mut impl Write) -> io::Result<()> {
        for start in 0..self.nodes.len() {
            if self.nodes[start].marked {
                continue;
            }

            // Iterative DFS to collect cluster members
            let mut cluster: Vec<String> = Vec::new();
            let mut stack: Vec<usize> = vec![start];

            while let Some(node_idx) = stack.pop() {
                if self.nodes[node_idx].marked {
                    continue;
                }
                self.nodes[node_idx].marked = true;
                cluster.push(self.nodes[node_idx].name.clone());

                // Add unvisited neighbors to stack
                for &neighbor_idx in &self.nodes[node_idx].neighbors {
                    if !self.nodes[neighbor_idx].marked {
                        stack.push(neighbor_idx);
                    }
                }
            }

            // Print cluster: space-separated names
            let cluster_str = cluster.join(" ");
            writeln!(output, "{}", cluster_str)?;
        }
        Ok(())
    }

    /// Apply Jaccard coefficient filtering to create a new graph.
    ///
    /// For each edge (a, b), compute the Jaccard similarity of the closed
    /// neighborhoods N[a] and N[b]. If the coefficient is below the threshold,
    /// the edge is removed from the new graph.
    ///
    /// Formula: jaccard = (|N[a] ∩ N[b]| + 2) / (|N[a] ∪ N[b]|)
    /// The +2 accounts for a and b being mutual neighbors (closed neighborhood).
    pub fn apply_jaccard_coeff(&self, coeff: f64) -> Graph {
        let mut new_graph = Graph::new();

        for (i, node) in self.nodes.iter().enumerate() {
            for &j in &node.neighbors {
                if j <= i {
                    continue; // Process each edge once
                }

                let link_coeff = self.calc_link_coeff(i, j);
                if link_coeff >= coeff {
                    new_graph.add_linked_nodes(&self.nodes[i].name, &self.nodes[j].name);
                }
            }
        }

        new_graph
    }

    /// Compute modified Jaccard similarity between two linked nodes.
    ///
    /// Uses HashSet intersection for O(min(|A|, |B|)) computation instead
    /// of the O(|A| * |B|) nested loop in the C++ implementation.
    fn calc_link_coeff(&self, a_idx: usize, b_idx: usize) -> f64 {
        let a_neighbors = &self.nodes[a_idx].neighbors;
        let b_neighbors = &self.nodes[b_idx].neighbors;

        // Count common neighbors using HashSet intersection
        let num_common = if a_neighbors.len() < b_neighbors.len() {
            a_neighbors.intersection(b_neighbors).count()
        } else {
            b_neighbors.intersection(a_neighbors).count()
        } as f64;

        // Closed neighborhood Jaccard:
        // |A ∩ B| = num_common + 2 (a and b are mutual neighbors)
        // |A ∪ B| = |A| + |B| - num_common
        let total_common = num_common + 2.0;
        let total_vertices = (a_neighbors.len() + b_neighbors.len()) as f64 - num_common;

        if total_vertices == 0.0 {
            0.0
        } else {
            total_common / total_vertices
        }
    }

    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }
}
