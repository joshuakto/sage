use crate::error::MaxLfqError;
use crate::matrix::IntensityMatrix;
use petgraph::algo::tarjan_scc;
use petgraph::graph::{NodeIndex, UnGraph};

pub struct NormalizationFactors {
    pub factors: Vec<f64>, // Log-scale normalization per sample
}

impl NormalizationFactors {
    pub fn compute(matrix: &IntensityMatrix, min_peptides: usize) -> Result<Self, MaxLfqError> {
        if matrix.n_samples == 0 {
            return Ok(Self {
                factors: Vec::new(),
            });
        }

        let graph = build_sample_graph(matrix, min_peptides);
        let components = tarjan_scc(&graph);
        let mut factors = vec![0.0; matrix.n_samples];

        for component in components {
            if component.len() == 1 {
                continue;
            }

            let component_factors = optimize_component(matrix, &component, &graph)?;
            for (node, factor) in component.iter().zip(component_factors) {
                factors[node.index()] = factor;
            }
        }

        Ok(Self { factors })
    }
}

fn build_sample_graph(matrix: &IntensityMatrix, min_peptides: usize) -> UnGraph<usize, usize> {
    let mut graph = UnGraph::new_undirected();

    let nodes: Vec<_> = (0..matrix.n_samples).map(|i| graph.add_node(i)).collect();

    for i in 0..matrix.n_samples {
        for j in (i + 1)..matrix.n_samples {
            let shared = matrix.count_shared_peptides(i, j);
            if shared >= min_peptides {
                graph.add_edge(nodes[i], nodes[j], shared);
            }
        }
    }

    graph
}

fn optimize_component(
    matrix: &IntensityMatrix,
    component: &[NodeIndex],
    _graph: &UnGraph<usize, usize>,
) -> Result<Vec<f64>, MaxLfqError> {
    if component.is_empty() {
        return Ok(Vec::new());
    }

    let mut averages = vec![0.0; component.len()];
    let mut counts = vec![0usize; component.len()];

    for row in matrix.matrix.outer_iterator() {
        let mut entries: Vec<(usize, f64)> = Vec::new();
        for (col, &value) in row.iter() {
            if let Some(pos) = component.iter().position(|node| node.index() == col) {
                entries.push((pos, value as f64));
            }
        }

        for (pos, value) in entries {
            averages[pos] += value;
            counts[pos] += 1;
        }
    }

    for (avg, count) in averages.iter_mut().zip(counts.iter()) {
        if *count > 0 {
            *avg /= *count as f64;
        }
    }

    Ok(averages)
}
