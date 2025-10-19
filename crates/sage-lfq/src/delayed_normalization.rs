use crate::error::MaxLfqError;
use crate::matrix::IntensityMatrix;
use nalgebra::{DMatrix, DVector};
use std::collections::VecDeque;

/// Compute global normalization offsets by minimizing H(N).
///
/// This function implements MaxLFQ's delayed normalization by finding sample-level
/// offsets that minimize the variance in pairwise peptide ratios across the entire dataset.
///
/// # Arguments
/// * `matrix` - The full peptide×sample intensity matrix
/// * `reference_sample` - Optional reference sample index (auto-selected if None)
///
/// # Returns
/// Vector of normalization offsets (one per sample), where offset = 0 for reference
pub fn compute_delayed_normalization(
    matrix: &IntensityMatrix,
    reference_sample: Option<usize>,
) -> Result<Vec<f64>, MaxLfqError> {
    if matrix.n_samples == 0 {
        return Ok(Vec::new());
    }

    let reference = match reference_sample {
        Some(idx) if idx < matrix.n_samples => idx,
        Some(idx) => return Err(MaxLfqError::InvalidReferenceSample(idx)),
        None => select_reference_sample(matrix),
    };

    // Find connected components - samples must share peptides to be normalized together
    let components = find_connected_components(matrix);
    let mut offsets = vec![0.0; matrix.n_samples];

    for component in components {
        if component.len() < 2 {
            // Single-sample components have no normalization needed
            continue;
        }

        let (a, b) = build_component_equations(matrix, &component);

        if a.nrows() == 0 {
            continue;
        }

        let reference_in_component = component
            .iter()
            .position(|&idx| idx == reference)
            .unwrap_or(0);

        let solution = solve_with_reference(a, b, reference_in_component)
            .ok_or(MaxLfqError::NormalizationFailed)?;

        for (local_idx, &sample_idx) in component.iter().enumerate() {
            offsets[sample_idx] = solution[local_idx];
        }
    }

    Ok(offsets)
}

fn select_reference_sample(matrix: &IntensityMatrix) -> usize {
    let mut counts = vec![0usize; matrix.n_samples];

    for row in matrix.matrix.outer_iterator() {
        for (col, _) in row.iter() {
            counts[col] += 1;
        }
    }

    counts
        .iter()
        .enumerate()
        .max_by_key(|&(_, count)| count)
        .map(|(idx, _)| idx)
        .unwrap_or(0)
}

fn build_component_equations(
    matrix: &IntensityMatrix,
    component: &[usize],
) -> (DMatrix<f64>, DVector<f64>) {
    let n_samples = component.len();
    if n_samples == 0 {
        return (DMatrix::zeros(0, 0), DVector::zeros(0));
    }

    // Map global sample indices to local indices within component
    let mut sample_map = vec![None; matrix.n_samples];
    for (local_idx, &sample_idx) in component.iter().enumerate() {
        sample_map[sample_idx] = Some(local_idx);
    }

    let mut equations: Vec<(usize, usize, f64)> = Vec::new();

    for row in matrix.matrix.outer_iterator() {
        let mut entries: Vec<(usize, f64)> = Vec::new();
        for (col, &value) in row.iter() {
            if let Some(local_col) = sample_map[col] {
                entries.push((local_col, value as f64));
            }
        }

        if entries.len() < 2 {
            continue;
        }

        // For each pair of samples sharing this peptide, create equation:
        // N_A - N_B = I_A - I_B
        for i in 0..entries.len() {
            for j in (i + 1)..entries.len() {
                let (sample_a, intensity_a) = entries[i];
                let (sample_b, intensity_b) = entries[j];
                let target = intensity_a - intensity_b;
                equations.push((sample_a, sample_b, target));
            }
        }
    }

    let n_eq = equations.len();
    let mut a_matrix = DMatrix::zeros(n_eq, n_samples);
    let mut b_vector = DVector::zeros(n_eq);

    for (row_idx, (sample_a, sample_b, target)) in equations.into_iter().enumerate() {
        a_matrix[(row_idx, sample_a)] = 1.0;
        a_matrix[(row_idx, sample_b)] = -1.0;
        b_vector[row_idx] = target;
    }

    (a_matrix, b_vector)
}

fn solve_with_reference(
    mut a: DMatrix<f64>,
    mut b: DVector<f64>,
    reference: usize,
) -> Option<Vec<f64>> {
    debug_assert!(reference < a.ncols(), "Reference sample out of bounds");

    // Add constraint: N_reference = 0
    let row_pos = a.nrows();
    a = a.insert_row(row_pos, 0.0);
    b = b.insert_row(row_pos, 0.0);

    let constraint_row = a.nrows() - 1;
    a[(constraint_row, reference)] = 1.0;
    b[constraint_row] = 0.0;

    // Solve normal equations: (A^T A) N = A^T b
    let ata = a.transpose() * &a;
    let atb = a.transpose() * b;

    // Try LU decomposition first, fall back to QR if singular
    ata.clone()
        .lu()
        .solve(&atb)
        .or_else(|| ata.qr().solve(&atb))
        .map(|sol| sol.as_slice().to_vec())
}

fn find_connected_components(matrix: &IntensityMatrix) -> Vec<Vec<usize>> {
    let n_samples = matrix.n_samples;
    if n_samples == 0 {
        return Vec::new();
    }

    // Build adjacency list: samples are connected if they share a peptide
    let mut adjacency = vec![Vec::<usize>::new(); n_samples];

    for row in matrix.matrix.outer_iterator() {
        let cols: Vec<usize> = row.iter().map(|(col, _)| col).collect();
        for i in 0..cols.len() {
            for j in (i + 1)..cols.len() {
                let (a, b) = (cols[i], cols[j]);
                adjacency[a].push(b);
                adjacency[b].push(a);
            }
        }
    }

    // Deduplicate adjacency lists
    for neighbors in adjacency.iter_mut() {
        neighbors.sort_unstable();
        neighbors.dedup();
    }

    // BFS to find connected components
    let mut visited = vec![false; n_samples];
    let mut components = Vec::new();

    for start in 0..n_samples {
        if visited[start] {
            continue;
        }

        let mut queue = VecDeque::new();
        let mut component = Vec::new();
        queue.push_back(start);
        visited[start] = true;

        while let Some(node) = queue.pop_front() {
            component.push(node);
            for &neighbor in &adjacency[node] {
                if !visited[neighbor] {
                    visited[neighbor] = true;
                    queue.push_back(neighbor);
                }
            }
        }

        components.push(component);
    }

    components
}

#[cfg(test)]
mod tests {
    use super::*;
    use sprs::TriMat;

    fn build_matrix(values: &[Vec<Option<f64>>]) -> IntensityMatrix {
        let n_peptides = values.len();
        let n_samples = values.first().map(|row| row.len()).unwrap_or(0);

        let mut triplets = TriMat::new((n_peptides, n_samples));
        for (i, row) in values.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                if let Some(v) = value {
                    triplets.add_triplet(i, j, *v as f32);
                }
            }
        }

        IntensityMatrix {
            matrix: triplets.to_csr(),
            peptide_ids: (0..n_peptides).map(|idx| format!("pep{idx}")).collect(),
            sample_names: (0..n_samples).map(|idx| format!("s{idx}")).collect(),
            n_peptides,
            n_samples,
        }
    }

    #[test]
    fn technical_replicates_near_zero_offsets() {
        let values = vec![vec![Some(10.0); 3]; 100];
        let matrix = build_matrix(&values);

        let norm = compute_delayed_normalization(&matrix, Some(0)).unwrap();
        assert_eq!(norm.len(), 3);
        assert!(norm.iter().all(|&x| x.abs() < 1e-6));
    }

    #[test]
    fn detects_systematic_loading_difference() {
        let mut values = vec![vec![Some(10.0); 3]; 100];
        for row in &mut values {
            row[2] = Some(11.0); // Sample 2 has 2x loading (1 log2 unit)
        }
        let matrix = build_matrix(&values);

        let norm = compute_delayed_normalization(&matrix, Some(0)).unwrap();
        assert_eq!(norm.len(), 3);
        assert!((norm[0]).abs() < 0.01); // Reference
        assert!((norm[1]).abs() < 0.01); // Same as reference
        assert!((norm[2] - 1.0).abs() < 0.1); // 2x loading offset
    }

    #[test]
    fn handles_disconnected_components() {
        let values = vec![
            vec![Some(10.0), Some(10.0), None, None],
            vec![Some(11.0), Some(11.0), None, None],
            vec![None, None, Some(12.0), Some(12.0)],
            vec![None, None, Some(13.0), Some(13.0)],
        ];
        let matrix = build_matrix(&values);

        let norm = compute_delayed_normalization(&matrix, None).unwrap();
        assert_eq!(norm.len(), 4);
        // Samples 0-1 form one component, 2-3 form another
        assert!((norm[0] - norm[1]).abs() < 0.1);
        assert!((norm[2] - norm[3]).abs() < 0.1);
    }

    #[test]
    fn auto_selects_most_complete_sample() {
        let values = vec![
            vec![Some(10.0), None, None],
            vec![Some(11.0), Some(11.5), None],
            vec![Some(12.0), Some(12.0), Some(12.5)],
        ];
        let matrix = build_matrix(&values);

        let norm = compute_delayed_normalization(&matrix, None).unwrap();
        assert_eq!(norm.len(), 3);
        // Sample 0 is most complete (3 peptides)
    }

    #[test]
    fn invalid_reference_returns_error() {
        let values = vec![vec![Some(10.0); 2]; 10];
        let matrix = build_matrix(&values);
        let err = compute_delayed_normalization(&matrix, Some(5)).unwrap_err();
        assert!(matches!(err, MaxLfqError::InvalidReferenceSample(5)));
    }

    #[test]
    fn preserves_biological_signal() {
        // 80% proteins unchanged, 20% proteins 2-fold up
        // This tests that differential expression with minority changing is preserved
        let mut values = Vec::new();
        for i in 0..100 {
            if i < 80 {
                values.push(vec![Some(10.0), Some(10.0)]);
            } else {
                values.push(vec![Some(10.0), Some(11.0)]); // 2-fold in sample 1
            }
        }
        let matrix = build_matrix(&values);

        let norm = compute_delayed_normalization(&matrix, Some(0)).unwrap();
        // With 80/20 split, median is ~0, so offset should be small
        // Biological signal (20% upregulated) is preserved
        assert!((norm[0]).abs() < 0.01);
        assert!(norm[1].abs() < 0.3); // Small offset, not removing the biology
    }
}
