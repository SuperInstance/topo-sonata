//! Persistent homology of harmonic spaces.
//!
//! Computes H₀, H₁, H₂ via boundary matrix reduction (Gaussian elimination).
//! Betti barcodes serve as tension profiles: long bars = persistent dissonance.

use crate::{BettiSequence, Filtration, PersistenceBarcode, SimplicialComplex};

/// Compute persistent homology barcodes for a filtered complex.
///
/// Uses the standard matrix reduction algorithm on the boundary matrix,
/// processing simplices in filtration order.
pub fn compute_persistence(filtration: &Filtration) -> Vec<PersistenceBarcode> {
    let simplices = &filtration.complex.simplices;
    let values = &filtration.values;
    let n = simplices.len();
    if n == 0 {
        return vec![
            PersistenceBarcode { dimension: 0, bars: vec![] },
            PersistenceBarcode { dimension: 1, bars: vec![] },
            PersistenceBarcode { dimension: 2, bars: vec![] },
        ];
    }

    // Create index array sorted by filtration value (stable for ties)
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        values[a]
            .partial_cmp(&values[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Build boundary matrix: columns = simplices in filtration order,
    // rows = simplices in filtration order.
    // boundary[col_index][row_index] = true if simplex `order[row_index]` is
    // a face of simplex `order[col_index]` with the appropriate sign parity.
    //
    // Actually we'll use the standard approach: for each simplex σ,
    // its boundary is the sum of its (d-1)-dimensional faces.
    // We represent the boundary matrix as a list of low entries.

    // Build a map from simplex (as sorted vec) to its original index.
    let mut simplex_to_idx: std::collections::BTreeMap<Vec<usize>, usize> =
        std::collections::BTreeMap::new();
    for (idx, s) in simplices.iter().enumerate() {
        simplex_to_idx.insert(s.clone(), idx);
    }

    // Map from original simplex index to position in filtration order
    let mut pos_in_order = vec![0usize; n];
    for (pos, &idx) in order.iter().enumerate() {
        pos_in_order[idx] = pos;
    }

    // Rebuild boundary
    let mut boundary: Vec<Vec<usize>> = vec![vec![]; n];

    for (col_pos, &simp_idx) in order.iter().enumerate() {
        let simplex = &simplices[simp_idx];
        if simplex.len() <= 1 {
            continue;
        }
        for i in 0..simplex.len() {
            let mut face: Vec<usize> = simplex.clone();
            face.remove(i);
            if let Some(&face_orig_idx) = simplex_to_idx.get(&face) {
                let row_pos = pos_in_order[face_orig_idx];
                boundary[col_pos].push(row_pos);
            }
        }
        boundary[col_pos].sort_unstable();
    }

    // Reduce the boundary matrix using Gaussian elimination (mod 2)
    // We maintain a "low" array: low[col] = row index of lowest non-zero entry
    // or None if the column is zero.
    let mut low: Vec<Option<usize>> = vec![None; n];
    // Also maintain the reduced matrix columns
    let mut reduced: Vec<Vec<usize>> = boundary.clone();

    for col in 0..n {
        // Find the lowest non-zero row
        let mut lowest = get_lowest(&reduced[col]);
        while let Some(low_row) = lowest {
            // Check if another column already has this as its low
            let mut found_collision = false;
            for prev_col in 0..col {
                if low[prev_col] == Some(low_row) {
                    // Add (XOR) prev_col into col
                    reduced[col] = symmetric_difference(&reduced[col], &reduced[prev_col]);
                    lowest = get_lowest(&reduced[col]);
                    found_collision = true;
                    break;
                }
            }
            if !found_collision {
                low[col] = Some(low_row);
                break;
            }
        }
    }

    // Extract barcodes from the reduced matrix
    // A pair (low[col]=row, col) means a (dim(row))-dimensional class born at
    // filtration value of order[row] and dies at filtration value of order[col].
    let mut births: Vec<Option<usize>> = vec![None; n]; // index in order

    for col in 0..n {
        if let Some(row) = low[col] {
            births[row] = Some(col);
        }
    }

    let max_dim = filtration.complex.dimension;
    let mut barcodes: Vec<PersistenceBarcode> = (0..=max_dim)
        .map(|d| PersistenceBarcode {
            dimension: d,
            bars: vec![],
        })
        .collect();

    for row in 0..n {
        let simp_idx = order[row];
        let dim = simplices[simp_idx].len().saturating_sub(1);
        if dim > max_dim {
            continue;
        }
        let birth_val = values[simp_idx];
        if let Some(col) = births[row] {
            let death_idx = order[col];
            let death_val = values[death_idx];
            if death_val > birth_val {
                barcodes[dim].bars.push((birth_val, death_val));
            }
            // if death_val == birth_val, the bar is too short to be meaningful
        } else {
            // Class persists to infinity — use f64::INFINITY
            barcodes[dim].bars.push((birth_val, f64::INFINITY));
        }
    }

    // Sort bars within each barcode
    for bc in &mut barcodes {
        bc.bars.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    }

    barcodes
}

/// Get the lowest (max) row index in a sparse column.
fn get_lowest(col: &[usize]) -> Option<usize> {
    if col.is_empty() {
        None
    } else {
        // The column is sorted, so the last entry is the lowest (max row index)
        Some(col[col.len() - 1])
    }
}

/// Symmetric difference of two sorted vectors (XOR, mod 2 addition).
fn symmetric_difference(a: &[usize], b: &[usize]) -> Vec<usize> {
    let mut result = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        if a[i] < b[j] {
            result.push(a[i]);
            i += 1;
        } else if a[i] > b[j] {
            result.push(b[j]);
            j += 1;
        } else {
            // Both present → cancel (mod 2)
            i += 1;
            j += 1;
        }
    }
    while i < a.len() {
        result.push(a[i]);
        i += 1;
    }
    while j < b.len() {
        result.push(b[j]);
        j += 1;
    }
    result
}

/// Compute Betti numbers from a persistence barcode.
///
/// βₖ = number of bars in dimension k that span to infinity.
pub fn betti_from_barcodes(barcodes: &[PersistenceBarcode]) -> BettiSequence {
    if barcodes.is_empty() {
        return BettiSequence { betti: vec![0] };
    }
    let max_dim = barcodes.iter().map(|bc| bc.dimension).max().unwrap_or(0);
    let mut betti = vec![0usize; max_dim + 1];
    for bc in barcodes {
        for &(birth, death) in &bc.bars {
            if death == f64::INFINITY {
                betti[bc.dimension] += 1;
            }
        }
    }
    BettiSequence { betti }
}

/// Convenience: compute Betti numbers for an unfiltered simplicial complex.
///
/// Assigns filtration value 0.0 to every simplex and computes persistent homology.
pub fn betti_numbers(complex: &SimplicialComplex) -> BettiSequence {
    let filtration = Filtration {
        complex: complex.clone(),
        values: vec![0.0; complex.simplices.len()],
    };
    let barcodes = compute_persistence(&filtration);
    betti_from_barcodes(&barcodes)
}

/// Compute Betti numbers at a specific filtration threshold.
///
/// Counts bars in the barcode whose interval contains `threshold`.
pub fn betti_at_threshold(barcodes: &[PersistenceBarcode], threshold: f64) -> BettiSequence {
    if barcodes.is_empty() {
        return BettiSequence { betti: vec![0] };
    }
    let max_dim = barcodes.iter().map(|bc| bc.dimension).max().unwrap_or(0);
    let mut betti = vec![0usize; max_dim + 1];
    for bc in barcodes {
        for &(birth, death) in &bc.bars {
            let alive = birth <= threshold && (death > threshold || death == f64::INFINITY);
            if alive {
                betti[bc.dimension] += 1;
            }
        }
    }
    BettiSequence { betti }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Chord, SimplicialComplex};

    /// Build a complex with given simplices directly.
    fn make_complex(simplices: Vec<Vec<usize>>) -> SimplicialComplex {
        let mut s = simplices;
        for face in &mut s {
            face.sort_unstable();
        }
        s.sort_unstable();
        let dimension = s.iter().map(|f| f.len().saturating_sub(1)).max().unwrap_or(0);
        SimplicialComplex { simplices: s, dimension }
    }

    #[test]
    fn test_betti_point() {
        // A single point: β₀=1, β₁=0
        let complex = make_complex(vec![vec![0]]);
        let betti = betti_numbers(&complex);
        assert_eq!(betti.betti, vec![1]);
    }

    #[test]
    fn test_betti_two_points() {
        // Two disconnected points: β₀=2
        let complex = make_complex(vec![vec![0], vec![1]]);
        let betti = betti_numbers(&complex);
        assert_eq!(betti.betti, vec![2]);
    }

    #[test]
    fn test_betti_edge() {
        // Two points connected by an edge: β₀=1, β₁=0
        let complex = make_complex(vec![vec![0], vec![1], vec![0, 1]]);
        let betti = betti_numbers(&complex);
        assert_eq!(betti.betti, vec![1, 0]);
    }

    #[test]
    fn test_betti_triangle_filled() {
        // A filled triangle: β₀=1, β₁=0
        let complex = make_complex(vec![
            vec![0], vec![1], vec![2],
            vec![0, 1], vec![0, 2], vec![1, 2],
            vec![0, 1, 2],
        ]);
        let betti = betti_numbers(&complex);
        assert_eq!(betti.betti, vec![1, 0, 0]);
    }

    #[test]
    fn test_betti_triangle_hollow() {
        // A hollow triangle (3 edges, no face): β₀=1, β₁=1
        let complex = make_complex(vec![
            vec![0], vec![1], vec![2],
            vec![0, 1], vec![0, 2], vec![1, 2],
        ]);
        let betti = betti_numbers(&complex);
        assert_eq!(betti.betti, vec![1, 1]);
    }

    #[test]
    fn test_betti_square_hollow() {
        // A hollow square (4 edges, no faces): β₀=1, β₁=1
        let complex = make_complex(vec![
            vec![0], vec![1], vec![2], vec![3],
            vec![0, 1], vec![1, 2], vec![2, 3], vec![0, 3],
        ]);
        let betti = betti_numbers(&complex);
        assert_eq!(betti.betti, vec![1, 1]);
    }

    #[test]
    fn test_betti_tetrahedron_filled() {
        // A solid tetrahedron: β₀=1, β₁=0, β₂=0
        let complex = make_complex(vec![
            vec![0], vec![1], vec![2], vec![3],
            vec![0, 1], vec![0, 2], vec![0, 3], vec![1, 2], vec![1, 3], vec![2, 3],
            vec![0, 1, 2], vec![0, 1, 3], vec![0, 2, 3], vec![1, 2, 3],
            vec![0, 1, 2, 3],
        ]);
        let betti = betti_numbers(&complex);
        assert_eq!(betti.betti, vec![1, 0, 0, 0]);
    }

    #[test]
    fn test_betti_sphere() {
        // Hollow tetrahedron (boundary of a 3-simplex): β₀=1, β₁=0, β₂=1
        let complex = make_complex(vec![
            vec![0], vec![1], vec![2], vec![3],
            vec![0, 1], vec![0, 2], vec![0, 3], vec![1, 2], vec![1, 3], vec![2, 3],
            vec![0, 1, 2], vec![0, 1, 3], vec![0, 2, 3], vec![1, 2, 3],
        ]);
        let betti = betti_numbers(&complex);
        assert_eq!(betti.betti, vec![1, 0, 1]);
    }

    #[test]
    fn test_symmetric_difference() {
        assert_eq!(symmetric_difference(&[1, 3, 5], &[2, 3, 4]), vec![1, 2, 4, 5]);
        assert_eq!(symmetric_difference(&[1, 2], &[1, 2]), vec![]);
        assert_eq!(symmetric_difference(&[], &[1, 2]), vec![1, 2]);
    }

    #[test]
    fn test_persistence_barcode_extraction() {
        // Edge filtered: vertex 0 at t=0, vertex 1 at t=0, edge at t=1
        let complex = make_complex(vec![vec![0], vec![1], vec![0, 1]]);
        let filtration = Filtration {
            complex,
            values: vec![0.0, 0.0, 1.0],
        };
        let barcodes = compute_persistence(&filtration);
        // H₀: one bar [0, ∞), one bar [0, 1)
        assert_eq!(barcodes[0].dimension, 0);
        assert!(barcodes[0].bars.len() >= 1);
    }

    #[test]
    fn test_betti_from_barcodes() {
        let barcodes = vec![
            PersistenceBarcode {
                dimension: 0,
                bars: vec![(0.0, f64::INFINITY), (0.0, 1.0)],
            },
            PersistenceBarcode {
                dimension: 1,
                bars: vec![(1.0, 2.0)],
            },
        ];
        let betti = betti_from_barcodes(&barcodes);
        assert_eq!(betti.betti, vec![1, 0]); // only infinite bars count
    }

    #[test]
    fn test_betti_at_threshold() {
        let barcodes = vec![
            PersistenceBarcode {
                dimension: 0,
                bars: vec![(0.0, f64::INFINITY), (0.0, 1.0)],
            },
            PersistenceBarcode {
                dimension: 1,
                bars: vec![(0.5, 2.0)],
            },
        ];
        let betti = betti_at_threshold(&barcodes, 0.5);
        // At t=0.5: H₀ has 2 alive bars, H₁ has 1 alive bar
        assert_eq!(betti.betti, vec![2, 1]);
    }

    #[test]
    fn test_full_pipeline_chord_complex() {
        use crate::simplicial_chord::progression_to_complex;
        let c = Chord { notes: vec![0, 4, 7] };
        let g = Chord { notes: vec![2, 7, 11] };
        let complex = progression_to_complex(&[c, g]);
        let betti = betti_numbers(&complex);
        // Two filled triangles sharing an edge should have β₀=1, β₁=0
        assert_eq!(betti.betti[0], 1);
    }
}
