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
            PersistenceBarcode {
                dimension: 0,
                bars: vec![],
            },
            PersistenceBarcode {
                dimension: 1,
                bars: vec![],
            },
            PersistenceBarcode {
                dimension: 2,
                bars: vec![],
            },
        ];
    }

    // Create index array sorted by: filtration value, then dimension (ascending),
    // then lexicographic. This ensures lower-dimensional simplices are processed
    // before higher-dimensional ones at the same filtration value.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        values[*a]
            .partial_cmp(&values[*b])
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                simplices[a]
                    .len()
                    .cmp(&simplices[b].len())
            })
            .then_with(|| simplices[a].cmp(simplices[b]))
    });

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

    // Boundary matrix (sparse, mod 2): boundary[col] = sorted list of row positions
    let mut reduced: Vec<Vec<usize>> = vec![vec![]; n];

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
                reduced[col_pos].push(row_pos);
            }
        }
        reduced[col_pos].sort_unstable();
    }

    // Reduce the boundary matrix using Gaussian elimination (mod 2).
    // low[col] = row index of lowest non-zero entry, or None if zero.
    let mut low: Vec<Option<usize>> = vec![None; n];

    for col in 0..n {
        let mut lowest = get_lowest(&reduced[col]);
        while let Some(low_row) = lowest {
            // Check if another column already has this as its low
            let mut found_collision = false;
            for prev_col in 0..col {
                if low[prev_col] == Some(low_row) {
                    // XOR prev_col into col
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

    // Extract barcodes from the reduced matrix.
    //
    // A column j with low[j] = Some(i) is a negative (destroyer) simplex:
    //   simplex order[j] (dim d) kills the H_{d-1} class created by simplex order[i].
    //
    // A column j with low[j] = None is a positive (creator) simplex.
    //   If it is NOT used as a low value by any other column, it creates an
    //   infinite H_d bar.
    let max_dim = filtration.complex.dimension;
    let mut barcodes: Vec<PersistenceBarcode> = (0..=max_dim)
        .map(|d| PersistenceBarcode {
            dimension: d,
            bars: vec![],
        })
        .collect();

    // Track which columns are used as low values
    let mut used_as_low = vec![false; n];
    for col in 0..n {
        if let Some(row) = low[col] {
            used_as_low[row] = true;
        }
    }

    for col in 0..n {
        let simp_idx = order[col];
        let dim = simplices[simp_idx].len().saturating_sub(1);

        if let Some(row) = low[col] {
            // Negative simplex: destroys a class
            let creator_idx = order[row];
            let creator_dim = simplices[creator_idx].len().saturating_sub(1);
            let birth_val = values[creator_idx];
            let death_val = values[simp_idx];
            if death_val > birth_val && creator_dim <= max_dim {
                barcodes[creator_dim].bars.push((birth_val, death_val));
            }
        } else if !used_as_low[col] {
            // Unpaired positive simplex → infinite bar
            if dim <= max_dim {
                barcodes[dim].bars.push((values[simp_idx], f64::INFINITY));
            }
        }
    }

    // Sort bars within each barcode
    for bc in &mut barcodes {
        bc.bars.sort_by(|a, b| {
            a.0.partial_cmp(&b.0)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    barcodes
}

/// Get the lowest (max) row index in a sparse column.
fn get_lowest(col: &[usize]) -> Option<usize> {
    if col.is_empty() {
        None
    } else {
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
        for &(_birth, death) in &bc.bars {
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
    let mut simplices = complex.simplices.clone();
    // Ensure sorted by dimension then lexicographic for consistent ordering
    simplices.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
    let sorted_complex = SimplicialComplex {
        simplices: simplices.clone(),
        dimension: complex.dimension,
    };
    let filtration = Filtration {
        complex: sorted_complex,
        values: vec![0.0; simplices.len()],
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

    fn make_complex(simplices: Vec<Vec<usize>>) -> SimplicialComplex {
        let mut s = simplices;
        for face in &mut s {
            face.sort_unstable();
        }
        s.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
        s.dedup();
        let dimension = s.iter().map(|f| f.len().saturating_sub(1)).max().unwrap_or(0);
        SimplicialComplex {
            simplices: s,
            dimension,
        }
    }

    #[test]
    fn test_betti_point() {
        let complex = make_complex(vec![vec![0]]);
        let betti = betti_numbers(&complex);
        assert_eq!(betti.betti, vec![1]);
    }

    #[test]
    fn test_betti_two_points() {
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
        let complex = make_complex(vec![
            vec![0],
            vec![1],
            vec![2],
            vec![0, 1],
            vec![0, 2],
            vec![1, 2],
            vec![0, 1, 2],
        ]);
        let betti = betti_numbers(&complex);
        assert_eq!(betti.betti, vec![1, 0, 0]);
    }

    #[test]
    fn test_betti_triangle_hollow() {
        let complex = make_complex(vec![
            vec![0],
            vec![1],
            vec![2],
            vec![0, 1],
            vec![0, 2],
            vec![1, 2],
        ]);
        let betti = betti_numbers(&complex);
        assert_eq!(betti.betti, vec![1, 1]);
    }

    #[test]
    fn test_betti_square_hollow() {
        let complex = make_complex(vec![
            vec![0],
            vec![1],
            vec![2],
            vec![3],
            vec![0, 1],
            vec![1, 2],
            vec![2, 3],
            vec![0, 3],
        ]);
        let betti = betti_numbers(&complex);
        assert_eq!(betti.betti, vec![1, 1]);
    }

    #[test]
    fn test_betti_tetrahedron_filled() {
        let complex = make_complex(vec![
            vec![0],
            vec![1],
            vec![2],
            vec![3],
            vec![0, 1],
            vec![0, 2],
            vec![0, 3],
            vec![1, 2],
            vec![1, 3],
            vec![2, 3],
            vec![0, 1, 2],
            vec![0, 1, 3],
            vec![0, 2, 3],
            vec![1, 2, 3],
            vec![0, 1, 2, 3],
        ]);
        let betti = betti_numbers(&complex);
        assert_eq!(betti.betti, vec![1, 0, 0, 0]);
    }

    #[test]
    fn test_betti_sphere() {
        // Hollow tetrahedron (boundary of a 3-simplex): β₀=1, β₁=0, β₂=1
        let complex = make_complex(vec![
            vec![0],
            vec![1],
            vec![2],
            vec![3],
            vec![0, 1],
            vec![0, 2],
            vec![0, 3],
            vec![1, 2],
            vec![1, 3],
            vec![2, 3],
            vec![0, 1, 2],
            vec![0, 1, 3],
            vec![0, 2, 3],
            vec![1, 2, 3],
        ]);
        let betti = betti_numbers(&complex);
        assert_eq!(betti.betti, vec![1, 0, 1]);
    }

    #[test]
    fn test_symmetric_difference() {
        assert_eq!(
            symmetric_difference(&[1, 3, 5], &[2, 3, 4]),
            vec![1, 2, 4, 5]
        );
        assert_eq!(symmetric_difference(&[1, 2], &[1, 2]), vec![]);
        assert_eq!(symmetric_difference(&[], &[1, 2]), vec![1, 2]);
    }

    #[test]
    fn test_persistence_barcode_extraction() {
        let complex = make_complex(vec![vec![0], vec![1], vec![0, 1]]);
        let filtration = Filtration {
            complex,
            values: vec![0.0, 0.0, 1.0],
        };
        let barcodes = compute_persistence(&filtration);
        assert_eq!(barcodes[0].dimension, 0);
        // H₀: one bar born at 0 that never dies, one born at 0 that dies at 1
        assert!(barcodes[0].bars.len() >= 1);
        // There should be an infinite bar
        assert!(barcodes[0].bars.iter().any(|&(b, d)| d == f64::INFINITY));
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
        assert_eq!(betti.betti, vec![2, 1]);
    }

    #[test]
    fn test_full_pipeline_chord_complex() {
        use crate::simplicial_chord::progression_to_complex;
        let c = Chord {
            notes: vec![0, 4, 7],
        };
        let g = Chord {
            notes: vec![2, 7, 11],
        };
        let complex = progression_to_complex(&[c, g]);
        let betti = betti_numbers(&complex);
        // Two filled triangles sharing vertex 7 → contractible-like: β₀=1
        assert_eq!(betti.betti[0], 1);
    }
}
