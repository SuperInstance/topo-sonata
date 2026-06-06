//! Chords as simplices and chord progressions as simplicial complexes.
//!
//! | Musical concept | Topological concept   |
//! |-----------------|-----------------------|
//! | Note (pitch class) | Vertex (0-simplex) |
//! | Interval (dyad)   | Edge (1-simplex)   |
//! | Triad             | Triangle (2-simplex)|
//! | Seventh chord     | Tetrahedron (3-simplex) |

use crate::{Chord, SimplicialComplex};

/// Normalise a chord: sort and deduplicate pitch classes mod 12.
pub fn normalise_chord(notes: &[u32]) -> Vec<u32> {
    let mut pcs: Vec<u32> = notes.iter().map(|&n| n % 12).collect();
    pcs.sort_unstable();
    pcs.dedup();
    pcs
}

/// Convert a chord into a simplex (sorted vertex list).
///
/// Pitch classes (0–11) serve directly as vertex indices.
pub fn chord_to_simplex(chord: &Chord) -> Vec<usize> {
    chord.notes.iter().map(|&pc| pc as usize).collect()
}

/// Voice-leading distance between two chords (minimal total semitone movement).
///
/// Brute-force all permutations of b matched to a.
/// For chords of unequal size we pad the smaller with its nearest note.
pub fn voice_leading_distance(a: &Chord, b: &Chord) -> f64 {
    let na = normalise_chord(&a.notes);
    let nb = normalise_chord(&b.notes);

    let n = na.len().max(nb.len());
    let padded_a = pad_chord(&na, n);
    let padded_b = pad_chord(&nb, n);

    let mut best = f64::INFINITY;
    let mut perm: Vec<usize> = (0..n).collect();
    loop {
        let dist: f64 = perm
            .iter()
            .enumerate()
            .map(|(i, &j)| pc_distance(padded_a[i], padded_b[j]))
            .sum();
        if dist < best {
            best = dist;
        }
        if !next_permutation(&mut perm) {
            break;
        }
    }
    best
}

/// Pad a chord to size `n` by repeating the last note.
fn pad_chord(notes: &[u32], n: usize) -> Vec<u32> {
    if notes.is_empty() {
        return vec![0; n];
    }
    let mut result = notes.to_vec();
    while result.len() < n {
        result.push(*result.last().unwrap());
    }
    result
}

/// Minimal distance between two pitch classes on the circle of 12.
fn pc_distance(a: u32, b: u32) -> f64 {
    let d = (a as i32 - b as i32).abs() % 12;
    let d = d.min(12 - d);
    d as f64
}

/// Next lexicographic permutation in place. Returns false if already last.
fn next_permutation(arr: &mut [usize]) -> bool {
    let n = arr.len();
    if n < 2 {
        return false;
    }
    let mut i = n - 2;
    while arr[i] >= arr[i + 1] {
        if i == 0 {
            return false;
        }
        i -= 1;
    }
    let mut j = n - 1;
    while arr[j] <= arr[i] {
        j -= 1;
    }
    arr.swap(i, j);
    arr[i + 1..].reverse();
    true
}

/// Build a simplicial complex from a chord progression.
///
/// Each chord becomes a maximal simplex. All subsets are included
/// to satisfy closure under faces. Pitch classes (0–11) are used
/// directly as vertex indices.
pub fn progression_to_complex(chords: &[Chord]) -> SimplicialComplex {
    let mut simplices: Vec<Vec<usize>> = Vec::new();

    for chord in chords {
        let vertices: Vec<usize> = chord.notes.iter().map(|&pc| pc as usize).collect();
        add_all_faces(&vertices, &mut simplices);
    }

    // Deduplicate
    simplices.sort_unstable();
    simplices.dedup();

    // Compute dimension
    let dimension = simplices
        .iter()
        .map(|s| s.len().saturating_sub(1))
        .max()
        .unwrap_or(0);

    SimplicialComplex {
        simplices,
        dimension,
    }
}

/// Add all faces (subsets) of a simplex to the list, including the simplex itself.
fn add_all_faces(vertices: &[usize], simplices: &mut Vec<Vec<usize>>) {
    let n = vertices.len();
    if n == 0 {
        return;
    }
    let count = 1u32 << n;
    for mask in 1..count {
        let mut face: Vec<usize> = Vec::new();
        for i in 0..n {
            if mask & (1 << i) != 0 {
                face.push(vertices[i]);
            }
        }
        face.sort_unstable();
        simplices.push(face);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Chord;

    #[test]
    fn test_normalise_chord() {
        assert_eq!(normalise_chord(&[0, 4, 7]), vec![0, 4, 7]);
        assert_eq!(normalise_chord(&[12, 16, 19]), vec![0, 4, 7]);
        assert_eq!(normalise_chord(&[7, 4, 0]), vec![0, 4, 7]);
        assert_eq!(normalise_chord(&[0, 0, 0]), vec![0]);
    }

    #[test]
    fn test_chord_to_simplex() {
        let chord = Chord {
            notes: vec![0, 4, 7],
        };
        let simplex = chord_to_simplex(&chord);
        assert_eq!(simplex, vec![0, 4, 7]);
    }

    #[test]
    fn test_voice_leading_c_to_am() {
        let c = Chord {
            notes: vec![0, 4, 7],
        };
        let am = Chord {
            notes: vec![0, 4, 9],
        };
        let dist = voice_leading_distance(&c, &am);
        assert_eq!(dist, 2.0);
    }

    #[test]
    fn test_voice_leading_c_to_g() {
        let c = Chord {
            notes: vec![0, 4, 7],
        };
        let g = Chord {
            notes: vec![2, 7, 11],
        };
        let dist = voice_leading_distance(&c, &g);
        assert_eq!(dist, 3.0);
    }

    #[test]
    fn test_pc_distance() {
        assert_eq!(pc_distance(0, 0), 0.0);
        assert_eq!(pc_distance(0, 1), 1.0);
        assert_eq!(pc_distance(0, 11), 1.0);
        assert_eq!(pc_distance(0, 6), 6.0);
    }

    #[test]
    fn test_progression_to_complex_single_chord() {
        let c_major = Chord {
            notes: vec![0, 4, 7],
        };
        let complex = progression_to_complex(&[c_major]);
        assert_eq!(complex.simplices.len(), 7);
        assert_eq!(complex.dimension, 2);
        assert!(complex.simplices.contains(&vec![0]));
        assert!(complex.simplices.contains(&vec![0, 4, 7]));
    }

    #[test]
    fn test_progression_to_complex_two_chords() {
        let c = Chord {
            notes: vec![0, 4, 7],
        };
        let g = Chord {
            notes: vec![2, 7, 11],
        };
        let complex = progression_to_complex(&[c, g]);
        // C major: 7 simplices, G major: 7 simplices, shared: {7} — deduped
        assert!(complex.simplices.len() >= 13);
        assert!(complex.simplices.contains(&vec![0, 4, 7]));
        assert!(complex.simplices.contains(&vec![2, 7, 11]));
    }

    #[test]
    fn test_next_permutation() {
        let mut arr = vec![0, 1, 2];
        assert!(next_permutation(&mut arr));
        assert_eq!(arr, vec![0, 2, 1]);
        assert!(next_permutation(&mut arr));
        assert_eq!(arr, vec![1, 0, 2]);
    }

    #[test]
    fn test_dyad_simplex_dimension() {
        let dyad = Chord {
            notes: vec![0, 7],
        };
        let complex = progression_to_complex(&[dyad]);
        assert_eq!(complex.dimension, 1);
        assert_eq!(complex.simplices.len(), 3); // {0}, {7}, {0,7}
    }
}
