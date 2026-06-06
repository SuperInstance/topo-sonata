//! Voice-leading spaces as Vietoris-Rips complexes.
//!
//! Points in the space are chords; two chords are connected by an edge if
//! their voice-leading distance is ≤ ε. The Vietoris-Rips complex VR(ε) has
//! a simplex for every clique in this graph.
//!
//! Minimal voice-leading corresponds to geodesics in the filtration space.

use crate::persistence_tension::betti_numbers;
use crate::simplicial_chord::voice_leading_distance;
use crate::{Chord, SimplicialComplex};

/// Build a Vietoris-Rips complex from a set of chords with radius ε.
///
/// Two chords are connected if their voice-leading distance ≤ ε.
/// A set of chords forms a simplex if all pairs are within ε.
pub fn vietoris_rips(chords: &[Chord], epsilon: f64) -> SimplicialComplex {
    let n = chords.len();
    if n == 0 {
        return SimplicialComplex {
            simplices: vec![],
            dimension: 0,
        };
    }

    // Compute pairwise distances
    let mut dist = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let d = voice_leading_distance(&chords[i], &chords[j]);
            dist[i][j] = d;
            dist[j][i] = d;
        }
    }

    // Build adjacency: connected if distance ≤ epsilon
    let adjacent = |i: usize, j: usize| -> bool { dist[i][j] <= epsilon };

    let mut simplices: Vec<Vec<usize>> = Vec::new();

    // Add all 0-simplices (individual chords)
    for i in 0..n {
        simplices.push(vec![i]);
    }

    // Add all 1-simplices (edges)
    for i in 0..n {
        for j in (i + 1)..n {
            if adjacent(i, j) {
                simplices.push(vec![i, j]);
            }
        }
    }

    // Add higher simplices: check all subsets for clique property
    for size in 3..=n {
        let mut subset: Vec<usize> = (0..size).collect();
        loop {
            // Check if this subset forms a clique
            let is_clique = (0..subset.len()).all(|i| {
                (i + 1..subset.len()).all(|j| adjacent(subset[i], subset[j]))
            });
            if is_clique {
                simplices.push(subset.clone());
            }
            // Advance to next combination
            if !next_combination(&mut subset, n) {
                break;
            }
        }
    }

    simplices.sort_unstable();
    let dimension = simplices.iter().map(|s| s.len().saturating_sub(1)).max().unwrap_or(0);

    SimplicialComplex { simplices, dimension }
}

/// Advance a combination to the next one. Returns false if already the last.
fn next_combination(comb: &mut Vec<usize>, n: usize) -> bool {
    let k = comb.len();
    let mut i = k as i32 - 1;
    while i >= 0 {
        if comb[i as usize] < n - k + (i as usize) {
            comb[i as usize] += 1;
            for j in (i as usize + 1)..k {
                comb[j] = comb[j - 1] + 1;
            }
            return true;
        }
        i -= 1;
    }
    false
}

/// Find the minimal voice-leading path between two chords through a space.
///
/// Given a set of intermediate chords, finds the shortest path (by voice-leading
/// distance) from `start` to `end` using chords from `space`.
///
/// Uses Dijkstra's algorithm.
pub fn minimal_voice_leading(
    start: &Chord,
    end: &Chord,
    space: &[Chord],
) -> Vec<(Chord, f64)> {
    if space.is_empty() {
        let d = voice_leading_distance(start, end);
        return vec![(start.clone(), 0.0), (end.clone(), d)];
    }

    // Build full vertex set: start + space + end
    let mut vertices = vec![start.clone()];
    vertices.extend_from_slice(space);
    vertices.push(end.clone());
    let n = vertices.len();
    let start_idx = 0;
    let end_idx = n - 1;

    // Dijkstra
    let mut dist = vec![f64::INFINITY; n];
    let mut prev = vec![None::<usize>; n];
    let mut visited = vec![false; n];
    dist[start_idx] = 0.0;

    for _ in 0..n {
        // Find unvisited vertex with minimum distance
        let u = (0..n)
            .filter(|&i| !visited[i])
            .min_by(|&a, &b| dist[a].partial_cmp(&dist[b]).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap();

        visited[u] = true;
        if u == end_idx {
            break;
        }

        for v in 0..n {
            if visited[v] {
                continue;
            }
            let d = voice_leading_distance(&vertices[u], &vertices[v]);
            let new_dist = dist[u] + d;
            if new_dist < dist[v] {
                dist[v] = new_dist;
                prev[v] = Some(u);
            }
        }
    }

    // Reconstruct path
    let mut path = Vec::new();
    let mut current = end_idx;
    while current != start_idx {
        path.push((vertices[current].clone(), dist[current]));
        current = prev[current].unwrap_or(start_idx);
        if current == start_idx {
            break;
        }
    }
    path.push((vertices[start_idx].clone(), 0.0));
    path.reverse();

    path
}

/// Compute the topological complexity of a voice-leading space.
///
/// Returns the Betti numbers of the Vietoris-Rips complex.
pub fn voice_leading_topology(chords: &[Chord], epsilon: f64) -> crate::BettiSequence {
    let complex = vietoris_rips(chords, epsilon);
    betti_numbers(&complex)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vietoris_rips_small_epsilon() {
        let c = Chord { notes: vec![0, 4, 7] };
        let dm = Chord { notes: vec![2, 5, 9] };
        // With small epsilon, only isolated points
        let vr = vietoris_rips(&[c.clone(), dm.clone()], 0.5);
        assert_eq!(vr.simplices.len(), 2); // just two points
        assert_eq!(vr.dimension, 0);
    }

    #[test]
    fn test_vietoris_rips_large_epsilon() {
        let c = Chord { notes: vec![0, 4, 7] };
        let am = Chord { notes: vec![0, 4, 9] };
        // Distance between C and Am is 2.0
        let vr = vietoris_rips(&[c, am], 3.0);
        assert!(vr.simplices.len() >= 3); // two points + edge
        assert_eq!(vr.dimension, 1);
    }

    #[test]
    fn test_vietoris_rips_triangle() {
        let c = Chord { notes: vec![0, 4, 7] };
        let am = Chord { notes: vec![0, 4, 9] };
        let dm = Chord { notes: vec![2, 5, 9] };
        // C↔Am: 2.0, C↔Dm: larger, Am↔Dm: some distance
        let vr = vietoris_rips(&[c, am, dm], 6.0);
        // With epsilon=6, all pairs should connect → triangle
        assert!(vr.simplices.iter().any(|s| s.len() == 3));
    }

    #[test]
    fn test_minimal_voice_leading_direct() {
        let c = Chord { notes: vec![0, 4, 7] };
        let g = Chord { notes: vec![2, 7, 11] };
        let path = minimal_voice_leading(&c, &g, &[]);
        assert_eq!(path.len(), 2);
        assert_eq!(path[0].0, c);
        assert_eq!(path[1].0, g);
        assert_eq!(path[1].1, 3.0); // voice-leading distance C→G
    }

    #[test]
    fn test_minimal_voice_leading_via_intermediate() {
        let c = Chord { notes: vec![0, 4, 7] };
        let am = Chord { notes: vec![0, 4, 9] };
        let g = Chord { notes: vec![2, 7, 11] };
        let path = minimal_voice_leading(&c, &g, &[am.clone()]);
        // Direct: C→G = 3.0
        // Via Am: C→Am = 2.0, Am→G = ?
        // Should pick the shortest route
        assert!(path.len() >= 2);
    }

    #[test]
    fn test_voice_leading_topology() {
        let c = Chord { notes: vec![0, 4, 7] };
        let am = Chord { notes: vec![0, 4, 9] };
        let dm = Chord { notes: vec![2, 5, 9] };
        let f = Chord { notes: vec![5, 9, 0] };
        let betti = voice_leading_topology(&[c, am, dm, f], 3.0);
        // With reasonable epsilon, should be connected → β₀ = 1
        assert_eq!(betti.betti[0], 1);
    }

    #[test]
    fn test_empty_vr() {
        let vr = vietoris_rips(&[], 1.0);
        assert_eq!(vr.simplices.len(), 0);
    }

    #[test]
    fn test_single_point_vr() {
        let c = Chord { notes: vec![0, 4, 7] };
        let vr = vietoris_rips(&[c], 1.0);
        assert_eq!(vr.simplices.len(), 1); // just the point itself
    }
}
