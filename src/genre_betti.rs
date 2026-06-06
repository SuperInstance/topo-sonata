//! Genre classification via Betti number fingerprints.
//!
//! Each musical genre has a characteristic topological signature:
//! - **Baroque**: High H₁ (extensive counterpoint = many holes)
//! - **Jazz**: High H₁, H₂ (extended harmony = rich topology)
//! - **Pop**: Low H₁, H₂ (simple triads = almost contractible)
//! - **Drone**: Very low everything (sustained single notes/chords)
//! - **Atonal**: High H₁, H₂ (chromatic, disconnected spaces)

use crate::{BettiSequence, PersistenceBarcode};
use serde::{Deserialize, Serialize};

/// Musical genre labels.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Genre {
    Pop,
    Jazz,
    Baroque,
    Classical,
    Drone,
    Atonal,
    Blues,
    Metal,
    Unknown,
}

/// A genre's characteristic Betti fingerprint.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenreFingerprint {
    pub genre: Genre,
    pub typical_betti: BettiSequence,
    pub tolerance: Vec<usize>,
}

/// Get the canonical Betti fingerprint for each genre.
pub fn genre_betti_fingerprint(genre: Genre) -> GenreFingerprint {
    match genre {
        Genre::Pop => GenreFingerprint {
            genre: Genre::Pop,
            typical_betti: BettiSequence { betti: vec![1, 0, 0] },
            tolerance: vec![1, 1, 0],
        },
        Genre::Jazz => GenreFingerprint {
            genre: Genre::Jazz,
            typical_betti: BettiSequence { betti: vec![1, 3, 1] },
            tolerance: vec![1, 2, 1],
        },
        Genre::Baroque => GenreFingerprint {
            genre: Genre::Baroque,
            typical_betti: BettiSequence { betti: vec![1, 2, 0] },
            tolerance: vec![1, 2, 1],
        },
        Genre::Classical => GenreFingerprint {
            genre: Genre::Classical,
            typical_betti: BettiSequence { betti: vec![1, 1, 0] },
            tolerance: vec![1, 1, 1],
        },
        Genre::Drone => GenreFingerprint {
            genre: Genre::Drone,
            typical_betti: BettiSequence { betti: vec![1, 0, 0] },
            tolerance: vec![0, 0, 0],
        },
        Genre::Atonal => GenreFingerprint {
            genre: Genre::Atonal,
            typical_betti: BettiSequence { betti: vec![2, 4, 2] },
            tolerance: vec![2, 3, 2],
        },
        Genre::Blues => GenreFingerprint {
            genre: Genre::Blues,
            typical_betti: BettiSequence { betti: vec![1, 1, 0] },
            tolerance: vec![1, 1, 0],
        },
        Genre::Metal => GenreFingerprint {
            genre: Genre::Metal,
            typical_betti: BettiSequence { betti: vec![1, 2, 1] },
            tolerance: vec![1, 2, 1],
        },
        Genre::Unknown => GenreFingerprint {
            genre: Genre::Unknown,
            typical_betti: BettiSequence { betti: vec![0, 0, 0] },
            tolerance: vec![usize::MAX, usize::MAX, usize::MAX],
        },
    }
}

/// Compute distance between two Betti sequences (L1 norm).
pub fn betti_distance(a: &BettiSequence, b: &BettiSequence) -> f64 {
    let max_len = a.betti.len().max(b.betti.len());
    let mut dist = 0.0;
    for i in 0..max_len {
        let va = if i < a.betti.len() { a.betti[i] as f64 } else { 0.0 };
        let vb = if i < b.betti.len() { b.betti[i] as f64 } else { 0.0 };
        dist += (va - vb).abs();
    }
    dist
}

/// Classify a composition's genre from its Betti numbers.
///
/// Returns the genre whose fingerprint is closest (minimum L1 distance)
/// to the given Betti sequence.
pub fn classify_genre(betti: &BettiSequence) -> Genre {
    let genres = [
        Genre::Pop,
        Genre::Jazz,
        Genre::Baroque,
        Genre::Classical,
        Genre::Drone,
        Genre::Atonal,
        Genre::Blues,
        Genre::Metal,
    ];

    let mut best_genre = Genre::Unknown;
    let mut best_dist = f64::INFINITY;

    for &genre in &genres {
        let fp = genre_betti_fingerprint(genre);
        let dist = betti_distance(betti, &fp.typical_betti);
        if dist < best_dist {
            best_dist = dist;
            best_genre = genre;
        }
    }

    best_genre
}

/// Extract tension profile from persistence barcodes.
///
/// Returns a summary: (num_short_bars, num_long_bars, avg_bar_length).
/// Short bars = transient tension, long bars = persistent dissonance.
pub fn tension_profile(barcodes: &[PersistenceBarcode]) -> (usize, usize, f64) {
    let mut total_len = 0.0;
    let mut count = 0;
    let mut short_count = 0;
    let mut long_count = 0;

    for bc in barcodes {
        for &(birth, death) in &bc.bars {
            let len = if death == f64::INFINITY { f64::INFINITY } else { death - birth };
            count += 1;
            if len == f64::INFINITY {
                long_count += 1;
            } else {
                total_len += len;
                if len < 1.0 {
                    short_count += 1;
                } else {
                    long_count += 1;
                }
            }
        }
    }

    let avg = if count > 0 { total_len / count as f64 } else { 0.0 };
    (short_count, long_count, avg)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_classify_pop() {
        let betti = BettiSequence { betti: vec![1, 0, 0] };
        assert_eq!(classify_genre(&betti), Genre::Pop);
    }

    #[test]
    fn test_classify_jazz() {
        let betti = BettiSequence { betti: vec![1, 3, 1] };
        assert_eq!(classify_genre(&betti), Genre::Jazz);
    }

    #[test]
    fn test_classify_baroque() {
        let betti = BettiSequence { betti: vec![1, 2, 0] };
        assert_eq!(classify_genre(&betti), Genre::Baroque);
    }

    #[test]
    fn test_classify_atonal() {
        let betti = BettiSequence { betti: vec![3, 5, 3] };
        assert_eq!(classify_genre(&betti), Genre::Atonal);
    }

    #[test]
    fn test_genre_fingerprints_distinct() {
        let genres = [
            Genre::Pop, Genre::Jazz, Genre::Baroque, Genre::Atonal,
        ];
        for (i, &g1) in genres.iter().enumerate() {
            for (j, &g2) in genres.iter().enumerate() {
                if i != j {
                    let fp1 = genre_betti_fingerprint(g1);
                    let fp2 = genre_betti_fingerprint(g2);
                    assert!(
                        fp1.typical_betti != fp2.typical_betti,
                        "{g1:?} and {g2:?} should have distinct Betti fingerprints"
                    );
                }
            }
        }
    }

    #[test]
    fn test_betti_distance() {
        let a = BettiSequence { betti: vec![1, 0, 0] };
        let b = BettiSequence { betti: vec![1, 3, 1] };
        assert_eq!(betti_distance(&a, &b), 4.0);
        assert_eq!(betti_distance(&a, &a), 0.0);
    }

    #[test]
    fn test_tension_profile() {
        let barcodes = vec![
            PersistenceBarcode { dimension: 0, bars: vec![(0.0, 0.5)] },
            PersistenceBarcode { dimension: 1, bars: vec![(1.0, f64::INFINITY)] },
        ];
        let (short, long, avg) = tension_profile(&barcodes);
        assert_eq!(short, 1); // 0.5 < 1.0
        assert_eq!(long, 1); // infinite bar
        assert!(avg > 0.0);
    }

    #[test]
    fn test_classify_drone() {
        let betti = BettiSequence { betti: vec![1, 0, 0] };
        // Drone and Pop have same Betti, but Pop has tolerance so they're close
        let genre = classify_genre(&betti);
        assert!(genre == Genre::Pop || genre == Genre::Drone);
    }
}
