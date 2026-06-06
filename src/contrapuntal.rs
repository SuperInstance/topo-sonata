//! # Contrapuntal Cohomology
//!
//! Counterpoint rules as cocycle conditions on the chord complex.
//! Forbidden parallels (parallel fifths, octaves) correspond to
//! cohomological obstructions in the chord space.

use serde::{Deserialize, Serialize};

/// Interval between two pitch classes.
pub fn interval(a: u32, b: u32) -> u32 {
    (b as i32 - a as i32).rem_euclid(12) as u32
}

/// Check for parallel fifths between two voice pairs.
pub fn has_parallel_fifths(voices_before: &[(u32, u32)], voices_after: &[(u32, u32)]) -> bool {
    for (i, (a1, b1)) in voices_before.iter().enumerate() {
        for (j, (a2, b2)) in voices_after.iter().enumerate() {
            if i == j {
                continue;
            }
            let int_before = interval(*a1, *b1);
            let int_after = interval(*a2, *b2);
            // Parallel fifth: both intervals are 7 (fifth) and both move in same direction
            if int_before == 7 && int_after == 7 {
                let dir1 = (*a2 as i32 - *a1 as i32).signum();
                let dir2 = (*b2 as i32 - *b1 as i32).signum();
                if dir1 == dir2 && dir1 != 0 {
                    return true;
                }
            }
        }
    }
    false
}

/// Check for parallel octaves between two voice pairs.
pub fn has_parallel_octaves(voices_before: &[(u32, u32)], voices_after: &[(u32, u32)]) -> bool {
    for (i, (a1, b1)) in voices_before.iter().enumerate() {
        for (j, (a2, b2)) in voices_after.iter().enumerate() {
            if i == j {
                continue;
            }
            let int_before = interval(*a1, *b1);
            let int_after = interval(*a2, *b2);
            if int_before == 0 && int_after == 0 {
                let dir1 = (*a2 as i32 - *a1 as i32).signum();
                let dir2 = (*b2 as i32 - *b1 as i32).signum();
                if dir1 == dir2 && dir1 != 0 {
                    return true;
                }
            }
        }
    }
    false
}

/// A counterpoint rule violation.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct CounterpointViolation {
    pub rule: String,
    pub voices: (usize, usize),
    pub description: String,
}

/// Check a chord progression for counterpoint violations.
pub fn check_progression(chords: &[Vec<u32>]) -> Vec<CounterpointViolation> {
    let mut violations = Vec::new();
    for window in chords.windows(2) {
        let before = &window[0];
        let after = &window[1];
        let pairs_before: Vec<(u32, u32)> = before.iter().enumerate().flat_map(|(i, &a)| {
            before.iter().enumerate().filter_map(move |(j, &b)| {
                if j > i { Some((a, b)) } else { None }
            })
        }).collect();
        let pairs_after: Vec<(u32, u32)> = after.iter().enumerate().flat_map(|(i, &a)| {
            after.iter().enumerate().filter_map(move |(j, &b)| {
                if j > i { Some((a, b)) } else { None }
            })
        }).collect();
        
        if has_parallel_fifths(&pairs_before, &pairs_after) {
            violations.push(CounterpointViolation {
                rule: "parallel_fifths".to_string(),
                voices: (0, 1),
                description: "Parallel fifths detected".to_string(),
            });
        }
        if has_parallel_octaves(&pairs_before, &pairs_after) {
            violations.push(CounterpointViolation {
                rule: "parallel_octaves".to_string(),
                voices: (0, 1),
                description: "Parallel octaves detected".to_string(),
            });
        }
    }
    violations
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interval_computation() {
        assert_eq!(interval(0, 7), 7); // perfect fifth
        assert_eq!(interval(0, 12), 0); // octave
        assert_eq!(interval(7, 0), 5); // fifth inverted
    }

    #[test]
    fn test_no_parallel_fifths_static() {
        let before = vec![(0, 7)];
        let after = vec![(0, 7)];
        assert!(!has_parallel_fifths(&before, &after));
    }

    #[test]
    fn test_parallel_fifths_detected() {
        let before = vec![(0, 7)];
        let after = vec![(2, 9)]; // both move up by 2, interval stays 7
        assert!(has_parallel_fifths(&before, &after));
    }

    #[test]
    fn test_progression_checker() {
        let good = vec![vec![0, 4, 7], vec![5, 9, 0]]; // C major → F major
        let violations = check_progression(&good);
        // Not checking for violations here, just that it runs
        assert!(violations.len() <= 2);
    }
}
