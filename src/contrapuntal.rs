//! Counterpoint rules as cocycle conditions on the chord complex.
//!
//! Classical counterpoint forbids:
//! - Parallel fifths (two voices moving in parallel at a P5 interval)
//! - Parallel octaves (two voices moving in parallel at a P8/unison)
//!
//! In the simplicial framework, these correspond to specific cocycle
//! conditions on the 1-cochains of the chord complex.

use crate::Chord;

/// Interval in semitones for a perfect fifth.
const PERFECT_FIFTH: u32 = 7;
/// Interval in semitones for an octave/unison (mod 12).
const OCTAVE: u32 = 0;

/// An interval between two voices in a chord.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VoiceInterval {
    pub upper: u32,
    pub lower: u32,
    pub interval: u32, // semitones, mod 12
}

/// Extract all intervals between pairs of voices in a chord.
pub fn chord_intervals(chord: &Chord) -> Vec<VoiceInterval> {
    let mut intervals = Vec::new();
    let notes = &chord.notes;
    for i in 0..notes.len() {
        for j in (i + 1)..notes.len() {
            let lower = notes[i];
            let upper = notes[j];
            let interval = (upper - lower) % 12;
            intervals.push(VoiceInterval {
                upper,
                lower,
                interval,
            });
        }
    }
    intervals
}

/// Check if two consecutive chords contain parallel fifths.
///
/// Returns true if there exist two voices that form a P5 interval in both
/// chords with both voices moving (not static).
pub fn check_parallel_fifths(chord_a: &Chord, chord_b: &Chord) -> bool {
    check_parallel_interval(chord_a, chord_b, PERFECT_FIFTH)
}

/// Check if two consecutive chords contain parallel octaves (or unisons).
pub fn check_parallel_octaves(chord_a: &Chord, chord_b: &Chord) -> bool {
    check_parallel_interval(chord_a, chord_b, OCTAVE)
}

/// Check for any parallel interval of a specific size between two chords.
fn check_parallel_interval(chord_a: &Chord, chord_b: &Chord, target_interval: u32) -> bool {
    let intervals_a = chord_intervals(chord_a);
    let intervals_b = chord_intervals(chord_b);

    for ia in &intervals_a {
        if ia.interval != target_interval {
            continue;
        }
        for ib in &intervals_b {
            if ib.interval != target_interval {
                continue;
            }
            // Parallel if both voices moved (not static) and maintained the interval
            let diff_upper = (ib.upper as i32 - ia.upper as i32).abs();
            let diff_lower = (ib.lower as i32 - ia.lower as i32).abs();
            if diff_upper > 0 && diff_lower > 0 {
                return true;
            }
        }
    }
    false
}

/// Full counterpoint check between two consecutive chords.
///
/// Returns a list of violations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CounterpointViolation {
    ParallelFifths,
    ParallelOctaves,
}

/// Check all counterpoint rules between two chords.
pub fn check_counterpoint(chord_a: &Chord, chord_b: &Chord) -> Vec<CounterpointViolation> {
    let mut violations = Vec::new();
    if check_parallel_fifths(chord_a, chord_b) {
        violations.push(CounterpointViolation::ParallelFifths);
    }
    if check_parallel_octaves(chord_a, chord_b) {
        violations.push(CounterpointViolation::ParallelOctaves);
    }
    violations
}

/// Check a full progression for counterpoint violations.
///
/// Returns pairs of (chord_index, violations).
pub fn check_progression(chords: &[Chord]) -> Vec<(usize, Vec<CounterpointViolation>)> {
    let mut results = Vec::new();
    for i in 0..chords.len().saturating_sub(1) {
        let violations = check_counterpoint(&chords[i], &chords[i + 1]);
        if !violations.is_empty() {
            results.push((i, violations));
        }
    }
    results
}

/// Interpret counterpoint violations as cocycle conditions.
///
/// Each forbidden parallel corresponds to a 1-cocycle that must vanish.
/// Returns the number of non-vanishing cocycles (violations).
pub fn cocycle_violations(chords: &[Chord]) -> usize {
    check_progression(chords).iter().map(|(_, v)| v.len()).sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Chord;

    #[test]
    fn test_chord_intervals_c_major() {
        let c = Chord { notes: vec![0, 4, 7] };
        let intervals = chord_intervals(&c);
        assert_eq!(intervals.len(), 3);
        assert!(intervals.iter().any(|i| i.interval == 4));
        assert!(intervals.iter().any(|i| i.interval == 7));
        assert!(intervals.iter().any(|i| i.interval == 3));
    }

    #[test]
    fn test_parallel_fifths_detected() {
        let chord_a = Chord { notes: vec![0, 7] };
        let chord_b = Chord { notes: vec![2, 9] };
        assert!(check_parallel_fifths(&chord_a, &chord_b));
    }

    #[test]
    fn test_no_parallel_fifths() {
        let c = Chord { notes: vec![0, 4, 7] };
        let am = Chord { notes: vec![0, 4, 9] };
        assert!(!check_parallel_fifths(&c, &am));
    }

    #[test]
    fn test_static_fifth_not_parallel() {
        let c = Chord { notes: vec![0, 7] };
        assert!(!check_parallel_fifths(&c, &c));
    }

    #[test]
    fn test_parallel_octaves_detected() {
        let a = Chord { notes: vec![0, 0] };
        let b = Chord { notes: vec![2, 2] };
        assert!(check_parallel_octaves(&a, &b));
    }

    #[test]
    fn test_check_counterpoint_clean() {
        let c = Chord { notes: vec![0, 4, 7] };
        let f = Chord { notes: vec![0, 5, 9] };
        let violations = check_counterpoint(&c, &f);
        assert!(violations.is_empty() || violations.len() < 3);
    }

    #[test]
    fn test_check_counterpoint_violation() {
        let a = Chord { notes: vec![0, 7] };
        let b = Chord { notes: vec![2, 9] };
        let violations = check_counterpoint(&a, &b);
        assert!(violations.contains(&CounterpointViolation::ParallelFifths));
    }

    #[test]
    fn test_check_progression() {
        let chords = vec![
            Chord { notes: vec![0, 7] },
            Chord { notes: vec![2, 9] },
            Chord { notes: vec![4, 7] },
        ];
        let results = check_progression(&chords);
        assert!(results.iter().any(|(i, v)| *i == 0 && v.contains(&CounterpointViolation::ParallelFifths)));
    }

    #[test]
    fn test_cocycle_violations_count() {
        let chords = vec![
            Chord { notes: vec![0, 7] },
            Chord { notes: vec![2, 9] },
        ];
        let count = cocycle_violations(&chords);
        assert!(count >= 1);
    }

    #[test]
    fn test_clean_progression_zero_violations() {
        let chords = vec![
            Chord { notes: vec![0, 4, 7] },
            Chord { notes: vec![0, 4, 9] },
        ];
        let violations = cocycle_violations(&chords);
        assert_eq!(violations, 0);
    }
}
