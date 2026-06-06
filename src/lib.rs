//! # topo-sonata
//!
//! Musical compositions as simplicial complexes where **chords are simplices**
//! and **persistent homology detects holes in harmonic space**.
//!
//! Unresolved tension is literally a topological hole. Genre classification
//! via Betti numbers.

pub mod simplicial_chord;
pub mod persistence_tension;
pub mod genre_betti;
pub mod filtration_cadence;
pub mod voice_lead;
pub mod contrapuntal;

use serde::{Deserialize, Serialize};

// ── Core types ──────────────────────────────────────────────────────────────

/// A chord represented as a collection of MIDI pitch classes (0–11).
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Chord {
    /// Sorted, deduplicated pitch-class notes (0–11).
    pub notes: Vec<u32>,
}

/// A finite abstract simplicial complex.
///
/// `simplices[i]` is the set of vertex indices forming simplex *i*,
/// sorted in increasing order. The complex is closed under subsets.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SimplicialComplex {
    /// All simplices, each represented as a sorted vec of vertex indices.
    pub simplices: Vec<Vec<usize>>,
    /// Maximum dimension of any simplex in the complex.
    pub dimension: usize,
}

/// A single bar in a persistence barcode.
///
/// The bar `[birth, death)` represents a homology class that appears at
/// filtration value `birth` and disappears at `death`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PersistenceBarcode {
    /// Homology dimension (0 = H₀, 1 = H₁, …).
    pub dimension: usize,
    /// List of (birth, death) intervals.
    pub bars: Vec<(f64, f64)>,
}

/// Betti numbers [β₀, β₁, β₂, …] summarising the topology of a complex.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct BettiSequence {
    pub betti: Vec<usize>,
}

/// A filtered simplicial complex — each simplex is assigned a filtration value.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Filtration {
    pub complex: SimplicialComplex,
    /// Filtration value for each simplex (same index as `complex.simplices`).
    pub values: Vec<f64>,
}

// ── Re-exports ──────────────────────────────────────────────────────────────

pub use simplicial_chord::{chord_to_simplex, progression_to_complex, voice_leading_distance};
pub use persistence_tension::{compute_persistence, betti_numbers, betti_from_barcodes};
pub use genre_betti::{classify_genre, Genre, genre_betti_fingerprint};
pub use filtration_cadence::{build_filtration, generate_cadence, cadence_kills_h1};
pub use voice_lead::{vietoris_rips, minimal_voice_leading};
pub use contrapuntal::{check_parallel_fifths, check_parallel_octaves, check_counterpoint};
