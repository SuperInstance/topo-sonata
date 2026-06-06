# topo-sonata

**Musical compositions as simplicial complexes — persistent homology detects harmonic holes.**

```
topo-sonata
├── Chords are simplices (vertices = notes, edges = intervals, triangles = triads)
├── Chord progressions are simplicial complexes filtered by voice-leading distance
├── Persistent homology detects "holes" in harmonic space
├── Unresolved tension = a topological hole
└── Genre classification via Betti numbers
```

## The Big Idea

Music has topology. A chord progression traces a path through harmonic space. When that path doesn't close — when there's unresolved dissonance — it creates a literal hole in the simplicial complex. Persistent homology finds these holes. The Betti numbers tell you how many there are, and the persistence barcodes tell you how stubborn they are.

**Long bar = persistent dissonance. A perfect authentic cadence kills H₁.**

## Modules

### `simplicial_chord` — Chords as Simplices

| Musical concept       | Topological concept    |
|-----------------------|------------------------|
| Note (pitch class)    | Vertex (0-simplex)     |
| Interval (dyad)       | Edge (1-simplex)       |
| Triad                 | Triangle (2-simplex)   |
| Seventh chord         | Tetrahedron (3-simplex)|

Chord progressions become simplicial complexes. Voice-leading distance provides the filtration.

### `persistence_tension` — Persistent Homology of Harmonic Space

Computes H₀, H₁, H₂ via boundary matrix reduction (Gaussian elimination from scratch — no external math deps). Betti barcodes serve as tension profiles:

- **H₀** = number of connected components (harmonic regions)
- **H₁** = number of harmonic loops (unresolved cycles of tension)
- **H₂** = enclosed harmonic volumes (deep structural dissonance)

### `genre_betti` — Genre as Equivalence Class of Betti Sequences

| Genre     | Typical β  | Why                                                    |
|-----------|------------|--------------------------------------------------------|
| Pop       | [1, 0, 0]  | Simple triads, almost contractible                      |
| Drone     | [1, 0, 0]  | Sustained single chord, no movement                     |
| Blues     | [1, 1, 0]  | Dominant 7ths create one persistent loop                |
| Classical | [1, 1, 0]  | Functional harmony with controlled tension              |
| Baroque   | [1, 2, 0]  | Extensive counterpoint = many independent lines          |
| Metal     | [1, 2, 1]  | Power chords + chromaticism = richer topology           |
| Jazz      | [1, 3, 1]  | Extended harmony (9ths, 13ths) = complex simplicial structure |
| Atonal    | [2, 4, 2]  | Chromatic, disconnected — topological chaos             |

Classification uses minimum L¹ distance to canonical Betti fingerprints.

### `filtration_cadence` — Cadences That Kill Homology Classes

The **perfect authentic cadence** is a topological operation: it adds simplices that fill the 1-dimensional holes, sending H₁ → 0. The V→I resolution isn't just a convention — it's the simplest chord change that provably eliminates all loops in the harmonic space.

### `voice_lead` — Voice-Leading as Vietoris-Rips Complexes

Points in the space are chords. Two chords are ε-connected if their voice-leading distance ≤ ε. The Vietoris-Rips complex VR(ε) has a simplex for every clique. Minimal voice-leading = geodesic in filtration space (Dijkstra on the VR complex).

### `contrapuntal` — Counterpoint as Cocycle Conditions

Forbidden parallels (P5, P8) are cocycle obstructions on the chord complex. Each violation is a non-vanishing cocycle. The total count of violations measures contrapuntal "incorrectness" topologically.

## Core Types

```rust
struct Chord { notes: Vec<u32> }                                    // MIDI pitch classes
struct SimplicialComplex { simplices: Vec<Vec<usize>>, dimension: usize }
struct PersistenceBarcode { dimension: usize, bars: Vec<(f64, f64)> }
struct BettiSequence { betti: Vec<usize> }                          // [H₀, H₁, H₂]
struct Filtration { complex: SimplicialComplex, values: Vec<f64> }
```

All public types derive `Serialize`/`Deserialize` via serde.

## Quick Start

```rust
use topo_sonata::*;

// Build a chord progression
let c_major  = Chord { notes: vec![0, 4, 7] };
let a_minor  = Chord { notes: vec![0, 4, 9] };

// Create a simplicial complex
let complex = progression_to_complex(&[c_major.clone(), a_minor.clone()]);

// Compute Betti numbers
let betti = betti_numbers(&complex);
println!("Betti: {:?}", betti.betti); // [1, ...]

// Classify genre
let genre = classify_genre(&betti);
println!("Genre: {:?}", genre);

// Build a filtration and compute persistence
let filt = build_filtration(&[c_major.clone(), a_minor.clone()]);
let barcodes = compute_persistence(&filt);
for bc in &barcodes {
    println!("H{}: {} bars", bc.dimension, bc.bars.len());
}

// Check counterpoint
let violations = check_counterpoint(&c_major, &a_minor);
println!("Violations: {:?}", violations);
```

## Topological Invariants

| Complex           | β₀ | β₁ | β₂ |
|-------------------|----|----|-----|
| Point             |  1 |  0 |  0 |
| Circle            |  1 |  1 |  0 |
| Sphere (S²)       |  1 |  0 |  1 |
| Filled triangle   |  1 |  0 |  0 |
| Hollow triangle   |  1 |  1 |  0 |
| Filled tetrahedron|  1 |  0 |  0 |
| Hollow tetrahedron|  1 |  0 |  1 |

All verified by the test suite.

## Architecture

```
Chord progression
       │
       ▼
  Simplicial Complex (pitch classes = vertices)
       │
       ▼
  Filtration (voice-leading distance)
       │
       ▼
  Boundary Matrix Reduction (Gaussian elimination, mod 2)
       │
       ▼
  Persistence Barcodes ──► Betti Numbers ──► Genre Classification
       │
       ▼
  Tension Profile (short bars = transient, long bars = persistent)
```

## Building & Testing

```bash
cargo build
cargo test
```

55 tests covering:
- Chord-as-simplex construction and normalisation
- Simplicial complex from chord progressions
- Betti numbers for canonical complexes (point, edge, triangle, sphere, tetrahedron)
- Persistent homology barcode extraction
- Genre Betti fingerprint distinctness and classification
- Cadence generation and H₁ reduction
- Voice-leading distance computation
- Filtration ordering
- Counterpoint violation detection (parallel fifths/octaves)
- Vietoris-Rips complex construction
- Full pipeline: chords → complex → filtration → barcodes

## Dependencies

- `serde` (with `derive`) — serialization of all public types

No external math libraries. Boundary matrix reduction implemented from scratch.

## License

MIT
