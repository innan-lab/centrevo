# Centrevo Implementation Guide: Clean Rewrite from Scratch

**Version:** 2.0 Clean Architecture  
**Date:** November 10, 2025  
**Purpose:** Complete rewrite guide with optimal data structures and performance patterns

---

## Table of Contents

1. [Core Principles](#core-principles)
2. [Project Structure](#project-structure)
3. [Data Types Foundation](#data-types-foundation)
4. [Module-by-Module Implementation](#module-by-module-implementation)
5. [Performance Patterns](#performance-patterns)
6. [Testing Strategy](#testing-strategy)
7. [Migration from V1](#migration-from-v1)

---

## Core Principles

### Design Philosophy
1. **Performance by Default** - Fast operations without tricks
2. **Zero-Cost Abstractions** - Type safety without runtime overhead
3. **Clear Ownership** - Explicit about who owns what
4. **Testable** - Every module independently testable
5. **Composable** - Build complex from simple pieces

### Data Structure Rules
```rust
// ✅ DO: Use these for their intended purposes
Vec<u8>         // Mutable sequences, active operations
Arc<[u8]>       // Immutable shared sequences, parallel reads
Arc<str>        // Shared identifiers and labels
Arc<[char]>     // Shared alphabets

// ⚠️ MAYBE: Use with care
Box<[u8]>       // Owned immutable sequences (if no sharing needed)
Cow<'a, [u8]>   // Function parameters that might mutate

// ❌ DON'T: Avoid these
Vec<char>       // Use Vec<u8> instead (4x memory)
String          // Use Arc<str> for shared data
ByteChunk       // Custom 2-bit encoding (too slow)
```

---

## Project Structure

### Directory Layout
```
centrevo/
├── Cargo.toml
├── README.md
├── IMPLEMENTATION_GUIDE.md
├── RECOMMENDATIONS.md
├── src/
│   ├── lib.rs                  # Public API
│   ├── prelude.rs              # Common imports
│   │
│   ├── base/                   # Foundation types
│   │   ├── mod.rs
│   │   ├── nucleotide.rs       # Base representation
│   │   ├── alphabet.rs         # Alphabet handling
│   │   └── sequence.rs         # Sequence types
│   │
│   ├── genome/                 # Genome structures
│   │   ├── mod.rs
│   │   ├── chromosome.rs       # Chromosome implementation
│   │   ├── haplotype.rs        # Haplotype collection
│   │   └── individual.rs       # Individual organism
│   │
│   ├── evolution/              # Evolutionary processes
│   │   ├── mod.rs
│   │   ├── mutation.rs         # Mutation operations
│   │   ├── recombination.rs    # Recombination logic
│   │   └── selection.rs        # Selection and fitness
│   │
│   ├── simulation/             # Simulation engine
│   │   ├── mod.rs
│   │   ├── engine.rs           # Main simulation loop
│   │   ├── population.rs       # Population management
│   │   └── parameters.rs       # Configuration
│   │
│   ├── utils/                  # Utilities
│   │   ├── mod.rs
│   │   ├── rng.rs              # Random number generation
│   │   ├── distributions.rs    # Probability distributions
│   │   └── pool.rs             # Object pooling
│   │
│   ├── storage/                # Persistence
│   │   ├── mod.rs
│   │   ├── database.rs         # SQLite recording
│   │   └── serialization.rs    # Import/export
│   │
│   ├── python/                 # Python bindings
│   │   ├── mod.rs
│   │   └── bindings.rs         # PyO3 wrappers
│   │
│   └── bin/
│       └── centrevo.rs         # CLI interface
│
├── benches/                    # Benchmarks
│   ├── sequence_ops.rs
│   ├── mutation.rs
│   └── simulation.rs
│
└── tests/                      # Integration tests
    ├── sequence_tests.rs
    ├── evolution_tests.rs
    └── simulation_tests.rs
```

---

## Data Types Foundation

### 1. Base Types (`src/base/nucleotide.rs`)

```rust
/// Represents a single nucleotide base.
/// Uses u8 internally (0=A, 1=C, 2=G, 3=T) for efficiency.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum Nucleotide {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
}

impl Nucleotide {
    /// Convert from u8 index (0-3)
    #[inline(always)]
    pub const fn from_index(idx: u8) -> Option<Self> {
        match idx {
            0 => Some(Self::A),
            1 => Some(Self::C),
            2 => Some(Self::G),
            3 => Some(Self::T),
            _ => None,
        }
    }

    /// Convert to u8 index (0-3)
    #[inline(always)]
    pub const fn to_index(self) -> u8 {
        self as u8
    }

    /// Convert from ASCII byte
    #[inline]
    pub const fn from_ascii(byte: u8) -> Option<Self> {
        match byte {
            b'A' | b'a' => Some(Self::A),
            b'C' | b'c' => Some(Self::C),
            b'G' | b'g' => Some(Self::G),
            b'T' | b't' => Some(Self::T),
            _ => None,
        }
    }

    /// Convert to ASCII byte (uppercase)
    #[inline(always)]
    pub const fn to_ascii(self) -> u8 {
        match self {
            Self::A => b'A',
            Self::C => b'C',
            Self::G => b'G',
            Self::T => b'T',
        }
    }

    /// Convert to char
    #[inline(always)]
    pub const fn to_char(self) -> char {
        self.to_ascii() as char
    }

    /// Get complement base
    #[inline(always)]
    pub const fn complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::T => Self::A,
            Self::C => Self::G,
            Self::G => Self::C,
        }
    }

    /// Check if purine (A or G)
    #[inline(always)]
    pub const fn is_purine(self) -> bool {
        matches!(self, Self::A | Self::G)
    }

    /// Check if pyrimidine (C or T)
    #[inline(always)]
    pub const fn is_pyrimidine(self) -> bool {
        matches!(self, Self::C | Self::T)
    }
}

impl TryFrom<u8> for Nucleotide {
    type Error = InvalidNucleotide;

    fn try_from(byte: u8) -> Result<Self, Self::Error> {
        Self::from_ascii(byte).ok_or(InvalidNucleotide(byte))
    }
}

impl From<Nucleotide> for u8 {
    #[inline(always)]
    fn from(nuc: Nucleotide) -> u8 {
        nuc.to_index()
    }
}

impl From<Nucleotide> for char {
    #[inline(always)]
    fn from(nuc: Nucleotide) -> char {
        nuc.to_char()
    }
}

#[derive(Debug, Clone, Copy)]
pub struct InvalidNucleotide(pub u8);

impl std::fmt::Display for InvalidNucleotide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Invalid nucleotide byte: {} ('{}')", self.0, self.0 as char)
    }
}

impl std::error::Error for InvalidNucleotide {}
```

### 2. Alphabet (`src/base/alphabet.rs`)

```rust
use std::sync::Arc;
use super::Nucleotide;

/// Shared, immutable alphabet.
/// Use Arc to share one instance across all chromosomes in a population.
#[derive(Debug, Clone)]
pub struct Alphabet {
    /// Character representation of bases
    chars: Arc<[char]>,
    /// Mapping from char to index for fast lookup
    char_to_index: Arc<std::collections::HashMap<char, u8>>,
}

impl Alphabet {
    /// Create a new alphabet from characters.
    /// The order determines the index mapping.
    pub fn new(chars: impl Into<Vec<char>>) -> Self {
        let chars: Vec<char> = chars.into();
        let char_to_index = chars
            .iter()
            .enumerate()
            .map(|(i, &c)| (c, i as u8))
            .collect();

        Self {
            chars: chars.into(),
            char_to_index: Arc::new(char_to_index),
        }
    }

    /// Standard DNA alphabet (A, C, G, T)
    pub fn dna() -> Self {
        Self::new(vec!['A', 'C', 'G', 'T'])
    }

    /// Get the number of bases in this alphabet
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.chars.len()
    }

    /// Check if empty (should never be)
    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.chars.is_empty()
    }

    /// Get character by index
    #[inline]
    pub fn get_char(&self, index: u8) -> Option<char> {
        self.chars.get(index as usize).copied()
    }

    /// Get index by character
    #[inline]
    pub fn get_index(&self, c: char) -> Option<u8> {
        self.char_to_index.get(&c).copied()
    }

    /// Get all characters as slice
    #[inline]
    pub fn chars(&self) -> &[char] {
        &self.chars
    }

    /// Check if character is in alphabet
    #[inline]
    pub fn contains(&self, c: char) -> bool {
        self.char_to_index.contains_key(&c)
    }
}

impl Default for Alphabet {
    fn default() -> Self {
        Self::dna()
    }
}

impl PartialEq for Alphabet {
    fn eq(&self, other: &Self) -> bool {
        // Fast path: check if they point to the same Arc
        Arc::ptr_eq(&self.chars, &other.chars)
            || self.chars == other.chars
    }
}

impl Eq for Alphabet {}
```

### 3. Sequence Types (`src/base/sequence.rs`)

```rust
use std::sync::Arc;
use super::{Nucleotide, Alphabet};

/// Mutable sequence - use for active operations
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sequence {
    /// Indices into alphabet (0-3 for standard DNA)
    data: Vec<u8>,
    /// Shared reference to alphabet
    alphabet: Alphabet,
}

impl Sequence {
    /// Create a new empty sequence
    pub fn new(alphabet: Alphabet) -> Self {
        Self {
            data: Vec::new(),
            alphabet,
        }
    }

    /// Create with capacity
    pub fn with_capacity(capacity: usize, alphabet: Alphabet) -> Self {
        Self {
            data: Vec::with_capacity(capacity),
            alphabet,
        }
    }

    /// Create from raw indices
    pub fn from_indices(indices: Vec<u8>, alphabet: Alphabet) -> Self {
        Self {
            data: indices,
            alphabet,
        }
    }

    /// Create from string
    pub fn from_str(s: &str, alphabet: Alphabet) -> Result<Self, InvalidSequence> {
        let data: Result<Vec<u8>, _> = s
            .chars()
            .map(|c| {
                alphabet
                    .get_index(c.to_ascii_uppercase())
                    .ok_or(InvalidSequence::InvalidChar(c))
            })
            .collect();

        Ok(Self {
            data: data?,
            alphabet,
        })
    }

    /// Get length
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Check if empty
    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Get base at position
    #[inline]
    pub fn get(&self, index: usize) -> Option<Nucleotide> {
        self.data.get(index).and_then(|&idx| Nucleotide::from_index(idx))
    }

    /// Set base at position
    #[inline]
    pub fn set(&mut self, index: usize, base: Nucleotide) -> Result<(), OutOfBounds> {
        self.data
            .get_mut(index)
            .map(|slot| *slot = base.to_index())
            .ok_or(OutOfBounds { index, len: self.len() })
    }

    /// Get raw indices
    #[inline]
    pub fn indices(&self) -> &[u8] {
        &self.data
    }

    /// Get mutable raw indices
    #[inline]
    pub fn indices_mut(&mut self) -> &mut [u8] {
        &mut self.data
    }

    /// Get alphabet
    #[inline]
    pub fn alphabet(&self) -> &Alphabet {
        &self.alphabet
    }

    /// Push a base
    #[inline]
    pub fn push(&mut self, base: Nucleotide) {
        self.data.push(base.to_index());
    }

    /// Pop a base
    #[inline]
    pub fn pop(&mut self) -> Option<Nucleotide> {
        self.data.pop().and_then(Nucleotide::from_index)
    }

    /// Insert a base
    #[inline]
    pub fn insert(&mut self, index: usize, base: Nucleotide) {
        self.data.insert(index, base.to_index());
    }

    /// Remove a base
    #[inline]
    pub fn remove(&mut self, index: usize) -> Nucleotide {
        Nucleotide::from_index(self.data.remove(index))
            .expect("Invalid base in sequence")
    }

    /// Convert to string
    pub fn to_string(&self) -> String {
        self.data
            .iter()
            .filter_map(|&idx| self.alphabet.get_char(idx))
            .collect()
    }

    /// Convert to immutable shared sequence
    pub fn into_shared(self) -> SharedSequence {
        SharedSequence {
            data: self.data.into(),
            alphabet: self.alphabet,
        }
    }

    /// Convert to immutable shared sequence (cloning data)
    pub fn to_shared(&self) -> SharedSequence {
        SharedSequence {
            data: self.data.clone().into(),
            alphabet: self.alphabet.clone(),
        }
    }
}

/// Immutable shared sequence - use for read-only operations
#[derive(Debug, Clone)]
pub struct SharedSequence {
    /// Shared immutable indices
    data: Arc<[u8]>,
    /// Shared reference to alphabet
    alphabet: Alphabet,
}

impl SharedSequence {
    /// Create from Arc
    pub fn new(data: Arc<[u8]>, alphabet: Alphabet) -> Self {
        Self { data, alphabet }
    }

    /// Get length
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Check if empty
    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Get base at position
    #[inline]
    pub fn get(&self, index: usize) -> Option<Nucleotide> {
        self.data.get(index).and_then(|&idx| Nucleotide::from_index(idx))
    }

    /// Get raw indices
    #[inline]
    pub fn indices(&self) -> &[u8] {
        &self.data
    }

    /// Get alphabet
    #[inline]
    pub fn alphabet(&self) -> &Alphabet {
        &self.alphabet
    }

    /// Convert to string
    pub fn to_string(&self) -> String {
        self.data
            .iter()
            .filter_map(|&idx| self.alphabet.get_char(idx))
            .collect()
    }

    /// Clone data into mutable sequence
    pub fn to_mutable(&self) -> Sequence {
        Sequence {
            data: self.data.to_vec(),
            alphabet: self.alphabet.clone(),
        }
    }

    /// Get strong reference count (for debugging)
    pub fn strong_count(&self) -> usize {
        Arc::strong_count(&self.data)
    }
}

impl PartialEq for SharedSequence {
    fn eq(&self, other: &Self) -> bool {
        // Fast path: check if same Arc
        Arc::ptr_eq(&self.data, &other.data)
            || (self.data == other.data && self.alphabet == other.alphabet)
    }
}

impl Eq for SharedSequence {}

// Errors
#[derive(Debug, Clone)]
pub enum InvalidSequence {
    InvalidChar(char),
    EmptySequence,
}

impl std::fmt::Display for InvalidSequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidChar(c) => write!(f, "Invalid character in sequence: '{}'", c),
            Self::EmptySequence => write!(f, "Empty sequence not allowed"),
        }
    }
}

impl std::error::Error for InvalidSequence {}

#[derive(Debug, Clone, Copy)]
pub struct OutOfBounds {
    pub index: usize,
    pub len: usize,
}

impl std::fmt::Display for OutOfBounds {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Index {} out of bounds (len = {})", self.index, self.len)
    }
}

impl std::error::Error for OutOfBounds {}
```

---

## Module-by-Module Implementation

### 4. Chromosome (`src/genome/chromosome.rs`)

```rust
use std::sync::Arc;
use crate::base::{Alphabet, Sequence, SharedSequence, Nucleotide};

/// A chromosome is a sequence organized into repeat units and higher-order repeats.
#[derive(Debug, Clone)]
pub struct Chromosome {
    /// Unique identifier (shared across clones)
    id: Arc<str>,
    /// The sequence data (mutable during simulation)
    sequence: Sequence,
    /// Repeat unit length in bases
    ru_length: usize,
    /// Number of repeat units per higher-order repeat
    rus_per_hor: usize,
}

impl Chromosome {
    /// Create a new chromosome
    pub fn new(
        id: impl Into<Arc<str>>,
        sequence: Sequence,
        ru_length: usize,
        rus_per_hor: usize,
    ) -> Self {
        Self {
            id: id.into(),
            sequence,
            ru_length,
            rus_per_hor,
        }
    }

    /// Create a uniform chromosome (all same base)
    pub fn uniform(
        id: impl Into<Arc<str>>,
        base: Nucleotide,
        length: usize,
        ru_length: usize,
        rus_per_hor: usize,
        alphabet: Alphabet,
    ) -> Self {
        let mut sequence = Sequence::with_capacity(length, alphabet);
        for _ in 0..length {
            sequence.push(base);
        }

        Self::new(id, sequence, ru_length, rus_per_hor)
    }

    /// Get chromosome ID
    #[inline]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Get sequence reference
    #[inline]
    pub fn sequence(&self) -> &Sequence {
        &self.sequence
    }

    /// Get mutable sequence reference
    #[inline]
    pub fn sequence_mut(&mut self) -> &mut Sequence {
        &mut self.sequence
    }

    /// Get length
    #[inline]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Check if empty
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Get repeat unit length
    #[inline]
    pub fn ru_length(&self) -> usize {
        self.ru_length
    }

    /// Get RUs per HOR
    #[inline]
    pub fn rus_per_hor(&self) -> usize {
        self.rus_per_hor
    }

    /// Get HOR length in bases
    #[inline]
    pub fn hor_length(&self) -> usize {
        self.ru_length * self.rus_per_hor
    }

    /// Get number of complete HORs
    #[inline]
    pub fn num_hors(&self) -> usize {
        self.len() / self.hor_length()
    }

    /// Get alphabet
    #[inline]
    pub fn alphabet(&self) -> &Alphabet {
        self.sequence.alphabet()
    }

    /// Calculate GC content
    pub fn gc_content(&self) -> f64 {
        let mut gc_count = 0;
        let mut total = 0;

        for &idx in self.sequence.indices() {
            if let Some(nuc) = Nucleotide::from_index(idx) {
                total += 1;
                if matches!(nuc, Nucleotide::G | Nucleotide::C) {
                    gc_count += 1;
                }
            }
        }

        if total == 0 {
            0.0
        } else {
            gc_count as f64 / total as f64
        }
    }

    /// Convert to string representation
    pub fn to_string(&self) -> String {
        self.sequence.to_string()
    }

    /// Convert to formatted string with RU and HOR delimiters
    pub fn to_formatted_string(&self, ru_delim: char, hor_delim: char) -> String {
        let chars: Vec<char> = self.sequence
            .indices()
            .iter()
            .filter_map(|&idx| self.alphabet().get_char(idx))
            .collect();

        let mut result = String::with_capacity(chars.len() * 2);
        let hor_len = self.hor_length();

        for (i, c) in chars.iter().enumerate() {
            if i > 0 {
                if i % hor_len == 0 {
                    result.push(hor_delim);
                } else if i % self.ru_length == 0 {
                    result.push(ru_delim);
                }
            }
            result.push(*c);
        }

        result
    }

    /// Create a shared immutable view
    pub fn to_shared(&self) -> SharedChromosome {
        SharedChromosome {
            id: self.id.clone(),
            sequence: self.sequence.to_shared(),
            ru_length: self.ru_length,
            rus_per_hor: self.rus_per_hor,
        }
    }
}

/// Immutable shared chromosome - use for parallel operations
#[derive(Debug, Clone)]
pub struct SharedChromosome {
    id: Arc<str>,
    sequence: SharedSequence,
    ru_length: usize,
    rus_per_hor: usize,
}

impl SharedChromosome {
    /// Get chromosome ID
    #[inline]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Get sequence reference
    #[inline]
    pub fn sequence(&self) -> &SharedSequence {
        &self.sequence
    }

    /// Get length
    #[inline]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Check if empty
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Get repeat unit length
    #[inline]
    pub fn ru_length(&self) -> usize {
        self.ru_length
    }

    /// Get RUs per HOR
    #[inline]
    pub fn rus_per_hor(&self) -> usize {
        self.rus_per_hor
    }

    /// Calculate GC content
    pub fn gc_content(&self) -> f64 {
        let mut gc_count = 0;
        let mut total = 0;

        for &idx in self.sequence.indices() {
            if let Some(nuc) = Nucleotide::from_index(idx) {
                total += 1;
                if matches!(nuc, Nucleotide::G | Nucleotide::C) {
                    gc_count += 1;
                }
            }
        }

        if total == 0 {
            0.0
        } else {
            gc_count as f64 / total as f64
        }
    }

    /// Clone data into mutable chromosome
    pub fn to_mutable(&self) -> Chromosome {
        Chromosome {
            id: self.id.clone(),
            sequence: self.sequence.to_mutable(),
            ru_length: self.ru_length,
            rus_per_hor: self.rus_per_hor,
        }
    }
}
```

### 5. Haplotype (`src/genome/haplotype.rs`)

```rust
use crate::genome::Chromosome;

/// A haplotype is a collection of chromosomes.
#[derive(Debug, Clone)]
pub struct Haplotype {
    /// Chromosomes in this haplotype
    chromosomes: Vec<Chromosome>,
}

impl Haplotype {
    /// Create a new empty haplotype
    pub fn new() -> Self {
        Self {
            chromosomes: Vec::new(),
        }
    }

    /// Create with specific capacity
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            chromosomes: Vec::with_capacity(capacity),
        }
    }

    /// Create from chromosomes
    pub fn from_chromosomes(chromosomes: Vec<Chromosome>) -> Self {
        Self { chromosomes }
    }

    /// Get number of chromosomes
    #[inline]
    pub fn len(&self) -> usize {
        self.chromosomes.len()
    }

    /// Check if empty
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.chromosomes.is_empty()
    }

    /// Get chromosome by index
    #[inline]
    pub fn get(&self, index: usize) -> Option<&Chromosome> {
        self.chromosomes.get(index)
    }

    /// Get mutable chromosome by index
    #[inline]
    pub fn get_mut(&mut self, index: usize) -> Option<&mut Chromosome> {
        self.chromosomes.get_mut(index)
    }

    /// Add a chromosome
    pub fn push(&mut self, chromosome: Chromosome) {
        self.chromosomes.push(chromosome);
    }

    /// Get all chromosomes
    #[inline]
    pub fn chromosomes(&self) -> &[Chromosome] {
        &self.chromosomes
    }

    /// Get mutable chromosomes
    #[inline]
    pub fn chromosomes_mut(&mut self) -> &mut [Chromosome] {
        &mut self.chromosomes
    }

    /// Iterate over chromosomes
    pub fn iter(&self) -> impl Iterator<Item = &Chromosome> {
        self.chromosomes.iter()
    }

    /// Iterate mutably over chromosomes
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut Chromosome> {
        self.chromosomes.iter_mut()
    }

    /// Total length across all chromosomes
    pub fn total_length(&self) -> usize {
        self.chromosomes.iter().map(|chr| chr.len()).sum()
    }
}

impl Default for Haplotype {
    fn default() -> Self {
        Self::new()
    }
}

impl IntoIterator for Haplotype {
    type Item = Chromosome;
    type IntoIter = std::vec::IntoIter<Chromosome>;

    fn into_iter(self) -> Self::IntoIter {
        self.chromosomes.into_iter()
    }
}
```

### 6. Individual (`src/genome/individual.rs`)

```rust
use std::sync::Arc;
use crate::genome::Haplotype;

/// An individual organism with diploid genome.
#[derive(Debug, Clone)]
pub struct Individual {
    /// Unique identifier
    id: Arc<str>,
    /// First haplotype
    haplotype1: Haplotype,
    /// Second haplotype
    haplotype2: Haplotype,
    /// Fitness value (cached)
    fitness: f64,
}

impl Individual {
    /// Create a new individual
    pub fn new(
        id: impl Into<Arc<str>>,
        haplotype1: Haplotype,
        haplotype2: Haplotype,
    ) -> Self {
        Self {
            id: id.into(),
            haplotype1,
            haplotype2,
            fitness: 0.0,
        }
    }

    /// Get individual ID
    #[inline]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Get first haplotype
    #[inline]
    pub fn haplotype1(&self) -> &Haplotype {
        &self.haplotype1
    }

    /// Get mutable first haplotype
    #[inline]
    pub fn haplotype1_mut(&mut self) -> &mut Haplotype {
        &mut self.haplotype1
    }

    /// Get second haplotype
    #[inline]
    pub fn haplotype2(&self) -> &Haplotype {
        &self.haplotype2
    }

    /// Get mutable second haplotype
    #[inline]
    pub fn haplotype2_mut(&mut self) -> &mut Haplotype {
        &mut self.haplotype2
    }

    /// Get fitness
    #[inline]
    pub fn fitness(&self) -> f64 {
        self.fitness
    }

    /// Set fitness
    #[inline]
    pub fn set_fitness(&mut self, fitness: f64) {
        self.fitness = fitness;
    }

    /// Get both haplotypes
    pub fn haplotypes(&self) -> (&Haplotype, &Haplotype) {
        (&self.haplotype1, &self.haplotype2)
    }

    /// Get both haplotypes mutably
    pub fn haplotypes_mut(&mut self) -> (&mut Haplotype, &mut Haplotype) {
        (&mut self.haplotype1, &mut self.haplotype2)
    }
}
```

---

## Performance Patterns

### Pattern 1: Zero-Copy Sharing

```rust
// ✅ DO: Share alphabet across entire population
let alphabet = Alphabet::dna();  // Created once

let mut chromosomes = Vec::new();
for i in 0..1000 {
    // All chromosomes share the SAME alphabet Arc
    let chr = Chromosome::uniform(
        format!("chr{}", i),
        Nucleotide::A,
        1000,
        20,
        5,
        alphabet.clone(),  // Just increments ref count!
    );
    chromosomes.push(chr);
}

// ❌ DON'T: Create new alphabet for each chromosome
for i in 0..1000 {
    let alphabet = Alphabet::dna();  // BAD: Creates 1000 copies!
    let chr = Chromosome::uniform(..., alphabet);
}
```

### Pattern 2: Parallel Operations with SharedSequence

```rust
use rayon::prelude::*;

// ✅ DO: Convert to shared before parallel operations
let chromosomes: Vec<Chromosome> = /* ... */;
let shared: Vec<SharedChromosome> = chromosomes
    .iter()
    .map(|chr| chr.to_shared())
    .collect();

// Parallel fitness calculation (no cloning needed!)
let fitnesses: Vec<f64> = shared
    .par_iter()
    .map(|chr| calculate_fitness(chr))
    .collect();

// ❌ DON'T: Clone mutable chromosomes in parallel
let fitnesses: Vec<f64> = chromosomes
    .par_iter()
    .map(|chr| {
        let cloned = chr.clone();  // BAD: Clones entire sequence!
        calculate_fitness(&cloned)
    })
    .collect();
```

### Pattern 3: Mutation with Object Pooling

```rust
use crate::utils::pool::SequencePool;

pub struct MutationEngine {
    pool: SequencePool,
}

impl MutationEngine {
    pub fn mutate(&mut self, chromosome: &mut Chromosome) {
        // Get buffer from pool (reused allocation)
        let mut buffer = self.pool.get(chromosome.len());
        
        // Do mutation work using buffer
        buffer.extend_from_slice(chromosome.sequence().indices());
        // ... mutation logic ...
        
        // Update chromosome
        *chromosome.sequence_mut() = Sequence::from_indices(
            buffer.clone(),
            chromosome.alphabet().clone(),
        );
        
        // Return buffer to pool for reuse
        self.pool.recycle(buffer);
    }
}
```

### Pattern 4: Efficient String Building

```rust
// ✅ DO: Pre-allocate and build once
pub fn format_population(individuals: &[Individual]) -> String {
    let mut result = String::with_capacity(individuals.len() * 100);
    for ind in individuals {
        use std::fmt::Write;
        writeln!(result, ">{}\\n{}", ind.id(), ind.haplotype1().chromosomes()[0].to_string()).unwrap();
    }
    result
}

// ❌ DON'T: Concatenate repeatedly
pub fn format_population_bad(individuals: &[Individual]) -> String {
    let mut result = String::new();
    for ind in individuals {
        result = result + ">" + ind.id() + "\\n";  // BAD: Creates new String each time!
    }
    result
}
```

### Pattern 5: SIMD Operations (Advanced)

```rust
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

/// Count nucleotides using SIMD
pub fn count_bases_simd(sequence: &[u8]) -> [usize; 4] {
    let mut counts = [0usize; 4];
    
    // Process 16 bytes at a time with SSE
    let chunks = sequence.chunks_exact(16);
    let remainder = chunks.remainder();
    
    for chunk in chunks {
        // SIMD operations here
        // ... implementation ...
    }
    
    // Handle remainder
    for &base in remainder {
        if base < 4 {
            counts[base as usize] += 1;
        }
    }
    
    counts
}

/// Fallback for non-x86 platforms
#[cfg(not(target_arch = "x86_64"))]
pub fn count_bases_simd(sequence: &[u8]) -> [usize; 4] {
    let mut counts = [0usize; 4];
    for &base in sequence {
        if base < 4 {
            counts[base as usize] += 1;
        }
    }
    counts
}
```

---

## Testing Strategy

### Unit Tests Structure

```rust
// In src/base/sequence.rs
#[cfg(test)]
mod tests {
    use super::*;

    fn test_alphabet() -> Alphabet {
        Alphabet::dna()
    }

    #[test]
    fn sequence_creation() {
        let seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.get(0), Some(Nucleotide::A));
    }

    #[test]
    fn sequence_mutation() {
        let mut seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        seq.set(1, Nucleotide::T).unwrap();
        assert_eq!(seq.to_string(), "ATGT");
    }

    #[test]
    fn shared_sequence_immutable() {
        let seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let shared = seq.to_shared();
        assert_eq!(shared.len(), 4);
        // Cannot call shared.set() - compile error!
    }

    #[test]
    fn arc_sharing() {
        let seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let shared1 = seq.to_shared();
        let shared2 = shared1.clone();
        
        // Both point to same data
        assert_eq!(Arc::strong_count(&shared1.data), 2);
    }
}
```

### Benchmark Template

```rust
// In benches/sequence_ops.rs
use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use centrevo::base::{Alphabet, Sequence, Nucleotide};

fn bench_sequence_creation(c: &mut Criterion) {
    let alphabet = Alphabet::dna();
    let sizes = [100, 1000, 10000];

    let mut group = c.benchmark_group("sequence_creation");
    for size in sizes {
        group.bench_with_input(
            BenchmarkId::new("from_str", size),
            &size,
            |b, &size| {
                let s: String = "ACGT".repeat(size / 4);
                b.iter(|| {
                    Sequence::from_str(&s, alphabet.clone())
                })
            },
        );
    }
    group.finish();
}

fn bench_clone_vs_share(c: &mut Criterion) {
    let alphabet = Alphabet::dna();
    let seq = Sequence::from_str(&"ACGT".repeat(1000), alphabet).unwrap();
    let shared = seq.to_shared();

    let mut group = c.benchmark_group("clone_vs_share");
    
    group.bench_function("clone_mutable", |b| {
        b.iter(|| black_box(seq.clone()))
    });

    group.bench_function("clone_shared", |b| {
        b.iter(|| black_box(shared.clone()))
    });

    group.finish();
}

criterion_group!(benches, bench_sequence_creation, bench_clone_vs_share);
criterion_main!(benches);
```

---

## Migration from V1

### Step-by-Step Migration

#### Phase 1: New Foundation (Week 1)
1. Create new directory structure
2. Implement base types (Nucleotide, Alphabet, Sequence)
3. Write comprehensive unit tests
4. Benchmark against V1

#### Phase 2: Genome Structures (Week 2)
1. Implement Chromosome with new types
2. Implement Haplotype
3. Implement Individual
4. Port existing tests

#### Phase 3: Evolution Logic (Week 3)
1. Implement mutation module
2. Implement recombination module
3. Implement fitness/selection module
4. Verify against V1 results

#### Phase 4: Simulation Engine (Week 4)
1. Implement simulation loop
2. Implement population management
3. Add parallelization with Rayon
4. Performance tuning

#### Phase 5: Integration (Week 5)
1. Implement CLI interface
2. Implement Python bindings
3. Add database recording
4. Documentation and examples

### Compatibility Layer (Optional)

```rust
// src/compat/v1.rs - Compatibility with V1 API

use crate::base::{Sequence, Alphabet};
use crate::genome::Chromosome;

/// Wrapper to maintain V1 API compatibility
pub struct V1Chromosome {
    inner: Chromosome,
}

impl V1Chromosome {
    pub fn new(chr_id: String, sequence: Vec<char>, /* ... */) -> Result<Self, String> {
        let alphabet = Alphabet::dna();
        let seq_str: String = sequence.into_iter().collect();
        let seq = Sequence::from_str(&seq_str, alphabet)
            .map_err(|e| format!("{}", e))?;
        
        Ok(Self {
            inner: Chromosome::new(chr_id, seq, /* ... */),
        })
    }

    // ... other V1 methods ...
}
```

---

## Quick Reference

### Memory Ownership Rules

| Type | Ownership | Mutability | Cost to Clone | Use Case |
|------|-----------|------------|---------------|----------|
| `Vec<u8>` | Owned | Mutable | O(n) | Active sequences |
| `Arc<[u8]>` | Shared | Immutable | O(1) | Parallel reads |
| `Arc<str>` | Shared | Immutable | O(1) | IDs, labels |
| `Arc<[char]>` | Shared | Immutable | O(1) | Shared alphabets |
| `&[u8]` | Borrowed | Immutable | N/A | Function params |
| `&mut [u8]` | Borrowed | Mutable | N/A | Function params |

### Common Operations Cheat Sheet

```rust
// Create sequence
let seq = Sequence::from_str("ACGT", alphabet)?;

// Clone cheaply
let shared = seq.to_shared();
let clone = shared.clone();  // Just increments ref count

// Mutate
let mut mut_seq = shared.to_mutable();
mut_seq.set(0, Nucleotide::T)?;

// Share alphabet
let alphabet = Alphabet::dna();
let seq1 = Sequence::from_str("ACG", alphabet.clone());
let seq2 = Sequence::from_str("TGC", alphabet.clone());

// Parallel iteration
use rayon::prelude::*;
chromosomes.par_iter()
    .map(|chr| chr.to_shared())
    .map(|shared| calculate_fitness(&shared))
    .collect()
```

---

## Performance Checklist

Before committing code, verify:

- [ ] Alphabets are shared via `Arc<[char]>`, not cloned
- [ ] IDs use `Arc<str>`, not `String`
- [ ] Parallel operations use `SharedSequence`, not cloned `Sequence`
- [ ] No unnecessary `to_string()` calls in hot loops
- [ ] String building uses pre-allocated capacity
- [ ] Temporary sequences use object pooling when possible
- [ ] Benchmark shows improvement over V1
- [ ] No allocations in inner loops (use `cargo flamegraph`)

---

## Next Steps

1. **Start Fresh**: Create new branch `v2-rewrite`
2. **Build Foundation**: Implement base types first
3. **Test Everything**: Write tests before moving to next module
4. **Benchmark Early**: Compare each module to V1
5. **Iterate**: Refine based on benchmarks

**Remember:**
- Simple is fast
- Share, don't clone
- Measure, don't guess
- Test, then optimize

Good luck with the rewrite! 🚀
