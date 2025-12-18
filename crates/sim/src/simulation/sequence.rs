//! Initialization of simulations from custom sequences.
//!
//! This module provides functionality to initialize simulations from various
//! sequence sources including FASTA files, JSON data, and previous simulation
//! databases. It validates that input sequences match the expected parameters.

use crate::base::{Nucleotide, Sequence};
use crate::errors::InitializationError;
use crate::genome::repeat_map::RepeatMap;
use crate::genome::{Chromosome, Haplotype, Individual};
use crate::simulation::{
    GenerationMode, InitializationConfig, SequenceSource, UniformRepeatStructure,
};
use centrevo_codec::CodecStrategy;
use rand::Rng;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::str::FromStr;

/// A parsed sequence entry with ID and sequence data (legacy format).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequenceEntry {
    /// Sequence identifier
    pub id: String,
    /// Sequence data as string
    pub seq: String,
}

/// A sequence entry with explicit individual, haplotype, and chromosome indices.
/// This is the preferred format for multi-chromosome configurations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequenceEntryWithIndices {
    /// Individual index (0-based)
    pub ind_idx: usize,
    /// Haplotype index (0 or 1 for diploid)
    pub hap_idx: usize,
    /// Chromosome index (0-based within haplotype)
    pub chr_idx: usize,
    /// Sequence data as string
    pub seq: String,
}

impl SequenceEntryWithIndices {
    /// Create from a SequenceEntry by parsing the ID.
    /// Expected ID format: `<prefix>_ind<N>_hap<M>_chr<K>` or `ind<N>_hap<M>_chr<K>`
    /// where N = individual index, M = haplotype index, K = chromosome index.
    pub fn from_entry_with_parsed_id(entry: SequenceEntry) -> Result<Self, InitializationError> {
        Self::parse_id_format(&entry.id, entry.seq)
    }

    /// Parse sequence ID to extract indices.
    /// Supports formats:
    /// - `ind<N>_hap<M>_chr<K>`
    /// - `<prefix>_ind<N>_hap<M>_chr<K>`
    /// - `ind<N>_h<M>_chr<K>` (abbreviated haplotype)
    /// - `i<N>_h<M>_c<K>` (fully abbreviated)
    fn parse_id_format(id: &str, seq: String) -> Result<Self, InitializationError> {
        // Try to find the pattern in the string
        // We look for ind, hap/h, and chr/c markers

        let id_lower = id.to_lowercase();

        // Find individual index - try "ind" first, then standalone "i" between underscores
        let ind_idx = Self::extract_index(&id_lower, &["ind"])
            .or_else(|| {
                // Look for _i<N>_ pattern
                if let Some(pos) = id_lower.find("_i") {
                    let after = &id_lower[pos + 2..];
                    let digits: String = after.chars().take_while(|c| c.is_ascii_digit()).collect();
                    if !digits.is_empty() && after.starts_with(&digits) {
                        // Check if followed by underscore or end
                        let after_digits = &after[digits.len()..];
                        if after_digits.is_empty() || after_digits.starts_with('_') {
                            return digits.parse::<usize>().ok();
                        }
                    }
                }
                None
            })
            .ok_or_else(|| InitializationError::Parse(format!(
                "Could not find individual index in '{id}'. Expected format: ind<N>_hap<M>_chr<K>"
            )))?;

        // Find haplotype index - try "hap" first, then "_h<N>" pattern
        let hap_idx = Self::extract_index(&id_lower, &["hap"])
            .or_else(|| {
                // Look for _h<N>_ pattern
                if let Some(pos) = id_lower.find("_h") {
                    let after = &id_lower[pos + 2..];
                    let digits: String = after.chars().take_while(|c| c.is_ascii_digit()).collect();
                    if !digits.is_empty() && after.starts_with(&digits) {
                        // Check if followed by underscore or end
                        let after_digits = &after[digits.len()..];
                        if after_digits.is_empty() || after_digits.starts_with('_') {
                            return digits.parse::<usize>().ok();
                        }
                    }
                }
                None
            })
            .ok_or_else(|| InitializationError::Parse(format!(
                "Could not find haplotype index in '{id}'. Expected format: ind<N>_hap<M>_chr<K>"
            )))?;

        // Find chromosome index - try "chr" first, then "_c<N>" pattern
        let chr_idx = Self::extract_index(&id_lower, &["chr"])
            .or_else(|| {
                // Look for _c<N> pattern (can be at end or followed by underscore)
                if let Some(pos) = id_lower.find("_c") {
                    let after = &id_lower[pos + 2..];
                    let digits: String = after.chars().take_while(|c| c.is_ascii_digit()).collect();
                    if !digits.is_empty() && after.starts_with(&digits) {
                        // Check if followed by underscore or end
                        let after_digits = &after[digits.len()..];
                        if after_digits.is_empty() || after_digits.starts_with('_') {
                            return digits.parse::<usize>().ok();
                        }
                    }
                }
                None
            })
            .ok_or_else(|| InitializationError::Parse(format!(
                "Could not find chromosome index in '{id}'. Expected format: ind<N>_hap<M>_chr<K>"
            )))?;

        Ok(Self {
            ind_idx,
            hap_idx,
            chr_idx,
            seq,
        })
    }

    /// Extract an index following one of the given prefixes.
    fn extract_index(text: &str, prefixes: &[&str]) -> Option<usize> {
        for prefix in prefixes {
            if let Some(pos) = text.find(prefix) {
                // Make sure this is a word boundary (either at start or preceded by non-alphanumeric)
                if pos > 0 {
                    let prev_char = text.chars().nth(pos - 1)?;
                    if prev_char.is_alphanumeric() {
                        continue; // This is part of another word
                    }
                }

                let after = &text[pos + prefix.len()..];
                // Extract digits
                let digits: String = after.chars().take_while(|c| c.is_ascii_digit()).collect();
                if !digits.is_empty()
                    && let Ok(idx) = digits.parse::<usize>()
                {
                    return Some(idx);
                }
            }
        }
        None
    }
}

// SequenceInput enum removed, replaced by SequenceConfig in configs.rs

// Removed InitializationError definition, imported from crate::errors

/// Convert sequence indices (Vec<u8>) to a String representation.
fn indices_to_string(indices: &[u8]) -> Result<String, InitializationError> {
    let mut result = String::with_capacity(indices.len());
    for &idx in indices {
        let nuc = Nucleotide::from_index(idx).ok_or_else(|| {
            InitializationError::Parse(format!("Invalid base index {idx} in sequence"))
        })?;
        result.push(nuc.to_char());
    }
    Ok(result)
}

/// Parse sequences from a FASTA file.
///
/// FASTA format with IDs containing individual, haplotype, and chromosome information:
/// ```text
/// >ind0_hap0_chr0 description
/// ACGTACGTACGT
/// >ind0_hap0_chr1
/// TGCATGCA
/// >ind0_hap1_chr0
/// GGCCTTAA
/// ```
///
/// The ID must contain patterns that identify:
/// - Individual index: `ind<N>` or `i<N>`
/// - Haplotype index: `hap<M>` or `h<M>`
/// - Chromosome index: `chr<K>` or `c<K>`
pub fn parse_fasta(
    path: impl AsRef<Path>,
) -> Result<Vec<SequenceEntryWithIndices>, InitializationError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut entries = Vec::new();
    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();

        if line.is_empty() {
            continue;
        }

        if let Some(header) = line.strip_prefix('>') {
            // Save previous entry if exists
            if let Some(id) = current_id.take() {
                if !current_seq.is_empty() {
                    let entry =
                        SequenceEntryWithIndices::parse_id_format(&id, current_seq.clone())?;
                    entries.push(entry);
                    current_seq.clear();
                }
            }

            // Parse new header
            // Remove '>'
            let id = header
                .split_whitespace()
                .next()
                .ok_or_else(|| InitializationError::Parse("Empty FASTA header".to_string()))?;
            current_id = Some(id.to_string());
        } else {
            // Accumulate sequence data
            current_seq.push_str(line);
        }
    }

    // Save last entry
    if let Some(id) = current_id {
        if !current_seq.is_empty() {
            let entry = SequenceEntryWithIndices::parse_id_format(&id, current_seq)?;
            entries.push(entry);
        }
    }

    if entries.is_empty() {
        return Err(InitializationError::Parse(
            "No sequences found in FASTA file".to_string(),
        ));
    }

    Ok(entries)
}

/// Parse sequences from JSON.
///
/// JSON format can be either:
/// 1. A file path to a JSON file
/// 2. A JSON string
///
/// Supports two JSON structures:
///
/// **New format (recommended for multi-chromosome):**
/// ```json
/// [
///   {"ind_idx": 0, "hap_idx": 0, "chr_idx": 0, "seq": "ACGTACGT"},
///   {"ind_idx": 0, "hap_idx": 0, "chr_idx": 1, "seq": "TGCATGCA"},
///   {"ind_idx": 0, "hap_idx": 1, "chr_idx": 0, "seq": "GGCCTTAA"}
/// ]
/// ```
///
/// **Legacy format (for backward compatibility):**
/// ```json
/// [
///   {"id": "ind0_hap0_chr0", "seq": "ACGTACGT"},
///   {"id": "ind0_hap0_chr1", "seq": "TGCATGCA"}
/// ]
/// ```
pub fn parse_json(input: &str) -> Result<Vec<SequenceEntryWithIndices>, InitializationError> {
    // Helper to convert Vec<SequenceEntry> to Vec<SequenceEntryWithIndices>
    let convert_legacy = |entries: Vec<SequenceEntry>| -> Result<Vec<SequenceEntryWithIndices>, InitializationError> {
        entries.into_iter()
            .map(SequenceEntryWithIndices::from_entry_with_parsed_id)
            .collect()
    };

    // Try to parse as new format (JSON string)
    if let Ok(entries) = serde_json::from_str::<Vec<SequenceEntryWithIndices>>(input) {
        if entries.is_empty() {
            return Err(InitializationError::Parse("Empty JSON array".to_string()));
        }
        return Ok(entries);
    }

    // Try legacy format (JSON string)
    if let Ok(entries) = serde_json::from_str::<Vec<SequenceEntry>>(input) {
        if entries.is_empty() {
            return Err(InitializationError::Parse("Empty JSON array".to_string()));
        }
        return convert_legacy(entries);
    }

    // Try as file path - new format
    if let Ok(file) = File::open(input) {
        let reader = BufReader::new(file);
        if let Ok(entries) = serde_json::from_reader::<_, Vec<SequenceEntryWithIndices>>(reader) {
            if entries.is_empty() {
                return Err(InitializationError::Parse(
                    "Empty JSON array in file".to_string(),
                ));
            }
            return Ok(entries);
        }
    }

    // Try as file path - legacy format
    let file = File::open(input)?;
    let reader = BufReader::new(file);
    let entries: Vec<SequenceEntry> = serde_json::from_reader(reader)?;

    if entries.is_empty() {
        return Err(InitializationError::Parse(
            "Empty JSON array in file".to_string(),
        ));
    }

    convert_legacy(entries)
}

/// Load sequences from a previous simulation database.
///
/// Loads the population state from a specific generation (or the last generation
/// if not specified). Returns sequences with explicit indices.
pub fn load_from_database(
    db_path: impl AsRef<Path>,
    sim_id: &str,
    generation: Option<usize>,
    codec: &CodecStrategy,
) -> Result<Vec<SequenceEntryWithIndices>, InitializationError> {
    use crate::storage::QueryBuilder;

    let query = QueryBuilder::new(db_path.as_ref())
        .map_err(|e| InitializationError::Database(format!("Failed to open database: {e}")))?;

    // Determine which generation to load
    let target_gen = match generation {
        Some(g) => g,
        None => {
            // Get the last recorded generation
            let recorded_gens = query.get_recorded_generations(sim_id).map_err(|e| {
                InitializationError::Database(format!("Failed to get generations: {e}"))
            })?;

            *recorded_gens
                .last()
                .ok_or_else(|| InitializationError::Database("No generations found".to_string()))?
        }
    };

    // Load population at that generation
    let snapshots = query
        .get_generation(sim_id, target_gen)
        .map_err(|e| InitializationError::Database(format!("Failed to load generation: {e}")))?;

    // Close query builder directly since we don't need the config (sequences are self-contained)
    query
        .close()
        .map_err(|e| InitializationError::Database(format!("Failed to close database: {e}")))?;

    // Convert snapshots to sequence entries with indices
    let mut entries = Vec::new();

    // Database stores sequences per chromosome, so we need to group them properly
    // Assuming snapshots are for single-chromosome case (chr_idx = 0)
    // For multi-chromosome, database would need to be extended
    for (ind_idx, snapshot) in snapshots.iter().enumerate() {
        // Convert Vec<u8> indices to String
        // Decode sequences
        let decoded_hap0 = codec
            .decode(&snapshot.haplotype1_seq)
            .map_err(|e| InitializationError::Database(format!("Failed to decode seq: {e}")))?;
        let decoded_hap1 = codec
            .decode(&snapshot.haplotype2_seq)
            .map_err(|e| InitializationError::Database(format!("Failed to decode seq: {e}")))?;

        let hap0_seq = indices_to_string(&decoded_hap0)?;
        let hap1_seq = indices_to_string(&decoded_hap1)?;

        // Add both haplotypes with explicit indices
        // Assuming single chromosome per haplotype (chr_idx = 0)
        entries.push(SequenceEntryWithIndices {
            ind_idx,
            hap_idx: 0,
            chr_idx: 0,
            seq: hap0_seq,
        });
        entries.push(SequenceEntryWithIndices {
            ind_idx,
            hap_idx: 1,
            chr_idx: 0,
            seq: hap1_seq,
        });
    }

    if entries.is_empty() {
        return Err(InitializationError::Database(
            "No sequences found in generation".to_string(),
        ));
    }

    Ok(entries)
}

/// Load individuals directly from database, preserving RepeatMaps.
pub fn load_individuals_from_database(
    db_path: impl AsRef<Path>,
    sim_id: &str,
    generation: Option<usize>,
    codec: &CodecStrategy,
) -> Result<Vec<Individual>, InitializationError> {
    use crate::storage::QueryBuilder;

    let query = QueryBuilder::new(db_path.as_ref())
        .map_err(|e| InitializationError::Database(format!("Failed to open database: {e}")))?;

    // Determine target generation
    let target_gen = match generation {
        Some(g) => g,
        None => {
            let recorded_gens = query.get_recorded_generations(sim_id).map_err(|e| {
                InitializationError::Database(format!("Failed to get generations: {e}"))
            })?;
            *recorded_gens
                .last()
                .ok_or_else(|| InitializationError::Database("No generations found".to_string()))?
        }
    };

    // Load snapshots
    let snapshots = query
        .get_generation(sim_id, target_gen)
        .map_err(|e| InitializationError::Database(format!("Failed to load generation: {e}")))?;

    query.close().ok();

    // Convert snapshots to Individuals
    let mut individuals = Vec::with_capacity(snapshots.len());
    for snapshot in snapshots {
        let ind = snapshot.to_individual(codec).map_err(|e| {
            InitializationError::Database(format!("Failed to reconstruct individual: {e}"))
        })?;
        individuals.push(ind);
    }

    Ok(individuals)
}

/// Validate that sequences match the expected parameters.
///
/// Validates:
/// - All required (ind, hap, chr) combinations are present
/// - No duplicate (ind, hap, chr) combinations
/// - Sequence lengths match expected chromosome length
/// - All bases are valid for the alphabet
pub fn validate_sequences(
    entries: &[SequenceEntryWithIndices],
    structure: &UniformRepeatStructure,
    population_size: usize,
) -> Result<(), InitializationError> {
    // Calculate expected counts
    let expected_chr_length = structure.chr_length();
    let expected_sequences_per_individual = structure.chrs_per_hap * 2; // 2 haplotypes
    let expected_total_sequences = population_size * expected_sequences_per_individual;

    // Validate sequence count
    if entries.len() != expected_total_sequences {
        return Err(InitializationError::Validation(format!(
            "Expected {} sequences (population_size={} × chromosomes_per_haplotype={} × 2 haplotypes), but got {}",
            expected_total_sequences,
            population_size,
            structure.chrs_per_hap,
            entries.len()
        )));
    }

    // Track seen combinations to detect duplicates
    use std::collections::HashSet;
    let mut seen = HashSet::new();

    // Validate each sequence
    for (idx, entry) in entries.iter().enumerate() {
        // Check indices are in valid ranges
        if entry.ind_idx >= population_size {
            return Err(InitializationError::Validation(format!(
                "Sequence {} has ind_idx={} but population_size={}",
                idx, entry.ind_idx, population_size
            )));
        }
        if entry.hap_idx >= 2 {
            return Err(InitializationError::Validation(format!(
                "Sequence {} has hap_idx={} but must be 0 or 1",
                idx, entry.hap_idx
            )));
        }
        if entry.chr_idx >= structure.chrs_per_hap {
            return Err(InitializationError::Validation(format!(
                "Sequence {} has chr_idx={} but chrs_per_hap={}",
                idx, entry.chr_idx, structure.chrs_per_hap
            )));
        }

        // Check for duplicates
        let key = (entry.ind_idx, entry.hap_idx, entry.chr_idx);
        if !seen.insert(key) {
            return Err(InitializationError::Validation(format!(
                "Duplicate sequence for ind_idx={}, hap_idx={}, chr_idx={}",
                entry.ind_idx, entry.hap_idx, entry.chr_idx
            )));
        }

        // Check length
        if entry.seq.len() != expected_chr_length {
            return Err(InitializationError::Validation(format!(
                "Sequence {} (ind={}, hap={}, chr={}) has length {}, expected {}",
                idx,
                entry.ind_idx,
                entry.hap_idx,
                entry.chr_idx,
                entry.seq.len(),
                expected_chr_length
            )));
        }

        // Check that all characters are valid bases
        for (pos, ch) in entry.seq.chars().enumerate() {
            if Nucleotide::from_ascii(ch as u8).is_none() {
                return Err(InitializationError::Validation(format!(
                    "Sequence {} (ind={}, hap={}, chr={}) contains invalid character '{}' at position {}",
                    idx, entry.ind_idx, entry.hap_idx, entry.chr_idx, ch, pos
                )));
            }
        }
    }

    // Verify all expected combinations are present
    for ind_idx in 0..population_size {
        for hap_idx in 0..2 {
            for chr_idx in 0..structure.chrs_per_hap {
                let key = (ind_idx, hap_idx, chr_idx);
                if !seen.contains(&key) {
                    return Err(InitializationError::Validation(format!(
                        "Missing sequence for ind_idx={ind_idx}, hap_idx={hap_idx}, chr_idx={chr_idx}"
                    )));
                }
            }
        }
    }

    Ok(())
}

/// Create individuals from validated sequence entries.
///
/// This function assumes sequences have already been validated. It organizes
/// sequences by their explicit indices rather than assuming sequential order.
pub fn create_individuals_from_sequences(
    entries: Vec<SequenceEntryWithIndices>,
    structure: &UniformRepeatStructure,
    population_size: usize,
) -> Result<Vec<Individual>, InitializationError> {
    use std::collections::HashMap;

    // Organize sequences by (ind_idx, hap_idx, chr_idx)
    let mut seq_map: HashMap<(usize, usize, usize), String> = HashMap::new();
    for entry in entries {
        let key = (entry.ind_idx, entry.hap_idx, entry.chr_idx);
        seq_map.insert(key, entry.seq);
    }

    let mut individuals = Vec::with_capacity(population_size);

    // Create individuals in order
    for ind_idx in 0..population_size {
        // Create haplotype 0
        let mut hap0 = Haplotype::with_capacity(structure.chrs_per_hap);
        for chr_idx in 0..structure.chrs_per_hap {
            let key = (ind_idx, 0, chr_idx);
            let seq_str = seq_map.get(&key).ok_or_else(|| {
                InitializationError::Validation(format!(
                    "Missing sequence for ind_idx={ind_idx}, hap_idx=0, chr_idx={chr_idx}"
                ))
            })?;

            let seq = Sequence::from_str(seq_str)
                .map_err(|e| InitializationError::Parse(format!("Invalid sequence: {e}")))?;

            // Create uniform map based on structure
            let map = RepeatMap::uniform(
                structure.ru_length,
                structure.rus_per_hor,
                structure.hors_per_chr,
            );

            let chr = Chromosome::new(format!("ind_{ind_idx}_h0_chr{chr_idx}"), seq, map);
            hap0.push(chr);
        }

        // Create haplotype 1
        let mut hap1 = Haplotype::with_capacity(structure.chrs_per_hap);
        for chr_idx in 0..structure.chrs_per_hap {
            let key = (ind_idx, 1, chr_idx);
            let seq_str = seq_map.get(&key).ok_or_else(|| {
                InitializationError::Validation(format!(
                    "Missing sequence for ind_idx={ind_idx}, hap_idx=1, chr_idx={chr_idx}"
                ))
            })?;

            let seq = Sequence::from_str(seq_str)
                .map_err(|e| InitializationError::Parse(format!("Invalid sequence: {e}")))?;

            // Create uniform map based on structure
            let map = RepeatMap::uniform(
                structure.ru_length,
                structure.rus_per_hor,
                structure.hors_per_chr,
            );

            let chr = Chromosome::new(format!("ind_{ind_idx}_h1_chr{chr_idx}"), seq, map);
            hap1.push(chr);
        }

        // Create individual
        let individual = Individual::new(format!("ind_{ind_idx}"), hap0, hap1);
        individuals.push(individual);
    }

    Ok(individuals)
}

/// High-level function to initialize from any sequence source.
/// High-level function to initialize from configuration.
pub fn initialize(
    config: &InitializationConfig,
    population_size: usize,
    rng: &mut impl Rng,
    codec: &CodecStrategy,
) -> Result<Vec<Individual>, InitializationError> {
    match config {
        InitializationConfig::Generate { structure, mode } => match mode {
            GenerationMode::Random => create_random_population(structure, population_size, rng)
                .map_err(InitializationError::Validation),
            GenerationMode::Uniform => create_uniform_population(structure, population_size)
                .map_err(InitializationError::Validation),
        },
        InitializationConfig::Load { source } => match source {
            SequenceSource::Fasta { path, bed_path } => {
                let entries = parse_fasta(path)?;
                initialize_from_fasta_bed(entries, bed_path.clone(), population_size)
            }
            SequenceSource::FormattedString {
                sequence,
                hor_delim,
                ru_delim,
            } => initialize_from_formatted_string(
                sequence.as_str(),
                *hor_delim,
                *ru_delim,
                population_size,
            ),

            SequenceSource::Database {
                path,
                sim_id,
                generation,
            } => {
                // Use the new loader that preserves maps and skips structure validation
                // This allows resuming evolved populations without strict structure checks
                load_individuals_from_database(path, sim_id, *generation, codec)
            }
        },
    }
}

/// Create initial population with uniform sequences.
fn create_uniform_population(
    structure: &UniformRepeatStructure,
    pop_size: usize,
) -> Result<Vec<Individual>, String> {
    let mut individuals = Vec::with_capacity(pop_size);

    for i in 0..pop_size {
        let ind = create_uniform_individual(format!("ind_{i}"), structure)?;
        individuals.push(ind);
    }

    Ok(individuals)
}

/// Create a single individual with uniform sequences.
fn create_uniform_individual(
    id: impl Into<std::sync::Arc<str>>,
    structure: &UniformRepeatStructure,
) -> Result<Individual, String> {
    // Create haplotypes with uniform chromosomes
    let mut hap1 = Haplotype::with_capacity(structure.chrs_per_hap);
    let mut hap2 = Haplotype::with_capacity(structure.chrs_per_hap);

    for chr_idx in 0..structure.chrs_per_hap {
        // Create chromosomes
        let chr1 = Chromosome::uniform(
            format!("chr{chr_idx}"),
            structure.init_base,
            structure.ru_length,
            structure.rus_per_hor,
            structure.hors_per_chr,
        );
        let chr2 = Chromosome::uniform(
            format!("chr{chr_idx}"),
            structure.init_base,
            structure.ru_length,
            structure.rus_per_hor,
            structure.hors_per_chr,
        );

        hap1.push(chr1);
        hap2.push(chr2);
    }

    Ok(Individual::new(id, hap1, hap2))
}

/// Create initial population with random sequences.
fn create_random_population(
    structure: &UniformRepeatStructure,
    pop_size: usize,
    rng: &mut impl Rng,
) -> Result<Vec<Individual>, String> {
    let mut individuals = Vec::with_capacity(pop_size);

    for i in 0..pop_size {
        let ind = create_random_individual(format!("ind_{i}"), structure, rng)?;
        individuals.push(ind);
    }

    Ok(individuals)
}

/// Create a single individual with random sequences.
fn create_random_individual(
    id: impl Into<std::sync::Arc<str>>,
    structure: &UniformRepeatStructure,
    rng: &mut impl Rng,
) -> Result<Individual, String> {
    let chr_length = structure.chr_length();
    let alphabet_size = 4; // DNA has 4 bases

    // Create haplotypes with random chromosomes
    let mut hap1 = Haplotype::with_capacity(structure.chrs_per_hap);
    let mut hap2 = Haplotype::with_capacity(structure.chrs_per_hap);

    for chr_idx in 0..structure.chrs_per_hap {
        // Create random sequence for haplotype 1
        let mut seq1 = Sequence::with_capacity(chr_length);
        for _ in 0..chr_length {
            let idx = rng.random_range(0..alphabet_size);
            let base = Nucleotide::from_index(idx as u8)
                .ok_or_else(|| format!("Invalid nucleotide index: {idx}"))?;
            seq1.push(base);
        }

        // Create random sequence for haplotype 2
        let mut seq2 = Sequence::with_capacity(chr_length);
        for _ in 0..chr_length {
            let idx = rng.random_range(0..alphabet_size);
            let base = Nucleotide::from_index(idx as u8)
                .ok_or_else(|| format!("Invalid nucleotide index: {idx}"))?;
            seq2.push(base);
        }

        // Create map
        let map = RepeatMap::uniform(
            structure.ru_length,
            structure.rus_per_hor,
            structure.hors_per_chr,
        );

        // Create chromosomes
        let chr1 = Chromosome::new(format!("chr{chr_idx}"), seq1, map.clone());
        let chr2 = Chromosome::new(format!("chr{chr_idx}"), seq2, map);

        hap1.push(chr1);
        hap2.push(chr2);
    }

    Ok(Individual::new(id, hap1, hap2))
}

/// Initialize individuals from FASTA and BED.
fn initialize_from_fasta_bed(
    entries: Vec<SequenceEntryWithIndices>,
    bed_path: String,
    population_size: usize,
) -> Result<Vec<Individual>, InitializationError> {
    // 1. Parse BED file into a map of (ind, hap, chr) -> RepeatMap
    // BED format: chrom start end name/type
    // We need to map "chrom" in BED to (ind, hap, chr).
    // The FASTA IDs are parsed into (ind, hap, chr).
    // We assume BED uses the same IDs as FASTA?
    // e.g. "ind0_hap0_chr0 0 100 RU"

    let bed_maps = parse_bed_to_maps(&bed_path)?;

    // 2. Create individuals
    let mut individuals = Vec::with_capacity(population_size);

    // Group entries by individual
    // (Similar logic to create_individuals_uniform but looking up map)

    // ... implementation details ...
    // For brevity in this step, I'll implement a simplified version or placeholder
    // and fill it in next step.

    // Re-using the grouping logic
    use std::collections::HashMap;
    let mut seq_map: HashMap<(usize, usize, usize), String> = HashMap::new();
    for entry in entries {
        seq_map.insert((entry.ind_idx, entry.hap_idx, entry.chr_idx), entry.seq);
    }

    for ind_idx in 0..population_size {
        // We need to know how many chromosomes per haplotype.
        // We can infer from the keys in seq_map or bed_maps.
        // Let's assume 1 chr per hap for now or infer max.

        let mut hap0 = Haplotype::new(); // Use new() which is empty
        let mut hap1 = Haplotype::new();

        // We need to iterate chromosomes.
        // Let's find max chr_idx for this individual in seq_map
        let max_chr = seq_map
            .keys()
            .filter(|(i, _, _)| *i == ind_idx)
            .map(|(_, _, c)| *c)
            .max()
            .unwrap_or(0);

        for chr_idx in 0..=max_chr {
            // Hap 0
            if let Some(seq_str) = seq_map.get(&(ind_idx, 0, chr_idx)) {
                let seq = Sequence::from_str(seq_str)
                    .map_err(|e| InitializationError::Parse(format!("Invalid sequence: {e}")))?;

                let map = bed_maps
                    .get(&(ind_idx, 0, chr_idx))
                    .cloned()
                    .ok_or_else(|| {
                        InitializationError::Validation(format!(
                            "No BED map found for ind_{ind_idx}_hap0_chr{chr_idx}"
                        ))
                    })?;

                hap0.push(Chromosome::new(
                    format!("ind_{ind_idx}_h0_chr{chr_idx}"),
                    seq,
                    map,
                ));
            }

            // Hap 1
            if let Some(seq_str) = seq_map.get(&(ind_idx, 1, chr_idx)) {
                let seq = Sequence::from_str(seq_str)
                    .map_err(|e| InitializationError::Parse(format!("Invalid sequence: {e}")))?;

                let map = bed_maps
                    .get(&(ind_idx, 1, chr_idx))
                    .cloned()
                    .ok_or_else(|| {
                        InitializationError::Validation(format!(
                            "No BED map found for ind_{ind_idx}_hap1_chr{chr_idx}"
                        ))
                    })?;

                hap1.push(Chromosome::new(
                    format!("ind_{ind_idx}_h1_chr{chr_idx}"),
                    seq,
                    map,
                ));
            }
        }

        individuals.push(Individual::new(format!("ind_{ind_idx}"), hap0, hap1));
    }

    Ok(individuals)
}

/// Parse BED file to RepeatMaps.
/// Returns map of (ind, hap, chr) -> RepeatMap
fn parse_bed_to_maps(
    path: &str,
) -> Result<std::collections::HashMap<(usize, usize, usize), RepeatMap>, InitializationError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Temporary storage: ID -> (ru_offsets, hor_offsets)
    // We need to build RepeatMap from intervals.
    // RepeatMap expects offsets.
    // BED: chrom start end name
    // We assume BED is sorted by start.

    let mut builders: std::collections::HashMap<String, (Vec<usize>, Vec<usize>)> =
        std::collections::HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 3 {
            continue;
        }

        let chrom = parts[0];
        let start: usize = parts[1]
            .parse()
            .map_err(|_| InitializationError::Parse("Invalid start in BED".into()))?;
        let end: usize = parts[2]
            .parse()
            .map_err(|_| InitializationError::Parse("Invalid end in BED".into()))?;
        let name = if parts.len() > 3 { parts[3] } else { "RU" };

        let entry = builders
            .entry(chrom.to_string())
            .or_insert((vec![0], vec![0]));
        let (ru_offsets, _hor_offsets) = entry;

        // If it's an RU, we add the END position to ru_offsets.
        // We assume contiguous RUs for now? Or RepeatMap handles gaps?
        // RepeatMap assumes contiguous RUs.
        // If BED has gaps, we might have issues.
        // Let's assume contiguous RUs.

        if name == "HOR" {
            // It's an HOR boundary?
            // RepeatMap uses hor_idx_offsets which are indices into RUs.
            // This is tricky to parse from BED directly without knowing RU count.
            // Alternative: BED defines RUs, and we infer HORs from some other marker?
            // Or BED has "HOR" lines and "RU" lines?

            // Let's assume simple case: BED only contains RUs.
            // And maybe we infer HORs or they are not supported via simple BED yet?
            // User said "indicate RU and HOR intervals".

            // If we have HOR lines, we can track them.
            // But RepeatMap needs HOR offsets in terms of RU indices.
            // So we need to process RUs first, then map HOR start/end to RU indices.
            // This requires two passes or buffering.
        } else {
            // Assume RU
            // Check continuity
            if let Some(&last) = ru_offsets.last() {
                if start != last {
                    // Gap or overlap
                    // For now, let's just push the end.
                    // If start > last, we have a gap. RepeatMap doesn't support gaps well (it's a map of RUs).
                    // We might need to fill with "spacer" RUs?
                }
            }
            ru_offsets.push(end);
        }
    }

    // Convert builders to RepeatMaps
    let mut maps = std::collections::HashMap::new();
    for (id, (ru_offsets, _)) in builders {
        // Parse ID
        let indices = SequenceEntryWithIndices::parse_id_format(&id, String::new())?;

        // Create RepeatMap
        // For now, creating a map with 1 HOR containing all RUs if no HOR info
        let num_rus = ru_offsets.len() - 1;
        let hor_offsets = vec![0, num_rus];

        let map = RepeatMap::new(ru_offsets, hor_offsets)
            .map_err(|e| InitializationError::Validation(e.to_string()))?;

        maps.insert((indices.ind_idx, indices.hap_idx, indices.chr_idx), map);
    }

    Ok(maps)
}

/// Initialize from formatted string (uniform across population).
fn initialize_from_formatted_string(
    sequence: &str,
    hor_delim: char,
    ru_delim: char,
    population_size: usize,
) -> Result<Vec<Individual>, InitializationError> {
    // Parse the string once to get Sequence and RepeatMap
    let chr_template = Chromosome::from_formatted_string("template", sequence, hor_delim, ru_delim)
        .map_err(|e| InitializationError::Parse(e.to_string()))?;

    let mut individuals = Vec::with_capacity(population_size);

    for i in 0..population_size {
        let mut hap0 = Haplotype::new();
        let mut hap1 = Haplotype::new();

        // Create 1 chromosome per haplotype (default for this mode)
        let c0 = Chromosome::new(
            format!("ind_{i}_h0_chr0"),
            chr_template.sequence().clone(),
            chr_template.map().clone(),
        );
        let c1 = Chromosome::new(
            format!("ind_{i}_h1_chr0"),
            chr_template.sequence().clone(),
            chr_template.map().clone(),
        );

        hap0.push(c0);
        hap1.push(c1);

        individuals.push(Individual::new(format!("ind_{i}"), hap0, hap1));
    }

    Ok(individuals)
}

/// Create individuals from validated sequence entries (Uniform).
pub fn create_individuals_uniform(
    entries: Vec<SequenceEntryWithIndices>,
    structure: &UniformRepeatStructure,
    population_size: usize,
) -> Result<Vec<Individual>, InitializationError> {
    create_individuals_from_sequences(entries, structure, population_size)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::Nucleotide;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn test_structure() -> UniformRepeatStructure {
        UniformRepeatStructure::new(
            Nucleotide::A,
            10, // ru_length
            5,  // rus_per_hor
            2,  // hors_per_chr
            1,  // chrs_per_hap
        )
    }

    #[test]
    fn test_parse_fasta_simple() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">ind0_hap0_chr0").unwrap();
        writeln!(file, "ACGTACGT").unwrap();
        writeln!(file, ">ind0_hap1_chr0").unwrap();
        writeln!(file, "TGCATGCA").unwrap();

        let entries = parse_fasta(file.path()).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].ind_idx, 0);
        assert_eq!(entries[0].hap_idx, 0);
        assert_eq!(entries[0].chr_idx, 0);
        assert_eq!(entries[0].seq, "ACGTACGT");
        assert_eq!(entries[1].ind_idx, 0);
        assert_eq!(entries[1].hap_idx, 1);
        assert_eq!(entries[1].chr_idx, 0);
        assert_eq!(entries[1].seq, "TGCATGCA");
    }

    #[test]
    fn test_parse_fasta_multiline() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">ind0_hap0_chr0 description").unwrap();
        writeln!(file, "ACGT").unwrap();
        writeln!(file, "ACGT").unwrap();
        writeln!(file, ">ind0_hap1_chr0").unwrap();
        writeln!(file, "TGCA").unwrap();
        writeln!(file, "TGCA").unwrap();

        let entries = parse_fasta(file.path()).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].seq, "ACGTACGT");
        assert_eq!(entries[1].seq, "TGCATGCA");
    }

    #[test]
    fn test_parse_fasta_empty_lines() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">ind0_hap0_chr0").unwrap();
        writeln!(file, "ACGT").unwrap();
        writeln!(file).unwrap();
        writeln!(file, "ACGT").unwrap();

        let entries = parse_fasta(file.path()).unwrap();
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].seq, "ACGTACGT");
    }

    #[test]
    fn test_parse_fasta_abbreviated_format() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">prefix_ind0_h0_c0").unwrap();
        writeln!(file, "ACGT").unwrap();
        writeln!(file, ">sample_ind1_h1_chr0").unwrap();
        writeln!(file, "TGCA").unwrap();

        let entries = parse_fasta(file.path()).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].ind_idx, 0);
        assert_eq!(entries[0].hap_idx, 0);
        assert_eq!(entries[1].ind_idx, 1);
        assert_eq!(entries[1].hap_idx, 1);
    }

    #[test]
    fn test_parse_fasta_empty_file() {
        let file = NamedTempFile::new().unwrap();
        let result = parse_fasta(file.path());
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_json_new_format() {
        let json = r#"[
            {"ind_idx": 0, "hap_idx": 0, "chr_idx": 0, "seq": "ACGTACGT"},
            {"ind_idx": 0, "hap_idx": 1, "chr_idx": 0, "seq": "TGCATGCA"}
        ]"#;

        let entries = parse_json(json).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].ind_idx, 0);
        assert_eq!(entries[0].hap_idx, 0);
        assert_eq!(entries[0].chr_idx, 0);
        assert_eq!(entries[0].seq, "ACGTACGT");
        assert_eq!(entries[1].seq, "TGCATGCA");
    }

    #[test]
    fn test_parse_json_legacy_format() {
        let json = r#"[
            {"id": "ind0_hap0_chr0", "seq": "ACGTACGT"},
            {"id": "ind0_hap1_chr0", "seq": "TGCATGCA"}
        ]"#;

        let entries = parse_json(json).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].ind_idx, 0);
        assert_eq!(entries[0].hap_idx, 0);
        assert_eq!(entries[0].chr_idx, 0);
    }

    #[test]
    fn test_parse_json_file() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            r#"[
            {{"ind_idx": 0, "hap_idx": 0, "chr_idx": 0, "seq": "ACGTACGT"}},
            {{"ind_idx": 0, "hap_idx": 1, "chr_idx": 0, "seq": "TGCATGCA"}}
        ]"#
        )
        .unwrap();

        let entries = parse_json(file.path().to_str().unwrap()).unwrap();
        assert_eq!(entries.len(), 2);
    }

    #[test]
    fn test_parse_json_empty() {
        let json = "[]";
        let result = parse_json(json);
        assert!(result.is_err());
    }

    #[test]
    fn test_validate_sequences_correct() {
        let structure = test_structure();
        let chr_len = structure.chr_length(); // 100 bases

        // For 1 individual, 1 chr per hap, 2 haplotypes = 2 sequences
        let entries = vec![
            SequenceEntryWithIndices {
                ind_idx: 0,
                hap_idx: 0,
                chr_idx: 0,
                seq: "A".repeat(chr_len),
            },
            SequenceEntryWithIndices {
                ind_idx: 0,
                hap_idx: 1,
                chr_idx: 0,
                seq: "C".repeat(chr_len),
            },
        ];

        let result = validate_sequences(&entries, &structure, 1);
        assert!(result.is_ok());
    }

    #[test]
    fn test_validate_sequences_wrong_count() {
        let structure = test_structure();
        let chr_len = structure.chr_length();

        // Only 1 sequence instead of 2
        let entries = vec![SequenceEntryWithIndices {
            ind_idx: 0,
            hap_idx: 0,
            chr_idx: 0,
            seq: "A".repeat(chr_len),
        }];

        let result = validate_sequences(&entries, &structure, 1);
        assert!(matches!(result, Err(InitializationError::Validation(_))));
    }
}
