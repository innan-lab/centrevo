//! Initialization of simulations from custom sequences.
//!
//! This module provides functionality to initialize simulations from various
//! sequence sources including FASTA files, JSON data, and previous simulation
//! databases. It validates that input sequences match the expected parameters.

use crate::base::{Alphabet, Sequence};
use crate::genome::{Chromosome, Haplotype, Individual};
use crate::simulation::RepeatStructure;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

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
                let digits: String = after.chars()
                    .take_while(|c| c.is_ascii_digit())
                    .collect();
                if !digits.is_empty() && let Ok(idx) = digits.parse::<usize>() {
                    return Some(idx);
                }
            }
        }
        None
    }
}

/// Input source for custom sequences.
#[derive(Debug, Clone)]
pub enum SequenceInput {
    /// FASTA file path
    Fasta(String),
    /// JSON file path or JSON string
    Json(String),
    /// Database path and simulation ID
    Database { path: String, sim_id: String, generation: Option<usize> },
}

/// Error types for sequence initialization.
#[derive(Debug)]
pub enum InitializationError {
    /// IO error
    Io(std::io::Error),
    /// Parse error
    Parse(String),
    /// Validation error
    Validation(String),
    /// Database error
    Database(String),
}

impl std::fmt::Display for InitializationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(e) => write!(f, "IO error: {e}"),
            Self::Parse(msg) => write!(f, "Parse error: {msg}"),
            Self::Validation(msg) => write!(f, "Validation error: {msg}"),
            Self::Database(msg) => write!(f, "Database error: {msg}"),
        }
    }
}

impl std::error::Error for InitializationError {}

impl From<std::io::Error> for InitializationError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}

impl From<serde_json::Error> for InitializationError {
    fn from(e: serde_json::Error) -> Self {
        Self::Parse(format!("JSON error: {e}"))
    }
}

/// Convert sequence indices (Vec<u8>) to a String representation.
fn indices_to_string(indices: &[u8], alphabet: &Alphabet) -> Result<String, InitializationError> {
    let mut result = String::with_capacity(indices.len());
    for &idx in indices {
        let ch = alphabet.get_char(idx)
            .ok_or_else(|| InitializationError::Parse(
                format!("Invalid base index {idx} in sequence")
            ))?;
        result.push(ch);
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
pub fn parse_fasta(path: impl AsRef<Path>) -> Result<Vec<SequenceEntryWithIndices>, InitializationError> {
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
                    let entry = SequenceEntryWithIndices::parse_id_format(&id, current_seq.clone())?;
                    entries.push(entry);
                    current_seq.clear();
                }
            }
            
            // Parse new header
            // Remove '>'
            let id = header.split_whitespace().next()
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
        return Err(InitializationError::Parse("No sequences found in FASTA file".to_string()));
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
                return Err(InitializationError::Parse("Empty JSON array in file".to_string()));
            }
            return Ok(entries);
        }
    }
    
    // Try as file path - legacy format
    let file = File::open(input)?;
    let reader = BufReader::new(file);
    let entries: Vec<SequenceEntry> = serde_json::from_reader(reader)?;
    
    if entries.is_empty() {
        return Err(InitializationError::Parse("Empty JSON array in file".to_string()));
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
) -> Result<Vec<SequenceEntryWithIndices>, InitializationError> {
    use crate::storage::QueryBuilder;
    
    let query = QueryBuilder::new(db_path.as_ref())
        .map_err(|e| InitializationError::Database(format!("Failed to open database: {e}")))?;
    
    // Determine which generation to load
    let target_gen = match generation {
        Some(g) => g,
        None => {
            // Get the last recorded generation
            let recorded_gens = query.get_recorded_generations(sim_id)
                .map_err(|e| InitializationError::Database(format!("Failed to get generations: {e}")))?;
            
            *recorded_gens.last()
                .ok_or_else(|| InitializationError::Database("No generations found".to_string()))?
        }
    };
    
    // Load population at that generation
    let snapshots = query.get_generation(sim_id, target_gen)
        .map_err(|e| InitializationError::Database(format!("Failed to load generation: {e}")))?;
    
    // Get config to convert indices to strings
    let config = query.get_full_config(sim_id)
        .map_err(|e| InitializationError::Database(format!("Failed to load config: {e}")))?;
    
    query.close()
        .map_err(|e| InitializationError::Database(format!("Failed to close database: {e}")))?;
    
    // Convert snapshots to sequence entries with indices
    let mut entries = Vec::new();
    
    // Database stores sequences per chromosome, so we need to group them properly
    // Assuming snapshots are for single-chromosome case (chr_idx = 0)
    // For multi-chromosome, database would need to be extended
    for (ind_idx, snapshot) in snapshots.iter().enumerate() {
        // Convert Vec<u8> indices to String
        let hap0_seq = indices_to_string(&snapshot.haplotype1_seq, &config.structure.alphabet)?;
        let hap1_seq = indices_to_string(&snapshot.haplotype2_seq, &config.structure.alphabet)?;
        
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
        return Err(InitializationError::Database("No sequences found in generation".to_string()));
    }
    
    Ok(entries)
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
    structure: &RepeatStructure,
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
                idx, entry.ind_idx, entry.hap_idx, entry.chr_idx,
                entry.seq.len(),
                expected_chr_length
            )));
        }
        
        // Check that all characters are valid bases
        for (pos, ch) in entry.seq.chars().enumerate() {
            if structure.alphabet.get_index(ch.to_ascii_uppercase()).is_none() {
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
    structure: &RepeatStructure,
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
            let seq_str = seq_map.get(&key)
                .ok_or_else(|| InitializationError::Validation(format!(
                    "Missing sequence for ind_idx={ind_idx}, hap_idx=0, chr_idx={chr_idx}"
                )))?;
            
            let seq = Sequence::from_str(seq_str, structure.alphabet.clone())
                .map_err(|e| InitializationError::Parse(format!("Invalid sequence: {e}")))?;
            
            let chr = Chromosome::new(
                format!("ind_{ind_idx}_h0_chr{chr_idx}"),
                seq,
                structure.ru_length,
                structure.rus_per_hor,
            );
            hap0.push(chr);
        }
        
        // Create haplotype 1
        let mut hap1 = Haplotype::with_capacity(structure.chrs_per_hap);
        for chr_idx in 0..structure.chrs_per_hap {
            let key = (ind_idx, 1, chr_idx);
            let seq_str = seq_map.get(&key)
                .ok_or_else(|| InitializationError::Validation(format!(
                    "Missing sequence for ind_idx={ind_idx}, hap_idx=1, chr_idx={chr_idx}"
                )))?;
            
            let seq = Sequence::from_str(seq_str, structure.alphabet.clone())
                .map_err(|e| InitializationError::Parse(format!("Invalid sequence: {e}")))?;
            
            let chr = Chromosome::new(
                format!("ind_{ind_idx}_h1_chr{chr_idx}"),
                seq,
                structure.ru_length,
                structure.rus_per_hor,
            );
            hap1.push(chr);
        }
        
        // Create individual
        let individual = Individual::new(format!("ind_{ind_idx}"), hap0, hap1);
        individuals.push(individual);
    }
    
    Ok(individuals)
}

/// High-level function to initialize from any sequence source.
pub fn initialize_from_source(
    source: SequenceInput,
    structure: &RepeatStructure,
    population_size: usize,
) -> Result<Vec<Individual>, InitializationError> {
    // Parse sequences based on source type
    let entries = match source {
        SequenceInput::Fasta(path) => parse_fasta(path)?,
        SequenceInput::Json(input) => parse_json(&input)?,
        SequenceInput::Database { path, sim_id, generation } => {
            load_from_database(path, &sim_id, generation)?
        }
    };
    
    // Validate sequences
    validate_sequences(&entries, structure, population_size)?;
    
    // Create individuals
    create_individuals_from_sequences(entries, structure, population_size)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::Nucleotide;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn test_structure() -> RepeatStructure {
        RepeatStructure::new(
            Alphabet::dna(),
            Nucleotide::A,
            10,  // ru_length
            5,   // rus_per_hor
            2,   // hors_per_chr
            1,   // chrs_per_hap
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
        writeln!(file, "").unwrap();
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
        writeln!(file, r#"[
            {{"ind_idx": 0, "hap_idx": 0, "chr_idx": 0, "seq": "ACGTACGT"}},
            {{"ind_idx": 0, "hap_idx": 1, "chr_idx": 0, "seq": "TGCATGCA"}}
        ]"#).unwrap();
        
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
        let entries = vec![
            SequenceEntryWithIndices {
                ind_idx: 0,
                hap_idx: 0,
                chr_idx: 0,
                seq: "A".repeat(chr_len),
            },
        ];
        
        let result = validate_sequences(&entries, &structure, 1);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Expected 2"));
    }

    #[test]
    fn test_validate_sequences_wrong_length() {
        let structure = test_structure();
        
        let entries = vec![
            SequenceEntryWithIndices {
                ind_idx: 0,
                hap_idx: 0,
                chr_idx: 0,
                seq: "A".repeat(50), // Wrong length
            },
            SequenceEntryWithIndices {
                ind_idx: 0,
                hap_idx: 1,
                chr_idx: 0,
                seq: "C".repeat(100),
            },
        ];
        
        let result = validate_sequences(&entries, &structure, 1);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("length 50, expected 100"));
    }

    #[test]
    fn test_validate_sequences_invalid_base() {
        let structure = test_structure();
        let chr_len = structure.chr_length();
        
        let mut invalid_seq = "A".repeat(chr_len - 1);
        invalid_seq.push('N'); // Invalid base
        
        let entries = vec![
            SequenceEntryWithIndices {
                ind_idx: 0,
                hap_idx: 0,
                chr_idx: 0,
                seq: invalid_seq,
            },
            SequenceEntryWithIndices {
                ind_idx: 0,
                hap_idx: 1,
                chr_idx: 0,
                seq: "C".repeat(chr_len),
            },
        ];
        
        let result = validate_sequences(&entries, &structure, 1);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("invalid character"));
    }

    #[test]
    fn test_validate_sequences_missing_combination() {
        let structure = test_structure();
        let chr_len = structure.chr_length();
        
        // Missing ind0_hap1_chr0 (only providing one sequence when we need 2)
        let entries = vec![
            SequenceEntryWithIndices {
                ind_idx: 0,
                hap_idx: 0,
                chr_idx: 0,
                seq: "A".repeat(chr_len),
            },
        ];
        
        let result = validate_sequences(&entries, &structure, 1);
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        // Could be either "Missing sequence" or "Expected 2" depending on which check fails first
        assert!(err_msg.contains("Missing sequence") || err_msg.contains("Expected 2"));
    }

    #[test]
    fn test_validate_sequences_duplicate() {
        let structure = test_structure();
        let chr_len = structure.chr_length();
        
        // Duplicate ind0_hap0_chr0
        let entries = vec![
            SequenceEntryWithIndices {
                ind_idx: 0,
                hap_idx: 0,
                chr_idx: 0,
                seq: "A".repeat(chr_len),
            },
            SequenceEntryWithIndices {
                ind_idx: 0,
                hap_idx: 0,
                chr_idx: 0,
                seq: "A".repeat(chr_len),
            },
        ];
        
        let result = validate_sequences(&entries, &structure, 1);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Duplicate"));
    }

    #[test]
    fn test_create_individuals_from_sequences() {
        let structure = test_structure();
        let chr_len = structure.chr_length();
        
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
        
        let individuals = create_individuals_from_sequences(entries, &structure, 1).unwrap();
        assert_eq!(individuals.len(), 1);
        assert_eq!(individuals[0].id(), "ind_0");
        assert_eq!(individuals[0].haplotype1().len(), 1);
        assert_eq!(individuals[0].haplotype2().len(), 1);
    }

    #[test]
    fn test_create_individuals_multiple() {
        let structure = test_structure();
        let chr_len = structure.chr_length();
        
        // 2 individuals × 1 chr per hap × 2 haplotypes = 4 sequences
        // Test that order doesn't matter
        let entries = vec![
            SequenceEntryWithIndices { ind_idx: 1, hap_idx: 0, chr_idx: 0, seq: "G".repeat(chr_len) },
            SequenceEntryWithIndices { ind_idx: 0, hap_idx: 1, chr_idx: 0, seq: "C".repeat(chr_len) },
            SequenceEntryWithIndices { ind_idx: 0, hap_idx: 0, chr_idx: 0, seq: "A".repeat(chr_len) },
            SequenceEntryWithIndices { ind_idx: 1, hap_idx: 1, chr_idx: 0, seq: "T".repeat(chr_len) },
        ];
        
        let individuals = create_individuals_from_sequences(entries, &structure, 2).unwrap();
        assert_eq!(individuals.len(), 2);
        assert_eq!(individuals[0].id(), "ind_0");
        assert_eq!(individuals[1].id(), "ind_1");
    }

    #[test]
    fn test_create_individuals_multi_chromosome() {
        // Structure with 2 chromosomes per haplotype
        let structure = RepeatStructure::new(
            Alphabet::dna(),
            Nucleotide::A,
            10,  // ru_length
            5,   // rus_per_hor
            2,   // hors_per_chr
            2,   // chrs_per_hap - THIS IS THE KEY
        );
        let chr_len = structure.chr_length();
        
        // 1 individual × 2 chr per hap × 2 haplotypes = 4 sequences
        let entries = vec![
            SequenceEntryWithIndices { ind_idx: 0, hap_idx: 0, chr_idx: 0, seq: "A".repeat(chr_len) },
            SequenceEntryWithIndices { ind_idx: 0, hap_idx: 0, chr_idx: 1, seq: "C".repeat(chr_len) },
            SequenceEntryWithIndices { ind_idx: 0, hap_idx: 1, chr_idx: 0, seq: "G".repeat(chr_len) },
            SequenceEntryWithIndices { ind_idx: 0, hap_idx: 1, chr_idx: 1, seq: "T".repeat(chr_len) },
        ];
        
        let individuals = create_individuals_from_sequences(entries, &structure, 1).unwrap();
        assert_eq!(individuals.len(), 1);
        assert_eq!(individuals[0].haplotype1().len(), 2);  // 2 chromosomes
        assert_eq!(individuals[0].haplotype2().len(), 2);  // 2 chromosomes
    }

    #[test]
    fn test_sequence_entry_serialization() {
        let entry = SequenceEntryWithIndices {
            ind_idx: 0,
            hap_idx: 1,
            chr_idx: 2,
            seq: "ACGT".to_string(),
        };
        
        let json = serde_json::to_string(&entry).unwrap();
        let deserialized: SequenceEntryWithIndices = serde_json::from_str(&json).unwrap();
        
        assert_eq!(entry.ind_idx, deserialized.ind_idx);
        assert_eq!(entry.hap_idx, deserialized.hap_idx);
        assert_eq!(entry.chr_idx, deserialized.chr_idx);
        assert_eq!(entry.seq, deserialized.seq);
    }

    #[test]
    fn test_legacy_sequence_entry_serialization() {
        let entry = SequenceEntry {
            id: "test".to_string(),
            seq: "ACGT".to_string(),
        };
        
        let json = serde_json::to_string(&entry).unwrap();
        let deserialized: SequenceEntry = serde_json::from_str(&json).unwrap();
        
        assert_eq!(entry.id, deserialized.id);
        assert_eq!(entry.seq, deserialized.seq);
    }

    #[test]
    fn test_initialization_error_display() {
        let err = InitializationError::Validation("test error".to_string());
        assert!(err.to_string().contains("Validation error"));
        assert!(err.to_string().contains("test error"));
    }

    #[test]
    fn test_id_parsing_variants() {
        // Test various ID formats
        let test_cases = vec![
            ("ind0_hap0_chr0", 0, 0, 0),
            ("prefix_ind1_hap1_chr2", 1, 1, 2),
            ("ind5_h3_chr4", 5, 3, 4),
            ("sample_ind10_h1_c7", 10, 1, 7),
        ];

        for (id, exp_ind, exp_hap, exp_chr) in test_cases {
            let entry = SequenceEntryWithIndices::parse_id_format(id, "ACGT".to_string()).unwrap();
            assert_eq!(entry.ind_idx, exp_ind, "Failed for ID: {}", id);
            assert_eq!(entry.hap_idx, exp_hap, "Failed for ID: {}", id);
            assert_eq!(entry.chr_idx, exp_chr, "Failed for ID: {}", id);
        }
    }

    #[test]
    fn test_id_parsing_invalid() {
        let invalid_ids = vec![
            "no_indices_here",
            "ind_only_ind0",
            "hap_only_hap0",
        ];

        for id in invalid_ids {
            let result = SequenceEntryWithIndices::parse_id_format(id, "ACGT".to_string());
            assert!(result.is_err(), "Should have failed for ID: {}", id);
        }
    }
}
