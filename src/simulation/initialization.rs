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

/// A parsed sequence entry with ID and sequence data.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequenceEntry {
    /// Sequence identifier
    pub id: String,
    /// Sequence data as string
    pub seq: String,
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
            Self::Io(e) => write!(f, "IO error: {}", e),
            Self::Parse(msg) => write!(f, "Parse error: {}", msg),
            Self::Validation(msg) => write!(f, "Validation error: {}", msg),
            Self::Database(msg) => write!(f, "Database error: {}", msg),
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
        Self::Parse(format!("JSON error: {}", e))
    }
}

/// Convert sequence indices (Vec<u8>) to a String representation.
fn indices_to_string(indices: &[u8], alphabet: &Alphabet) -> Result<String, InitializationError> {
    let mut result = String::with_capacity(indices.len());
    for &idx in indices {
        let ch = alphabet.get_char(idx)
            .ok_or_else(|| InitializationError::Parse(
                format!("Invalid base index {} in sequence", idx)
            ))?;
        result.push(ch);
    }
    Ok(result)
}

/// Parse sequences from a FASTA file.
///
/// FASTA format:
/// ```text
/// >sequence_id description
/// ACGTACGTACGT
/// >another_id
/// TGCATGCA
/// ```
pub fn parse_fasta(path: impl AsRef<Path>) -> Result<Vec<SequenceEntry>, InitializationError> {
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
        
        if line.starts_with('>') {
            // Save previous entry if exists
            if let Some(id) = current_id.take() {
                if !current_seq.is_empty() {
                    entries.push(SequenceEntry {
                        id,
                        seq: current_seq.clone(),
                    });
                    current_seq.clear();
                }
            }
            
            // Parse new header
            let header = &line[1..]; // Remove '>'
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
            entries.push(SequenceEntry {
                id,
                seq: current_seq,
            });
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
/// Expected JSON structure:
/// ```json
/// [
///   {"id": "seq1", "seq": "ACGTACGT"},
///   {"id": "seq2", "seq": "TGCATGCA"}
/// ]
/// ```
pub fn parse_json(input: &str) -> Result<Vec<SequenceEntry>, InitializationError> {
    // Try to parse as JSON string first
    if let Ok(entries) = serde_json::from_str::<Vec<SequenceEntry>>(input) {
        if entries.is_empty() {
            return Err(InitializationError::Parse("Empty JSON array".to_string()));
        }
        return Ok(entries);
    }
    
    // Try as file path
    let file = File::open(input)?;
    let reader = BufReader::new(file);
    let entries: Vec<SequenceEntry> = serde_json::from_reader(reader)?;
    
    if entries.is_empty() {
        return Err(InitializationError::Parse("Empty JSON array in file".to_string()));
    }
    
    Ok(entries)
}

/// Load sequences from a previous simulation database.
///
/// Loads the population state from a specific generation (or the last generation
/// if not specified).
pub fn load_from_database(
    db_path: impl AsRef<Path>,
    sim_id: &str,
    generation: Option<usize>,
) -> Result<Vec<SequenceEntry>, InitializationError> {
    use crate::storage::QueryBuilder;
    
    let query = QueryBuilder::new(db_path.as_ref())
        .map_err(|e| InitializationError::Database(format!("Failed to open database: {}", e)))?;
    
    // Determine which generation to load
    let target_gen = match generation {
        Some(g) => g,
        None => {
            // Get the last recorded generation
            let recorded_gens = query.get_recorded_generations(sim_id)
                .map_err(|e| InitializationError::Database(format!("Failed to get generations: {}", e)))?;
            
            *recorded_gens.last()
                .ok_or_else(|| InitializationError::Database("No generations found".to_string()))?
        }
    };
    
    // Load population at that generation
    let snapshots = query.get_generation(sim_id, target_gen)
        .map_err(|e| InitializationError::Database(format!("Failed to load generation: {}", e)))?;
    
    // Get config to convert indices to strings
    let config = query.get_full_config(sim_id)
        .map_err(|e| InitializationError::Database(format!("Failed to load config: {}", e)))?;
    
    query.close()
        .map_err(|e| InitializationError::Database(format!("Failed to close database: {}", e)))?;
    
    // Convert snapshots to sequence entries
    let mut entries = Vec::new();
    
    for snapshot in snapshots {
        // Convert Vec<u8> indices to String
        let hap1_seq = indices_to_string(&snapshot.haplotype1_seq, &config.structure.alphabet)?;
        let hap2_seq = indices_to_string(&snapshot.haplotype2_seq, &config.structure.alphabet)?;
        
        // Add both haplotypes as separate sequences
        entries.push(SequenceEntry {
            id: format!("{}_hap1", snapshot.individual_id),
            seq: hap1_seq,
        });
        entries.push(SequenceEntry {
            id: format!("{}_hap2", snapshot.individual_id),
            seq: hap2_seq,
        });
    }
    
    if entries.is_empty() {
        return Err(InitializationError::Database("No sequences found in generation".to_string()));
    }
    
    Ok(entries)
}

/// Validate that sequences match the expected parameters.
pub fn validate_sequences(
    entries: &[SequenceEntry],
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
    
    // Validate each sequence
    for (idx, entry) in entries.iter().enumerate() {
        // Check length
        if entry.seq.len() != expected_chr_length {
            return Err(InitializationError::Validation(format!(
                "Sequence '{}' (index {}) has length {}, expected {}",
                entry.id,
                idx,
                entry.seq.len(),
                expected_chr_length
            )));
        }
        
        // Check that all characters are valid bases
        for (pos, ch) in entry.seq.chars().enumerate() {
            if structure.alphabet.get_index(ch.to_ascii_uppercase()).is_none() {
                return Err(InitializationError::Validation(format!(
                    "Sequence '{}' contains invalid character '{}' at position {}",
                    entry.id,
                    ch,
                    pos
                )));
            }
        }
    }
    
    Ok(())
}

/// Create individuals from validated sequence entries.
///
/// This function assumes sequences have already been validated and are in the
/// correct order. Sequences are grouped by individual and haplotype.
pub fn create_individuals_from_sequences(
    entries: Vec<SequenceEntry>,
    structure: &RepeatStructure,
    population_size: usize,
) -> Result<Vec<Individual>, InitializationError> {
    let seqs_per_individual = structure.chrs_per_hap * 2; // 2 haplotypes
    
    if entries.len() != population_size * seqs_per_individual {
        return Err(InitializationError::Validation(
            "Sequence count does not match expected population structure".to_string()
        ));
    }
    
    let mut individuals = Vec::with_capacity(population_size);
    
    // Process sequences in groups
    for ind_idx in 0..population_size {
        let start_idx = ind_idx * seqs_per_individual;
        let end_idx = start_idx + seqs_per_individual;
        let ind_seqs = &entries[start_idx..end_idx];
        
        // Split into two haplotypes
        let hap1_seqs = &ind_seqs[0..structure.chrs_per_hap];
        let hap2_seqs = &ind_seqs[structure.chrs_per_hap..];
        
        // Create haplotype 1
        let mut hap1 = Haplotype::with_capacity(structure.chrs_per_hap);
        for (chr_idx, seq_entry) in hap1_seqs.iter().enumerate() {
            let seq = Sequence::from_str(&seq_entry.seq, structure.alphabet.clone())
                .map_err(|e| InitializationError::Parse(format!("Invalid sequence: {}", e)))?;
            
            let chr = Chromosome::new(
                format!("ind_{}_h1_chr{}", ind_idx, chr_idx),
                seq,
                structure.ru_length,
                structure.rus_per_hor,
            );
            hap1.push(chr);
        }
        
        // Create haplotype 2
        let mut hap2 = Haplotype::with_capacity(structure.chrs_per_hap);
        for (chr_idx, seq_entry) in hap2_seqs.iter().enumerate() {
            let seq = Sequence::from_str(&seq_entry.seq, structure.alphabet.clone())
                .map_err(|e| InitializationError::Parse(format!("Invalid sequence: {}", e)))?;
            
            let chr = Chromosome::new(
                format!("ind_{}_h2_chr{}", ind_idx, chr_idx),
                seq,
                structure.ru_length,
                structure.rus_per_hor,
            );
            hap2.push(chr);
        }
        
        // Create individual
        let individual = Individual::new(format!("ind_{}", ind_idx), hap1, hap2);
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
        writeln!(file, ">seq1").unwrap();
        writeln!(file, "ACGTACGT").unwrap();
        writeln!(file, ">seq2").unwrap();
        writeln!(file, "TGCATGCA").unwrap();
        
        let entries = parse_fasta(file.path()).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].id, "seq1");
        assert_eq!(entries[0].seq, "ACGTACGT");
        assert_eq!(entries[1].id, "seq2");
        assert_eq!(entries[1].seq, "TGCATGCA");
    }

    #[test]
    fn test_parse_fasta_multiline() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">seq1 description").unwrap();
        writeln!(file, "ACGT").unwrap();
        writeln!(file, "ACGT").unwrap();
        writeln!(file, ">seq2").unwrap();
        writeln!(file, "TGCA").unwrap();
        writeln!(file, "TGCA").unwrap();
        
        let entries = parse_fasta(file.path()).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].id, "seq1");
        assert_eq!(entries[0].seq, "ACGTACGT");
        assert_eq!(entries[1].id, "seq2");
        assert_eq!(entries[1].seq, "TGCATGCA");
    }

    #[test]
    fn test_parse_fasta_empty_lines() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">seq1").unwrap();
        writeln!(file, "ACGT").unwrap();
        writeln!(file, "").unwrap();
        writeln!(file, "ACGT").unwrap();
        
        let entries = parse_fasta(file.path()).unwrap();
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].seq, "ACGTACGT");
    }

    #[test]
    fn test_parse_fasta_empty_file() {
        let file = NamedTempFile::new().unwrap();
        let result = parse_fasta(file.path());
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_json_string() {
        let json = r#"[
            {"id": "seq1", "seq": "ACGTACGT"},
            {"id": "seq2", "seq": "TGCATGCA"}
        ]"#;
        
        let entries = parse_json(json).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].id, "seq1");
        assert_eq!(entries[0].seq, "ACGTACGT");
        assert_eq!(entries[1].id, "seq2");
        assert_eq!(entries[1].seq, "TGCATGCA");
    }

    #[test]
    fn test_parse_json_file() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, r#"[
            {{"id": "seq1", "seq": "ACGTACGT"}},
            {{"id": "seq2", "seq": "TGCATGCA"}}
        ]"#).unwrap();
        
        let entries = parse_json(file.path().to_str().unwrap()).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].id, "seq1");
        assert_eq!(entries[1].id, "seq2");
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
            SequenceEntry {
                id: "seq1".to_string(),
                seq: "A".repeat(chr_len),
            },
            SequenceEntry {
                id: "seq2".to_string(),
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
            SequenceEntry {
                id: "seq1".to_string(),
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
            SequenceEntry {
                id: "seq1".to_string(),
                seq: "A".repeat(50), // Wrong length
            },
            SequenceEntry {
                id: "seq2".to_string(),
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
            SequenceEntry {
                id: "seq1".to_string(),
                seq: invalid_seq,
            },
            SequenceEntry {
                id: "seq2".to_string(),
                seq: "C".repeat(chr_len),
            },
        ];
        
        let result = validate_sequences(&entries, &structure, 1);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("invalid character"));
    }

    #[test]
    fn test_create_individuals_from_sequences() {
        let structure = test_structure();
        let chr_len = structure.chr_length();
        
        let entries = vec![
            SequenceEntry {
                id: "ind0_h1_chr0".to_string(),
                seq: "A".repeat(chr_len),
            },
            SequenceEntry {
                id: "ind0_h2_chr0".to_string(),
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
        let entries = vec![
            SequenceEntry { id: "ind0_h1".to_string(), seq: "A".repeat(chr_len) },
            SequenceEntry { id: "ind0_h2".to_string(), seq: "C".repeat(chr_len) },
            SequenceEntry { id: "ind1_h1".to_string(), seq: "G".repeat(chr_len) },
            SequenceEntry { id: "ind1_h2".to_string(), seq: "T".repeat(chr_len) },
        ];
        
        let individuals = create_individuals_from_sequences(entries, &structure, 2).unwrap();
        assert_eq!(individuals.len(), 2);
        assert_eq!(individuals[0].id(), "ind_0");
        assert_eq!(individuals[1].id(), "ind_1");
    }

    #[test]
    fn test_sequence_entry_serialization() {
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
}
