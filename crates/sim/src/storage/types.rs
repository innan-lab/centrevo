//! Shared types for storage.

use crate::base::FitnessValue;
use crate::genome::Individual;
use crate::simulation::{
    FitnessConfig, MutationConfig, Population, RecombinationConfig, SimulationConfig,
    UniformRepeatStructure,
};
use centrevo_codec::CodecStrategy;
use serde::{Deserialize, Serialize};

/// Complete simulation configuration for resumability.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationSnapshot {
    pub structure: UniformRepeatStructure,
    pub mutation: MutationConfig,
    pub recombination: RecombinationConfig,
    pub fitness: FitnessConfig,
    pub config: SimulationConfig,
}

/// Recording strategy for when to persist simulation state.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum RecordingStrategy {
    /// Record every N generations.
    EveryN(usize),

    /// Record at specific generations.
    Specific(Vec<usize>),

    /// Record all generations.
    All,

    /// No recording.
    None,
}

impl RecordingStrategy {
    /// Check if generation should be recorded
    pub fn should_record(&self, generation: usize) -> bool {
        match self {
            Self::EveryN(n) => generation.is_multiple_of(*n),
            Self::Specific(gens) => gens.contains(&generation),
            Self::All => true,
            Self::None => false,
        }
    }
}

/// Snapshot of an individual for database storage.
#[derive(Debug, Clone)]
pub struct IndividualSnapshot {
    pub individual_id: String,
    pub haplotype1_chr_id: String,
    pub haplotype1_seq: Vec<u8>,         // Encoded
    pub haplotype1_map: Option<Vec<u8>>, // Encoded serialized map
    pub haplotype2_chr_id: String,
    pub haplotype2_seq: Vec<u8>,         // Encoded
    pub haplotype2_map: Option<Vec<u8>>, // Encoded serialized map
    pub haplotype1_fitness: Option<f64>,
    pub haplotype2_fitness: Option<f64>,
    pub fitness: Option<f64>,
}

impl IndividualSnapshot {
    /// Create a snapshot from an individual.
    /// Encodes sequence using the provided strategy.
    /// Encodes structure using UnpackedRS (always).
    pub fn from_individual(ind: &Individual, codec: &CodecStrategy) -> Self {
        let h1 = ind.haplotype1();
        let h2 = ind.haplotype2();

        // Helper to process chromosome
        let process_chr =
            |chr_opt: Option<&crate::genome::Chromosome>| -> (String, Vec<u8>, Option<Vec<u8>>) {
                if let Some(chr) = chr_opt {
                    // 1. Sequence -> Configured Codec
                    // Convert Nucleotide -> u8 indices
                    let raw_seq: Vec<u8> = chr
                        .sequence()
                        .as_slice()
                        .iter()
                        .map(|n| n.to_index())
                        .collect();

                    let encoded_seq = codec.encode(&raw_seq).expect("Failed to encode sequence");

                    // 2. Map -> Serialize -> UnpackedRS (Fixed)
                    let map_bytes = bincode::serialize(chr.map()).unwrap_or_default();
                    let encoded_map = CodecStrategy::UnpackedRS.encode(&map_bytes).ok(); // Option

                    (chr.id().to_string(), encoded_seq, encoded_map)
                } else {
                    (String::new(), Vec::new(), None)
                }
            };

        let (h1_chr_id, h1_seq, h1_map) = process_chr(h1.get(0));
        let (h2_chr_id, h2_seq, h2_map) = process_chr(h2.get(0));

        Self {
            individual_id: ind.id().to_string(),
            haplotype1_chr_id: h1_chr_id,
            haplotype1_seq: h1_seq,
            haplotype1_map: h1_map,
            haplotype2_chr_id: h2_chr_id,
            haplotype2_seq: h2_seq,
            haplotype2_map: h2_map,
            haplotype1_fitness: h1.cached_fitness().map(|f| *f),
            haplotype2_fitness: h2.cached_fitness().map(|f| *f),
            fitness: ind.cached_fitness().map(|f| *f),
        }
    }

    /// Reconstruct an Individual from a snapshot.
    /// Uses the provided codec for sequence, and UnpackedRS for map.
    pub fn to_individual(&self, codec: &CodecStrategy) -> Result<Individual, String> {
        use crate::base::{Nucleotide, Sequence};
        use crate::genome::repeat_map::RepeatMap;
        use crate::genome::{Chromosome, Haplotype};

        // Helper to reconstruct chromosome
        let reconstruct_chr =
            |id: &str, seq_data: &[u8], map_data: Option<&[u8]>| -> Result<Chromosome, String> {
                // 1. Decode Sequence
                let raw_seq = codec
                    .decode(seq_data)
                    .map_err(|e| format!("Seq Decode: {e}"))?;

                let nucs: Vec<Nucleotide> = raw_seq
                    .iter()
                    .map(|&i| Nucleotide::from_index(i).unwrap_or(Nucleotide::A))
                    .collect();
                let seq = Sequence::from_nucleotides(nucs);

                // 2. Decode Map
                let map = if let Some(bytes) = map_data {
                    let raw_map_bytes = CodecStrategy::UnpackedRS
                        .decode(bytes)
                        .map_err(|e| format!("Map Decode: {e}"))?;

                    bincode::deserialize::<RepeatMap>(&raw_map_bytes)
                        .map_err(|e| format!("Map Deser: {e}"))?
                } else {
                    return Err("Missing RepeatMap in snapshot".to_string());
                };

                Ok(Chromosome::new(id.to_string(), seq, map))
            };

        let chr1 = reconstruct_chr(
            &self.haplotype1_chr_id,
            &self.haplotype1_seq,
            self.haplotype1_map.as_deref(),
        )?;
        let chr2 = reconstruct_chr(
            &self.haplotype2_chr_id,
            &self.haplotype2_seq,
            self.haplotype2_map.as_deref(),
        )?;

        // Create haplotypes
        let mut hap1 = Haplotype::new();
        hap1.push(chr1);
        if let Some(f) = self.haplotype1_fitness {
            hap1.set_cached_fitness(FitnessValue::new(f));
        }

        let mut hap2 = Haplotype::new();
        hap2.push(chr2);
        if let Some(f) = self.haplotype2_fitness {
            hap2.set_cached_fitness(FitnessValue::new(f));
        }

        // Create individual
        let mut individual = Individual::new(self.individual_id.as_str(), hap1, hap2);
        if let Some(f) = self.fitness {
            individual.set_cached_fitness(FitnessValue::new(f));
        }

        Ok(individual)
    }
}

/// Aggregated fitness statistics for a generation.
#[derive(Debug, Clone, Copy)]
pub struct FitnessStats {
    pub mean: f64,
    pub min: f64,
    pub max: f64,
    pub std: f64,
}

impl FitnessStats {
    /// Calculate fitness statistics from a population.
    pub fn from_population(pop: &Population) -> Self {
        let fitnesses: Vec<f64> = pop
            .individuals()
            .iter()
            .filter_map(|ind| ind.cached_fitness())
            .map(|f| *f)
            .collect();

        if fitnesses.is_empty() {
            return Self {
                mean: 0.0,
                min: 0.0,
                max: 0.0,
                std: 0.0,
            };
        }

        let sum: f64 = fitnesses.iter().sum();
        let mean = sum / fitnesses.len() as f64;

        let min = fitnesses
            .iter()
            .copied()
            .min_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(0.0);

        let max = fitnesses
            .iter()
            .copied()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(0.0);

        // Calculate standard deviation
        let variance: f64 =
            fitnesses.iter().map(|&f| (f - mean).powi(2)).sum::<f64>() / fitnesses.len() as f64;
        let std = variance.sqrt();

        Self {
            mean,
            min,
            max,
            std,
        }
    }
}
