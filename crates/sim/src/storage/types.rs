use crate::base::FitnessValue;
use crate::genome::Individual;
use crate::simulation::Population;
use centrevo_codec::CodecStrategy;
use serde::{Deserialize, Serialize};

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
    pub haplotype1_seq: Vec<u8>,         // Encoded
    pub haplotype1_map: Option<Vec<u8>>, // Encoded serialized map
    pub haplotype2_seq: Vec<u8>,         // Encoded
    pub haplotype2_map: Option<Vec<u8>>, // Encoded serialized map
    pub fitness: Option<f64>,
}

impl IndividualSnapshot {
    /// Create a snapshot from an individual.
    pub fn from_individual(ind: &Individual, codec: &CodecStrategy) -> Self {
        let h1 = ind.haplotype1();
        let h2 = ind.haplotype2();

        // Helper to process chromosome (assumes index 0 for now as per current limitation)
        let process_chr =
            |chr_opt: Option<&crate::genome::Chromosome>| -> (Vec<u8>, Option<Vec<u8>>) {
                if let Some(chr) = chr_opt {
                    // 1. Sequence -> Configured Codec
                    let raw_seq: Vec<u8> = chr
                        .sequence()
                        .as_slice()
                        .iter()
                        .map(|n| n.to_index())
                        .collect();

                    let encoded_seq = codec.encode(&raw_seq).expect("Failed to encode sequence");

                    // 2. Map -> Serialize -> UnpackedRS (Fixed)
                    let map_bytes = bincode::serialize(chr.map()).unwrap_or_default();
                    let encoded_map = CodecStrategy::UnpackedRS.encode(&map_bytes).ok();

                    (encoded_seq, encoded_map)
                } else {
                    (Vec::new(), None)
                }
            };

        let (h1_seq, h1_map) = process_chr(h1.get(0));
        let (h2_seq, h2_map) = process_chr(h2.get(0));

        Self {
            haplotype1_seq: h1_seq,
            haplotype1_map: h1_map,
            haplotype2_seq: h2_seq,
            haplotype2_map: h2_map,
            fitness: ind.cached_fitness().map(|f| f64::from(*f)),
        }
    }

    /// Reconstruct an Individual from a snapshot.
    pub fn to_individual(&self, id: &str, codec: &CodecStrategy) -> Result<Individual, String> {
        use crate::base::{Nucleotide, Sequence};
        use crate::genome::repeat_map::RepeatMap;
        use crate::genome::{Chromosome, Haplotype};

        // Helper to reconstruct chromosome
        let reconstruct_chr = |chr_suffix: &str,
                               seq_data: &[u8],
                               map_data: Option<&[u8]>|
         -> Result<Chromosome, String> {
            // 1. Decode Sequence
            if seq_data.is_empty() {
                // Handle empty chromosome case if needed, though usually means logic error if expected
                // For now, assume if check handles it upstream or returns empty chromosome
                return Ok(Chromosome::new(
                    format!("{}_{}", id, chr_suffix),
                    Sequence::new(),
                    RepeatMap::new(vec![0], vec![0])
                        .unwrap_or_else(|_| RepeatMap::uniform(0, 0, 0)),
                ));
            }

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
                RepeatMap::new(vec![0], vec![0]).unwrap_or_else(|_| RepeatMap::uniform(0, 0, 0))
            };

            Ok(Chromosome::new(format!("{}_{}", id, chr_suffix), seq, map))
        };

        let chr1 = reconstruct_chr(
            "h1_chr1",
            &self.haplotype1_seq,
            self.haplotype1_map.as_deref(),
        )?;
        let chr2 = reconstruct_chr(
            "h2_chr1",
            &self.haplotype2_seq,
            self.haplotype2_map.as_deref(),
        )?;

        let mut hap1 = Haplotype::new();
        hap1.push(chr1);
        let mut hap2 = Haplotype::new();
        hap2.push(chr2);

        let mut individual = Individual::new(id, hap1, hap2);
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
