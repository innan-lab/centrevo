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
    pub fn from_individual(
        ind: &Individual,
        codec: &CodecStrategy,
        arena: &crate::base::GenomeArena,
    ) -> Self {
        let h1 = ind.haplotype1();
        let h2 = ind.haplotype2();

        // Helper to process chromosome (assumes index 0 for now as per current limitation)
        let process_chr =
            |chr_opt: Option<&crate::genome::Chromosome>| -> (Vec<u8>, Option<Vec<u8>>) {
                if let Some(chr) = chr_opt {
                    // 1. Sequence -> Configured Codec
                    let raw_seq: Vec<u8> =
                        chr.sequence(arena).iter().map(|n| n.to_index()).collect();

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
            fitness: ind.cached_fitness().map(|f| f.ln().into()),
        }
    }

    /// Convert to Individual (restoring full state).
    ///
    /// The sequences are decoded and allocated in the provided `GenomeArena`.
    pub fn to_individual(
        &self,
        id: &str,
        codec: &CodecStrategy,
        arena: &mut crate::base::GenomeArena,
    ) -> Result<Individual, String> {
        use crate::base::FitnessValue;
        use crate::genome::{Chromosome, Haplotype};

        let hap1_seq = codec
            .decode(&self.haplotype1_seq)
            .map_err(|e| format!("Failed to decode hap1: {e}"))?;
        let hap2_seq = codec
            .decode(&self.haplotype2_seq)
            .map_err(|e| format!("Failed to decode hap2: {e}"))?;

        // 3. Reconstruct Chromosomes (using decoded sequences)
        // Note: Codec output is Vec<u8> (indices). Convert to Vec<Nucleotide> first.
        let hap1_seq: crate::base::Sequence = crate::base::Sequence::from_indices(hap1_seq)
            .map_err(|e| format!("Invalid nucleotides in hap1: {e}"))?;
        let hap2_seq: crate::base::Sequence = crate::base::Sequence::from_indices(hap2_seq)
            .map_err(|e| format!("Invalid nucleotides in hap2: {e}"))?;

        // Recover fitness from log-space (if needed, but Individual constructor doesn't take fitness)
        // Fitness is cached in the Individual.
        // We set it after creation.

        // Assuming single chromosome per haplotype for now (as stored in DB)
        // We allocate the sequences in the arena
        let slice1 = arena.alloc(hap1_seq.as_slice());
        let slice2 = arena.alloc(hap2_seq.as_slice());

        // Deserialize RepeatMap from bytes if present, otherwise assume simple default or error?
        // Note: RepeatMap::uniform... but here we want to restore from DB.
        // Assuming the Option<Vec<u8>> contains bincode encoded RepeatMap.
        let map1 = if let Some(bytes) = &self.haplotype1_map {
            let decoded = CodecStrategy::UnpackedRS
                .decode(bytes)
                .map_err(|e| format!("Failed to decode hap1 map: {e}"))?;
            bincode::deserialize::<crate::genome::RepeatMap>(&decoded)
                .map_err(|e| format!("Failed to deserialize hap1 map: {e}"))?
        } else {
            // Fallback or error? For now fallback to simple default until map storage is strictly required
            // Or maybe we can't create valid Chromosome without it?
            // Assuming default uniform for legacy data?
            crate::genome::RepeatMap::uniform(10, 5, 2)
        };

        // Similar for hap2
        let map2 = if let Some(bytes) = &self.haplotype2_map {
            let decoded = CodecStrategy::UnpackedRS
                .decode(bytes)
                .map_err(|e| format!("Failed to decode hap2 map: {e}"))?;
            bincode::deserialize::<crate::genome::RepeatMap>(&decoded)
                .map_err(|e| format!("Failed to deserialize hap2 map: {e}"))?
        } else {
            crate::genome::RepeatMap::uniform(10, 5, 2)
        };

        let chr1 = Chromosome::new(format!("{id}_h1_c0"), slice1, map1);

        let chr2 = Chromosome::new(format!("{id}_h2_c0"), slice2, map2);

        let mut hap1 = Haplotype::new();
        hap1.push(chr1);

        let mut hap2 = Haplotype::new();
        hap2.push(chr2);

        let mut ind = Individual::new(id, hap1, hap2);

        // Restore fitness
        // fitness is recorded in log_fitness
        // 10^log_fitness = fitness_value
        // If fitness was LETHAL (0), log is -inf.
        // We handle this carefully.
        if let Some(log_f) = self.fitness {
            // Assuming 'fitness' field now stores log10 fitness
            if log_f.is_infinite() && log_f.is_sign_negative() {
                ind.set_cached_fitness(FitnessValue::LETHAL_FITNESS);
            } else {
                // Convert back from log10
                let val = log_f.exp();
                ind.set_cached_fitness(FitnessValue::new(val));
            }
        }

        Ok(ind)
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
    /// Returns stats in log-scale for consistency with stored fitness.
    pub fn from_population(pop: &Population) -> Self {
        let log_fitnesses: Vec<f64> = pop
            .individuals()
            .iter()
            .filter_map(|ind| ind.cached_fitness())
            .map(|f| f.ln().into())
            .collect();

        if log_fitnesses.is_empty() {
            return Self {
                mean: 0.0,
                min: 0.0,
                max: 0.0,
                std: 0.0,
            };
        }

        let sum: f64 = log_fitnesses.iter().sum();
        let mean = sum / log_fitnesses.len() as f64;

        let min = log_fitnesses
            .iter()
            .copied()
            .min_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(0.0);

        let max = log_fitnesses
            .iter()
            .copied()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(0.0);

        // Calculate standard deviation in log-space
        let variance: f64 = log_fitnesses
            .iter()
            .map(|&f| (f - mean).powi(2))
            .sum::<f64>()
            / log_fitnesses.len() as f64;
        let std = variance.sqrt();

        Self {
            mean,
            min,
            max,
            std,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::{FitnessValue, Nucleotide};
    use crate::genome::{Chromosome, Haplotype, Individual};

    #[test]
    fn test_individual_snapshot_log_fitness_roundtrip() {
        let mut arena = crate::base::GenomeArena::new();
        let h1 = Haplotype::from_chromosomes(vec![Chromosome::uniform(
            "c1",
            Nucleotide::A,
            10,
            1,
            1,
            &mut arena,
        )]);
        let h2 = Haplotype::from_chromosomes(vec![Chromosome::uniform(
            "c1",
            Nucleotide::T,
            10,
            1,
            1,
            &mut arena,
        )]);
        let mut ind = Individual::new("test", h1, h2);

        // Linear fitness 0.5 -> Log fitness ln(0.5)
        let linear_fitness = FitnessValue::new(0.5);
        ind.set_cached_fitness(linear_fitness);

        let codec = CodecStrategy::default();
        let snapshot = IndividualSnapshot::from_individual(&ind, &codec, &arena);

        // Verify stored fitness is in log-space
        let expected_log_fitness: f64 = linear_fitness.ln().into();
        assert!((snapshot.fitness.unwrap() - expected_log_fitness).abs() < 1e-12);

        // Verify round-trip back to linear-space
        let recovered_ind = snapshot.to_individual("test", &codec, &mut arena).unwrap();
        let recovered_fitness: f64 = recovered_ind.cached_fitness().unwrap().into();
        assert!((recovered_fitness - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_fitness_stats_log_space() {
        use std::f64::consts::LN_2;
        let h1 = Haplotype::new();
        let h2 = Haplotype::new();

        let mut ind1 = Individual::new("1", h1.clone(), h2.clone());
        ind1.set_cached_fitness(FitnessValue::new(0.5)); // ln(0.5) = -LN_2

        let mut ind2 = Individual::new("2", h1, h2);
        ind2.set_cached_fitness(FitnessValue::new(0.25)); // ln(0.25) = -2 * LN_2

        let pop = Population::new("test", vec![ind1, ind2]);
        let stats = FitnessStats::from_population(&pop);

        let log_f1 = -LN_2;
        let log_f2 = -2.0 * LN_2;

        assert!((stats.min - log_f2).abs() < 1e-12);
        assert!((stats.max - log_f1).abs() < 1e-12);
        assert!((stats.mean - (log_f1 + log_f2) / 2.0).abs() < 1e-12);
    }
}
