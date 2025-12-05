//! Mutation operations for sequences.
//!
//! This module simulates realistic DNA mutation processes at the molecular level.
//! It provides two main types of mutations:
//!
//! ## Point Substitutions (Base Replacements)
//! When a DNA replication error or chemical damage causes one nucleotide to be
//! replaced with another (e.g., A → G). These are modeled using substitution models
//! (e.g., JC69, K2P) that specify different rates for different types of changes.
//!
//! ## Insertions and Deletions (Indels)
//! When stretches of DNA are added or removed, typically due to slippage during
//! replication or mechanical errors. These change the sequence length.
//!
//! The module uses stochastic (random) simulation based on evolutionary biology
//! principles. Mutation rates are probabilities per generation, allowing realistic
//! simulation of genetic drift and molecular evolution under specific mutation
//! regimes (e.g., conservation vs. rapidly evolving regions).

use crate::base::{Nucleotide, Sequence};
pub use crate::errors::MutationError;
use rand::Rng;
use rand_distr::{Distribution, Poisson};
use serde::{Deserialize, Serialize};

/// Substitution model for nucleotide mutations.
///
/// This represents a Markov chain model of DNA sequence evolution. It specifies
/// how frequently each nucleotide changes to another over evolutionary time.
///
/// The model uses a symmetric rate matrix (a biological assumption that reflects
/// reversible evolution - the process is reversible in time). For each nucleotide,
/// you can transition to any of the other 3 nucleotides, with biologically
/// motivated rates. For example:
/// - In JC69 (Jukes-Cantor): All substitutions are equally likely
/// - In K2P: Transitions (A↔G, C↔T) are more common than transversions (others)
///
/// The matrix is symmetric (rate A→C equals C→A) and has zero diagonal
/// (no nucleotide can "mutate" to itself).
///
/// Nucleotides are indexed as: A=0, C=1, G=2, T=3
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SubstitutionModel {
    /// Rate for A <-> C substitutions
    rate_a_c: f64,
    /// Rate for A <-> G substitutions
    rate_a_g: f64,
    /// Rate for A <-> T substitutions
    rate_a_t: f64,
    /// Rate for C <-> G substitutions
    rate_c_g: f64,
    /// Rate for C <-> T substitutions
    rate_c_t: f64,
    /// Rate for G <-> T substitutions
    rate_g_t: f64,

    /// Total mutation rate for each base (cached for performance)
    /// Index corresponds to nucleotide index: [0]=A, [1]=C, [2]=G, [3]=T
    total_rates: [f64; 4],
}

impl SubstitutionModel {
    /// Create a new substitution model with the given rate matrix.
    ///
    /// The rate matrix specifies mutation probabilities. Each entry matrix[i][j]
    /// is the per-base, per-generation probability that nucleotide i mutates to j.
    /// Rates should be small (e.g., 0.001 to 0.01) for realistic evolutionary timescales.
    ///
    /// The total mutation rate for each base (sum across all possible targets)
    /// must not exceed 1.0, since probabilities cannot exceed 1. In practice,
    /// total rates are much lower (e.g., 0.01-0.10) except in very fast-evolving
    /// sequences or extreme simulation conditions.
    ///
    /// # Arguments
    /// * `matrix` - 4x4 rate matrix where matrix[i][j] is the rate from nucleotide i to j.
    ///   Rows/columns are ordered as [A, C, G, T] (indices 0-3).
    ///   Matrix must be symmetric (reversible evolution) and have zero diagonal
    ///   (biology: nucleotides don't "mutate" to themselves).
    ///
    /// # Errors
    /// Returns an error if:
    /// - Any rate is negative or greater than 1.0
    /// - Matrix is not symmetric (violates reversibility assumption)
    /// - Diagonal is not zero (biology violation)
    /// - Total mutation rate for any base exceeds 1.0 (probability violation)
    pub fn new(matrix: [[f64; 4]; 4]) -> Result<Self, MutationError> {
        // Validate diagonal is zero
        for (i, row) in matrix.iter().enumerate() {
            if row[i] != 0.0 {
                return Err(MutationError::InvalidMutationRate(row[i]));
            }
        }

        // Validate symmetry and rates are in valid range
        for i in 0..4 {
            for j in (i + 1)..4 {
                check_matrix_values(matrix[i][j], matrix[j][i])?;
            }
        }

        // Extract unique rates (upper triangle)
        // Matrix indices: A=0, C=1, G=2, T=3
        let rate_a_c = matrix[0][1]; // A->C
        let rate_a_g = matrix[0][2]; // A->G
        let rate_a_t = matrix[0][3]; // A->T
        let rate_c_g = matrix[1][2]; // C->G
        let rate_c_t = matrix[1][3]; // C->T
        let rate_g_t = matrix[2][3]; // G->T

        // Compute total mutation rate for each base
        let mut total_rates = [0.0; 4];
        // A (index 0): sum of A->C, A->G, A->T
        total_rates[0] = rate_a_c + rate_a_g + rate_a_t;
        // C (index 1): sum of C->A, C->G, C->T
        total_rates[1] = rate_a_c + rate_c_g + rate_c_t;
        // G (index 2): sum of G->A, G->C, G->T
        total_rates[2] = rate_a_g + rate_c_g + rate_g_t;
        // T (index 3): sum of T->A, T->C, T->G
        total_rates[3] = rate_a_t + rate_c_t + rate_g_t;

        // Validate that total rates don't exceed 1.0
        for &total in total_rates.iter() {
            if total > 1.0 {
                return Err(MutationError::InvalidMutationRate(total));
            }
        }

        Ok(Self {
            rate_a_c,
            rate_a_g,
            rate_a_t,
            rate_c_g,
            rate_c_t,
            rate_g_t,
            total_rates,
        })
    }

    /// Create a JC69 (Jukes-Cantor 1969) substitution model for DNA sequences.
    ///
    /// JC69 is the simplest and most symmetric substitution model. It assumes all
    /// mutations are equally likely - an A is equally likely to become a C, G, or T.
    /// While unrealistic (transitions are typically more common than transversions
    /// in real DNA), JC69 is useful for baseline comparisons and neutral evolution
    /// simulations.
    ///
    /// Given a total mutation rate mu, the model divides it equally:
    /// - Each of the 3 possible substitutions has rate mu/3
    /// - Total rate per base = mu
    /// - Biologically: mu ≈ 0.001-0.01 for typical sequences
    ///
    /// # Arguments
    /// * `mu` - Total mutation rate per base per generation (must be between 0.0 and 1.0).
    ///   Typical values: 0.0 (no evolution), 0.01 (1% per base per generation),
    ///   0.1 (highly mutable region)
    pub fn jc69(mu: f64) -> Result<Self, MutationError> {
        if !(0.0..=1.0).contains(&mu) {
            return Err(MutationError::InvalidMutationRate(mu));
        }
        let rate = mu / 3.0;

        // Create symmetric matrix with equal off-diagonal rates
        // Rows/columns: [A, C, G, T]
        let matrix = [
            [0.0, rate, rate, rate], // A -> [A, C, G, T]
            [rate, 0.0, rate, rate], // C -> [A, C, G, T]
            [rate, rate, 0.0, rate], // G -> [A, C, G, T]
            [rate, rate, rate, 0.0], // T -> [A, C, G, T]
        ];

        Self::new(matrix)
    }

    /// Create a uniform substitution model.
    ///
    /// All off-diagonal transitions occur with equal probability.
    /// This is equivalent to JC69.
    pub fn uniform(mu: f64) -> Result<Self, MutationError> {
        Self::jc69(mu)
    }

    /// Get the rate matrix (upper triangle).
    /// Returns array: [A->C, A->G, A->T, C->G, C->T, G->T]
    #[inline]
    pub fn rates(&self) -> [f64; 6] {
        [
            self.rate_a_c,
            self.rate_a_g,
            self.rate_a_t,
            self.rate_c_g,
            self.rate_c_t,
            self.rate_g_t,
        ]
    }

    /// Get the total mutation rate for a specific base.
    #[inline]
    pub fn total_rate(&self, base: Nucleotide) -> f64 {
        self.total_rates[base.to_index() as usize]
    }

    /// Get the rate for a specific substitution.
    ///
    /// # Arguments
    /// * `from` - Source nucleotide
    /// * `to` - Target nucleotide
    ///
    /// # Returns
    /// The substitution rate, or 0.0 if from == to
    #[inline]
    pub fn rate(&self, from: Nucleotide, to: Nucleotide) -> f64 {
        // Use individual fields for clarity
        match (from, to) {
            // Self-mutations have zero rate
            (Nucleotide::A, Nucleotide::A)
            | (Nucleotide::C, Nucleotide::C)
            | (Nucleotide::G, Nucleotide::G)
            | (Nucleotide::T, Nucleotide::T) => 0.0,

            // Symmetric substitutions
            (Nucleotide::A, Nucleotide::C) | (Nucleotide::C, Nucleotide::A) => self.rate_a_c,
            (Nucleotide::A, Nucleotide::G) | (Nucleotide::G, Nucleotide::A) => self.rate_a_g,
            (Nucleotide::A, Nucleotide::T) | (Nucleotide::T, Nucleotide::A) => self.rate_a_t,
            (Nucleotide::C, Nucleotide::G) | (Nucleotide::G, Nucleotide::C) => self.rate_c_g,
            (Nucleotide::C, Nucleotide::T) | (Nucleotide::T, Nucleotide::C) => self.rate_c_t,
            (Nucleotide::G, Nucleotide::T) | (Nucleotide::T, Nucleotide::G) => self.rate_g_t,
        }
    }

    /// Mutate a single base according to the substitution model.
    ///
    /// This implements the stochastic process: for a given base, we randomly decide
    /// whether a mutation occurs (based on the total mutation rate), and if so,
    /// randomly select which nucleotide it mutates to (weighted by substitution rates).
    ///
    /// Biologically: this represents a single potential mutation event that might
    /// occur during DNA replication or repair.
    ///
    /// # Returns
    /// Either the original base (no mutation) or a new base (if mutation occurred).
    #[inline]
    pub fn mutate_base<R: Rng + ?Sized>(&self, base: Nucleotide, rng: &mut R) -> Nucleotide {
        let base_idx = base.to_index() as usize;
        let total_rate = self.total_rates[base_idx];

        // Check if mutation occurs
        if rng.random::<f64>() >= total_rate {
            return base; // No mutation
        }

        // Mutation occurs - select target base proportional to rates
        let r = rng.random::<f64>() * total_rate;

        // Get the three possible target bases and their rates
        let (targets, rates) = self.get_target_bases_and_rates(base);

        // Select target based on cumulative probabilities
        let mut cumulative = 0.0;
        for i in 0..3 {
            cumulative += rates[i];
            if r < cumulative {
                return targets[i];
            }
        }

        // Fallback (should not reach here due to floating point precision)
        targets[2]
    }

    /// Get the three possible target bases and their rates for a given source base.
    #[inline]
    fn get_target_bases_and_rates(&self, base: Nucleotide) -> ([Nucleotide; 3], [f64; 3]) {
        match base {
            Nucleotide::A => (
                [Nucleotide::C, Nucleotide::G, Nucleotide::T],
                [self.rate_a_c, self.rate_a_g, self.rate_a_t],
            ),
            Nucleotide::C => (
                [Nucleotide::A, Nucleotide::G, Nucleotide::T],
                [self.rate_a_c, self.rate_c_g, self.rate_c_t],
            ),
            Nucleotide::G => (
                [Nucleotide::A, Nucleotide::C, Nucleotide::T],
                [self.rate_a_g, self.rate_c_g, self.rate_g_t],
            ),
            Nucleotide::T => (
                [Nucleotide::A, Nucleotide::C, Nucleotide::G],
                [self.rate_a_t, self.rate_c_t, self.rate_g_t],
            ),
        }
    }

    /// Mutate a sequence in place according to the substitution model.
    ///
    /// This simulates the standard molecular evolution process: go through each base,
    /// decide independently if it mutates (based on the substitution model's rates),
    /// and if so, randomly select its new nucleotide.
    ///
    /// Biologically: this represents one generation of sequence evolution, where
    /// mutations accumulate across the sequence according to the model's parameters.
    /// Expected number of mutations = sequence_length × average_mutation_rate.
    ///
    /// For example, a 1000 bp sequence with mu=0.01 would have ~10 mutations per generation.
    ///
    /// Performance note: this method examines every base (O(N) complexity), making it
    /// efficient for high mutation rates but wasteful when mutation rates are very low
    /// (in which case see mutate_sequence_sparse).
    ///
    /// # Arguments
    /// * `sequence` - The sequence to mutate (modified in place).
    /// * `rng` - Random number generator.
    ///
    /// # Returns
    /// The number of mutations that actually occurred.
    pub fn mutate_sequence<R: Rng + ?Sized>(&self, sequence: &mut Sequence, rng: &mut R) -> usize {
        let len = sequence.len();
        if len == 0 {
            return 0;
        }

        let mut mutation_count = 0;
        let indices = sequence.as_mut_slice();

        // Bulk generate random floats for both mutation decisions and target selection
        let mut random_floats = vec![0.0f64; len * 2];
        rng.fill(&mut random_floats[..]);

        for i in 0..len {
            let base = indices[i];
            let base_idx = base.to_index() as usize;
            let total_rate = self.total_rates[base_idx];

            // Check if mutation occurs using pre-generated random float
            if random_floats[i] < total_rate {
                // Mutation occurs - select target base
                let r = random_floats[len + i] * total_rate;

                let (targets, rates) = self.get_target_bases_and_rates(base);

                // Select target based on cumulative probabilities
                let mut cumulative = 0.0;
                let mut new_base = targets[2]; // Default fallback
                for j in 0..3 {
                    cumulative += rates[j];
                    if r < cumulative {
                        new_base = targets[j];
                        break;
                    }
                }

                indices[i] = new_base;
                mutation_count += 1;
            }
        }

        mutation_count
    }

    /// Mutate a sequence using a sparse sampling approach for better performance.
    ///
    /// This is an optimization for very low mutation rates (< 1%). Instead of
    /// checking every base (most of which don't mutate), we use inverse transform
    /// sampling to skip directly to the next likely mutation site. This is a
    /// form of rejection sampling (also called "thinning").
    ///
    /// The process:
    /// 1. Sample the distance to the next potential mutation using Geometric distribution
    /// 2. Check if that position actually mutates based on its specific rate
    /// 3. Repeat until we pass the end of the sequence
    ///
    /// Why this works:
    /// - With max_rate (the highest mutation rate across all bases), we know a mutation
    ///   event is rare
    /// - We can efficiently find gaps between potential mutations
    /// - When we land on a base, we accept/reject based on individual rates
    /// - The overall probability is mathematically identical to checking every base
    ///
    /// Biologically: same outcome as mutate_sequence, but much faster for sparse mutations.
    /// This is crucial for large genomes with conserved regions.
    ///
    /// # Arguments
    /// * `sequence` - The sequence to mutate (modified in place).
    /// * `rng` - Random number generator.
    ///
    /// # Returns
    /// The number of mutations that occurred.
    ///
    /// # Performance
    /// - Low mutation rates (< 1%): O(K) where K = number of mutations (10-1000× faster)
    /// - Medium mutation rates: O(N) worst case, but still efficient
    /// - High mutation rates (> 10%): Automatically falls back to standard O(N) method
    #[inline]
    pub fn mutate_sequence_sparse<R: Rng + ?Sized>(
        &self,
        sequence: &mut Sequence,
        rng: &mut R,
    ) -> usize {
        let len = sequence.len();
        if len == 0 {
            return 0;
        }

        // Calculate max mutation rate across all bases
        let max_rate = self.total_rates.iter().fold(0.0f64, |a, &b| a.max(b));

        if max_rate <= 0.0 {
            return 0;
        }

        // Special case: if max mutation rate is very high, use standard approach
        // The crossover point is roughly where we expect to visit every base anyway
        if max_rate > 0.1 {
            return self.mutate_sequence(sequence, rng);
        }

        let mut mutation_count = 0;
        let indices = sequence.as_mut_slice();
        let mut current_pos = 0;

        // Geometric distribution: how many non-mutating bases until the next mutation?
        // P(X=k) = (1-p)^k * p, where p = max_rate
        // We use inverse transform sampling: if u ~ Uniform(0,1), then
        // X = floor(log(u) / log(1-p)) gives a geometric random variate
        // This efficiently computes the skip distance without iterating each base
        let log_1_minus_p = (1.0 - max_rate).ln();

        loop {
            // Step 1: Jump to next candidate position
            // Use inverse transform sampling to skip non-mutating bases
            let u: f64 = rng.random();
            let skip = (u.ln() / log_1_minus_p).floor() as usize;
            current_pos += skip;

            // Have we passed the end of the sequence?
            if current_pos >= len {
                break;
            }

            // Step 2: At this candidate position, check if mutation actually occurs
            // (Rejection sampling step)
            let base = indices[current_pos];
            let base_idx = base.to_index() as usize;
            let total_rate = self.total_rates[base_idx];

            // Acceptance: did we land on a position that actually mutates?
            // We sampled using max_rate, but this position's actual rate might be lower.
            // Accept with probability = (actual_rate / max_rate)
            // This ensures the final distribution is identical to checking every base.
            let acceptance_prob = total_rate / max_rate;

            if rng.random_bool(acceptance_prob) {
                // Mutation occurs: select new nucleotide
                indices[current_pos] = self.mutate_base_direct(base, rng);
                mutation_count += 1;
            }

            // Move to next position to continue sampling
            current_pos += 1;
        }

        mutation_count
    }

    /// Mutate a base directly (used internally, assumes mutation should occur)
    #[inline]
    fn mutate_base_direct<R: Rng + ?Sized>(&self, base: Nucleotide, rng: &mut R) -> Nucleotide {
        let base_idx = base.to_index() as usize;
        let total_rate = self.total_rates[base_idx];

        // Select target base proportional to rates
        let r = rng.random::<f64>() * total_rate;

        let (targets, rates) = self.get_target_bases_and_rates(base);

        let mut cumulative = 0.0;
        for i in 0..3 {
            cumulative += rates[i];
            if r < cumulative {
                return targets[i];
            }
        }

        targets[2]
    }
}

/// Model for insertion and deletion (indel) mutations.
///
/// Indels are structural mutations that change sequence length. They arise from:
/// - Replication slippage: DNA polymerase stutters on repetitive sequences
/// - Mechanical damage: breaks in the DNA strand
/// - Repair errors: incorrect rejoining of DNA breaks
///
/// This model is separate from point mutations (substitutions) because:
/// 1. They have different molecular causes
/// 2. They're often context-dependent (more common in homopolymer runs)
/// 3. They have different fitness effects (often deleterious, especially in coding regions)
///
/// The model uses independent rates for insertions and deletions, and models
/// indel lengths using a Geometric distribution (which reflects their observed
/// biology: small indels are more common than large ones).
///
/// # Indel Length Distribution
/// Indel lengths follow a Geometric distribution with parameter p. This gives
/// a long tail (occasional large indels) but concentrates probability on small sizes.
/// - Mean indel length = 1/p
/// - Example: p=0.5 means average indels are ~2bp long
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IndelModel {
    /// Rate of insertion per base per generation
    insertion_rate: f64,
    /// Rate of deletion per base per generation
    deletion_rate: f64,
    /// Parameter p for Geometric distribution of indel lengths (0.0 < p <= 1.0).
    /// P(X=k) = (1-p)^(k-1) * p
    /// Mean length = 1/p
    length_p: f64,
}

impl IndelModel {
    /// Create a new indel model.
    ///
    /// # Arguments
    /// * `insertion_rate` - Rate of insertion per base per generation (0.0 to 1.0)
    /// * `deletion_rate` - Rate of deletion per base per generation (0.0 to 1.0)
    /// * `length_p` - Parameter for Geometric distribution of lengths (0.0 < p <= 1.0).
    ///   Higher p means shorter indels. Mean length is 1/p.
    ///
    /// # Errors
    /// Returns an error if rates are invalid or length_p is out of range.
    pub fn new(
        insertion_rate: f64,
        deletion_rate: f64,
        length_p: f64,
    ) -> Result<Self, MutationError> {
        if !(0.0..=1.0).contains(&insertion_rate) {
            return Err(MutationError::InvalidMutationRate(insertion_rate));
        }
        if !(0.0..=1.0).contains(&deletion_rate) {
            return Err(MutationError::InvalidMutationRate(deletion_rate));
        }
        if !(0.0..=1.0).contains(&length_p) || length_p == 0.0 {
            return Err(MutationError::InvalidMutationRate(length_p));
        }

        Ok(Self {
            insertion_rate,
            deletion_rate,
            length_p,
        })
    }

    /// Apply indels to a sequence.
    ///
    /// This simulates one generation of indel mutations. For each base, we determine
    /// how many insertion and deletion events occur (using Poisson distribution,
    /// which models rare independent events). Then we apply them, being careful to
    /// maintain valid sequence indices.
    ///
    /// Algorithm:
    /// 1. Sample number of insertion events from Poisson(length × insertion_rate)
    /// 2. Sample number of deletion events from Poisson(length × deletion_rate)
    /// 3. For each event, randomly select position and indel length (from Geometric)
    /// 4. Sort events by position (descending) and apply to avoid index corruption
    /// 5. Handle edge cases (position beyond sequence length, etc.)
    ///
    /// Biologically: this represents indel mutability in a genomic region during
    /// DNA replication. The Poisson model is appropriate for independent, rare events.
    ///
    /// # Returns
    /// The number of indel events (insertions + deletions) applied.
    pub fn apply_indels<R: Rng + ?Sized>(&self, sequence: &mut Sequence, rng: &mut R) -> usize {
        let len = sequence.len();
        if len == 0 {
            return 0;
        }

        let mut events = 0;

        // Calculate expected numbers of insertion and deletion events
        // For a sequence of length L with rate r, we expect L×r events
        // The actual count is stochastic (follows a Poisson distribution)
        // This reflects the molecular reality: some bases mutate, some don't
        let expected_insertions = len as f64 * self.insertion_rate;
        let expected_deletions = len as f64 * self.deletion_rate;

        let num_insertions = if expected_insertions > 0.0 {
            match Poisson::new(expected_insertions) {
                Ok(p) => p.sample(rng) as usize,
                Err(_) => 0,
            }
        } else {
            0
        };

        let num_deletions = if expected_deletions > 0.0 {
            match Poisson::new(expected_deletions) {
                Ok(p) => p.sample(rng) as usize,
                Err(_) => 0,
            }
        } else {
            0
        };

        if num_insertions == 0 && num_deletions == 0 {
            return 0;
        }

        // Generate all events first, then apply them in reverse position order
        // Why reverse order? If we delete positions from low to high, we'd corrupt indices.
        // By going high to low, we avoid index shifting issues during application.
        // For insertions at the same position, we apply them before deletions (arbitrary choice).

        #[derive(Debug)]
        enum IndelEvent {
            Insertion { pos: usize, len: usize },
            Deletion { pos: usize, len: usize },
        }

        let mut event_list = Vec::with_capacity(num_insertions + num_deletions);
        let geo = rand_distr::Geometric::new(self.length_p).unwrap(); // p validated in new()

        // Generate insertion events
        for _ in 0..num_insertions {
            // Position: insertions can happen before position 0 or after the last base
            let pos = rng.random_range(0..=len);
            // Length: sample from Geometric distribution (add 1 because Geometric starts at 0)
            // This gives realistic indel length distribution (small indels most common)
            let indel_len = (geo.sample(rng) + 1) as usize;
            event_list.push(IndelEvent::Insertion { pos, len: indel_len });
        }

        // Generate deletion events
        for _ in 0..num_deletions {
            if len == 0 {
                break; // Can't delete from empty sequence
            }
            // Position: can only delete from existing bases (0 to len-1)
            let pos = rng.random_range(0..len);
            // Length: again from Geometric (realistic distribution)
            let indel_len = (geo.sample(rng) + 1) as usize;
            event_list.push(IndelEvent::Deletion { pos, len: indel_len });
        }

        // Sort events by position descending (highest position first)
        // This prevents index shifting issues when modifying the sequence
        // For same position, deletions are processed before insertions (prevents index conflicts)
        event_list.sort_by(|a, b| {
            let pos_a = match a {
                IndelEvent::Insertion { pos, .. } => *pos,
                IndelEvent::Deletion { pos, .. } => *pos,
            };
            let pos_b = match b {
                IndelEvent::Insertion { pos, .. } => *pos,
                IndelEvent::Deletion { pos, .. } => *pos,
            };
            pos_b.cmp(&pos_a)
        });

        // Apply events in the order they were sorted (descending position)
        for event in event_list {
            match event {
                IndelEvent::Insertion { pos, len } => {
                    // Insert 'len' random nucleotides starting at position 'pos'
                    // Generate random bases - biologically, inserted sequences are random
                    for _ in 0..len {
                        let base = Nucleotide::from_index(rng.random_range(0..4)).unwrap();
                        // Clamp position in case previous deletions shifted the sequence
                        let current_len = sequence.len();
                        let insert_pos = pos.min(current_len);
                        sequence.insert(insert_pos, base);
                    }
                    events += 1;
                }
                IndelEvent::Deletion { pos, len } => {
                    // Delete 'len' bases starting at position 'pos'
                    let current_len = sequence.len();
                    if current_len == 0 {
                        continue; // Sequence was already emptied
                    }

                    // Clamp position (can't delete past the sequence)
                    let start_pos = pos.min(current_len - 1);
                    // Don't delete more bases than exist after this position
                    let delete_len = len.min(current_len - start_pos);

                    // Remove bases one at a time from the same position
                    // (each remove shifts subsequent bases left)
                    for _ in 0..delete_len {
                        sequence.remove(start_pos);
                    }
                    events += 1;
                }
            }
        }

        events
    }
}

/// Sample k positions from [0, n) without replacement using reservoir sampling.
///
/// This is an efficient algorithm for sampling without replacement when k << n.
/// Time complexity: O(k), Space complexity: O(k)
///
/// # Arguments
/// * `n` - The range to sample from [0, n)
/// * `k` - The number of samples to take (must be <= n)
/// * `rng` - Random number generator
///
/// # Returns
/// A vector of k unique positions in [0, n)
#[inline]
fn check_matrix_values(rate_ij: f64, rate_ji: f64) -> Result<(), MutationError> {
    // Check symmetry
    if (rate_ij - rate_ji).abs() > 1e-10 {
        return Err(MutationError::InvalidMutationRate(rate_ij));
    }

    // Check valid range
    if !(0.0..=1.0).contains(&rate_ij) {
        return Err(MutationError::InvalidMutationRate(rate_ij));
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;
    use std::str::FromStr;

    #[test]
    fn test_substitution_model_new() {
        let rate = 0.01;
        let matrix = [
            [0.0, rate, rate, rate],
            [rate, 0.0, rate, rate],
            [rate, rate, 0.0, rate],
            [rate, rate, rate, 0.0],
        ];
        let model = SubstitutionModel::new(matrix).unwrap();
        let expected_rates = [rate, rate, rate, rate, rate, rate];
        assert_eq!(model.rates(), expected_rates);
    }

    #[test]
    fn test_substitution_model_invalid_rate() {
        // Negative rate
        let invalid_matrix = [
            [0.0, -0.1, 0.01, 0.01],
            [-0.1, 0.0, 0.01, 0.01],
            [0.01, 0.01, 0.0, 0.01],
            [0.01, 0.01, 0.01, 0.0],
        ];
        assert!(SubstitutionModel::new(invalid_matrix).is_err());

        // Rate > 1.0
        let invalid_matrix2 = [
            [0.0, 1.5, 0.01, 0.01],
            [1.5, 0.0, 0.01, 0.01],
            [0.01, 0.01, 0.0, 0.01],
            [0.01, 0.01, 0.01, 0.0],
        ];
        assert!(SubstitutionModel::new(invalid_matrix2).is_err());

        // Total rate exceeds 1.0 for a base
        let invalid_matrix3 = [
            [0.0, 0.4, 0.4, 0.4],
            [0.4, 0.0, 0.01, 0.01],
            [0.4, 0.01, 0.0, 0.01],
            [0.4, 0.01, 0.01, 0.0],
        ];
        assert!(SubstitutionModel::new(invalid_matrix3).is_err());

        // Non-zero diagonal
        let invalid_matrix4 = [
            [0.5, 0.01, 0.01, 0.01],
            [0.01, 0.0, 0.01, 0.01],
            [0.01, 0.01, 0.0, 0.01],
            [0.01, 0.01, 0.01, 0.0],
        ];
        assert!(SubstitutionModel::new(invalid_matrix4).is_err());

        // Non-symmetric matrix
        let invalid_matrix5 = [
            [0.0, 0.01, 0.01, 0.01],
            [0.02, 0.0, 0.01, 0.01], // Different from [0][1]
            [0.01, 0.01, 0.0, 0.01],
            [0.01, 0.01, 0.01, 0.0],
        ];
        assert!(SubstitutionModel::new(invalid_matrix5).is_err());
    }

    #[test]
    fn test_substitution_model_jc69() {
        let model = SubstitutionModel::jc69(0.03).unwrap();
        // JC69 with mu=0.03 gives rate=0.01 for each substitution
        let expected = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01];
        assert_eq!(model.rates(), expected);

        // Total rate for each base should be 0.03
        assert!((model.total_rate(Nucleotide::A) - 0.03).abs() < 1e-10);
        assert!((model.total_rate(Nucleotide::C) - 0.03).abs() < 1e-10);
        assert!((model.total_rate(Nucleotide::G) - 0.03).abs() < 1e-10);
        assert!((model.total_rate(Nucleotide::T) - 0.03).abs() < 1e-10);
    }

    #[test]
    fn test_substitution_model_zero_rate() {
        let model = SubstitutionModel::jc69(0.0).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

        // With zero mutation rate, base should never change
        for _ in 0..100 {
            let base = Nucleotide::A;
            let mutated = model.mutate_base(base, &mut rng);
            assert_eq!(mutated, base);
        }
    }

    #[test]
    fn test_substitution_model_mutate_base() {
        let model = SubstitutionModel::jc69(1.0).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

        // With mutation rate 1.0, base should always change
        let mut mutations = 0;
        for _ in 0..100 {
            let base = Nucleotide::A;
            let mutated = model.mutate_base(base, &mut rng);
            if mutated != base {
                mutations += 1;
            }
        }

        // Should have many mutations (close to 100)
        assert!(mutations > 90);
    }

    #[test]
    fn test_substitution_model_mutate_base_distribution() {
        // Use JC69 with rate 1.0 (mu=1.0 means total rate 1.0, each transition 1/3)
        let model = SubstitutionModel::jc69(1.0).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

        let mut counts = [0; 4]; // Indices 0-3 for A,C,G,T
        for _ in 0..1000 {
            let mutated = model.mutate_base(Nucleotide::A, &mut rng);
            counts[mutated.to_index() as usize] += 1;
        }

        // A (index 0) should never appear (always mutates)
        assert_eq!(counts[0], 0);

        // C, G, T (indices 1-3) should be roughly equally distributed
        for &count in &counts[1..4] {
            assert!(count > 250);
            assert!(count < 450);
        }
    }

    #[test]
    fn test_mutate_sequence_zero_rate() {
        let model = SubstitutionModel::jc69(0.0).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

        let mut seq = Sequence::from_str("ACGTACGT").unwrap();
        let original = seq.to_string();

        let count = model.mutate_sequence(&mut seq, &mut rng);

        assert_eq!(count, 0);
        assert_eq!(seq.to_string(), original);
    }

    #[test]
    fn test_mutate_sequence_low_rate() {
        let model = SubstitutionModel::jc69(0.1).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

        let mut seq = Sequence::from_str("ACGTACGTACGTACGT").unwrap();
        let count = model.mutate_sequence(&mut seq, &mut rng);

        // With low rate, should have few mutations
        assert!(count < 5);
    }

    #[test]
    fn test_mutate_sequence_high_rate() {
        let model = SubstitutionModel::jc69(0.9).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

        let mut seq = Sequence::from_str("ACGTACGTACGTACGT").unwrap();
        let count = model.mutate_sequence(&mut seq, &mut rng);

        // With high rate, should have many mutations
        assert!(count > 10);
    }

    #[test]
    fn test_mutate_sequence_empty() {
        let model = SubstitutionModel::jc69(0.5).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

        let mut seq = Sequence::new();
        let count = model.mutate_sequence(&mut seq, &mut rng);

        assert_eq!(count, 0);
        assert!(seq.is_empty());
    }

    #[test]
    fn test_mutate_sequence_deterministic() {
        let model = SubstitutionModel::jc69(0.1).unwrap();

        let mut seq1 = Sequence::from_str("ACGTACGTACGTACGT").unwrap();
        let mut seq2 = Sequence::from_str("ACGTACGTACGTACGT").unwrap();

        let mut rng1 = Xoshiro256PlusPlus::seed_from_u64(123);
        let mut rng2 = Xoshiro256PlusPlus::seed_from_u64(123);

        let count1 = model.mutate_sequence(&mut seq1, &mut rng1);
        let count2 = model.mutate_sequence(&mut seq2, &mut rng2);

        // Same seed should produce same results
        assert_eq!(count1, count2);
        assert_eq!(seq1.to_string(), seq2.to_string());
    }

    #[test]
    fn test_mutation_error_display() {
        let err = MutationError::InvalidMutationRate(1.5);
        let msg = format!("{err}");
        assert!(msg.contains("Invalid mutation rate"));
        assert!(msg.contains("1.5"));
    }

    #[test]
    fn test_substitution_model_clone() {
        let model1 = SubstitutionModel::jc69(0.03).unwrap();
        let model2 = model1.clone();

        assert_eq!(model1.rates(), model2.rates());
        assert_eq!(
            model1.total_rate(Nucleotide::A),
            model2.total_rate(Nucleotide::A)
        );
    }

    #[test]
    fn test_substitution_model_rate_query() {
        let model = SubstitutionModel::jc69(0.03).unwrap();

        // All non-diagonal rates should be 0.01
        assert!((model.rate(Nucleotide::A, Nucleotide::C) - 0.01).abs() < 1e-10);
        assert!((model.rate(Nucleotide::A, Nucleotide::G) - 0.01).abs() < 1e-10);
        assert!((model.rate(Nucleotide::C, Nucleotide::T) - 0.01).abs() < 1e-10);

        // Diagonal rates should be 0
        assert_eq!(model.rate(Nucleotide::A, Nucleotide::A), 0.0);
        assert_eq!(model.rate(Nucleotide::C, Nucleotide::C), 0.0);

        // Symmetric rates
        assert_eq!(
            model.rate(Nucleotide::A, Nucleotide::C),
            model.rate(Nucleotide::C, Nucleotide::A)
        );
        assert_eq!(
            model.rate(Nucleotide::G, Nucleotide::T),
            model.rate(Nucleotide::T, Nucleotide::G)
        );
    }

    #[test]
    fn test_substitution_model_asymmetric() {
        // Create a transition/transversion model where transitions are more common
        // Transitions: A<->G, C<->T (purines <-> purines, pyrimidines <-> pyrimidines)
        // Transversions: all others
        // Matrix rows/columns: [A, C, G, T]
        let matrix = [
            [0.0, 0.005, 0.02, 0.005], // A -> [A, C, G, T]
            [0.005, 0.0, 0.005, 0.02], // C -> [A, C, G, T]
            [0.02, 0.005, 0.0, 0.005], // G -> [A, C, G, T]
            [0.005, 0.02, 0.005, 0.0], // T -> [A, C, G, T]
        ];

        let model = SubstitutionModel::new(matrix).unwrap();

        // Check transitions are higher than transversions
        assert!(
            model.rate(Nucleotide::A, Nucleotide::G) > model.rate(Nucleotide::A, Nucleotide::C)
        );
        assert!(
            model.rate(Nucleotide::C, Nucleotide::T) > model.rate(Nucleotide::C, Nucleotide::G)
        );

        // Check symmetry
        assert_eq!(
            model.rate(Nucleotide::A, Nucleotide::G),
            model.rate(Nucleotide::G, Nucleotide::A)
        );
    }

    // Tests for Sparse-based mutation

    #[test]
    fn test_mutate_sequence_sparse_zero_rate() {
        let model = SubstitutionModel::jc69(0.0).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

        let mut seq = Sequence::from_str("ACGTACGT").unwrap();
        let original = seq.to_string();

        let count = model.mutate_sequence_sparse(&mut seq, &mut rng);

        assert_eq!(count, 0);
        assert_eq!(seq.to_string(), original);
    }

    #[test]
    fn test_mutate_sequence_sparse_low_rate() {
        let model = SubstitutionModel::jc69(0.001).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

        let mut seq = Sequence::from_str("ACGT".repeat(250).as_str()).unwrap();
        let count = model.mutate_sequence_sparse(&mut seq, &mut rng);

        // With low rate on 1000bp, expect 0-5 mutations
        assert!(count < 10);
    }

    #[test]
    fn test_mutate_sequence_sparse_medium_rate() {
        let model = SubstitutionModel::jc69(0.01).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

        let mut seq = Sequence::from_str("ACGT".repeat(250).as_str()).unwrap();
        let count = model.mutate_sequence_sparse(&mut seq, &mut rng);

        // With medium rate on 1000bp, expect 5-20 mutations
        assert!(count > 0);
        assert!(count < 30);
    }

    #[test]
    fn test_mutate_sequence_sparse_high_rate_fallback() {
        // High rate should fallback to standard method
        let model = SubstitutionModel::jc69(0.6).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

        let mut seq = Sequence::from_str("ACGTACGTACGTACGT").unwrap();
        let count = model.mutate_sequence_sparse(&mut seq, &mut rng);

        // Should have many mutations
        assert!(count > 5);
    }

    #[test]
    fn test_mutate_sequence_sparse_empty() {
        let model = SubstitutionModel::jc69(0.01).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

        let mut seq = Sequence::new();
        let count = model.mutate_sequence_sparse(&mut seq, &mut rng);

        assert_eq!(count, 0);
        assert!(seq.is_empty());
    }

    #[test]
    fn test_mutate_sequence_sparse_deterministic() {
        let model = SubstitutionModel::jc69(0.01).unwrap();

        let mut seq1 = Sequence::from_str("ACGT".repeat(250).as_str()).unwrap();
        let mut seq2 = Sequence::from_str("ACGT".repeat(250).as_str()).unwrap();

        let mut rng1 = Xoshiro256PlusPlus::seed_from_u64(123);
        let mut rng2 = Xoshiro256PlusPlus::seed_from_u64(123);

        let count1 = model.mutate_sequence_sparse(&mut seq1, &mut rng1);
        let count2 = model.mutate_sequence_sparse(&mut seq2, &mut rng2);

        // Same seed should produce same results
        assert_eq!(count1, count2);
        assert_eq!(seq1.to_string(), seq2.to_string());
    }

    #[test]
    fn test_mutate_sequence_sparse_statistical_equivalence() {
        // Test that Sparse method produces statistically similar results to standard method
        let model = SubstitutionModel::jc69(0.01).unwrap();
        let seq_str = "ACGT".repeat(250); // 1000bp

        let mut standard_mutations = Vec::new();
        let mut sparse_mutations = Vec::new();

        // Run multiple trials
        for seed in 0..100 {
            let mut seq_standard = Sequence::from_str(&seq_str).unwrap();
            let mut seq_sparse = Sequence::from_str(&seq_str).unwrap();

            let mut rng_standard = Xoshiro256PlusPlus::seed_from_u64(seed);
            let mut rng_sparse = Xoshiro256PlusPlus::seed_from_u64(seed + 10000); // Different seed

            standard_mutations.push(model.mutate_sequence(&mut seq_standard, &mut rng_standard));
            sparse_mutations.push(model.mutate_sequence_sparse(&mut seq_sparse, &mut rng_sparse));
        }

        // Calculate means
        let mean_standard: f64 = standard_mutations.iter().sum::<usize>() as f64 / 100.0;
        let mean_sparse: f64 = sparse_mutations.iter().sum::<usize>() as f64 / 100.0;

        // Expected: 1000 * 0.01 = 10 mutations
        // Both methods should have similar means (within 20% of each other)
        let diff_ratio = (mean_standard - mean_sparse).abs() / mean_standard.max(mean_sparse);
        assert!(
            diff_ratio < 0.2,
            "Mean difference too large: standard={mean_standard}, sparse={mean_sparse}",
        );
    }

    #[test]
    fn test_indel_model_new() {
        let model = IndelModel::new(0.1, 0.1, 0.5);
        assert!(model.is_ok());
    }

    #[test]
    fn test_indel_model_invalid() {
        assert!(IndelModel::new(-0.1, 0.1, 0.5).is_err());
        assert!(IndelModel::new(0.1, 1.5, 0.5).is_err());
        assert!(IndelModel::new(0.1, 0.1, 0.0).is_err()); // p must be > 0
        assert!(IndelModel::new(0.1, 0.1, 1.1).is_err());
    }

    #[test]
    fn test_indel_model_apply_no_events() {
        let model = IndelModel::new(0.0, 0.0, 0.5).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let mut seq = Sequence::from_str("ACGT").unwrap();
        let original = seq.to_string();

        let events = model.apply_indels(&mut seq, &mut rng);
        assert_eq!(events, 0);
        assert_eq!(seq.to_string(), original);
    }

    #[test]
    fn test_indel_model_insertions() {
        // High insertion rate
        let model = IndelModel::new(1.0, 0.0, 0.5).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let mut seq = Sequence::from_str("AAAA").unwrap();

        let events = model.apply_indels(&mut seq, &mut rng);

        assert!(events > 0);
        assert!(seq.len() > 4);
    }

    #[test]
    fn test_indel_model_deletions() {
        // High deletion rate
        let model = IndelModel::new(0.0, 0.5, 0.5).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let mut seq = Sequence::from_str("ACGTACGTACGT").unwrap(); // Length 12

        let events = model.apply_indels(&mut seq, &mut rng);

        assert!(events > 0);
        assert!(seq.len() < 12);
    }

    #[test]
    fn test_indel_model_mixed() {
        let model = IndelModel::new(0.1, 0.1, 0.5).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let mut seq = Sequence::from_str("ACGT".repeat(100).as_str()).unwrap();

        let events = model.apply_indels(&mut seq, &mut rng);

        // Should have some events
        assert!(events > 0);
    }
}
