use super::repeat_map::RepeatMap;
use crate::base::{Nucleotide, Sequence, SharedSequence};
use crate::errors::{ChromosomeError, InvalidNucleotide};
use std::sync::Arc;

/// A chromosome: a named sequence organized into repeat units (RUs) and
/// higher-order repeats (HORs).
///
/// A `Chromosome` owns a mutable `Sequence` and provides helpers that expose
/// the repeat structure via `RepeatMap`.
/// For parallel/read-only operations prefer `SharedChromosome`
/// which shares sequence storage.
///
/// # Overview
///
/// `Chromosome` is a lightweight container for a sequence of
/// `Nucleotide`s together with a `RepeatMap` that describes the logical
/// boundaries of repeat units (RUs) and higher-order repeats (HORs). It
/// intentionally separates the raw sequence storage from the repeat
/// structure so consumers can modify sequences in-place while preserving
/// metadata.
///
/// # Examples
///
/// Create a short uniform chromosome, inspect its sequence and formatted
/// representation:
///
/// ```rust
/// # use centrevo_sim::genome::Chromosome;
/// # use centrevo_sim::base::Nucleotide;
/// let chr = Chromosome::uniform("chr1", Nucleotide::A, 2, 2, 1); // AA|AA
/// assert_eq!(chr.len(), 4);
/// assert_eq!(chr.to_string(), "AAAA");
/// assert_eq!(chr.to_formatted_string('-', '='), "AA-AA");
/// ```
#[derive(Debug, Clone)]
pub struct Chromosome {
    /// Type identifier. E.g., "chr1", "chr2", etc.
    /// Unique in a haplotype, but not unique in a population
    id: Arc<str>,
    /// The sequence data (mutable during simulation)
    sequence: Sequence,
    /// The repeat structure map
    map: RepeatMap,
}

impl Chromosome {
    /// Create a new chromosome
    pub fn new(id: impl Into<Arc<str>>, sequence: Sequence, map: RepeatMap) -> Self {
        Self {
            id: id.into(),
            sequence,
            map,
        }
    }

    /// Create a uniform chromosome where every base is the same `base`.
    ///
    /// This is a convenience constructor often used in tests and for creating
    /// deterministic initial populations.
    ///
    /// # Arguments
    /// - `id`: Chromosome identifier
    /// - `base`: Nucleotide base to fill the sequence with
    /// - `ru_length`: Length of each repeat unit in bases
    /// - `rus_per_hor`: Number of repeat units per higher-order repeat
    /// - `num_hors`: Number of higher-order repeats in the chromosome
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use centrevo_sim::genome::Chromosome;
    /// # use centrevo_sim::base::Nucleotide;
    /// // Two RUs per HOR (2 * 2 * 1 = 4 bases, all `A`)
    /// let chr = Chromosome::uniform("chr1", Nucleotide::A, 2, 2, 1);
    /// assert_eq!(chr.len(), 4);
    /// assert_eq!(chr.to_string(), "AAAA");
    /// ```
    pub fn uniform(
        id: impl Into<Arc<str>>,
        base: Nucleotide,
        ru_length: usize,
        rus_per_hor: usize,
        num_hors: usize,
    ) -> Self {
        let length = ru_length * rus_per_hor * num_hors;
        let mut sequence = Sequence::with_capacity(length);
        for _ in 0..length {
            sequence.push(base);
        }
        let map = RepeatMap::uniform(ru_length, rus_per_hor, num_hors);

        Self::new(id, sequence, map)
    }

    /// Create a chromosome from a nested vector structure `Vec<Vec<Vec<u8>>>`.
    ///
    /// The structure is:
    /// - Outer Vec `Vec<_>`: List of Higher-Order Repeats (HORs)
    /// - Middle Vec `Vec<Vec<_>>`: List of Repeat Units (RUs) within an HOR
    /// - Inner Vec `Vec<u8>`: Sequence of bytes (nucleotides) for an RU
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use centrevo_sim::genome::Chromosome;
    /// let nested = vec![vec![b"AC".to_vec(), b"GT".to_vec()]];
    /// let chr = Chromosome::from_nested_vec("chr1", nested).unwrap();
    /// assert_eq!(chr.to_string(), "ACGT");
    /// ```
    pub fn from_nested_vec(
        id: impl Into<Arc<str>>,
        nested: Vec<Vec<Vec<u8>>>,
    ) -> Result<Self, ChromosomeError> {
        let mut sequence_vec = Vec::new();
        let mut ru_offsets = vec![0];
        let mut hor_idx_offsets = vec![0];
        let mut current_ru_index = 0;

        for hor in nested {
            for ru in hor {
                for byte in ru {
                    let nuc = Nucleotide::from_ascii(byte).ok_or(InvalidNucleotide(byte))?;
                    sequence_vec.push(nuc);
                }
                ru_offsets.push(sequence_vec.len());
                current_ru_index += 1;
            }
            hor_idx_offsets.push(current_ru_index);
        }

        let sequence = Sequence::from_nucleotides(sequence_vec);
        let map =
            RepeatMap::new(ru_offsets, hor_idx_offsets).map_err(ChromosomeError::RepeatMapError)?;

        Ok(Self::new(id, sequence, map))
    }

    /// Create a chromosome from a nested string structure.
    ///
    /// The structure is `Vec<Vec<String>>`:
    /// - Outer Vec: List of Higher-Order Repeats (HORs)
    /// - Inner Vec: List of Repeat Units (RUs) as strings
    pub fn from_nested_strings(
        id: impl Into<Arc<str>>,
        nested: Vec<Vec<String>>,
    ) -> Result<Self, ChromosomeError> {
        let nested_bytes: Vec<Vec<Vec<u8>>> = nested
            .into_iter()
            .map(|hor| hor.into_iter().map(|ru| ru.into_bytes()).collect())
            .collect();

        Self::from_nested_vec(id, nested_bytes)
    }

    /// Create a chromosome from a formatted string with delimiters.
    ///
    /// `hor_delim` separates HORs.
    /// `ru_delim` separates RUs within an HOR.
    ///
    /// Example: "ACGT-ACGT=ACGT-ACGT" with ru_delim='-' and hor_delim='='
    /// parses to 2 HORs, each with 2 RUs.
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use centrevo_sim::genome::Chromosome;
    /// let input = "AC-GT=AA-TT";
    /// let chr = Chromosome::from_formatted_string("chr1", input, '=', '-').unwrap();
    /// assert_eq!(chr.to_string(), "ACGTAATT");
    /// ```
    pub fn from_formatted_string(
        id: impl Into<Arc<str>>,
        sequence: &str,
        hor_delim: char,
        ru_delim: char,
    ) -> Result<Self, ChromosomeError> {
        let hor_strings: Vec<&str> = sequence.split(hor_delim).collect();
        let mut nested_bytes = Vec::new();

        for hor_str in hor_strings {
            if hor_str.is_empty() {
                continue;
            }
            let ru_strings: Vec<&str> = hor_str.split(ru_delim).collect();
            let mut hor_vec = Vec::new();
            for ru_str in ru_strings {
                if ru_str.is_empty() {
                    continue;
                }
                hor_vec.push(ru_str.as_bytes().to_vec());
            }
            if !hor_vec.is_empty() {
                nested_bytes.push(hor_vec);
            }
        }

        Self::from_nested_vec(id, nested_bytes)
    }

    /// Return the chromosome identifier.
    ///
    /// The ID is stored in an `Arc<str>` so cloning the `Chromosome` cheaply
    /// shares the identifier.
    #[inline]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Borrow the underlying mutable `Sequence` for read-only operations.
    #[inline]
    pub fn sequence(&self) -> &Sequence {
        &self.sequence
    }

    /// Borrow the underlying `Sequence` mutably to apply in-place modifications
    /// such as mutation or recombination.
    #[inline]
    pub fn sequence_mut(&mut self) -> &mut Sequence {
        &mut self.sequence
    }

    /// Return the length of the chromosome in bases.
    #[inline]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Return true if the chromosome contains no bases.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Get the repeat map.
    #[inline]
    pub fn map(&self) -> &RepeatMap {
        &self.map
    }

    /// Return the number of complete HORs contained in this chromosome.
    #[inline]
    pub fn num_hors(&self) -> usize {
        self.map.num_hors()
    }

    /// Calculate GC content of the chromosome as a proportion in [0.0, 1.0].
    pub fn gc_content(&self) -> f64 {
        let mut gc_count = 0;
        let total = self.sequence.len();

        if total == 0 {
            return 0.0;
        }

        for &nuc in self.sequence.as_slice() {
            if matches!(nuc, Nucleotide::G | Nucleotide::C) {
                gc_count += 1;
            }
        }

        gc_count as f64 / total as f64
    }

    /// Convert the chromosome sequence to a formatted string inserting
    /// delimiters between repeat units and higher-order repeats.
    ///
    /// `ru_delim` is used between repeat units, `hor_delim` between HORs.
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use centrevo_sim::genome::Chromosome;
    /// # use centrevo_sim::base::Nucleotide;
    /// let chr = Chromosome::uniform("chr1", Nucleotide::A, 2, 2, 1);
    /// assert_eq!(chr.to_formatted_string('-', '='), "AA-AA");
    /// ```
    pub fn to_formatted_string(&self, ru_delim: char, hor_delim: char) -> String {
        let chars: Vec<char> = self
            .sequence
            .as_slice()
            .iter()
            .map(|&nuc| nuc.to_char())
            .collect();

        let mut result = String::with_capacity(chars.len() * 2);

        // We iterate through the sequence and insert delimiters based on the map.

        let num_rus = self.map.num_rus();
        for r in 0..num_rus {
            let (start, end) = self.map.get_ru_interval(r).unwrap();
            let seq_slice = &chars[start..end];
            for c in seq_slice {
                result.push(*c);
            }

            if r < num_rus - 1 {
                // Check if r+1 is a start of an HOR.
                // We need to check if `r+1` is in `hor_offsets`.
                // Since `hor_offsets` is sorted, we can check efficiently or just check if `r+1` matches any.
                // But `hor_offsets` contains RU indices.
                // If `r+1` is in `hor_offsets` (excluding 0), then it's an HOR boundary.

                // We can't easily check "contains" on Vec without iteration.
                // But we can iterate HORs.

                // Optimization: Pre-calculate boundary set? Or just iterate.
                // Since this is for display, performance is not critical.

                let mut is_hor_boundary = false;
                // hor_offsets[0] is 0. hor_offsets[1] is start of HOR 1 (end of HOR 0).
                // If r+1 == hor_offsets[k], then we are at boundary between HOR k-1 and HOR k.

                // We can iterate hor_offsets starting from 1.
                // But we don't have direct access to hor_offsets field here (it's private in RepeatMap).
                // We should add a method `is_hor_boundary(ru_index)` to RepeatMap.
                // For now, I'll assume I can add it.

                // Wait, I can't modify RepeatMap in this tool call.
                // I'll use `get_hor_interval` to check.
                // If `get_hor_interval(h).start == start of RU r+1`, then it's a boundary.

                for h in 1..self.map.num_hors() {
                    let (h_start, _) = self.map.get_hor_interval(h).unwrap();
                    let (next_ru_start, _) = self.map.get_ru_interval(r + 1).unwrap();
                    if h_start == next_ru_start {
                        is_hor_boundary = true;
                        break;
                    }
                }

                if is_hor_boundary {
                    result.push(hor_delim);
                } else {
                    result.push(ru_delim);
                }
            }
        }

        result
    }

    /// Create a `SharedChromosome` which provides an immutable, cheaply
    /// clonable view of this chromosome suitable for parallel read-only use.
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use centrevo_sim::genome::Chromosome;
    /// # use centrevo_sim::base::Nucleotide;
    /// let chr = Chromosome::uniform("chr1", Nucleotide::A, 2, 2, 1);
    /// let shared = chr.to_shared();
    /// assert_eq!(shared.len(), chr.len());
    /// let owned = shared.to_mutable();
    /// assert_eq!(owned.to_string(), chr.to_string());
    /// ```
    pub fn to_shared(&self) -> SharedChromosome {
        SharedChromosome {
            id: self.id.clone(),
            sequence: self.sequence.to_shared(),
            map: self.map.clone(),
        }
    }

    /// Get the sequence slice corresponding to the RU at `index`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use centrevo_sim::genome::Chromosome;
    /// # use centrevo_sim::base::Nucleotide;
    /// let chr = Chromosome::from_formatted_string("chr1", "AC-GT=AA-TT", '=', '-').unwrap();
    /// let slice = chr.get_ru_slice(1).unwrap();
    /// let s: String = slice.iter().map(|n| n.to_char()).collect();
    /// assert_eq!(s, "GT");
    /// ```
    pub fn get_ru_slice(&self, index: usize) -> Option<&[Nucleotide]> {
        let (start, end) = self.map.get_ru_interval(index)?;
        Some(&self.sequence.as_slice()[start..end])
    }

    /// Get the sequence slice corresponding to the HOR at `index`.
    /// 
    /// # Examples
    /// 
    /// ```rust
    /// # use centrevo_sim::genome::Chromosome;
    /// # use centrevo_sim::base::Nucleotide;
    /// let chr = Chromosome::from_formatted_string("chr1", "AC-GT=AA-TT", '=', '-').unwrap();
    /// let slice = chr.get_hor_slice(0).unwrap();
    /// let s: String = slice.iter().map(|n| n.to_char()).collect();
    /// assert_eq!(s, "ACGT");
    /// ```
    pub fn get_hor_slice(&self, index: usize) -> Option<&[Nucleotide]> {
        let (start, end) = self.map.get_hor_interval(index)?;
        Some(&self.sequence.as_slice()[start..end])
    }

    /// Calculate the Jaccard similarity of K-mers between an RU in this chromosome
    /// and an RU in another chromosome.
    ///
    /// # Arguments
    /// * `ru_index_self` - Index of RU in this chromosome
    /// * `other` - The other chromosome
    /// * `ru_index_other` - Index of RU in the other chromosome
    /// * `k` - K-mer size (e.g., 5 or 7)
    ///
    /// # Returns
    /// Similarity score in [0.0, 1.0]. Returns 0.0 if indices are invalid.
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use centrevo_sim::genome::Chromosome;
    /// # use centrevo_sim::base::Nucleotide;
    /// // small example with exact match
    /// let a = Chromosome::from_formatted_string("a", "ACGT=ACGT", '=', '-').unwrap();
    /// let b = Chromosome::from_formatted_string("b", "ACGT=ACGT", '=', '-').unwrap();
    /// // compare RU 0 of both (k=2 => "AC","CG","GT")
    /// assert_eq!(a.calculate_similarity(0, &b, 0, 2), 1.0);
    /// ```
    pub fn calculate_similarity(
        &self,
        ru_index_self: usize,
        other: &Self,
        ru_index_other: usize,
        k: usize,
    ) -> f64 {
        let slice1 = match self.get_ru_slice(ru_index_self) {
            Some(s) => s,
            None => return 0.0,
        };
        let slice2 = match other.get_ru_slice(ru_index_other) {
            Some(s) => s,
            None => return 0.0,
        };

        // Optimization: Exact match check
        if slice1 == slice2 {
            return 1.0;
        }

        if slice1.len() < k || slice2.len() < k {
            return 0.0;
        }

        // Simple K-mer counting using a HashSet is standard, but for performance
        // with small alphabet (ACGT) and small K, we could use a fixed-size array or bit packing.
        // However, for clarity and "good enough" performance given the window size,
        // we'll use a HashSet of u64 hashes or just the bytes if K is small.
        // Since K is small (e.g. 5-7), we can pack into u16 or u32.
        // A = 00, C = 01, G = 10, T = 11.
        // K=8 fits in u16.

        // Let's use a simple HashSet approach first for correctness.
        // We can optimize to bit-packing if needed.

        use std::collections::HashSet;

        let mut kmers1 = HashSet::with_capacity(slice1.len() - k + 1);
        for window in slice1.windows(k) {
            kmers1.insert(window);
        }

        let mut kmers2 = HashSet::with_capacity(slice2.len() - k + 1);
        for window in slice2.windows(k) {
            kmers2.insert(window);
        }

        let intersection = kmers1.intersection(&kmers2).count();
        let union = kmers1.len() + kmers2.len() - intersection;

        if union == 0 {
            0.0
        } else {
            intersection as f64 / union as f64
        }
    }

    /// Perform crossover between two chromosomes at potentially different positions.
    ///
    /// `pos1` is the break point in `self` and `pos2` is the break point in `other`.
    /// Chromosomes involved do not need to be the same total length; the crossover
    /// operation concatenates the left segment from one and the right segment from the other.
    /// Returns two new chromosomes as offspring.
    ///
    /// # Arguments
    /// * `other` - The other chromosome to crossover with
    /// * `pos1` - Crossover position in `self`
    /// * `pos2` - Crossover position in `other`
    /// 
    /// # Returns
    /// A Result containing a tuple of two new `Chromosome`s if successful,
    /// or an error string if the positions are out of bounds.
    /// 
    /// # Examples
    ///
    /// ```rust
    /// # use centrevo_sim::genome::Chromosome;
    /// # use centrevo_sim::base::Nucleotide;
    /// let a = Chromosome::uniform("a", Nucleotide::A, 1, 4, 1); // "AAAA"
    /// let b = Chromosome::uniform("b", Nucleotide::T, 1, 4, 1); // "TTTT"
    /// let (a1, b1) = a.crossover(&b, 2, 1).unwrap();
    /// // a1 = a[..2] + b[1..] = "AA" + "TTT" = "AATTT"
    /// assert_eq!(a1.to_string(), "AATTT");
    /// // b1 = b[..1] + a[2..] = "T" + "AA" = "TAA"
    /// assert_eq!(b1.to_string(), "TAA");
    /// ```
    pub fn crossover(
        &self,
        other: &Self,
        pos1: usize,
        pos2: usize,
    ) -> Result<(Self, Self), String> {
        if pos1 > self.len() {
            return Err(format!(
                "Position {pos1} out of bounds for self (len {})",
                self.len()
            ));
        }
        if pos2 > other.len() {
            return Err(format!(
                "Position {pos2} out of bounds for other (len {})",
                other.len()
            ));
        }

        // Create new sequences
        // Seq1 New = Self[..pos1] + Other[pos2..]
        let mut seq1_vec = Vec::with_capacity(pos1 + (other.len() - pos2));
        seq1_vec.extend_from_slice(&self.sequence.as_slice()[..pos1]);
        seq1_vec.extend_from_slice(&other.sequence.as_slice()[pos2..]);
        let seq1_new = Sequence::from_nucleotides(seq1_vec);

        // Seq2 New = Other[..pos2] + Self[pos1..]
        let mut seq2_vec = Vec::with_capacity(pos2 + (self.len() - pos1));
        seq2_vec.extend_from_slice(&other.sequence.as_slice()[..pos2]);
        seq2_vec.extend_from_slice(&self.sequence.as_slice()[pos1..]);
        let seq2_new = Sequence::from_nucleotides(seq2_vec);

        // Map crossover
        let (map1_left, map1_right) = self.map.split_at(pos1).map_err(|e| e.to_string())?;
        let (map2_left, map2_right) = other.map.split_at(pos2).map_err(|e| e.to_string())?;

        let map1_new = map1_left.merge(&map2_right).map_err(|e| e.to_string())?;
        let map2_new = map2_left.merge(&map1_right).map_err(|e| e.to_string())?;

        Ok((
            Self::new(self.id.clone(), seq1_new, map1_new),
            Self::new(other.id.clone(), seq2_new, map2_new),
        ))
    }

    /// Perform gene conversion operation.
    ///
    /// Replaces a tract in `self` (recipient) starting at `recipient_start` of length `length`
    /// with a tract from `other` (donor) starting at `donor_start` of length `length`.
    ///
    /// Chromosomes involved do not need to be the same total length; the donor tract
    /// is copied into the recipient and the resulting chromosome length may change
    /// accordingly.
    /// 
    /// # Arguments
    /// * `other` - Donor chromosome
    /// * `recipient_start` - Start position in recipient (self)
    /// * `donor_start` - Start position in donor (other)
    /// * `length` - Length of tract to replace
    /// 
    /// # Returns
    /// A Result containing a new `Chromosome` if successful,
    /// or an error string if the specified ranges are out of bounds.
    ///
    /// # Examples
    ///
    /// Simple example with equal total lengths:
    ///
    /// ```rust
    /// # use centrevo_sim::genome::Chromosome;
    /// # use centrevo_sim::base::Nucleotide;
    /// let rec = Chromosome::uniform("r", Nucleotide::A, 1, 4, 1); // "AAAA"
    /// let donor = Chromosome::uniform("d", Nucleotide::T, 1, 4, 1); // "TTTT"
    /// let out = rec.gene_conversion(&donor, 1, 1, 2).unwrap();
    /// // Replace rec[1..3] with donor[1..3] => "A" + "TT" + "A" = "ATTA"
    /// assert_eq!(out.to_string(), "ATTA");
    /// ```
    ///
    /// Different total lengths are supported. Donor tracts are copied into the
    /// recipient and may change the resulting chromosome length.
    ///
    /// ```rust
    /// # use centrevo_sim::genome::Chromosome;
    /// # use centrevo_sim::base::Nucleotide;
    /// // Recipient length 4, donor length 6
    /// let rec = Chromosome::uniform("r", Nucleotide::A, 1, 4, 1); // "AAAA"
    /// let donor = Chromosome::uniform("d", Nucleotide::T, 1, 6, 1); // "TTTTTT"
    /// let out = rec.gene_conversion(&donor, 1, 3, 2).unwrap();
    /// // Replace rec[1..3] with donor[3..5] => "A" + "TT" + "A" = "ATTA"
    /// assert_eq!(out.to_string(), "ATTA");
    /// ```
    pub fn gene_conversion(
        &self,
        other: &Self,
        recipient_start: usize,
        donor_start: usize,
        length: usize,
    ) -> Result<Self, String> {
        let end1 = recipient_start + length;
        let end2 = donor_start + length;

        if end1 > self.len() {
            return Err(format!(
                "Recipient range {recipient_start}..{end1} out of bounds (len {})",
                self.len()
            ));
        }
        if end2 > other.len() {
            return Err(format!(
                "Donor range {donor_start}..{end2} out of bounds (len {})",
                other.len()
            ));
        }

        // 1. Split Recipient: Left | Tract | Right
        // Split at recipient_start -> (Left, MidRight)
        let (map_left, map_mid_right) = self.map.split_at(recipient_start).map_err(|e| e.to_string())?;
        // Split MidRight at length -> (OldTract, Right)
        let (_map_old_tract, map_right) =
            map_mid_right.split_at(length).map_err(|e| e.to_string())?;

        // 2. Split Donor: DLeft | NewTract | DRight
        // Split at donor_start -> (DLeft, DMidRight)
        let (_map_d_left, map_d_mid_right) =
            other.map.split_at(donor_start).map_err(|e| e.to_string())?;
        // Split DMidRight at length -> (NewTract, DRight)
        let (map_new_tract, _map_d_right) = map_d_mid_right
            .split_at(length)
            .map_err(|e| e.to_string())?;

        // 3. Construct New Map: Left + NewTract + Right
        let map_temp = map_left.merge(&map_new_tract).map_err(|e| e.to_string())?;
        let map_final = map_temp.merge(&map_right).map_err(|e| e.to_string())?;

        // 4. Construct New Sequence
        // Left part
        // Compute final length for capacity: recipient left + donor tract + recipient right
            let final_len = recipient_start + (end2 - donor_start) + (self.len() - end1);
        let mut new_seq_vec = Vec::with_capacity(final_len);
        new_seq_vec.extend_from_slice(&self.sequence.as_slice()[..recipient_start]);
        // New Tract
        new_seq_vec.extend_from_slice(&other.sequence.as_slice()[donor_start..end2]);
        // Right part
        new_seq_vec.extend_from_slice(&self.sequence.as_slice()[end1..]);

        let new_seq = Sequence::from_nucleotides(new_seq_vec);

        Ok(Self::new(self.id.clone(), new_seq, map_final))
    }


}

impl std::fmt::Display for Chromosome {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.sequence)
    }
}

/// Immutable, shared view of a chromosome.
///
/// `SharedChromosome` contains a reference-counted view of the sequence data
/// and the repeat-structure metadata. It is intended for read-only access and
/// cheap cloning across threads. Converting a `SharedChromosome` back to an
/// owned `Chromosome` performs minimal copying where necessary and keeps the
/// shared, read-only operations efficient.
///
/// # Examples
///
/// ```rust
/// # use centrevo_sim::genome::Chromosome;
/// # use centrevo_sim::base::Nucleotide;
/// let chr = Chromosome::uniform("chr1", Nucleotide::A, 2, 2, 1);
/// let shared = chr.to_shared();
/// let owned = shared.to_mutable();
/// assert_eq!(owned.to_string(), chr.to_string());
/// ```
#[derive(Debug, Clone)]
pub struct SharedChromosome {
    id: Arc<str>,
    sequence: SharedSequence,
    map: RepeatMap,
}

impl SharedChromosome {
    /// Return the chromosome identifier.
    #[inline]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Borrow the shared, immutable sequence.
    #[inline]
    pub fn sequence(&self) -> &SharedSequence {
        &self.sequence
    }

    /// Return the length of the shared chromosome in bases.
    #[inline]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Return true if the shared chromosome contains no bases.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Get the repeat map.
    #[inline]
    pub fn map(&self) -> &RepeatMap {
        &self.map
    }

    /// Calculate GC content for the shared chromosome.
    pub fn gc_content(&self) -> f64 {
        let mut gc_count = 0;
        let total = self.sequence.len();

        if total == 0 {
            return 0.0;
        }

        for &nuc in self.sequence.as_slice() {
            if matches!(nuc, Nucleotide::G | Nucleotide::C) {
                gc_count += 1;
            }
        }

        gc_count as f64 / total as f64
    }

    /// Convert the shared chromosome into an owned `Chromosome` with
    /// mutable sequence data (clones the underlying indices if necessary).
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use centrevo_sim::genome::Chromosome;
    /// # use centrevo_sim::base::Nucleotide;
    /// let chr = Chromosome::uniform("chr1", Nucleotide::A, 2, 2, 1);
    /// let shared = chr.to_shared();
    /// let owned = shared.to_mutable();
    /// assert_eq!(owned.to_string(), chr.to_string());
    /// ```
    pub fn to_mutable(&self) -> Chromosome {
        Chromosome {
            id: self.id.clone(),
            sequence: self.sequence.to_mutable(),
            map: self.map.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    fn test_sequence() -> Sequence {
        Sequence::from_str("ACGTACGTACGTACGT").unwrap()
    }

    fn test_map(ru_len: usize, rus_per_hor: usize, num_hors: usize) -> RepeatMap {
        RepeatMap::uniform(ru_len, rus_per_hor, num_hors)
    }

    // ===== Chromosome Tests =====

    #[test]
    fn test_chromosome_new() {
        let seq = test_sequence();
        // 16 bp. RU=4. RUs=4. HOR=2 RUs. HORs=2.
        let map = test_map(4, 2, 2);
        let chr = Chromosome::new("chr1", seq.clone(), map);

        assert_eq!(chr.id(), "chr1");
        assert_eq!(chr.len(), 16);
        assert_eq!(chr.map().num_rus(), 4);
        assert_eq!(chr.map().num_hors(), 2);
    }

    #[test]
    fn test_chromosome_uniform() {
        // ru_length=10, rus_per_hor=10, num_hors=1 -> total length = 10*10*1 = 100
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 10, 10, 1);

        assert_eq!(chr.len(), 100);
        assert_eq!(chr.to_string(), "A".repeat(100));
        assert_eq!(chr.map().num_rus(), 10);
        assert_eq!(chr.map().num_hors(), 1);
    }

    #[test]
    fn test_chromosome_id() {
        let chr = Chromosome::uniform("test_id", Nucleotide::A, 10, 5, 2);
        assert_eq!(chr.id(), "test_id");
    }

    #[test]
    fn test_chromosome_id_shared() {
        let chr1 = Chromosome::uniform("chr1", Nucleotide::A, 10, 5, 2);
        let chr2 = chr1.clone();

        // Both should share the same ID Arc
        assert_eq!(chr1.id(), chr2.id());
    }

    #[test]
    fn test_chromosome_sequence() {
        let seq = test_sequence();
        let map = test_map(4, 2, 2);
        let chr = Chromosome::new("chr1", seq.clone(), map);

        assert_eq!(chr.sequence().to_string(), "ACGTACGTACGTACGT");
    }

    #[test]
    fn test_chromosome_sequence_mut() {
        let seq = test_sequence();
        let map = test_map(4, 2, 2);
        let mut chr = Chromosome::new("chr1", seq, map);

        chr.sequence_mut().set(0, Nucleotide::T).unwrap();
        assert_eq!(chr.sequence().to_string().chars().next(), Some('T'));
    }

    #[test]
    fn test_chromosome_len() {
        // ru_length=10, rus_per_hor=10, num_hors=1 -> total length = 100
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 10, 10, 1);
        assert_eq!(chr.len(), 100);
    }

    #[test]
    fn test_chromosome_is_empty() {
        let empty_seq = Sequence::new();
        let map = test_map(4, 2, 0); // 0 HORs
        let chr = Chromosome::new("chr1", empty_seq, map);
        assert!(chr.is_empty());

        let non_empty = Chromosome::uniform("chr2", Nucleotide::A, 10, 5, 2);
        assert!(!non_empty.is_empty());
    }

    #[test]
    fn test_chromosome_num_hors() {
        // ru_length=10, rus_per_hor=10, num_hors=5 -> total length = 500
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 10, 10, 5);
        assert_eq!(chr.num_hors(), 5);
    }

    #[test]
    fn test_chromosome_gc_content_all_gc() {
        let seq = Sequence::from_str("GCGCGCGC").unwrap();
        let map = test_map(2, 2, 2);
        let chr = Chromosome::new("chr1", seq, map);
        assert_eq!(chr.gc_content(), 1.0);
    }

    #[test]
    fn test_chromosome_gc_content_all_at() {
        let seq = Sequence::from_str("ATATATATAT").unwrap();
        // test_map(2, 2, 2) -> length 8.
        // We need map matching sequence.
        let map = RepeatMap::uniform(2, 2, 2); // Length 8.
        // Sequence is 10.
        // Chromosome::new doesn't validate map vs sequence length yet (it should).
        // But for GC content it doesn't matter.
        let chr = Chromosome::new("chr1", seq, map);
        assert_eq!(chr.gc_content(), 0.0);
    }

    #[test]
    fn test_chromosome_gc_content_half() {
        let seq = Sequence::from_str("ACGT").unwrap();
        let map = test_map(2, 2, 1);
        let chr = Chromosome::new("chr1", seq, map);
        assert_eq!(chr.gc_content(), 0.5);
    }

    #[test]
    fn test_chromosome_gc_content_empty() {
        let seq = Sequence::new();
        let map = test_map(2, 2, 0);
        let chr = Chromosome::new("chr1", seq, map);
        assert_eq!(chr.gc_content(), 0.0);
    }

    #[test]
    fn test_chromosome_to_string() {
        let seq = Sequence::from_str("ACGT").unwrap();
        let map = test_map(2, 2, 1);
        let chr = Chromosome::new("chr1", seq, map);
        assert_eq!(chr.to_string(), "ACGT");
    }

    #[test]
    fn test_chromosome_to_formatted_string_basic() {
        let seq = Sequence::from_str("ACGTACGTACGTACGT").unwrap();
        // RU=4. HOR=2 RUs (8bp). Total 16bp.
        let map = test_map(4, 2, 2);
        let chr = Chromosome::new("chr1", seq, map);

        let formatted = chr.to_formatted_string('|', '#');
        // RU length = 4, HOR length = 8
        // Expected: ACGT|ACGT#ACGT|ACGT
        // My implementation:
        // RU 0: ACGT. Next is RU 1. RU 1 is NOT HOR start (HOR starts at 0, 2).
        // Delim: |
        // RU 1: ACGT. Next is RU 2. RU 2 IS HOR start (HOR 1 starts at RU 2).
        // Delim: #
        // RU 2: ACGT. Next is RU 3. RU 3 is NOT HOR start.
        // Delim: |
        // RU 3: ACGT. Last RU. No delim.
        assert_eq!(formatted, "ACGT|ACGT#ACGT|ACGT");
    }

    #[test]
    fn test_chromosome_to_formatted_string_single_ru() {
        let seq = Sequence::from_str("ACGT").unwrap();
        let map = test_map(4, 1, 1);
        let chr = Chromosome::new("chr1", seq, map);

        let formatted = chr.to_formatted_string('|', '#');
        assert_eq!(formatted, "ACGT"); // No delimiters for single RU
    }

    #[test]
    fn test_chromosome_to_formatted_string_custom_delimiters() {
        let seq = Sequence::from_str("ACGTACGT").unwrap();
        let map = test_map(2, 2, 2); // RU=2. HOR=2 RUs (4bp). Total 8bp.
        let chr = Chromosome::new("chr1", seq, map);

        let formatted = chr.to_formatted_string('-', '=');
        // RU 0: AC. Next RU 1. Not HOR start. Delim -
        // RU 1: GT. Next RU 2. HOR start. Delim =
        // RU 2: AC. Next RU 3. Not HOR start. Delim -
        // RU 3: GT.
        // Expected: AC-GT=AC-GT
        assert_eq!(formatted, "AC-GT=AC-GT");
    }

    #[test]
    fn test_chromosome_clone() {
        let chr1 = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5);
        let chr2 = chr1.clone();

        assert_eq!(chr1.id(), chr2.id());
        assert_eq!(chr1.len(), chr2.len());
        assert_eq!(chr1.map(), chr2.map());
        assert_eq!(chr1.to_string(), chr2.to_string());
    }

    // ===== SharedChromosome Tests =====

    #[test]
    fn test_chromosome_to_shared() {
        // ru_length=10, rus_per_hor=10, num_hors=1 -> total length = 100
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 10, 10, 1);
        let shared = chr.to_shared();

        assert_eq!(shared.id(), "chr1");
        assert_eq!(shared.len(), 100);
        assert_eq!(shared.map().num_rus(), 10);
    }

    #[test]
    fn test_shared_chromosome_id() {
        let chr = Chromosome::uniform("test_id", Nucleotide::A, 10, 5, 2);
        let shared = chr.to_shared();
        assert_eq!(shared.id(), "test_id");
    }

    #[test]
    fn test_shared_chromosome_sequence() {
        let seq = test_sequence();
        let map = test_map(4, 2, 2);
        let chr = Chromosome::new("chr1", seq, map);
        let shared = chr.to_shared();

        assert_eq!(shared.sequence().len(), 16);
    }

    #[test]
    fn test_shared_chromosome_len() {
        // ru_length=10, rus_per_hor=10, num_hors=1 -> total length = 100
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 10, 10, 1);
        let shared = chr.to_shared();
        assert_eq!(shared.len(), 100);
    }

    #[test]
    fn test_shared_chromosome_is_empty() {
        let empty_seq = Sequence::new();
        let map = test_map(4, 2, 0);
        let chr = Chromosome::new("chr1", empty_seq, map);
        let shared = chr.to_shared();
        assert!(shared.is_empty());
    }

    #[test]
    fn test_shared_chromosome_gc_content() {
        let seq = Sequence::from_str("GCGCGCGC").unwrap();
        let map = test_map(2, 2, 2);
        let chr = Chromosome::new("chr1", seq, map);
        let shared = chr.to_shared();
        assert_eq!(shared.gc_content(), 1.0);
    }

    #[test]
    fn test_shared_chromosome_clone_is_cheap() {
        // ru_length=10, rus_per_hor=10, num_hors=10 -> total length = 1000
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 10, 10, 10);
        let shared1 = chr.to_shared();
        let shared2 = shared1.clone();

        // Both should share the same sequence data
        assert_eq!(shared1.sequence().strong_count(), 2);
        assert_eq!(shared2.sequence().strong_count(), 2);
    }

    #[test]
    fn test_shared_chromosome_to_mutable() {
        // ru_length=10, rus_per_hor=10, num_hors=1 -> total length = 100
        let chr1 = Chromosome::uniform("chr1", Nucleotide::A, 10, 10, 1);
        let shared = chr1.to_shared();
        let chr2 = shared.to_mutable();

        assert_eq!(chr2.id(), "chr1");
        assert_eq!(chr2.len(), 100);
        assert_eq!(chr2.map(), chr1.map());
        assert_eq!(chr2.to_string(), chr1.to_string());
    }

    #[test]
    fn test_roundtrip_mutable_to_shared_to_mutable() {
        // ru_length=10, rus_per_hor=10, num_hors=1 -> total length = 100
        let chr1 = Chromosome::uniform("chr1", Nucleotide::A, 10, 10, 1);
        let original_str = chr1.to_string();

        let shared = chr1.to_shared();
        let chr2 = shared.to_mutable();

        assert_eq!(chr2.to_string(), original_str);
        assert_eq!(chr2.id(), "chr1");
    }

    #[test]
    fn test_shared_chromosome_immutability() {
        // ru_length=10, rus_per_hor=10, num_hors=1 -> total length = 100
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 10, 10, 1);
        let shared = chr.to_shared();

        // Clone it
        let _cloned = shared.clone();

        // Original shared should be unchanged (immutable)
        assert_eq!(shared.len(), 100);
    }

    // ===== Integration Tests =====

    #[test]
    fn test_chromosome_complex_structure() {
        // Create a chromosome with realistic parameters
        // RU = 171 bp, 12 RUs per HOR = 2052 bp HOR, 10 HORs = 20520 bp total
        let chr = Chromosome::uniform(
            "chr1",
            Nucleotide::A,
            171, // ru_length
            12,  // rus_per_hor
            10,  // num_hors
        );

        assert_eq!(chr.num_hors(), 10);
    }

    #[test]
    fn test_chromosome_mutation_scenario() {
        // ru_length=10, rus_per_hor=10, num_hors=1 -> total length = 100
        let mut chr = Chromosome::uniform("chr1", Nucleotide::A, 10, 10, 1);

        // Mutate some bases
        chr.sequence_mut().set(0, Nucleotide::T).unwrap();
        chr.sequence_mut().set(50, Nucleotide::G).unwrap();

        // Check that mutations were applied
        assert_eq!(chr.sequence().get(0), Some(Nucleotide::T));
        assert_eq!(chr.sequence().get(50), Some(Nucleotide::G));
    }

    #[test]
    fn test_gene_conversion_different_chromosome_lengths() {
        // Recipient length 4, donor length 6
        let rec = Chromosome::uniform("r", Nucleotide::A, 1, 4, 1); // "AAAA"
        let donor = Chromosome::uniform("d", Nucleotide::T, 1, 6, 1); // "TTTTTT"

        // Replace rec[1..3] with donor[3..5] => "A" + "TT" + "A" = "ATTA"
        let out = rec.gene_conversion(&donor, 1, 3, 2).unwrap();
        assert_eq!(out.to_string(), "ATTA");
    }

    #[test]
    fn test_chromosome_different_ru_configurations() {
        // All create chromosomes with 4000 bp total length but different structures
        // chr1: ru_length=10, rus_per_hor=10, num_hors=40 -> 10*10*40 = 4000
        let chr1 = Chromosome::uniform("chr1", Nucleotide::A, 10, 10, 40);
        // chr2: ru_length=5, rus_per_hor=10, num_hors=80 -> 5*10*80 = 4000
        let chr2 = Chromosome::uniform("chr2", Nucleotide::A, 5, 10, 80);
        // chr3: ru_length=20, rus_per_hor=2, num_hors=100 -> 20*2*100 = 4000
        let chr3 = Chromosome::uniform("chr3", Nucleotide::A, 20, 2, 100);

        // All same length but different HOR structures
        assert_eq!(chr1.len(), chr2.len());
        assert_eq!(chr2.len(), chr3.len());
        assert_eq!(chr1.len(), 4000);

        // But different num_hors
        assert_eq!(chr1.num_hors(), 40);
        assert_eq!(chr2.num_hors(), 80);
        assert_eq!(chr3.num_hors(), 100);
    }

    #[test]
    fn test_chromosome_large() {
        // Test with larger, more realistic sizes
        // RU = 171 bp, 12 RUs per HOR = 2052 bp per HOR, 487 HORs ~ 1 Mbp
        let chr = Chromosome::uniform(
            "chr1",
            Nucleotide::A,
            171, // ru_length
            12,  // rus_per_hor
            487, // num_hors -> 171*12*487 = 999,324 bp
        );

        assert_eq!(chr.len(), 999_324);
        assert_eq!(chr.num_hors(), 487);
    }

    #[test]
    fn test_from_nested_vec() {
        // HOR 1: RU1="AC", RU2="GT"
        // HOR 2: RU3="AA", RU4="TT"
        let nested = vec![
            vec![b"AC".to_vec(), b"GT".to_vec()],
            vec![b"AA".to_vec(), b"TT".to_vec()],
        ];

        let chr = Chromosome::from_nested_vec("chr1", nested).unwrap();

        assert_eq!(chr.len(), 8);
        assert_eq!(chr.num_hors(), 2);
        assert_eq!(chr.map().num_rus(), 4);
        assert_eq!(chr.to_string(), "ACGTAATT");

        // Check structure
        // HOR 0: RUs 0, 1.
        // HOR 1: RUs 2, 3.
        let (start, end) = chr.map().get_hor_interval(0).unwrap();
        assert_eq!(start, 0);
        assert_eq!(end, 4); // "ACGT"

        let (start, end) = chr.map().get_hor_interval(1).unwrap();
        assert_eq!(start, 4);
        assert_eq!(end, 8); // "AATT"
    }

    #[test]
    fn test_from_nested_strings() {
        let nested = vec![
            vec!["AC".to_string(), "GT".to_string()],
            vec!["AA".to_string(), "TT".to_string()],
        ];

        let chr = Chromosome::from_nested_strings("chr1", nested).unwrap();
        assert_eq!(chr.to_string(), "ACGTAATT");
        assert_eq!(chr.num_hors(), 2);
    }

    #[test]
    fn test_from_string_with_delimiters() {
        let input = "AC-GT=AA-TT";
        let chr = Chromosome::from_formatted_string("chr1", input, '=', '-').unwrap();

        assert_eq!(chr.to_string(), "ACGTAATT");
        assert_eq!(chr.num_hors(), 2);
        assert_eq!(chr.map().num_rus(), 4);
    }

    #[test]
    fn test_from_string_with_delimiters_empty_parts() {
        // "AC-GT=" -> 1 HOR with 2 RUs. Trailing delimiter creates empty string which should be ignored.
        let input = "AC-GT=";
        let chr = Chromosome::from_formatted_string("chr1", input, '=', '-').unwrap();

        assert_eq!(chr.to_string(), "ACGT");
        assert_eq!(chr.num_hors(), 1);
        assert_eq!(chr.map().num_rus(), 2);
    }

    #[test]
    fn test_from_nested_vec_invalid_nucleotide() {
        let nested = vec![
            vec![vec![b'A', b'X']], // 'X' is invalid
        ];
        let result = Chromosome::from_nested_vec("chr1", nested);
        assert!(matches!(result, Err(ChromosomeError::InvalidNucleotide(_))));
    }

    #[test]
    fn test_from_nested_vec_empty() {
        let nested: Vec<Vec<Vec<u8>>> = vec![];
        let chr = Chromosome::from_nested_vec("chr1", nested).unwrap();
        assert!(chr.is_empty());
        assert_eq!(chr.num_hors(), 0);
        assert_eq!(chr.map().num_rus(), 0);
    }

    #[test]
    fn test_from_nested_vec_empty_hor() {
        let nested = vec![
            vec![], // Empty HOR
            vec![b"AC".to_vec()],
        ];
        // This creates an HOR with 0 RUs.
        // RepeatMap allows this?
        // hor_offsets: [0, 0, 1].
        // HOR 0: RUs 0..0 (empty).
        // HOR 1: RUs 0..1.
        let chr = Chromosome::from_nested_vec("chr1", nested).unwrap();
        assert_eq!(chr.num_hors(), 2);
        assert_eq!(chr.map().num_rus(), 1);
        assert_eq!(chr.to_string(), "AC");
    }

    #[test]
    fn test_from_nested_vec_empty_ru() {
        let nested = vec![
            vec![vec![], b"AC".to_vec()], // First RU is empty
        ];
        // RU 0: length 0.
        // RU 1: length 2.
        let chr = Chromosome::from_nested_vec("chr1", nested).unwrap();
        assert_eq!(chr.map().num_rus(), 2);
        assert_eq!(chr.to_string(), "AC");

        let (start, end) = chr.map().get_ru_interval(0).unwrap();
        assert_eq!(start, 0);
        assert_eq!(end, 0);

        let (start, end) = chr.map().get_ru_interval(1).unwrap();
        assert_eq!(start, 0);
        assert_eq!(end, 2);
    }

    #[test]
    fn test_from_string_with_delimiters_consecutive_delimiters() {
        // "AC--GT" -> Empty RU between AC and GT?
        // The implementation splits by delimiter.
        // "AC--GT".split('-') -> ["AC", "", "GT"].
        // The implementation checks `if ru_str.is_empty() { continue; }`.
        // So empty RUs are skipped.
        let input = "AC--GT";
        let chr = Chromosome::from_formatted_string("chr1", input, '=', '-').unwrap();
        assert_eq!(chr.to_string(), "ACGT");
        assert_eq!(chr.map().num_rus(), 2); // AC, GT. The empty one is skipped.
    }

    #[test]
    fn test_from_string_with_delimiters_empty_string() {
        let input = "";
        let chr = Chromosome::from_formatted_string("chr1", input, '=', '-').unwrap();
        assert!(chr.is_empty());
    }
}
