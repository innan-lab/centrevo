use super::repeat_map::RepeatMap;
use crate::base::{GenomeArena, Nucleotide, SequenceSlice};
use crate::errors::{ChromosomeError, InvalidNucleotide};
use std::sync::Arc;

/// A chromosome: a named sequence organized into repeat units (RUs) and
/// higher-order repeats (HORs).
///
/// This is a lightweight handle to data stored in a `GenomeArena`.
#[derive(Debug, Clone)]
pub struct Chromosome {
    id: Arc<str>,
    data: SequenceSlice,
    map: RepeatMap,
}

impl Chromosome {
    /// Create a new chromosome handle.
    ///
    /// # Examples
    ///
    /// ```
    /// use centrevo_sim::genome::{Chromosome, RepeatMap};
    /// use centrevo_sim::base::{Nucleotide, GenomeArena};
    ///
    /// let mut arena = GenomeArena::new();
    /// let seq: Vec<Nucleotide> = vec![Nucleotide::A; 10];
    /// let map = RepeatMap::uniform(1, 10, 1);
    /// let data = arena.alloc(&seq);
    ///
    /// let chr = Chromosome::new("chr1", data, map);
    /// assert_eq!(chr.id(), "chr1");
    /// assert_eq!(chr.len(), 10);
    /// ```
    pub fn new(id: impl Into<Arc<str>>, data: SequenceSlice, map: RepeatMap) -> Self {
        Self {
            id: id.into(),
            data,
            map,
        }
    }

    /// Create a uniform chromosome where every base is the same `base`.
    pub fn uniform(
        id: impl Into<Arc<str>>,
        base: Nucleotide,
        ru_length: usize,
        rus_per_hor: usize,
        num_hors: usize,
        arena: &mut GenomeArena,
    ) -> Self {
        let length = ru_length * rus_per_hor * num_hors;
        let mut sequence_vec = Vec::with_capacity(length);
        for _ in 0..length {
            sequence_vec.push(base);
        }
        let data = arena.alloc(&sequence_vec);
        let map = RepeatMap::uniform(ru_length, rus_per_hor, num_hors);

        Self::new(id, data, map)
    }

    /// Create from nested vec.
    pub fn from_nested_vec(
        id: impl Into<Arc<str>>,
        nested: Vec<Vec<Vec<u8>>>,
        arena: &mut GenomeArena,
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

        let data = arena.alloc(&sequence_vec);
        let map =
            RepeatMap::new(ru_offsets, hor_idx_offsets).map_err(ChromosomeError::RepeatMapError)?;

        Ok(Self::new(id, data, map))
    }

    /// Create from nested strings.
    pub fn from_nested_strings(
        id: impl Into<Arc<str>>,
        nested: Vec<Vec<String>>,
        arena: &mut GenomeArena,
    ) -> Result<Self, ChromosomeError> {
        let nested_bytes: Vec<Vec<Vec<u8>>> = nested
            .into_iter()
            .map(|hor| hor.into_iter().map(|ru| ru.into_bytes()).collect())
            .collect();

        Self::from_nested_vec(id, nested_bytes, arena)
    }

    /// Create from formatted string.
    ///
    /// The formatted string allows visualizing the repeat structure:
    /// - `hor_delim`: Separates higher-order repeats (HORs)
    /// - `ru_delim`: Separates repeat units (RUs) within an HOR
    ///
    /// # Examples
    ///
    /// ```
    /// use centrevo_sim::genome::Chromosome;
    /// use centrevo_sim::base::GenomeArena;
    ///
    /// let mut arena = GenomeArena::new();
    /// // 2 HORs, separated by '=', each with RUs separated by '-'
    /// let input = "AC-GT=AA-TT";
    /// let chr = Chromosome::from_formatted_string("chr1", input, '=', '-', &mut arena).unwrap();
    ///
    /// assert_eq!(chr.len(), 8);
    /// assert_eq!(chr.map().num_hors(), 2);
    /// assert_eq!(chr.to_formatted_string("-", "=", &arena), "AC-GT=AA-TT");
    /// ```
    pub fn from_formatted_string(
        id: impl Into<Arc<str>>,
        sequence: &str,
        hor_delim: char,
        ru_delim: char,
        arena: &mut GenomeArena,
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

        Self::from_nested_vec(id, nested_bytes, arena)
    }

    /// Return the chromosome identifier.
    #[inline]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Get the nucleotide sequence data.
    pub fn sequence<'a>(&self, arena: &'a GenomeArena) -> &'a [Nucleotide] {
        arena.get(self.data)
    }

    /// Get mutable nucleotide sequence data.
    pub fn sequence_mut<'a>(&self, arena: &'a mut GenomeArena) -> &'a mut [Nucleotide] {
        arena.get_mut(self.data)
    }

    /// Return the length of the chromosome in bases.
    #[inline]
    pub fn len(&self) -> usize {
        self.data.len as usize
    }

    /// Return true if the chromosome contains no bases.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.data.len == 0
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
    pub fn gc_content(&self, arena: &GenomeArena) -> f64 {
        let mut gc_count = 0;
        let total = self.len();

        if total == 0 {
            return 0.0;
        }

        for &nuc in self.sequence(arena) {
            if matches!(nuc, Nucleotide::G | Nucleotide::C) {
                gc_count += 1;
            }
        }

        gc_count as f64 / total as f64
    }

    /// Convert the chromosome sequence to a formatted string.
    pub fn to_formatted_string(
        &self,
        ru_delim: &str,
        hor_delim: &str,
        arena: &crate::base::GenomeArena,
    ) -> String {
        let chars: Vec<char> = self
            .sequence(arena)
            .iter()
            .map(|&nuc| nuc.to_char())
            .collect();

        let mut result = String::with_capacity(chars.len() * 2);

        let num_rus = self.map.num_rus();
        for r in 0..num_rus {
            let (start, end) = self.map.get_ru_interval(r).unwrap();
            let seq_slice = &chars[start..end];
            for c in seq_slice {
                result.push(*c);
            }

            if r < num_rus - 1 {
                let mut is_hor_boundary = false;
                for h in 1..self.map.num_hors() {
                    let (h_start, _) = self.map.get_hor_interval(h).unwrap();
                    let (next_ru_start, _) = self.map.get_ru_interval(r + 1).unwrap();
                    if h_start == next_ru_start {
                        is_hor_boundary = true;
                        break;
                    }
                }

                if is_hor_boundary {
                    result.push_str(hor_delim);
                } else {
                    result.push_str(ru_delim);
                }
            }
        }

        result
    }

    /// Get the sequence slice corresponding to the RU at `index`.
    pub fn get_ru_slice<'a>(
        &self,
        index: usize,
        arena: &'a GenomeArena,
    ) -> Option<&'a [Nucleotide]> {
        let (start, end) = self.map.get_ru_interval(index)?;
        Some(&self.sequence(arena)[start..end])
    }

    /// Get the sequence slice corresponding to the HOR at `index`.
    pub fn get_hor_slice<'a>(
        &self,
        index: usize,
        arena: &'a GenomeArena,
    ) -> Option<&'a [Nucleotide]> {
        let (start, end) = self.map.get_hor_interval(index)?;
        Some(&self.sequence(arena)[start..end])
    }

    /// Calculate the Jaccard similarity.
    pub fn calculate_similarity(
        &self,
        ru_index_self: usize,
        other: &Self,
        ru_index_other: usize,
        k: usize,
        arena: &GenomeArena,
    ) -> f64 {
        let slice1 = match self.get_ru_slice(ru_index_self, arena) {
            Some(s) => s,
            None => return 0.0,
        };
        let slice2 = match other.get_ru_slice(ru_index_other, arena) {
            Some(s) => s,
            None => return 0.0,
        };

        if slice1 == slice2 {
            return 1.0;
        }

        if slice1.len() < k || slice2.len() < k {
            return 0.0;
        }

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

    /// Perform crossover between two chromosomes at the given position.
    ///
    /// # Examples
    ///
    /// ```
    /// use centrevo_sim::genome::{Chromosome, RepeatMap};
    /// use centrevo_sim::base::{Nucleotide, GenomeArena};
    ///
    /// let mut arena = GenomeArena::new();
    /// let seq_a: Vec<Nucleotide> = vec![Nucleotide::A; 4];
    /// let seq_t: Vec<Nucleotide> = vec![Nucleotide::T; 4];
    /// let map = RepeatMap::uniform(1, 4, 1);
    /// let data_a = arena.alloc(&seq_a);
    /// let data_t = arena.alloc(&seq_t);
    /// let chr_a = Chromosome::new("A", data_a, map.clone());
    /// let chr_t = Chromosome::new("T", data_t, map.clone());
    ///
    /// let (off1, off2) = chr_a.crossover(&chr_t, 2, 2, &mut arena).unwrap();
    /// assert_eq!(off1.to_formatted_string("", "", &arena), "AATT");
    /// assert_eq!(off2.to_formatted_string("", "", &arena), "TTAA");
    /// ```
    pub fn crossover(
        &self,
        other: &Self,
        pos1: usize,
        pos2: usize,
        arena: &mut GenomeArena,
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
        seq1_vec.extend_from_slice(&self.sequence(arena)[..pos1]);
        seq1_vec.extend_from_slice(&other.sequence(arena)[pos2..]);
        let data1 = arena.alloc(&seq1_vec);

        // Seq2 New = Other[..pos2] + Self[pos1..]
        let mut seq2_vec = Vec::with_capacity(pos2 + (self.len() - pos1));
        seq2_vec.extend_from_slice(&other.sequence(arena)[..pos2]);
        seq2_vec.extend_from_slice(&self.sequence(arena)[pos1..]);
        let data2 = arena.alloc(&seq2_vec);

        // Map crossover
        let (map1_left, map1_right) = self.map.split_at(pos1).map_err(|e| e.to_string())?;
        let (map2_left, map2_right) = other.map.split_at(pos2).map_err(|e| e.to_string())?;

        let map1_new = map1_left.merge(&map2_right).map_err(|e| e.to_string())?;
        let map2_new = map2_left.merge(&map1_right).map_err(|e| e.to_string())?;

        Ok((
            Self::new(self.id.clone(), data1, map1_new),
            Self::new(other.id.clone(), data2, map2_new),
        ))
    }

    /// Perform gene conversion.
    ///
    /// # Examples
    ///
    /// ```
    /// use centrevo_sim::genome::{Chromosome, RepeatMap};
    /// use centrevo_sim::base::{Nucleotide, GenomeArena};
    ///
    /// let mut arena = GenomeArena::new();
    /// let seq_r: Vec<Nucleotide> = vec![Nucleotide::A; 5]; // Recipient must be mutable in concept, but here we generate new
    /// let seq_d: Vec<Nucleotide> = vec![Nucleotide::T; 5];
    /// let map = RepeatMap::uniform(1, 5, 1);
    /// let data_r = arena.alloc(&seq_r);
    /// let data_d = arena.alloc(&seq_d);
    /// let chr_r = Chromosome::new("R", data_r, map.clone());
    /// let chr_d = Chromosome::new("D", data_d, map.clone());
    ///
    /// // Copy 2 bases ("TT") from donor at 1 to recipient at 1
    /// let res = chr_r.gene_conversion(&chr_d, 1, 1, 2, &mut arena).unwrap();
    /// assert_eq!(res.to_formatted_string("", "", &arena), "ATTAA");
    /// ```
    pub fn gene_conversion(
        &self,
        other: &Self,
        recipient_start: usize,
        donor_start: usize,
        length: usize,
        arena: &mut GenomeArena,
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

        let (map_left, map_mid_right) = self
            .map
            .split_at(recipient_start)
            .map_err(|e| e.to_string())?;
        let (_map_old_tract, map_right) =
            map_mid_right.split_at(length).map_err(|e| e.to_string())?;

        let (_map_d_left, map_d_mid_right) =
            other.map.split_at(donor_start).map_err(|e| e.to_string())?;
        let (map_new_tract, _map_d_right) = map_d_mid_right
            .split_at(length)
            .map_err(|e| e.to_string())?;

        let map_temp = map_left.merge(&map_new_tract).map_err(|e| e.to_string())?;
        let map_final = map_temp.merge(&map_right).map_err(|e| e.to_string())?;

        let final_len = recipient_start + (end2 - donor_start) + (self.len() - end1);
        let mut new_seq_vec = Vec::with_capacity(final_len);
        new_seq_vec.extend_from_slice(&self.sequence(arena)[..recipient_start]);
        new_seq_vec.extend_from_slice(&other.sequence(arena)[donor_start..end2]);
        new_seq_vec.extend_from_slice(&self.sequence(arena)[end1..]);

        let data = arena.alloc(&new_seq_vec);

        Ok(Self::new(self.id.clone(), data, map_final))
    }

    /// Offset the page IDs of the sequence slice.
    /// Used when merging arenas.
    pub fn offset_page_ids(&mut self, offset: u32) {
        self.data.offset_page_id(offset);
    }

    /// Copy the sequence data from a source arena to a destination arena,
    /// and update the slice to point to the new location.
    pub fn localize(
        &mut self,
        source_arena: &crate::base::GenomeArena,
        dest_arena: &mut crate::base::GenomeArena,
    ) {
        let seq_data = self.sequence(source_arena);
        self.data = dest_arena.alloc(seq_data);
    }

    /// Get the underlying sequence slice handle.
    pub fn get_slice(&self) -> crate::base::SequenceSlice {
        self.data
    }

    /// Set the underlying sequence slice handle.
    pub fn set_slice(&mut self, slice: crate::base::SequenceSlice) {
        self.data = slice;
    }
}

impl std::fmt::Display for Chromosome {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Chromosome({})", self.id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_map(ru_len: usize, rus_per_hor: usize, num_hors: usize) -> RepeatMap {
        RepeatMap::uniform(ru_len, rus_per_hor, num_hors)
    }

    #[test]
    fn test_chromosome_new() {
        let mut arena = GenomeArena::new();
        // 16 bp. RU=4. RUs=4. HOR=2 RUs. HORs=2.
        // "ACGT" * 4
        let seq_vec = b"ACGTACGTACGTACGT".to_vec();
        let seq_slice = arena.alloc(
            &seq_vec
                .iter()
                .map(|&b| Nucleotide::from_ascii(b).unwrap())
                .collect::<Vec<_>>(),
        );

        let map = test_map(4, 2, 2);
        let chr = Chromosome::new("chr1", seq_slice, map);

        assert_eq!(chr.id(), "chr1");
        assert_eq!(chr.len(), 16);
        assert_eq!(chr.map().num_rus(), 4);
        assert_eq!(chr.map().num_hors(), 2);
    }

    #[test]
    fn test_uniform() {
        let mut arena = GenomeArena::new();
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 2, 2, 1, &mut arena);
        assert_eq!(chr.len(), 4);
        assert_eq!(chr.to_formatted_string("-", "=", &arena), "AA-AA");
    }

    #[test]
    fn test_from_formatted_string() {
        let mut arena = GenomeArena::new();
        let input = "AC-GT=AA-TT";
        let chr = Chromosome::from_formatted_string("chr1", input, '=', '-', &mut arena).unwrap();
        assert_eq!(chr.len(), 8);
        assert_eq!(chr.to_formatted_string("-", "=", &arena), "AC-GT=AA-TT");
    }

    #[test]
    fn test_crossover() {
        let mut arena = GenomeArena::new();
        let a = Chromosome::uniform("a", Nucleotide::A, 1, 4, 1, &mut arena); // "AAAA"
        let b = Chromosome::uniform("b", Nucleotide::T, 1, 4, 1, &mut arena); // "TTTT"

        let (a1, b1) = a.crossover(&b, 2, 1, &mut arena).unwrap();
        // a1 = a[..2] + b[1..] = "AA" + "TTT" = "AATTT"
        assert_eq!(a1.to_formatted_string("", "", &arena), "AATTT");
        // b1 = b[..1] + a[2..] = "T" + "AA" = "TAA"
        assert_eq!(b1.to_formatted_string("", "", &arena), "TAA");
    }

    #[test]
    fn test_gene_conversion() {
        let mut arena = GenomeArena::new();
        let rec = Chromosome::uniform("r", Nucleotide::A, 1, 4, 1, &mut arena); // "AAAA"
        let donor = Chromosome::uniform("d", Nucleotide::T, 1, 4, 1, &mut arena); // "TTTT"

        // Replace rec[1..3] (2 bases) with donor[1..3]
        let out = rec.gene_conversion(&donor, 1, 1, 2, &mut arena).unwrap();
        // "A" + "TT" + "A" = "ATTA"
        assert_eq!(out.to_formatted_string("", "", &arena), "ATTA");
    }
}
