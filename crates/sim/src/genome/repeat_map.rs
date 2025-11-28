use serde::{Deserialize, Serialize};
pub use crate::errors::RepeatMapError;

/// Maps the structure of a sequence into Repeat Units (RUs) and Higher-Order Repeats (HORs).
///
/// This structure maintains the boundaries of RUs and HORs relative to a flat sequence.
/// It allows for non-uniform lengths of RUs and HORs.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RepeatMap {
    /// Start indices (offsets) of every Repeat Unit in the sequence.
    ///
    /// - Must always start with 0.
    /// - Must always end with `sequence.len()`.
    /// - Length is `num_rus + 1`.
    /// - `ru_length(i) = ru_offsets[i+1] - ru_offsets[i]`
    ///
    /// # Example:
    /// For a sequence with RUs of lengths [100, 150, 120],
    /// the `ru_offsets` would be [0, 100, 250, 370].
    /// This indicates:
    /// - RU 0: 0 to 100, length of 100
    /// - RU 1: 100 to 250, length of 150
    /// - RU 2: 250 to 370, length of 120
    /// - Total sequence length: 370
    ru_offsets: Vec<usize>,

    /// Indices into `ru_offsets` that mark the start of each Higher-Order Repeat.
    ///
    /// - Must always start with 0.
    /// - Must always end with `num_rus` (i.e., `ru_offsets.len() - 1`).
    /// - Length is `num_hors + 1`.
    /// - `rus_in_hor(j) = hor_offsets[j+1] - hor_offsets[j]`
    ///
    /// # Example:
    /// For a sequence with HORs containing [2, 1] RUs respectively,
    /// the `hor_offsets` would be [0, 2, 3].
    /// This indicates:
    /// - HOR 0: RUs 0 to 2 (i.e., RUs 0 and 1)
    /// - HOR 1: RUs 2 to 3 (i.e., RU 2)
    /// - Total number of RUs: 3
    ///
    /// Combining with the previous example, this means:
    /// - HOR 0 spans RUs 0 and 1, covering sequence indices 0 to 250.
    /// - HOR 1 spans RU 2, covering sequence indices 250 to 370.
    hor_idx_offsets: Vec<usize>,
}

impl RepeatMap {
    /// Create a new RepeatMap from raw offsets.
    ///
    /// # Arguments
    /// * `ru_offsets` - Start indices of RUs. Must start with 0. Last element is total sequence length.
    /// * `hor_offsets` - Indices into `ru_offsets` for HOR starts. Must start with 0. Last element is number of RUs.
    pub fn new(ru_offsets: Vec<usize>, hor_idx_offsets: Vec<usize>) -> Result<Self, RepeatMapError> {
        // Validate RU offsets
        if ru_offsets.is_empty() || ru_offsets[0] != 0 {
            return Err(RepeatMapError::InvalidOffsets);
        }
        for i in 0..ru_offsets.len() - 1 {
            if ru_offsets[i] > ru_offsets[i+1] {
                return Err(RepeatMapError::InvalidOffsets);
            }
        }

        // Validate HOR offsets
        let num_rus = ru_offsets.len() - 1;
        if hor_idx_offsets.is_empty() || hor_idx_offsets[0] != 0 {
            return Err(RepeatMapError::InvalidHorOffsets);
        }
        if *hor_idx_offsets.last().unwrap() != num_rus {
            return Err(RepeatMapError::InvalidHorOffsets);
        }
        for i in 0..hor_idx_offsets.len() - 1 {
            if hor_idx_offsets[i] > hor_idx_offsets[i+1] {
                return Err(RepeatMapError::InvalidHorOffsets);
            }
        }

        Ok(Self {
            ru_offsets,
            hor_idx_offsets,
        })
    }

    /// Create a uniform RepeatMap.
    ///
    /// # Arguments
    /// * `ru_length` - Length of each RU.
    /// * `rus_per_hor` - Number of RUs per HOR.
    /// * `num_hors` - Total number of HORs.
    pub fn uniform(ru_length: usize, rus_per_hor: usize, num_hors: usize) -> Self {
        let num_rus = rus_per_hor * num_hors;
        let _total_len = ru_length * num_rus;

        // Always sorted offsets
        let ru_offsets: Vec<usize> = (0..=num_rus).map(|i| i * ru_length).collect();
        let hor_offsets: Vec<usize> = (0..=num_hors).map(|i| i * rus_per_hor).collect();

        Self {
            ru_offsets,
            hor_idx_offsets: hor_offsets,
        }
    }

    /// Get the total length of the sequence mapped by this structure.
    pub fn total_length(&self) -> usize {
        *self.ru_offsets.last().unwrap()
    }

    /// Get the number of Repeat Units.
    pub fn num_rus(&self) -> usize {
        self.ru_offsets.len() - 1
    }

    /// Get the number of Higher-Order Repeats.
    pub fn num_hors(&self) -> usize {
        self.hor_idx_offsets.len() - 1
    }

    /// Get the start and end indices (range) of the RU at `index`.
    pub fn get_ru_interval(&self, index: usize) -> Option<(usize, usize)> {
        if index >= self.num_rus() {
            None
        } else {
            Some((self.ru_offsets[index], self.ru_offsets[index + 1]))
        }
    }

    /// Get the start and end indices (range) of the HOR at `index`.
    /// Returns range in sequence coordinates.
    pub fn get_hor_interval(&self, index: usize) -> Option<(usize, usize)> {
        if index >= self.num_hors() {
            None
        } else {
            let start_ru_idx = self.hor_idx_offsets[index];
            let end_ru_idx = self.hor_idx_offsets[index + 1];
            Some((self.ru_offsets[start_ru_idx], self.ru_offsets[end_ru_idx]))
        }
    }

    /// Find which RU contains the given sequence position.
    /// Returns the RU index.
    pub fn find_ru_index(&self, position: usize) -> Option<usize> {
        if position >= self.total_length() {
            return None;
        }

        // Binary search for the position
        match self.ru_offsets.binary_search(&position) {
            Ok(idx) => Some(idx), // Exact match (start of RU)
            Err(idx) => Some(idx - 1), // Inside RU
        }
    }

    /// Find which HOR contains the given RU index.
    /// Returns the HOR index.
    pub fn find_hor_index(&self, ru_index: usize) -> Option<usize> {
        if ru_index >= self.num_rus() {
            return None;
        }
        match self.hor_idx_offsets.binary_search(&ru_index) {
            Ok(idx) => Some(idx), // Exact match (start of HOR)
            Err(idx) => Some(idx - 1), // Inside HOR
        }
    }

    /// Find the RU and HOR indices for a given sequence position.
    /// Returns (ru_index, hor_index).
    ///
    /// This is a convenience method combining `find_ru_index` and `find_hor_index`.
    pub fn find_indices(&self, position: usize) -> Option<(usize, usize)> {
        let ru_idx = self.find_ru_index(position)?;
        let hor_idx = self.find_hor_index(ru_idx)?;
        Some((ru_idx, hor_idx))
    }

    /// Split the map at a given sequence position.
    /// Position is 0-based and can be at RU boundaries or inside RUs.
    ///
    /// If the position is inside an RU, that RU is split into two.
    /// Returns (left_map, right_map).
    pub fn split_at(&self, position: usize) -> Result<(Self, Self), RepeatMapError> {
        if position > self.total_length() {
            return Err(RepeatMapError::IndexOutOfBounds);
        }

        // Handle edge cases
        if position == 0 {
            return Ok((
                Self::new(vec![0], vec![0])?, // Empty map
                self.clone()
            ));
        }

        if position == self.total_length() {
            return Ok((
                self.clone(),
                Self::new(vec![0], vec![0])? // Empty map
            ));
        }

        // Find which RU contains or starts at the split position
        let split_ru_idx = self.find_ru_index(position).ok_or(RepeatMapError::IndexOutOfBounds)?;

        // Check if we are splitting exactly at an RU boundary
        let is_boundary = self.ru_offsets[split_ru_idx] == position;

        // --- Construct Left Map RU offsets ---
        let left_ru_offsets = if is_boundary {
            // Split at boundary: include all RUs before the split point
            // Example: split at 100 with offsets [0, 100, 200]
            // split_ru_idx=1, we want [0, 100]
            self.ru_offsets[0..=split_ru_idx].to_vec()
        } else {
            // Split inside RU: include RU start positions up to and including split_ru_idx,
            // but replace the end boundary with the split position
            // Example: split at 50 with offsets [0, 100, 200]
            // split_ru_idx=0, we take [0, 100] and change to [0, 50]
            let mut offsets = self.ru_offsets[0..=split_ru_idx+1].to_vec();
            *offsets.last_mut().unwrap() = position;
            offsets
        };

        // --- Construct Right Map RU offsets ---
        // Both boundary and interior splits use the same logic:
        // - Boundary: RUs from split_ru_idx+1 onward (split_ru_idx is fully in left)
        // - Interior: RUs from split_ru_idx+1 onward (split_ru_idx is partially in both)
        // In both cases, we rebase by subtracting the split position
        let mut right_ru_offsets = vec![0];
        for &offset in &self.ru_offsets[split_ru_idx+1..] {
            right_ru_offsets.push(offset - position);
        }

        // --- Construct Left Map HOR offsets ---
        let left_ru_count = left_ru_offsets.len() - 1;
        let mut left_hor_offsets = vec![0];

        // Include all HOR boundaries that fall within the left map
        // When splitting at a boundary: include HOR boundaries < split_ru_idx
        // When splitting inside an RU: include HOR boundaries <= split_ru_idx
        for &h_off in &self.hor_idx_offsets[1..] {
            let should_include = if is_boundary {
                h_off < split_ru_idx
            } else {
                h_off <= split_ru_idx
            };

            if should_include {
                left_hor_offsets.push(h_off);
            } else {
                break;
            }
        }

        // Ensure the last HOR offset equals the number of RUs in left map
        if *left_hor_offsets.last().unwrap() != left_ru_count {
            left_hor_offsets.push(left_ru_count);
        }

        // --- Construct Right Map HOR offsets ---
        let right_ru_count = right_ru_offsets.len() - 1;
        let mut right_hor_offsets = vec![0];

        // Map HOR boundaries that fall in the right map
        // When splitting at a boundary: include HOR boundaries > split_ru_idx, rebased
        // When splitting inside an RU: include HOR boundaries > split_ru_idx, rebased
        // (Both cases are the same because split_ru_idx marks the last RU in left side)
        for &h_off in &self.hor_idx_offsets[1..] {
            if h_off > split_ru_idx {
                let rebased_idx = h_off - split_ru_idx;
                if rebased_idx <= right_ru_count {
                    right_hor_offsets.push(rebased_idx);
                }
            }
        }

        // Ensure the last HOR offset equals the number of RUs in right map
        if *right_hor_offsets.last().unwrap() != right_ru_count {
            right_hor_offsets.push(right_ru_count);
        }

        Ok((
            Self::new(left_ru_offsets, left_hor_offsets)?,
            Self::new(right_ru_offsets, right_hor_offsets)?
        ))
    }

    /// Merge another map onto the end of this one.
    pub fn merge(&self, other: &Self) -> Result<Self, RepeatMapError> {
        let current_len = self.total_length();
        let current_ru_count = self.num_rus();

        // 1. Merge RU offsets
        let mut new_ru_offsets = self.ru_offsets.clone();
        // Skip the first 0 of the other map
        for &offset in &other.ru_offsets[1..] {
            new_ru_offsets.push(offset + current_len);
        }

        // 2. Merge HOR offsets
        let mut new_hor_offsets = self.hor_idx_offsets.clone();
        // Skip the first 0 of the other map
        for &offset in &other.hor_idx_offsets[1..] {
            new_hor_offsets.push(offset + current_ru_count);
        }

        Self::new(new_ru_offsets, new_hor_offsets)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_uniform_creation() {
        let map = RepeatMap::uniform(10, 2, 2);
        // RUs: 4 total (2 per HOR * 2 HORs). Length 10 each. Total 40.
        assert_eq!(map.total_length(), 40);
        assert_eq!(map.num_rus(), 4);
        assert_eq!(map.num_hors(), 2);

        assert_eq!(map.ru_offsets, vec![0, 10, 20, 30, 40]);
        assert_eq!(map.hor_idx_offsets, vec![0, 2, 4]);
    }

    #[test]
    fn test_split_at_boundary() {
        let map = RepeatMap::uniform(10, 2, 2); // [0, 10, 20, 30, 40]
        // Split at 20 (End of RU 1, Start of RU 2. Also End of HOR 0).
        let (left, right) = map.split_at(20).unwrap();

        // Left: [0, 10, 20] -> 2 RUs. HORs: [0, 2] -> 1 HOR.
        assert_eq!(left.ru_offsets, vec![0, 10, 20]);
        assert_eq!(left.hor_idx_offsets, vec![0, 2]);

        // Right: [0, 10, 20] (Rebased from 20, 30, 40). 2 RUs. HORs: [0, 2].
        assert_eq!(right.ru_offsets, vec![0, 10, 20]);
        assert_eq!(right.hor_idx_offsets, vec![0, 2]);
    }

    #[test]
    fn test_split_inside_ru() {
        let map = RepeatMap::uniform(10, 1, 2); // [0, 10, 20]. HORs [0, 1, 2].
        // Split at 5 (Inside RU 0).
        let (left, right) = map.split_at(5).unwrap();

        // Left: [0, 5]. 1 RU. HORs: [0, 1].
        assert_eq!(left.ru_offsets, vec![0, 5]);
        assert_eq!(left.hor_idx_offsets, vec![0, 1]);

        // Right: [0, 5, 15]. (Rebased from 5, 10, 20).
        // Original RU 0 (5-10) becomes Right RU 0 (0-5).
        // Original RU 1 (10-20) becomes Right RU 1 (5-15).
        assert_eq!(right.ru_offsets, vec![0, 5, 15]);

        // Right HORs:
        // Original HORs: [0, 1, 2].
        // Split inside RU 0.
        // Right starts with partial RU 0.
        // Original boundary 1 (start of RU 1) is > 0.
        // New index = 1 - 0 = 1.
        // Original boundary 2 (end) is > 0.
        // New index = 2 - 0 = 2.
        // Result: [0, 1, 2].
        assert_eq!(right.hor_idx_offsets, vec![0, 1, 2]);
    }

    #[test]
    fn test_merge() {
        let map1 = RepeatMap::uniform(10, 1, 1); // [0, 10], [0, 1]
        let map2 = RepeatMap::uniform(20, 1, 1); // [0, 20], [0, 1]

        let merged = map1.merge(&map2).unwrap();

        // RUs: 10, 20. Offsets: [0, 10, 30].
        assert_eq!(merged.ru_offsets, vec![0, 10, 30]);

        // HORs: 1 from map1, 1 from map2. Total 2.
        // Offsets: [0, 1, 2].
        assert_eq!(merged.hor_idx_offsets, vec![0, 1, 2]);
    }

    #[test]
    fn test_split_at_edges() {
        let map = RepeatMap::uniform(10, 2, 2); // [0, 10, 20, 30, 40]

        // Split at 0 - left is empty
        let (left, right) = map.split_at(0).unwrap();
        assert_eq!(left.ru_offsets, vec![0]);
        assert_eq!(left.hor_idx_offsets, vec![0]);
        assert_eq!(right.ru_offsets, vec![0, 10, 20, 30, 40]);

        // Split at end - right is empty
        let (left, right) = map.split_at(40).unwrap();
        assert_eq!(left.ru_offsets, vec![0, 10, 20, 30, 40]);
        assert_eq!(right.ru_offsets, vec![0]);
        assert_eq!(right.hor_idx_offsets, vec![0]);
    }

    #[test]
    fn test_split_multiple_hors() {
        let map = RepeatMap::uniform(10, 2, 3); // 3 HORs, 2 RUs each
        // RUs: [0, 10, 20, 30, 40, 50, 60]
        // HORs: [0, 2, 4, 6]

        // Split inside second HOR at RU 3
        let (left, right) = map.split_at(35).unwrap();

        // Left: RU 0, 1, 2, partial 3
        assert_eq!(left.ru_offsets, vec![0, 10, 20, 30, 35]);
        // HORs: HOR 0 complete (RUs 0-1), HOR 1 complete (RUs 2-3 but 3 is partial)
        assert_eq!(left.hor_idx_offsets, vec![0, 2, 4]);

        // Right: partial RU 3, RU 4, 5
        assert_eq!(right.ru_offsets, vec![0, 5, 15, 25]);
        // HORs: partial HOR 1 (partial RU 3), HOR 2 (RUs 4-5)
        assert_eq!(right.hor_idx_offsets, vec![0, 1, 3]);
    }

    #[test]
    fn test_split_at_hor_boundary() {
        let map = RepeatMap::uniform(10, 2, 3); // 3 HORs, 2 RUs each
        // RUs: [0, 10, 20, 30, 40, 50, 60]
        // HORs: [0, 2, 4, 6]

        // Split exactly at HOR boundary (position 20 = start of HOR 1)
        let (left, right) = map.split_at(20).unwrap();

        // Left: HOR 0 (RUs 0-1)
        assert_eq!(left.ru_offsets, vec![0, 10, 20]);
        assert_eq!(left.hor_idx_offsets, vec![0, 2]);

        // Right: HOR 1, HOR 2
        assert_eq!(right.ru_offsets, vec![0, 10, 20, 30, 40]);
        assert_eq!(right.hor_idx_offsets, vec![0, 2, 4]);
    }

    #[test]
    fn test_split_roundtrip() {
        let map = RepeatMap::uniform(10, 2, 3);
        let original_len = map.total_length();

        // Split at boundary and merge should preserve structure
        let (left, right) = map.split_at(20).unwrap();
        let merged = left.merge(&right).unwrap();

        assert_eq!(merged.total_length(), original_len);
        assert_eq!(merged.num_rus(), map.num_rus());
        // When splitting at boundary, we preserve the HOR count
        assert_eq!(merged.num_hors(), map.num_hors());

        // Split inside RU creates an extra RU and HOR
        let (left2, right2) = map.split_at(25).unwrap();
        let merged2 = left2.merge(&right2).unwrap();

        assert_eq!(merged2.total_length(), original_len);
        // Splitting inside an RU increases the total RU count by 1
        assert_eq!(merged2.num_rus(), map.num_rus() + 1);
        // HORs may also increase when split inside
        assert!(merged2.num_hors() >= map.num_hors());
    }

    #[test]
    fn test_merge_split_converse() {
        // Test that merge correctly reconstructs sequence positions after split
        let map = RepeatMap::uniform(10, 2, 3); // [0, 10, 20, 30, 40, 50, 60]

        // Boundary split
        let (left, right) = map.split_at(20).unwrap();
        let merged = left.merge(&right).unwrap();

        // Verify all RU intervals match
        for i in 0..map.num_rus() {
            assert_eq!(map.get_ru_interval(i), merged.get_ru_interval(i));
        }

        // Verify all HOR intervals match
        for i in 0..map.num_hors() {
            assert_eq!(map.get_hor_interval(i), merged.get_hor_interval(i));
        }

        // Interior split creates a new boundary, but sequence coverage should be continuous
        let (left2, right2) = map.split_at(25).unwrap();
        let merged2 = left2.merge(&right2).unwrap();

        // Total length preserved
        assert_eq!(merged2.total_length(), map.total_length());

        // Verify split point is now a boundary
        assert_eq!(left2.total_length(), 25);
        assert_eq!(merged2.ru_offsets[left2.num_rus()], 25);

        // Verify continuity: left's last RU ends where right's first RU begins (after shift)
        assert_eq!(left2.total_length(), 25);
        let right_first_ru_len = right2.ru_offsets[1] - right2.ru_offsets[0];
        // Original RU 2 was [20, 30], split at 25
        // Left part: [20, 25] (length 5)
        // Right part: [0, 5] (length 5) in right map
        assert_eq!(right_first_ru_len, 5);
        assert_eq!(merged2.get_ru_interval(2), Some((20, 25))); // Left part of split RU
        assert_eq!(merged2.get_ru_interval(3), Some((25, 30))); // Right part of split RU
    }
}
