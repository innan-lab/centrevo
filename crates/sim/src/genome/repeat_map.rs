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
    ru_offsets: Vec<usize>,

    /// Indices into `ru_offsets` that mark the start of each Higher-Order Repeat.
    /// 
    /// - Must always start with 0.
    /// - Must always end with `num_rus` (i.e., `ru_offsets.len() - 1`).
    /// - Length is `num_hors + 1`.
    /// - `rus_in_hor(j) = hor_offsets[j+1] - hor_offsets[j]`
    hor_offsets: Vec<usize>,
}

impl RepeatMap {
    /// Create a new RepeatMap from raw offsets.
    ///
    /// # Arguments
    /// * `ru_offsets` - Start indices of RUs. Must start with 0. Last element is total sequence length.
    /// * `hor_offsets` - Indices into `ru_offsets` for HOR starts. Must start with 0. Last element is number of RUs.
    pub fn new(ru_offsets: Vec<usize>, hor_offsets: Vec<usize>) -> Result<Self, RepeatMapError> {
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
        if hor_offsets.is_empty() || hor_offsets[0] != 0 {
            return Err(RepeatMapError::InvalidHorOffsets);
        }
        if *hor_offsets.last().unwrap() != num_rus {
            return Err(RepeatMapError::InvalidHorOffsets);
        }
        for i in 0..hor_offsets.len() - 1 {
            if hor_offsets[i] > hor_offsets[i+1] {
                return Err(RepeatMapError::InvalidHorOffsets);
            }
        }

        Ok(Self {
            ru_offsets,
            hor_offsets,
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

        let ru_offsets: Vec<usize> = (0..=num_rus).map(|i| i * ru_length).collect();
        let hor_offsets: Vec<usize> = (0..=num_hors).map(|i| i * rus_per_hor).collect();

        Self {
            ru_offsets,
            hor_offsets,
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
        self.hor_offsets.len() - 1
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
            let start_ru_idx = self.hor_offsets[index];
            let end_ru_idx = self.hor_offsets[index + 1];
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

    /// Split the map at a given sequence position.
    /// 
    /// If the position is inside an RU, that RU is split into two.
    /// Returns (left_map, right_map).
    pub fn split_at(&self, position: usize) -> Result<(Self, Self), RepeatMapError> {
        if position > self.total_length() {
            return Err(RepeatMapError::IndexOutOfBounds);
        }

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

        // 1. Find split point in RUs
        let split_ru_idx = self.find_ru_index(position).ok_or(RepeatMapError::IndexOutOfBounds)?;
        
        // Check if we are splitting exactly at an RU boundary
        let is_boundary = self.ru_offsets[split_ru_idx] == position;
        
        // --- Construct Left Map ---
        let mut left_ru_offsets = self.ru_offsets[0..=split_ru_idx].to_vec();
        if !is_boundary {
            // If not at boundary, we include the split point as the end of the last RU
            left_ru_offsets.push(position);
        }
        
        // Fix left map last element if needed (it should be `position`)
        if *left_ru_offsets.last().unwrap() != position {
             // This happens if find_ru_index returned the RU *starting* at position
             // But we want the split to happen *before* it?
             // Let's re-evaluate find_ru_index behavior for exact matches.
             // If pos=100, offsets=[0, 100, 200].
             // binary_search(100) -> Ok(1). find_ru_index -> 1.
             // We want left to be [0, 100].
             // slice [0..=1] is [0, 100]. Correct.
             //
             // If pos=50. offsets=[0, 100].
             // binary_search(50) -> Err(1). find_ru_index -> 0.
             // slice [0..=0] is [0].
             // We push 50. -> [0, 50]. Correct.
        }

        // --- Construct Right Map ---
        // We need to rebase offsets starting from 0.
        let mut right_ru_offsets = Vec::new();
        right_ru_offsets.push(0);
        
        if !is_boundary {
            // If we split inside RU 0 (0-100) at 50.
            // Right part of RU 0 is 50-100. Length 50.
            // Right map starts with 0. Next offset should be 50.
            // Original next offset was 100. 100 - 50 = 50.
            // So we take offsets from split_ru_idx + 1 onwards.
            for &offset in &self.ru_offsets[split_ru_idx+1..] {
                right_ru_offsets.push(offset - position);
            }
        } else {
            // Split at 100. RU 1 starts at 100.
            // Right map should start with RU 1.
            // Offsets [100, 200]. Rebased: [0, 100].
            // split_ru_idx was 1.
            // We take from 1 onwards?
            // If split at 100, find_ru_index(100) = 1.
            // We want offsets[1..].
            // But we already pushed 0. So we want offsets[2..] rebased?
            // ru_offsets[1] is 100. 100-100=0. Already pushed.
            // ru_offsets[2] is 200. 200-100=100.
            // So we want split_ru_idx + 1.
            for &offset in &self.ru_offsets[split_ru_idx+1..] {
                right_ru_offsets.push(offset - position);
            }
        }

        // --- Handle HORs ---
        // This is tricky. We need to find which HORs belong to left and right.
        // And handle the split RU (which belongs to both? or split creates new RUs?)
        // If we split an RU, we increase the total number of RUs by 1 (one partial on left, one partial on right).
        
        // Let's reconstruct HOR offsets based on RU counts.
        // Left RUs count: left_ru_offsets.len() - 1.
        // Right RUs count: right_ru_offsets.len() - 1.
        
        // We need to find where the split RU falls in HORs.
        // split_ru_idx is the index of the RU being split (or the one starting at split point).
        
        // Find HOR containing split_ru_idx
        // hor_offsets: [0, 5, 10].
        // If split_ru_idx = 2. Inside HOR 0.
        // Left HORs: HOR 0 (partial).
        // Right HORs: HOR 0 (partial) + HOR 1.
        
        // We will simplify: The HOR structure is preserved as much as possible.
        // If an HOR is split, it becomes two HORs (one on left, one on right).
        
        // Left HOR offsets:
        let mut left_hor_offsets = Vec::new();
        left_hor_offsets.push(0);
        
        // Right HOR offsets:
        let mut right_hor_offsets = Vec::new();
        right_hor_offsets.push(0);
        
        // Iterate through original HORs
        
        // Let's rebuild HORs for Left Map
        // We want to include all HOR boundaries that are <= split_ru_idx_for_hor?
        // Actually, simpler:
        // The left map has `left_ru_count` RUs.
        // We just need to map the original HOR boundaries to these RUs.
        // If an HOR boundary was at RU 5, and we split at RU 2 (inside), 
        // RU 2 becomes RU 2a (left) and RU 2b (right).
        // Left map has RUs 0, 1, 2a. (3 RUs).
        // Original HOR boundary at 0 is kept.
        // Original HOR boundary at 5 is > 3. Not in left map.
        // But we must close the last HOR in left map.
        
        // Strategy:
        // 1. Filter `hor_offsets` for those <= split_ru_idx.
        // 2. Add them to `left_hor_offsets`.
        // 3. If the last added offset is not equal to `left_ru_count`, add `left_ru_count`.
        
        let left_ru_count = left_ru_offsets.len() - 1;
        for &h_off in &self.hor_offsets {
            if h_off <= split_ru_idx {
                if *left_hor_offsets.last().unwrap() != h_off {
                    left_hor_offsets.push(h_off);
                }
            } else {
                break;
            }
        }
        if *left_hor_offsets.last().unwrap() != left_ru_count {
            left_hor_offsets.push(left_ru_count);
        }
        
        // Strategy for Right Map:
        // 1. We need to map original HOR boundaries to the new right RUs.
        // 2. The first HOR in right map starts at 0.
        // 3. Subsequent HORs correspond to original boundaries > split_ru_idx.
        // 4. We need to adjust indices.
        //    If split was inside RU 2. Right starts with RU 2b.
        //    Original RU 3 becomes Right RU 1.
        //    Original RU index `k` becomes `k - split_ru_idx`.
        
        let right_ru_count = right_ru_offsets.len() - 1;
        
        // Find the first HOR boundary that is > split_ru_idx.
        // But wait, the partial HOR at the start of Right Map needs to be accounted for.
        // It starts at 0.
        
        for &h_off in &self.hor_offsets {
            if h_off > split_ru_idx {
                // This is a boundary in the right part.
                // Calculate new index.
                // If split at RU 2 (inside). split_ru_idx=0? No, 2.
                // Right RUs: 2b, 3, 4...
                // RU 2 is index 0 in right map.
                // RU 3 is index 1.
                // RU k is index k - 2.
                let new_idx = h_off - split_ru_idx;
                right_hor_offsets.push(new_idx);
            }
        }
        
        // Ensure it ends with right_ru_count
        if *right_hor_offsets.last().unwrap() != right_ru_count {
            right_hor_offsets.push(right_ru_count);
        }
        
        // Special case correction for boundary split
        // If split at 100 (RU 1 start). split_ru_idx=1.
        // Left RUs: 0. (Count 1).
        // Right RUs: 1, 2... (Count N-1).
        // Right RU 0 corresponds to Original RU 1.
        // Offset mapping: new = old - 1.
        // My logic `h_off - split_ru_idx` works (h_off - 1).
        
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
        let mut new_hor_offsets = self.hor_offsets.clone();
        // Skip the first 0 of the other map
        for &offset in &other.hor_offsets[1..] {
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
        assert_eq!(map.hor_offsets, vec![0, 2, 4]);
    }

    #[test]
    fn test_split_at_boundary() {
        let map = RepeatMap::uniform(10, 2, 2); // [0, 10, 20, 30, 40]
        // Split at 20 (End of RU 1, Start of RU 2. Also End of HOR 0).
        let (left, right) = map.split_at(20).unwrap();
        
        // Left: [0, 10, 20] -> 2 RUs. HORs: [0, 2] -> 1 HOR.
        assert_eq!(left.ru_offsets, vec![0, 10, 20]);
        assert_eq!(left.hor_offsets, vec![0, 2]);
        
        // Right: [0, 10, 20] (Rebased from 20, 30, 40). 2 RUs. HORs: [0, 2].
        assert_eq!(right.ru_offsets, vec![0, 10, 20]);
        assert_eq!(right.hor_offsets, vec![0, 2]);
    }

    #[test]
    fn test_split_inside_ru() {
        let map = RepeatMap::uniform(10, 1, 2); // [0, 10, 20]. HORs [0, 1, 2].
        // Split at 5 (Inside RU 0).
        let (left, right) = map.split_at(5).unwrap();
        
        // Left: [0, 5]. 1 RU. HORs: [0, 1].
        assert_eq!(left.ru_offsets, vec![0, 5]);
        assert_eq!(left.hor_offsets, vec![0, 1]);
        
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
        assert_eq!(right.hor_offsets, vec![0, 1, 2]);
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
        assert_eq!(merged.hor_offsets, vec![0, 1, 2]);
    }
}
