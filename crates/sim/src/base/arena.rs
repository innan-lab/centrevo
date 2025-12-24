use super::Nucleotide;
use serde::{Deserialize, Serialize};

/// Threshold multiplier for compaction.
/// If total allocated bytes > live bytes * FRAGMENTATION_THRESHOLD, we compact.
const FRAGMENTATION_THRESHOLD: f32 = 2.0;
/// Minimum bytes to trigger compaction (avoid compacting tiny arenas).
const MIN_COMPACTION_BYTES: usize = 1024 * 1024; // 1 MB

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct SequenceSlice {
    pub page_id: u32,
    pub start: u32,
    pub len: u32,
}

impl SequenceSlice {
    pub fn new(page_id: u32, start: u32, len: u32) -> Self {
        Self {
            page_id,
            start,
            len,
        }
    }

    /// Offset the page_id by a given amount.
    /// Used when merging arenas.
    pub fn offset_page_id(&mut self, offset: u32) {
        self.page_id += offset;
    }
}

/// A paged arena storage for nucleotide sequences.
///
/// This replaces individual heap allocations for each sequence with a few large
/// contiguous buffers (pages). This improves cache locality and reduces allocation overhead.
#[derive(Debug, Default, Serialize, Deserialize)]
pub struct GenomeArena {
    /// Pages of contiguous nucleotide storage.
    pages: Vec<Vec<Nucleotide>>,
    /// Total number of bytes allocated across all pages.
    total_allocated: usize,
}

impl GenomeArena {
    pub fn new() -> Self {
        Self {
            pages: vec![Vec::new()], // Start with one page
            total_allocated: 0,
        }
    }

    /// Access the nucleotide data for a given slice.
    ///
    /// # Panics
    /// Panics if the page_id or range is invalid.
    #[inline]
    pub fn get(&self, slice: SequenceSlice) -> &[Nucleotide] {
        &self.pages[slice.page_id as usize]
            [slice.start as usize..(slice.start + slice.len) as usize]
    }

    /// Mutably access the nucleotide data for a given slice.
    ///
    /// # Panics
    /// Panics if the page_id or range is invalid.
    #[inline]
    pub fn get_mut(&mut self, slice: SequenceSlice) -> &mut [Nucleotide] {
        &mut self.pages[slice.page_id as usize]
            [slice.start as usize..(slice.start + slice.len) as usize]
    }

    /// Allocates a new sequence in the current active page.
    /// If the page is too full? For now we just append unlimitedly to the active page
    /// until compaction forces a new page structure.
    ///
    /// Returns the handle (SequenceSlice) to the allocated data.
    pub fn alloc(&mut self, data: &[Nucleotide]) -> SequenceSlice {
        let current_page_idx = self.pages.len() - 1;
        let page = &mut self.pages[current_page_idx];

        let start = page.len() as u32;
        let len = data.len() as u32;

        page.extend_from_slice(data);
        self.total_allocated += data.len();

        SequenceSlice {
            page_id: current_page_idx as u32,
            start,
            len,
        }
    }

    /// Checks if compaction is needed based on fragmentation.
    pub fn needs_compaction(&self, live_bytes: usize) -> bool {
        if self.total_allocated < MIN_COMPACTION_BYTES {
            return false;
        }
        (self.total_allocated as f32) > (live_bytes as f32 * FRAGMENTATION_THRESHOLD)
    }

    /// Creates a NEW compacted arena containing only the data referenced by `live_slices`.
    /// Returns the new arena and a mapping of { old_slice -> new_slice }.
    ///
    /// Note: This function doesn't update the handles in place because they might be
    /// scattered across complex structs throughout the application. Instead, it returns
    /// a new purified Arena. The caller is responsible for mapping handles (or we can implement
    /// a visitor pattern later).
    ///
    /// Actually, to make this work smoothly with the "Handles" design, we probably want
    /// to take a closure that updates the world, OR just return the new Arena and let
    /// the simulation logic handle the swap.
    ///
    /// For this first pass implementation, let's provide a way to "Purify" a list of slices.
    pub fn compact(&self, live_slices: &[SequenceSlice]) -> (Self, Vec<SequenceSlice>) {
        let mut new_arena = Self::new();
        let mut new_slices = Vec::with_capacity(live_slices.len());

        for slice in live_slices {
            let data = self.get(*slice);
            new_slices.push(new_arena.alloc(data));
        }

        (new_arena, new_slices)
    }

    /// Merges another arena into this one.
    /// Returns the page_id offset to apply to slices from the other arena.
    pub fn merge(&mut self, other: GenomeArena) -> u32 {
        let page_offset = self.pages.len() as u32;
        self.pages.extend(other.pages);
        self.total_allocated += other.total_allocated;
        page_offset
    }
}
