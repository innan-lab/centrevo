# Refactoring Plan: Non-Uniform Repeat Structures

## 1. Overview
The goal is to transition `Centrevo` from a rigid, uniform repeat structure (fixed RU length, fixed RUs per HOR) to a flexible, non-uniform structure. This allows for the simulation of realistic evolutionary events like unequal crossover and indel mutations that alter the lengths of individual Repeat Units (RUs) and Higher-Order Repeats (HORs).

We will maintain the high-performance flat `Sequence` (vector of nucleotides) and introduce a metadata structure (`RepeatMap`) to track the boundaries of RUs and HORs.

## 2. Data Structures

### 2.1. The `RepeatMap` Struct
This new struct will live in `src/genome/repeat_map.rs` (new file). It acts as an index for the flat sequence.

```rust
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
```

### 2.2. Updated `Chromosome` Struct
The `Chromosome` struct in `src/genome/chromosome.rs` will be updated to replace the scalar length fields with the map.

```rust
pub struct Chromosome {
    id: Arc<str>,
    sequence: Sequence,

    // REPLACED: ru_length: usize,
    // REPLACED: rus_per_hor: usize,

    // NEW:
    map: RepeatMap,
}
```

## 3. Flexibility Analysis

This structure fully supports the user's requirements:

1.  **Non-uniform RU lengths**:
    *   Since `ru_offsets` stores explicit start points, the distance between `ru_offsets[i]` and `ru_offsets[i+1]` can be arbitrary.
    *   *Example*: `ru_offsets = [0, 100, 205, 300]` -> RU lengths are 100, 105, 95.

2.  **Non-uniform RUs per HOR**:
    *   Since `hor_offsets` stores indices into the RU list, the number of RUs in a HOR is flexible.
    *   *Example*: `hor_offsets = [0, 5, 12]` -> HOR 1 has 5 RUs (indices 0-4), HOR 2 has 7 RUs (indices 5-11).

3.  **Non-uniform HOR lengths**:
    *   Derived naturally from the sum of RU lengths within the HOR span.

## 4. Implementation Steps

### Step 1: Create `RepeatMap` Logic
*   Implement `src/genome/repeat_map.rs`.
*   Add methods:
    *   `new(ru_offsets, hor_offsets)`: With validation.
    *   `uniform(total_len, ru_len, rus_per_hor)`: For backward compatibility/initialization.
    *   `split_at(seq_index)`: Crucial for recombination. Returns two `RepeatMap`s.
    *   `merge(other)`: Joins two maps (adjusting offsets).
    *   `shift(amount)`: Handles insertions/deletions.

### Step 2: Refactor `Chromosome`
*   Update struct definition.
*   Update `Chromosome::uniform` to use `RepeatMap::uniform`.
*   Update accessors:
    *   `ru_length()` -> `average_ru_length()` or remove.
    *   `num_hors()` -> delegates to map.
*   Update `to_formatted_string` to iterate through the map.

### Step 3: Update Recombination Logic (`src/evolution/recombination.rs`)
*   **Crossover**:
    *   Currently assumes uniform length.
    *   New logic:
        1.  Perform sequence crossover (existing).
        2.  Find crossover point in `map`.
        3.  Split `map` at that point.
        4.  Merge the left map of Parent A with the right map of Parent B.
        5.  *Note*: If crossover is inside an RU, we need to decide whether to split the RU (add a new boundary) or keep the boundary (resulting in a hybrid RU length). **Decision: Split the RU.** This preserves the "flat sequence" truth.

### Step 4: Update Python Bindings
*   The Python API expects `ru_length` properties.
*   Update `PyChromosome` to expose `num_hors`, `num_rus`.
*   Maybe expose `get_ru_len(index)` for inspection.

### Step 5: Serialization
*   Ensure `RepeatMap` derives `Serialize` and `Deserialize` so checkpoints still work.

## 5. Example Scenarios

### Scenario A: Initialization
User requests: `ru_len=10`, `rus_per_hor=2`, `hors=2`. Total len=40.
*   `ru_offsets`: `[0, 10, 20, 30, 40]`
*   `hor_offsets`: `[0, 2, 4]`

### Scenario B: Deletion
Deletion of 5bp at index 15 (inside RU 1).
*   Sequence shrinks to 35bp.
*   `ru_offsets` update:
    *   Indices <= 15: Unchanged (`0, 10`).
    *   Indices > 15: Decrement by 5 (`20->15`, `30->25`, `40->35`).
*   Result `ru_offsets`: `[0, 10, 15, 25, 35]`.
*   RU 1 length becomes 5 (was 10).

### Scenario C: Crossover
Parent A (40bp) crosses with Parent B (40bp) at index 20.
*   Parent A Map Left (0-20): `ru_offsets=[0, 10, 20]`, `hor_offsets=[0, 2]` (Partial HOR).
*   Parent B Map Right (20-40): `ru_offsets=[20, 30, 40]`.
    *   *Rebase B*: Subtract 20 -> `[0, 10, 20]`.
*   New Map: Merge A-Left and Rebased-B-Right.
    *   `ru_offsets`: `[0, 10, 20]` + `[30, 40]` (shifted by 20) -> `[0, 10, 20, 30, 40]`.
