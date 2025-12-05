# Fitness Function Recommendations for HOR-Only Centromere Simulations

This document summarizes **existing** fitness functions in `selection.rs` and proposes **additional** functions tailored to simulations where:

- Only the **active centromere HOR array** (as a sequence of repeat units, RUs) is modeled.
- Flanking/pericentromeric context is **not** simulated.
- Fitness computations must be **very fast** (ideally 1–2 linear passes with simple integer / floating-point operations), even when **RUs mutate every generation**.

Each entry is described with:

- **Biological idea**: What biological signal or constraint it is approximating.
- **Implementation sketch**: How to compute it efficiently in this codebase.

---

## 1. Existing Fitness Functions

### 1.1 `GCContentFitness`

**Biological idea**

- Models selection on **global GC content** of the simulated region.
- Useful when GC content correlates with features like **DNA stability**, **chromatin state**, or **mutational biases** (e.g., GC-biased gene conversion).
- In a HOR-only centromere simulation, it can approximate a preference for arrays that maintain a particular **base composition** that supports proper chromatin structure.

**Implementation sketch**

- Let the haplotype sequence have length $L$ (in bases), with $G$ and $C$ counts.
- Compute GC fraction $g = (G + C) / L$.
- Define an optimal GC fraction $g^*$ and a **concentration/strength** parameter $\kappa$.
- Use a simple bell-shaped function around $g^*$, e.g. a Beta-like or quadratic penalty:

  - Example 1 (quadratic):
    - $f = \exp(-\kappa (g - g^*)^2)$
  - Example 2 (Beta-like shape, as in current code):
    - Convert $(g^*, \kappa)$ to $(\alpha, \beta)$ and compute a scaled density-like score.

- Complexity: **one linear pass** over bases, a few float ops.
- Implemented in `selection.rs` as a `HaplotypeFitness`:
  - Counts bases.
  - Computes GC fraction.
  - Applies a smooth penalty around the optimum.

---

### 1.2 `LengthFitness`

**Biological idea**

- Models selection on **array size / copy number**, via the total length of the sequence (or RU count, if interpreted that way).
- There is often an **optimal centromere size**: too short → insufficient kinetochore, too long → replication and segregation problems.
- In HOR-only simulations, this can directly target the **number of HOR copies** or **total base length** of the active centromere region.

**Implementation sketch**

- Let $L$ be the sequence length (or RU count).
- Choose an **optimal length** $L^*$ and a **spread** parameter $\sigma$ in log-space.
- Convert to log-space $x = \ln L$, $x^* = \ln L^*$, and penalize deviations:

  - $f = \exp\left( -\frac{(x - x^*)^2}{2\sigma^2} \right)$

- Complexity: constant time after retrieving $L$.
- Implemented in `selection.rs` as a `HaplotypeFitness`:
  - Reads length from haplotype.
  - Applies a log-normal-like fitness curve around the optimum.

---

### 1.3 `SequenceSimilarityFitness`

**Biological idea**

- Models selection on **sequence similarity** between two haplotypes (e.g., homologous chromosomes in a diploid individual).
- In a centromere context, large sequence divergence between homologs can impair **pairing**, **kinetochore assembly**, or chromosome segregation.
- Captures generic **compatibility**: more similar homologs typically yield higher fitness.

**Implementation sketch**

- Consider two sequences $A$ and $B$.
- Compute a distance based on **Hamming distance** plus **length difference**:
  - Let $d_{\text{ham}}$ be the Hamming distance over the aligned prefix.
  - Let $d_{\text{len}} = |\text{len}(A) - \text{len}(B)|$.
  - Define total distance $d = d_{\text{ham}} + d_{\text{len}}$.
- Convert to similarity $s \in [0, 1]$, e.g.:

  - $s = 1 - d / d_{\max}$ (clamped), where $d_{\max}$ is the maximum possible distance for normalization.

- Apply a shape parameter $\gamma$ to adjust selection strength:

  - $f = s^{\gamma}$

- Complexity: **linear** in shorter sequence length.
- Implemented in `selection.rs` as an `IndividualFitness`:
  - Iterates over both sequences in lockstep.
  - Counts mismatches and handles length differences.
  - Raises normalized similarity to a power.

---

### 1.4 `LengthSimilarityFitness`

**Biological idea**

- Models selection for **matching lengths** between homologous centromeres.
- Length mismatches between homologs can lead to **dosage imbalance**, pairing issues, or segregation defects.
- This is a coarse but very cheap proxy for structural compatibility.

**Implementation sketch**

- Let $L_1$ and $L_2$ be the lengths (or RU counts) of two homologous sequences.
- WLOG assume $L_1 \ge L_2$.
- Define a normalized length ratio $r = (L_1 - L_2) / L_1 \in [0,1]$.
- Fitness penalizes increasing $r$ using a shape parameter $\eta$:

  - $f = (1 - r)^{\eta}$

- Complexity: constant time given the two lengths.
- Implemented in `selection.rs` as an `IndividualFitness`:
  - Computes the normalized length difference.
  - Applies a power-law penalty around perfect equality.

---

## 2. New, HOR-Only and Ultra-Fast Fitness Ideas

All of the following assume that we only model the **active HOR array** and want extremely fast evaluation **under continuous mutation** of RUs:

- Operations are at the level of **RUs** (repeat units) or **HOR copies**.
- Prefer **one-pass** algorithms over the RU sequence.
- Avoid full motif scanning / large sliding windows.
- Avoid any per-generation, per-RU comparison to a large set of canonical sequences.

To keep the simulation practical when RUs can change every generation, we assume a two-level strategy:

1. **Local RU features / types** are derived using only *local, cheap rules* on each RU (e.g., a few positions, a short motif, simple GC or length bins).
2. These features are **updated only when an RU actually mutates**, at the end of the mutation routine that already touches its bases.

All new fitness functions below operate **only on these cached RU features** (small integers/booleans, plus RU counts/lengths), not on raw bases. This keeps fitness evaluation strictly linear-time with a very low constant factor.

### 2.1 RU-Type Composition Fitness

**Biological idea**

- Many centromeric HOR arrays consist of a limited set of **RU types** (e.g., monomers A/B/C, with or without specific binding motifs).
- Active centromeres often show a characteristic **mix** of those RU types (e.g., a particular fraction of CENP-B–box-containing RUs).
- This fitness function favors arrays whose **RU-type composition** matches a desired profile.

**Implementation sketch (mutation-aware, cached features)**

- Represent each RU by a small integer or enum `ru_type` that encodes **coarse, local features** rather than full canonical matching. For example, `ru_type` may combine:
  - `is_binding` (presence/absence of a short binding motif, detected by checking a few positions or a short pattern),
  - `gc_bin` (GC fraction binned into 2–4 levels, computed while applying mutations),
  - `length_bin` (RU length binned into short/normal/long).
- When an RU mutates (substitution or indel), recompute these features **in the mutation routine** that already scans its bases, and update its cached `ru_type`. RUs that do not mutate in a generation keep their existing `ru_type`.
- For a haplotype with $N$ RUs, fitness evaluation then does **one pass over cached `ru_type` values**:
  - Maintain a small histogram `counts[i]` of occurrences of each `ru_type`.
- Convert counts to proportions $p_i = \text{counts}_i / N$.
- Define a target composition $p_i^*$ and a strength parameter $\alpha$.
- Compute a simple deviation score and convert to fitness, for example:

  - Quadratic deviation:
    - $D = \sum_i (p_i - p_i^*)^2$
    - $f = \exp(-\alpha D)$

- Complexity: $O(N_\text{RUs})$ with tiny constants (integer increments, then a few float ops); there are **no base-level comparisons** during fitness evaluation.
- Integration:
  - New `HaplotypeFitness` struct (e.g., `RUTypeCompositionFitness`).
  - Stores target proportions and strength; on evaluation, walks RUs once, reads cached `ru_type`, and computes $f$.

---

### 2.2 HOR Pattern Regularity Fitness

**Biological idea**

- Centromeric HORs are often **highly regular**: a canonical pattern of RUs repeats many times.
- Disruptions to this pattern (insertions, deletions, rearranged RU types) reflect structural instability or degeneration.
- This fitness favors arrays that closely follow a fixed **canonical HOR pattern** or periodicity.

**Implementation sketch (type-based, no per-generation canonical matching)**

- Work entirely on cached **RU types** rather than raw sequences.
- Assume we know a **canonical pattern of RU types** as a short list:
  - `canonical_types = [t_0, t_1, ..., t_{k-1}]` with period $k`.
  - These types can be defined from an initial, well-characterized HOR, then treated as fixed labels.
- For a haplotype with cached `ru_type[i]` for each RU, in **one pass**:
  - For each position $i`, compare `ru_type[i]` to `canonical_types[i % k]`.
  - Count mismatches $m`.
- Define a mismatch fraction $u = m / N` and a strength parameter $\beta$.
- Convert to fitness, e.g.:

  - $f = \exp(-\beta u)$, or
  - $f = (1 - u)^{\beta}$ (clamped to \([0,1]\)).

- Complexity: $O(N_\text{RUs})`, just integer comparisons and a final exponential/power; **no base sequences are inspected during fitness evaluation**.
- Integration:
  - New `HaplotypeFitness` struct (e.g., `HORRegularityFitness`).
  - Stores `canonical_types` and strength; compares cached `ru_type` to the expected pattern on evaluation.

**Variant without explicit pattern**

- If only the **period length** $k$ is known, but not the exact pattern of types:
  - Compare cached `ru_type[i]` to `ru_type[i + k]` wherever defined; count matches and mismatches.
  - Let $M$ be the number of such comparisons, and $c$ the number of matches.
  - Fitness based on match fraction $c/M$ (higher = more regular repetition).

---

### 2.3 Local Homogeneity / RU-Clustering Fitness (Adjacency-Based)

**Biological idea**

- HOR arrays are often composed of **locally homogeneous blocks** of similar or identical RUs (e.g., homogenized HOR copies).
- Excessive fragmentation (alternating many different RU types) may indicate **structural instability**.
- Rather than full sliding-window statistics, we can use a cheap proxy based on **adjacent RUs**.

**Implementation sketch (feature-based adjacency)**

- For a haplotype with cached `ru_type[i]` for each RU:
  - For each $i$ in $0..N-2$, check whether `ru_type[i] == ru_type[i+1]`.
  - Let $c$ be the count of adjacent matches and $M = N - 1$ total adjacencies.
- Define homogeneity fraction $p = c / M$ and strength parameter $\lambda$.
- Fitness rewards larger $p$, for example:

  - $f = p^{\lambda}$ (with $p$ clamped away from 0), or
  - $f = \exp(\lambda p)$ with appropriate scaling.

- Complexity: $O(N_\text{RUs})`, one comparison per adjacency on small integers.
- Integration:
  - New `HaplotypeFitness` struct (e.g., `RUAdjacencyHomogeneityFitness`).
  - Evaluates with a single pass over the RU array, reading cached `ru_type` values only.

---

### 2.4 Binding-Competent RU Fraction Fitness

**Biological idea**

- Some RU types are **binding-competent** (e.g., they contain a motif recognized by centromeric proteins), others are not.
- There may be an optimal **fraction** of binding-competent RUs for robust kinetochore function.
- This fitness favors arrays with a binding-competent RU fraction near the optimum.

**Implementation sketch (locally detected binding, cached)**

- Assume each RU has a cached boolean property `is_binding` (e.g., presence of a CENP-B–like motif), determined by **local inspection** of that RU’s sequence:
  - When an RU mutates, re-check only that RU’s bases for the short binding motif and update `is_binding`.
  - Unmutated RUs retain their previous `is_binding` value.
- For a sequence of $N$ RUs at fitness evaluation time:
  - Count `b` = number of RUs with `is_binding == true`.
  - Compute fraction $q = b / N$.
- Define desired fraction $q^*$ and strength parameter $\alpha$.
- Use a simple bell-shaped curve, e.g.:

  - $f = \exp(-\alpha (q - q^*)^2)$

- Complexity: $O(N_\text{RUs})$ with a single boolean check per RU; no motif scanning occurs during fitness.
- Integration:
  - New `HaplotypeFitness` struct (e.g., `BindingRUFractionFitness`).
  - Stores $q^*$ and $\alpha$, reads cached `is_binding` flags, and evaluates $f$.

---

### 2.5 RU-Copy Number / Length Fitness (RU-Count Interpretation)

**Biological idea**

- Centromeres are often under selection for a specific **copy number** of HORs.
- Extremely short arrays may fail to form a functional kinetochore; extremely long arrays may be detrimental.
- This is essentially a re-interpretation of `LengthFitness` in terms of **RU count** rather than exact base pairs.

**Implementation sketch**

- Let $N$ be the number of RUs in the array.
- Choose an optimal copy number $N^*$ and log-space spread $\sigma$.
- As in `LengthFitness`, define:

  - $x = \ln N$, $x^* = \ln N^*$
  - $f = \exp\left( -\frac{(x - x^*)^2}{2\sigma^2} \right)$

- Complexity: constant time given $N$.
- Integration:
  - Either reuse `LengthFitness` by treating the length as `N` instead of base pairs.
  - Or add a thin wrapper struct that explicitly uses RU count if that’s more natural in your API.

---

## 3. Combining Fitness Components

For HOR-only, performance-sensitive simulations, a realistic yet cheap total fitness could be:

- `F_total = F_length × F_HOR_regular × F_RU_composition × F_binding_balance`

Where each factor is one of the functions above:

- `F_length`:
  - `LengthFitness` interpreted on RU count or base length.
- `F_HOR_regular`:
  - `HORRegularityFitness` based on canonical HOR pattern mismatches or period-based matches.
- `F_RU_composition`:
  - `RUTypeCompositionFitness` (histogram of RU types).
- `F_binding_balance`:
  - `BindingRUFractionFitness` (fraction of binding-competent RUs).

All of these:

- Run in **$O(N_\text{RUs})$ or better**.
- Use only **simple integer counting** and a small number of floats / exponentials.
- Are compatible with the existing `HaplotypeFitness` / `IndividualFitness` design in `selection.rs`.

You can implement them either as **separate fitness structs** combined multiplicatively, or as a **single composite** fitness that internally computes several of these signals in one pass over the RU sequence for maximum speed.
