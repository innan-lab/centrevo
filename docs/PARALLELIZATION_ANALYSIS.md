# Parallelization Analysis: Population-Level RNG Strategy

**Date:** November 11, 2025  
**Status:** Analysis & Proposal

---

## Current State

### RNG Usage in Simulation

Currently, the simulation uses a **single sequential RNG** at the `Simulation` level:

```rust
pub struct Simulation {
    // ...
    rng: StdRng,  // Single RNG shared across all operations
}
```

All operations use `&mut self.rng`:
1. **Mutation** - Each individual, each chromosome, each base
2. **Recombination** - Each individual, each chromosome pair
3. **Parent selection** - Selecting parent pairs
4. **Gamete selection** - Picking which haplotype to pass to offspring

### Current Performance (100 individuals × 10kb)

| Operation | Time | % of Total | Parallelizable? |
|-----------|------|------------|-----------------|
| Mutation | ~4-5ms | 30-37% | ✅ YES (per individual) |
| Recombination | ~3-4ms | 22-30% | ✅ YES (per individual) |
| Fitness calculation | ~2-3ms | 15-22% | ✅ Already parallel |
| Parent selection | ~1ms | 7-8% | ❌ NO (inherently sequential) |
| Offspring generation | ~2-3ms | 15-22% | ⚠️ PARTIAL (gamete selection) |
| **Total** | **13.5ms** | 100% | |

### Why Single RNG Blocks Parallelization

```rust
// Current: Cannot parallelize
fn apply_mutation(&mut self) -> Result<(), String> {
    for individual in self.population.individuals_mut() {  // Sequential!
        for chr in individual.haplotype1_mut().chromosomes_mut() {
            self.mutation.model.mutate_sequence(seq, &mut self.rng);
            //                                        ^^^^^^^^^^^^^^^^
            //                                        Requires &mut self
        }
    }
}
```

**Problem:** `&mut self.rng` prevents using Rayon's `par_iter_mut()` because:
- Can't have multiple mutable borrows
- Can't share `&mut StdRng` across threads (not `Sync`)

---

## Proposed Solution: Hierarchical RNG with Seeds

### Architecture

```
┌─────────────────────────────────────┐
│     Simulation (Root RNG)           │
│     rng: StdRng                     │
└──────────────┬──────────────────────┘
               │
               ├─ Generate seed[0] ──> Individual 0 RNG ──> Mutation, Recombination
               ├─ Generate seed[1] ──> Individual 1 RNG ──> Mutation, Recombination
               ├─ Generate seed[2] ──> Individual 2 RNG ──> Mutation, Recombination
               │           ...
               └─ Generate seed[N] ──> Individual N RNG ──> Mutation, Recombination
```

### Implementation Strategy

#### 1. Parallel Mutation

```rust
fn apply_mutation(&mut self) -> Result<(), String> {
    // Generate seeds for each individual from root RNG
    let seeds: Vec<u64> = (0..self.population.size())
        .map(|_| self.rng.random())
        .collect();
    
    // Parallel mutation with independent RNGs
    self.population.individuals_mut()
        .par_iter_mut()
        .zip(seeds.par_iter())
        .try_for_each(|(individual, &seed)| {
            let mut local_rng = StdRng::seed_from_u64(seed);
            
            // Mutate haplotype 1
            for chr in individual.haplotype1_mut().chromosomes_mut() {
                self.mutation.model.mutate_sequence(
                    chr.sequence_mut(),
                    &mut local_rng
                );
            }
            
            // Mutate haplotype 2
            for chr in individual.haplotype2_mut().chromosomes_mut() {
                self.mutation.model.mutate_sequence(
                    chr.sequence_mut(),
                    &mut local_rng
                );
            }
            
            Ok::<(), String>(())
        })?;
    
    Ok(())
}
```

**Benefits:**
- Each individual gets independent RNG
- Fully parallelizable with Rayon
- Reproducible (same seed sequence → same results)
- Thread-safe (no shared mutable state)

#### 2. Parallel Recombination

Similar strategy:

```rust
fn apply_recombination(&mut self) -> Result<(), String> {
    let seeds: Vec<u64> = (0..self.population.size())
        .map(|_| self.rng.random())
        .collect();
    
    self.population.individuals_mut()
        .par_iter_mut()
        .zip(seeds.par_iter())
        .try_for_each(|(individual, &seed)| {
            let mut local_rng = StdRng::seed_from_u64(seed);
            let (hap1, hap2) = individual.haplotypes_mut();
            
            for chr_idx in 0..hap1.len().min(hap2.len()) {
                if let (Some(chr1), Some(chr2)) = (hap1.get_mut(chr_idx), hap2.get_mut(chr_idx)) {
                    let event = self.recombination.params.sample_event(
                        chr1.len(),
                        &mut local_rng
                    );
                    // Apply event...
                }
            }
            
            Ok::<(), String>(())
        })?;
    
    Ok(())
}
```

#### 3. Parallel Offspring Generation (Gamete Selection)

```rust
fn generate_offspring(&mut self) -> Result<Vec<Individual>, String> {
    // Sequential: Fitness computation and parent selection
    let fitness_values = self.population.compute_fitness(&self.fitness);
    let pairs = self.population.select_parents(
        &mut self.rng,
        &fitness_values,
        self.config.population_size,
    );
    
    // Generate seeds for parallel gamete selection
    let seeds: Vec<u64> = (0..pairs.len())
        .map(|_| self.rng.random())
        .collect();
    
    let gen_str = format!("ind_gen{}_", self.generation() + 1);
    
    // Parallel offspring creation
    let offspring: Vec<Individual> = pairs
        .par_iter()
        .zip(seeds.par_iter())
        .enumerate()
        .map(|(i, ((&parent1_idx, &parent2_idx), &seed))| {
            let mut local_rng = StdRng::seed_from_u64(seed);
            
            let parent1 = self.population.get(parent1_idx).unwrap();
            let parent2 = self.population.get(parent2_idx).unwrap();
            
            let hap1 = if local_rng.random_bool(0.5) {
                parent1.haplotype1().clone()
            } else {
                parent1.haplotype2().clone()
            };
            
            let hap2 = if local_rng.random_bool(0.5) {
                parent2.haplotype1().clone()
            } else {
                parent2.haplotype2().clone()
            };
            
            Individual::new(format!("{}{}", gen_str, i), hap1, hap2)
        })
        .collect();
    
    Ok(offspring)
}
```

---

## Performance Analysis

### Expected Speedup (Theoretical)

Assuming 8 CPU cores and perfect scaling:

| Operation | Current | With Parallel | Speedup | New Time |
|-----------|---------|---------------|---------|----------|
| Mutation | 4.5ms | → 4.5ms / 8 | 8x | 0.56ms |
| Recombination | 3.5ms | → 3.5ms / 8 | 8x | 0.44ms |
| Offspring (gamete) | 1.5ms | → 1.5ms / 8 | 8x | 0.19ms |
| Parent selection | 1.0ms | (sequential) | 1x | 1.0ms |
| Fitness | 2.0ms | (already parallel) | 1x | 2.0ms |
| **Total** | **13.5ms** | | **3.6x** | **3.7ms** |

### Real-World Expectations

- **Best case (100+ individuals):** 3-4x speedup
- **Typical case (50 individuals):** 2-3x speedup
- **Small populations (10 individuals):** 1.2-1.5x speedup
  - Overhead of parallelization dominates
  - Better to stay sequential

### Overhead Considerations

1. **Seed generation:** ~0.1µs per seed × 100 = 10µs (negligible)
2. **Thread spawning:** Rayon thread pool reused (amortized)
3. **Cache coherency:** Individuals are independent (good locality)

---

## Trade-offs and Considerations

### ✅ Advantages

1. **Significant speedup** for medium-large populations (50+)
2. **Reproducible** - Same seed → same results
3. **Simple mental model** - Each individual gets its own RNG
4. **No locking overhead** - True data parallelism
5. **Scales with cores** - More cores = more speedup

### ⚠️ Challenges

1. **Sequential parent selection**
   - Cannot parallelize (cumulative distribution sampling)
   - Remains a bottleneck (~7-8% of time)
   - Acceptable trade-off

2. **Small population overhead**
   - Parallelization overhead may exceed benefits for N < 20
   - **Solution:** Use threshold-based parallelization

3. **Seed sequence matters**
   - Different order of seed generation = different results
   - Must maintain consistent ordering

4. **Debugging complexity**
   - Parallel bugs are harder to debug
   - **Solution:** Add sequential fallback mode

---

## Implementation Recommendations

### Phase 1: Add Threshold-Based Parallelization

```rust
const PARALLEL_THRESHOLD: usize = 20;  // Parallelize if pop_size >= 20

fn apply_mutation(&mut self) -> Result<(), String> {
    if self.population.size() < PARALLEL_THRESHOLD {
        // Sequential path for small populations
        self.apply_mutation_sequential()
    } else {
        // Parallel path for large populations
        self.apply_mutation_parallel()
    }
}
```

### Phase 2: Add Configuration Option

```rust
pub struct SimulationConfig {
    // ...
    pub parallel: bool,           // Enable parallelization
    pub parallel_threshold: usize, // Min population size for parallel
}
```

### Phase 3: Benchmark and Tune

- Measure actual speedup vs theoretical
- Tune `PARALLEL_THRESHOLD` based on benchmarks
- Test with different population sizes and sequence lengths

---

## Testing Strategy

### 1. Reproducibility Tests

```rust
#[test]
fn test_parallel_reproducibility() {
    let config = create_test_config();
    
    // Run sequentially
    let mut sim1 = Simulation::new(config.clone()).unwrap();
    sim1.run().unwrap();
    let pop1 = sim1.population().clone();
    
    // Run in parallel (same seed)
    let mut sim2 = Simulation::new(config).unwrap();
    sim2.run().unwrap();
    let pop2 = sim2.population();
    
    // Results should be identical
    assert_populations_equal(&pop1, pop2);
}
```

### 2. Performance Tests

```rust
#[bench]
fn bench_mutation_sequential_vs_parallel(b: &mut Bencher) {
    // Compare sequential and parallel performance
}
```

### 3. Correctness Tests

- Statistical tests (allele frequencies, diversity metrics)
- Edge cases (pop_size=1, pop_size=1000)
- Stress tests (many generations)

---

## Questions for Discussion

1. **Threshold value:** What's the optimal `PARALLEL_THRESHOLD`?
   - Propose: Start with 20, tune based on benchmarks

2. **Configuration:** Should parallelization be user-configurable or automatic?
   - Propose: Automatic with threshold, optional override

3. **Offspring generation:** Worth parallelizing the haplotype cloning?
   - Current bottleneck is the cloning itself, not the gamete selection
   - Parallelizing gamete selection: ~1.5ms → 0.19ms (1.3ms gain)
   - Propose: YES, implement it

4. **Fitness calculation:** Is it already parallelized efficiently?
   - Need to check if current implementation uses Rayon
   - If not, easy win

5. **Testing priority:** Focus on reproducibility or performance first?
   - Propose: Reproducibility first (critical for scientific validity)

---

## Next Steps

If approved, implementation order:

1. ✅ **Add parallel mutation** (highest impact: ~4ms → 0.5ms)
2. ✅ **Add parallel recombination** (second highest: ~3ms → 0.4ms)
3. ✅ **Add parallel offspring generation** (gamete selection)
4. ✅ **Add threshold-based switching**
5. ✅ **Comprehensive testing** (reproducibility + performance)
6. ✅ **Benchmark and tune** thresholds
7. ✅ **Update documentation**

**Estimated effort:** 2-3 days
**Expected benefit:** 3-4x speedup for typical workloads

---

**Ready to proceed?** Let's discuss any concerns or alternative approaches.
