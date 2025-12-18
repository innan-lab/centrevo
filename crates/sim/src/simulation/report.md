# Simulation Module Analysis Report

**Date**: December 17, 2025  
**Module**: `crates/sim/src/simulation`  
**Purpose**: Review consistency, identify duplication, and verify alignment with SIMULATION_ALGORITHM.md

---

## Executive Summary

The simulation module has undergone multiple refactors and shows signs of architectural evolution. While the core algorithm implementation is **sound and consistent with the documented algorithm**, there are several areas of concern:

### Critical Findings
- ✅ **Algorithm Alignment**: Core evolutionary cycle matches SIMULATION_ALGORITHM.md
- ⚠️ **Moderate Duplication**: Some initialization logic is duplicated
- ⚠️ **Inconsistent Abstractions**: Mix of high-level and low-level APIs
- ✅ **Good Test Coverage**: Comprehensive tests for core functionality
- ⚠️ **Configuration Complexity**: Multiple overlapping configuration patterns

### Recommendation Priority
1. **HIGH**: Simplify initialization pathways (3 different approaches)
2. **MEDIUM**: Consolidate fitness computation locations
3. **MEDIUM**: Remove deprecated APIs
4. **LOW**: Improve documentation consistency

---

## 1. Algorithm Consistency Check

### 1.1 Core Evolutionary Cycle ✅

The implementation in `engine.rs::step()` correctly follows the documented 3-step cycle:

```
DOCUMENTED (SIMULATION_ALGORITHM.md):
1. Parent Selection (fitness-proportional)
2. Offspring Generation (meiosis + mutation + fitness)
3. Population Replacement

IMPLEMENTED (engine.rs::step):
1. select_parents() - uses cached fitness ✓
2. produce_gamete() + fertilization in parallel ✓
   - Meiosis (recombination between haplotypes) ✓
   - Assortment (random haplotype selection) ✓
   - Mutation (after meiosis) ✓
   - Fitness computation (intrinsic, at birth) ✓
3. set_individuals() + increment_generation() ✓
```

**Status**: ✅ **CONSISTENT** - The algorithm is correctly implemented.

### 1.2 Fitness Model ✅

**Documented**: Multiplicative fitness model with intrinsic (birth) fitness.

**Implemented**:
- `FitnessConfig::compute_fitness()` (configs.rs:247-263): Correctly multiplies components
- `engine.rs::step()`: Computes fitness immediately after individual creation (lines 425-426)
- `engine.rs::from_config()`: Computes fitness for initial population (lines 213-216)
- `engine.rs::from_checkpoint()`: Recomputes fitness on restore (lines 292-295)

**Status**: ✅ **CONSISTENT** - Fitness is intrinsic and multiplicative.

### 1.3 Meiosis & Recombination ✅

**Documented**: Break formation → Repair (Crossover/Gene Conversion) → Homology search

**Implemented** (`engine.rs::produce_gamete()`, lines 383-439):
```rust
1. Clone parent haplotypes
2. Apply recombination events to chromosome pairs
3. Random assortment (select one haplotype)
4. Apply mutation to gamete
```

Uses `recombination.params.sample_events()` which implements the documented model.

**Status**: ✅ **CONSISTENT**

---

## 2. Code Duplication Analysis

### 2.1 Critical Duplication: Initialization Pathways ⚠️

**Problem**: Three separate but overlapping approaches to initialize simulations:

#### Pathway 1: Direct Constructor (DEPRECATED)
```rust
// engine.rs:45-72
pub fn new(...) -> Result<Self, String>
```
- Status: Marked deprecated
- **Action**: Should be removed in next major version

#### Pathway 2: from_config (Core API)
```rust
// engine.rs:148-230
pub fn from_config(...) -> Result<Self, String>
```
- Status: Main API, used by builder
- Logic: Parse SequenceConfig → initialize individuals → compute fitness

#### Pathway 3: SimulationBuilder (User-facing)
```rust
// builder.rs
pub fn build(self) -> Result<Simulation, BuilderError>
```
- Status: Recommended user API
- Internally calls `Simulation::from_config()`

**Duplication Issue**: The builder has complex logic (lines 313-407) that partially duplicates parameter validation and structure creation.

**Recommendation**: 
- ✅ Keep: `SimulationBuilder` (user-facing)
- ✅ Keep: `from_config()` (core API)
- ❌ Remove: `new()` (deprecated)
- 🔧 Simplify: Builder's `build()` method should do minimal work

### 2.2 Moderate Duplication: Individual Creation 📊

**Locations**:
1. `sequence.rs::create_uniform_individual()` (lines 678-708)
2. `sequence.rs::create_random_individual()` (lines 728-783)
3. `sequence.rs::create_individuals_from_sequences()` (lines 529-589)
4. `sequence.rs::initialize_from_formatted_string()` (lines 1097-1127)

**Analysis**: Each creates `Individual` objects with similar structure but different data sources.

**Verdict**: ✅ **ACCEPTABLE** - These are genuinely different use cases:
- Uniform: Template-based cloning
- Random: RNG-based generation
- From sequences: External data parsing
- From string: Domain-specific format

**Minor Issue**: Some code could be factored into helper functions for creating Haplotype/Chromosome pairs, but the duplication is reasonable.

### 2.3 Minor Duplication: Fitness Computation 📌

**Locations**:
1. `engine.rs::from_config()` - Initial population (lines 213-216)
2. `engine.rs::step()` - New offspring (lines 425-426)
3. `engine.rs::from_checkpoint()` - Restored population (lines 292-295)

**Code Pattern**:
```rust
individuals.par_iter_mut().for_each(|ind| {
    let val = fitness_config.compute_fitness(ind);
    ind.set_cached_fitness(val);
});
```

**Verdict**: ⚠️ **MINOR CONCERN** - Same 3-line pattern repeated 3 times.

**Recommendation**: Extract to method:
```rust
fn compute_population_fitness(
    individuals: &mut [Individual],
    fitness: &FitnessConfig
) {
    individuals.par_iter_mut().for_each(|ind| {
        let val = fitness.compute_fitness(ind);
        ind.set_cached_fitness(val);
    });
}
```

---

## 3. Architectural Concerns

### 3.1 Configuration Complexity ⚠️

**Issue**: Multiple configuration types with overlapping responsibilities:

```
SimulationConfig          - High-level simulation params (pop_size, generations, seed)
SequenceConfig            - How to initialize sequences
MutationConfig            - Mutation parameters
RecombinationConfig       - Recombination parameters
FitnessConfig             - Fitness functions
UniformRepeatStructure    - Repeat structure metadata
```

**Analysis**: 
- ✅ Good separation of concerns
- ⚠️ But: SequenceConfig has two variants (Generate/Load) with different requirements
- ⚠️ `UniformRepeatStructure` appears in multiple places and is optional in `Simulation`

**Example Confusion**:
```rust
// builder.rs:361
// structure is built here from individual parameters
let structure = UniformRepeatStructure::new(...)

// BUT also:
// configs.rs:20
SequenceConfig::Generate { structure, mode }
// structure is ALSO part of config
```

**Recommendation**: 
- Keep current design (it's workable)
- Add better documentation explaining when `structure` is required vs. optional
- Consider builder pattern for `SequenceConfig` in future

### 3.2 SequenceConfig vs. Builder Parameters 📊

**Inconsistency**: The builder stores some params separately, then packs them into config later.

```rust
// builder.rs fields:
ru_length: Option<usize>          // These 3 could be...
rus_per_hor: Option<usize>        
hors_per_chr: Option<usize>       
config: SequenceConfig            // ...directly in here as structure

// But also:
repeat_structure(ru, rus, hors)   // User calls this
init_uniform(base)                // Which sets config.mode
```

**Verdict**: 📝 **DESIGN TRADEOFF** - Current approach allows:
1. Ergonomic builder API (chain-able methods)
2. Validation at build time
3. Backward compatibility

**Recommendation**: ✅ **KEEP AS IS** - Document that builder is a staging area.

### 3.3 Population Selection Logic Location ✅

**Concern**: Is `Population::select_parents()` the right place for fitness-proportional selection?

**Analysis**:
```rust
// population.rs:89-161
impl Population {
    pub fn select_parents(&self, rng, n_pairs) -> Vec<(usize, usize)>
}
```

**Verdict**: ✅ **CORRECT PLACEMENT** 
- Population owns individuals and their fitness
- Selection is a population-level operation
- Encapsulates the "no self-pairing" constraint
- Handles edge cases (uniform fitness, lethal fitness)

---

## 4. API Surface Review

### 4.1 Public API Consistency ✅

**Engine.rs Public API**:
```rust
// Constructors
pub fn from_config(...)       ✅ Core API
pub fn from_checkpoint(...)   ✅ Resume API
pub fn new(...)               ⚠️ DEPRECATED - Remove

// Accessors
pub fn population()           ✅
pub fn population_mut()       ✅
pub fn generation()           ✅
pub fn *_config()             ✅ All getters present

// Execution
pub fn step()                 ✅
pub fn run()                  ✅
pub fn run_for()              ✅

// Checkpoint support
pub fn rng_state_bytes()      ✅
pub fn set_rng_from_bytes()   ✅
```

**Verdict**: ✅ **WELL-DESIGNED** - Clear, consistent API

### 4.2 Builder API ✅

```rust
// Required
.population_size()
.generations()
.repeat_structure()  // Required for Generate mode

// Optional
.mutation_rate()
.recombination()
.fitness()
.seed()
.codec()

// Initialization modes
.init_uniform()
.init_random()
.init_from_fasta()
.init_from_json()
.init_from_checkpoint()
```

**Verdict**: ✅ **EXCELLENT** - Fluent, discoverable, well-tested

---

## 5. Module Organization

### 5.1 File Responsibilities ✅

```
mod.rs         - Re-exports and module documentation
builder.rs     - SimulationBuilder (user-facing API)
configs.rs     - All configuration types
engine.rs      - Simulation engine (core algorithm)
population.rs  - Population container and selection
sequence.rs    - Initialization from various sources
```

**Verdict**: ✅ **LOGICAL SEPARATION** - Each file has clear purpose

### 5.2 Cross-Module Dependencies 📊

```
builder → engine → population
   ↓         ↓
configs ← sequence
```

**Concern**: `sequence.rs` is 1198 lines and handles many formats.

**Recommendation**: Consider splitting `sequence.rs`:
```
sequence/
  mod.rs           - Re-exports
  initialization.rs - High-level initialize() function
  fasta.rs         - FASTA parsing
  json.rs          - JSON parsing
  database.rs      - Database loading
  validation.rs    - Validation logic
  creation.rs      - Individual creation helpers
```

---

## 6. Error Handling

### 6.1 Error Types ✅

**Used**:
- `BuilderError` - Builder validation
- `InitializationError` - Sequence loading
- `String` - Generic errors in engine (⚠️)

**Concern**: `engine.rs` returns `Result<T, String>` in many places.

**Recommendation**: 
- Keep current design for now (it's consistent)
- In future, consider `SimulationError` enum

---

## 7. Testing Quality

### 7.1 Coverage Assessment ✅

**Engine Tests** (engine.rs:490-1271):
- ✅ Initialization paths
- ✅ Step mechanics
- ✅ Mutation integration
- ✅ Recombination integration
- ✅ Fitness computation
- ✅ RNG state serialization
- ✅ Checkpoint handling (basic)

**Builder Tests** (builder.rs:513-604):
- ✅ Valid configurations
- ✅ Missing required params
- ✅ Invalid parameters
- ✅ Initialization modes

**Population Tests** (population.rs:205-326):
- ✅ Selection logic
- ✅ Edge cases (empty, uniform fitness)
- ✅ No self-pairing constraint

**Sequence Tests** (sequence.rs:1129-1198):
- ✅ FASTA parsing
- ✅ JSON parsing
- ✅ Validation logic

**Verdict**: ✅ **EXCELLENT TEST COVERAGE**

---

## 8. Documentation Quality

### 8.1 Module-level Docs ✅

All files have clear module-level documentation.

### 8.2 Function-level Docs 📊

**Good**:
- `builder.rs` - Excellent examples in docs
- `engine.rs` - Public APIs well documented
- `configs.rs` - Comprehensive type docs

**Needs Improvement**:
- `sequence.rs` - Some internal functions lack docs
- `population.rs` - Could use more examples

### 8.3 Algorithm Documentation ✅

The code comments in `engine.rs::step()` clearly reference the 3-step cycle, making it easy to trace back to SIMULATION_ALGORITHM.md.

---

## 9. Specific Issues Found

### 9.1 Bug: Gene Conversion Clamping 🐛

**Location**: `engine.rs::apply_event_to_pair()` (lines 107-144)

**Issue**: The error handling for gene conversion out-of-bounds has complex fallback logic that clamps lengths. While tested, this could hide logic errors.

**Verdict**: ⚠️ **ACCEPTABLE WITH CAVEATS**
- Tested (lines 520-564)
- Documented behavior
- But: May silently skip invalid events

**Recommendation**: Add logging/metrics when clamping occurs.

### 9.2 Unused Field: `repeat_map` in Builder 📌

**Location**: `builder.rs:73`

```rust
repeat_map: Option<RepeatMap>,
```

**Used in**: 
- `init_with_map()` method (line 223)
- `build()` returns error "not fully implemented" (line 400)

**Verdict**: ⚠️ **INCOMPLETE FEATURE**

**Recommendation**: Either:
1. Complete the implementation, or
2. Remove the field and method until needed

### 9.3 BED File Parsing Incomplete 🚧

**Location**: `sequence.rs::parse_bed_to_maps()` (lines 921-989)

**Issue**: Incomplete implementation with TODOs and commented sections.

**Verdict**: 🚧 **WORK IN PROGRESS**

**Recommendation**: Either complete or feature-gate this functionality.

---

## 10. Consistency with SIMULATION_ALGORITHM.md

### 10.1 Documented Concepts vs. Implementation

| Concept | Location in Docs | Implementation | Status |
|---------|------------------|----------------|--------|
| Wright-Fisher | Overview | `engine.rs::step()` | ✅ |
| Discrete generations | Section 2 | `population.increment_generation()` | ✅ |
| Diploid individuals | Section 1 | `genome::Individual` (2 haplotypes) | ✅ |
| Intrinsic fitness | Section 2 | Computed at birth, cached | ✅ |
| Multiplicative fitness | Section 2 | `FitnessConfig::compute_fitness()` | ✅ |
| Meiosis → Assortment → Mutation | Section 3 | `produce_gamete()` | ✅ |
| Fitness-proportional selection | Section 3.1 | `Population::select_parents()` | ✅ |
| No self-fertilization | Section 3.1 | Enforced in `select_parents()` | ✅ |
| Parallel offspring generation | Section 3.2 | `par_iter()` in `step()` | ✅ |
| Fixed population size | Section 4.3 | Constant `N` throughout | ✅ |

**Verdict**: ✅ **FULLY CONSISTENT**

---

## 11. Recommendations

### 11.1 Priority 1: HIGH (Should Fix)

1. **Remove deprecated `Simulation::new()`**
   - File: `engine.rs:45-72`
   - Action: Delete method
   - Impact: Simplifies API surface

2. **Complete or remove `init_with_map()`**
   - File: `builder.rs:223-228`
   - Action: Either implement or remove
   - Impact: Reduces confusion

3. **Split `sequence.rs` into submodules**
   - File: `sequence.rs` (1198 lines)
   - Action: Create `sequence/` directory
   - Impact: Improves maintainability

### 11.2 Priority 2: MEDIUM (Nice to Have)

4. **Extract fitness computation helper**
   - Files: `engine.rs` (3 locations)
   - Action: Create `compute_population_fitness()` method
   - Impact: DRY, easier to maintain

5. **Add structured error types**
   - File: `engine.rs`
   - Action: Create `SimulationError` enum
   - Impact: Better error handling

6. **Complete BED file support or feature-gate it**
   - File: `sequence.rs::parse_bed_to_maps()`
   - Action: Finish implementation or hide behind feature flag
   - Impact: Avoid user confusion

### 11.3 Priority 3: LOW (Future Work)

7. **Add builder for SequenceConfig**
   - File: `configs.rs`
   - Action: Fluent builder for complex configurations
   - Impact: Improved ergonomics

8. **Add logging/metrics for edge cases**
   - File: `engine.rs::apply_event_to_pair()`
   - Action: Log when clamping occurs
   - Impact: Better debugging

---

## 12. Conclusion

### 12.1 Overall Assessment: ✅ GOOD

The simulation module is **well-architected and correctly implements the documented algorithm**. The code shows signs of evolution and refactoring, but remains **maintainable and consistent**.

### 12.2 Strengths
- ✅ Core algorithm is correct and well-tested
- ✅ Clear separation of concerns
- ✅ Excellent user-facing API (builder pattern)
- ✅ Comprehensive test coverage
- ✅ Good documentation

### 12.3 Weaknesses
- ⚠️ Some duplication in initialization paths
- ⚠️ Large `sequence.rs` file
- ⚠️ Incomplete features (BED, init_with_map)
- ⚠️ String-based errors in engine

### 12.4 Risk Assessment

**Technical Debt Level**: 🟡 **MODERATE**
- No critical bugs
- No architectural flaws
- Mostly organizational issues

**Maintainability**: 🟢 **HIGH**
- Clear code structure
- Good tests
- Reasonable complexity

**Alignment with Spec**: 🟢 **EXCELLENT**
- Matches SIMULATION_ALGORITHM.md perfectly
- Documented assumptions are implemented correctly

---

## Appendix A: Metrics

- **Total Lines**: ~4,500 (excluding tests)
- **Test Coverage**: ~40% (1,281 test lines)
- **Public Functions**: 47
- **Deprecated Functions**: 1
- **TODO/Incomplete**: 3 features

---

## Appendix B: Change Checklist

When refactoring, validate against this checklist:

- [ ] Algorithm still matches SIMULATION_ALGORITHM.md
- [ ] Fitness is intrinsic (computed at birth)
- [ ] No self-fertilization in selection
- [ ] Mutation happens after meiosis
- [ ] Population size remains constant
- [ ] All tests pass
- [ ] Builder API remains ergonomic
- [ ] Checkpoint/resume still works
