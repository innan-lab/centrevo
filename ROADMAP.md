# Centrevo Development Roadmap

**Project Status:** v0.1.0 - Core functionality implemented, ready for enhancement and production hardening

**Last Updated:** November 10, 2025

---

## 📊 Current Project Status

### ✅ Completed (Core v0.1.0)

**Solid Foundation**
- ✅ Core data structures (Base, Genome, Evolution modules) - ~7,091 LOC
- ✅ Comprehensive test coverage (247 tests passing, 0 failures)
- ✅ Full simulation engine with mutation, recombination, and selection
- ✅ SQLite-based persistence with recording strategies
- ✅ CLI interface with 4 commands (init, list, info, generations)
- ✅ Python bindings (PyO3) for all core functionality
- ✅ Extensive benchmark suite (Criterion-based)
- ✅ Parallel computation support (Rayon)
- ✅ Comprehensive documentation (CLI.md, PYTHON.md, storage README)

**Technical Achievements**
- Memory-efficient shared sequences (Arc-based)
- Professional error handling throughout
- Serializable configurations (JSON)
- Well-documented API with examples
- Performance benchmarking infrastructure

### ⚠️ Known Issues (Non-Critical)

**Code Quality (Clippy Warnings)**
- 88 clippy warnings (mostly deprecation notices in benchmarks)
  - 70+ warnings: `black_box` deprecation in benchmarks (use `std::hint::black_box`)
  - 3 warnings: `to_string()` methods should implement `Display` trait
  - Minor: manual clamp, collapsible if, etc.

**Missing Components**
- No main README.md file
- No LICENSE or CONTRIBUTING.md
- Utils module present but not implemented
- Limited analysis/visualization tools

---

## 🚀 Development Phases

---

## Phase 1: Quick Wins & Polish (1-2 weeks)
**Goal:** Fix code quality issues and add essential documentation

### 1.1 Code Quality Fixes (Priority: HIGH, Effort: LOW)

**Fix Clippy Warnings**
- [ ] Replace `criterion::black_box` with `std::hint::black_box` in all benchmarks
  - Files: `benches/*.rs` (4 files)
  - Impact: Eliminates 70+ deprecation warnings
  - Effort: 30 minutes (find-replace operation)

- [ ] Implement `Display` trait for `Sequence`, `SharedSequence`, and `Chromosome`
  - Files: `src/base/sequence.rs`, `src/genome/chromosome.rs`
  - Impact: Better Rust idioms, removes 3 warnings
  - Effort: 1 hour
  - Example:
    ```rust
    impl Display for Sequence {
        fn fmt(&self, f: &mut Formatter) -> fmt::Result {
            write!(f, "{}", self.iter().map(|n| n.to_char()).collect::<String>())
        }
    }
    ```

- [ ] Fix remaining clippy warnings (clamp, collapsible_if, is_multiple_of)
  - Files: `src/evolution/selection.rs`, `src/storage/*.rs`
  - Impact: Cleaner code
  - Effort: 30 minutes

**Expected Outcome:** Clean `cargo clippy` run with 0 warnings

### 1.2 Essential Documentation (Priority: HIGH, Effort: LOW)

- [ ] Create comprehensive **README.md**
  - Project overview and key features
  - Quick start guide
  - Installation instructions
  - Links to detailed docs
  - Example usage
  - Contributing guidelines
  - Effort: 2-3 hours

- [ ] Add **LICENSE** file
  - Choose appropriate license (MIT, Apache-2.0, GPL-3.0?)
  - Add to Cargo.toml
  - Effort: 15 minutes

- [ ] Create **CONTRIBUTING.md**
  - Code style guidelines
  - Testing requirements
  - Pull request process
  - Development setup
  - Effort: 1 hour

- [ ] Add **CHANGELOG.md**
  - Document v0.1.0 features
  - Template for future releases
  - Effort: 30 minutes

**Expected Outcome:** Professional, GitHub-ready project structure

### 1.3 CLI Enhancements (Priority: MEDIUM, Effort: LOW)

- [ ] Add `centrevo run` command to actually execute simulations
  - Currently only initializes; doesn't run generations
  - Add progress bar (indicatif crate)
  - Support for resuming interrupted simulations
  - Effort: 4-6 hours

- [ ] Add `centrevo export` command
  - Export data to CSV/JSON formats
  - Generate summary statistics
  - Effort: 3-4 hours

- [ ] Add `centrevo validate` command
  - Check database integrity
  - Verify simulation completeness
  - Effort: 2 hours

**Expected Outcome:** Full-featured CLI for production use

---

## Phase 2: Essential Features (2-4 weeks)
**Goal:** Add critical missing functionality for research use

### 2.1 Analysis & Visualization (Priority: HIGH, Effort: MEDIUM)

- [ ] Implement **analysis module** in Rust
  - Diversity metrics (nucleotide diversity π, Tajima's D)
  - Linkage disequilibrium
  - Sequence alignment utilities
  - Haplotype networks
  - Files: Create `src/analysis/` module
  - Effort: 1-2 weeks

- [ ] Add **plotting capabilities** (Python side)
  - Fitness trajectories
  - Allele frequency dynamics
  - GC content distribution
  - Pairwise distance matrices
  - Integration with matplotlib/seaborn
  - Files: Create `python/centrevo/plotting.py`
  - Effort: 1 week

- [ ] Create **example notebooks**
  - Jupyter notebooks demonstrating workflows
  - Tutorial for common analyses
  - Files: Create `notebooks/` directory
  - Effort: 3-4 days

**Expected Outcome:** Complete analysis pipeline from simulation to publication-ready figures

### 2.2 Performance Optimization (Priority: MEDIUM, Effort: MEDIUM)

- [ ] Add **compression** for sequence storage
  - Implement gzip/zstd compression for database BLOBs
  - Reduce disk usage by 60-80%
  - Files: `src/storage/database.rs`, `src/storage/recorder.rs`
  - Effort: 3-4 days

- [ ] Optimize **hot paths** identified in benchmarks
  - Profile mutation operations
  - Optimize fitness calculations
  - Consider SIMD for sequence operations
  - Effort: 1 week

- [ ] Add **streaming** support for large datasets
  - Iterator-based population queries
  - Avoid loading entire generation into memory
  - Files: `src/storage/query.rs`
  - Effort: 3-4 days

**Expected Outcome:** 2-3x speedup for large simulations, 70% reduction in storage

### 2.3 Testing & Quality (Priority: HIGH, Effort: LOW-MEDIUM)

- [ ] Add **integration tests**
  - End-to-end simulation tests
  - CLI integration tests
  - Python binding tests
  - Files: Create `tests/` directory
  - Effort: 1 week

- [ ] Set up **CI/CD pipeline**
  - GitHub Actions workflow
  - Automated testing on push
  - Benchmark regression detection
  - Clippy and rustfmt checks
  - Files: `.github/workflows/ci.yml`
  - Effort: 1-2 days

- [ ] Add **code coverage** tracking
  - Integrate with codecov.io or coveralls
  - Target: >80% coverage
  - Effort: 1 day

**Expected Outcome:** Robust, automated quality assurance

---

## Phase 3: Advanced Features (1-2 months)
**Goal:** Add sophisticated research capabilities

### 3.1 Advanced Evolution Models (Priority: MEDIUM, Effort: HIGH)

- [ ] **Complex selection models**
  - Frequency-dependent selection
  - Epistatic interactions
  - Environmental fluctuations
  - Files: Extend `src/evolution/selection.rs`
  - Effort: 2-3 weeks

- [ ] **Variable mutation rates**
  - Position-specific mutation rates
  - Hypermutable regions
  - Context-dependent mutation
  - Files: Extend `src/evolution/mutation.rs`
  - Effort: 1-2 weeks

- [ ] **Gene conversion models**
  - Biased gene conversion
  - Non-allelic homologous recombination
  - Files: Extend `src/evolution/recombination.rs`
  - Effort: 1 week

**Expected Outcome:** Realistic models for complex evolutionary scenarios

### 3.2 Population Genetics Features (Priority: MEDIUM, Effort: MEDIUM)

- [ ] **Demographic models**
  - Population size changes
  - Bottlenecks and expansions
  - Migration between populations
  - Files: Create `src/demography/` module
  - Effort: 2-3 weeks

- [ ] **Coalescent simulations**
  - Backward-time simulation
  - Integration with forward simulation
  - Files: Create `src/coalescent/` module
  - Effort: 3-4 weeks

- [ ] **Pedigree tracking**
  - Track parent-offspring relationships
  - Lineage reconstruction
  - Files: Extend `src/genome/individual.rs`, `src/storage/`
  - Effort: 1-2 weeks

**Expected Outcome:** Complete population genetics simulation suite

### 3.3 Advanced Storage & Retrieval (Priority: MEDIUM, Effort: MEDIUM)

- [ ] **Checkpoint/resume functionality**
  - Save and restore simulation state
  - Resume interrupted simulations
  - Files: `src/storage/checkpoint.rs`
  - Effort: 1 week

- [ ] **Distributed storage support**
  - PostgreSQL backend option
  - Cloud storage integration (S3, GCS)
  - Files: Create `src/storage/backends/`
  - Effort: 2-3 weeks

- [ ] **Query optimization**
  - Custom indices for common queries
  - Materialized views for statistics
  - Query builder improvements
  - Files: `src/storage/query.rs`
  - Effort: 1 week

**Expected Outcome:** Scalable storage for large-scale research projects

---

## Phase 4: Ecosystem & Adoption (Ongoing)
**Goal:** Build community and expand use cases

### 4.1 Documentation & Education (Priority: HIGH, Effort: MEDIUM)

- [ ] **Tutorial series**
  - Written tutorials for common scenarios
  - Video walkthroughs
  - Interactive web demos
  - Effort: Ongoing

- [ ] **API documentation improvements**
  - More examples in docstrings
  - Theory background in module docs
  - Mathematical formulations
  - Files: Throughout `src/`
  - Effort: 2-3 weeks

- [ ] **Publication & citation**
  - Write methods paper
  - Submit to bioinformatics journal
  - Create Zenodo DOI
  - Effort: 2-3 months

**Expected Outcome:** Well-documented, citable research tool

### 4.2 Interoperability (Priority: MEDIUM, Effort: MEDIUM)

- [ ] **Standard format support**
  - VCF file import/export
  - FASTA/FASTQ support
  - GFF/GTF annotation files
  - Files: Create `src/io/` module
  - Effort: 2-3 weeks

- [ ] **Integration with existing tools**
  - PopGen tools (ms, SLiM, msprime)
  - Phylogenetic software (BEAST, RAxML)
  - Files: Create `src/converters/`
  - Effort: 3-4 weeks

- [ ] **R bindings**
  - extendr for R interface
  - R package for analysis
  - Files: Create `R/` directory
  - Effort: 2-3 weeks

**Expected Outcome:** Seamless integration with existing research workflows

### 4.3 User Interface (Priority: LOW, Effort: HIGH)

- [ ] **Web interface**
  - Browser-based simulation runner
  - Real-time visualization
  - Parameter exploration
  - Tech: WebAssembly + React/Vue
  - Effort: 1-2 months

- [ ] **GUI application**
  - Desktop app (Tauri or egui)
  - Drag-and-drop configuration
  - Interactive plots
  - Effort: 1-2 months

**Expected Outcome:** Accessible to researchers without programming experience

---

## Phase 5: Research Extensions (Future)
**Goal:** Support cutting-edge research questions

### 5.1 Machine Learning Integration

- [ ] Parameter inference using ML
- [ ] Surrogate models for expensive simulations
- [ ] Automated model selection
- Effort: 2-3 months

### 5.2 Specialized Models

- [ ] Centromere-specific evolution (drive, meiotic drive)
- [ ] Satellite DNA dynamics
- [ ] Chromosome segregation modeling
- Effort: 3-6 months

### 5.3 Scalability

- [ ] GPU acceleration (CUDA/HIP)
- [ ] Distributed computing (MPI)
- [ ] Cloud-native deployment
- Effort: 2-4 months

---

## 📈 Success Metrics

### Technical Metrics
- **Code Quality:** 0 clippy warnings, >80% test coverage
- **Performance:** <1 second per generation for 1000 individuals
- **Storage:** <100 MB per 1000 generations (with compression)
- **Reliability:** >99% test pass rate, 0 critical bugs

### Adoption Metrics
- **Users:** 10+ research groups using Centrevo
- **Citations:** Published methods paper with citations
- **Contributions:** 5+ external contributors
- **Stars:** 100+ GitHub stars

### Research Impact
- **Publications:** Results used in 5+ peer-reviewed papers
- **Features:** All major population genetics scenarios supported
- **Validation:** Results match theoretical predictions

---

## 🎯 Prioritized Next Steps (Immediate)

### Week 1-2: Code Quality Sprint
1. **Day 1-2:** Fix all clippy warnings (black_box, Display, etc.)
2. **Day 3-4:** Create comprehensive README.md
3. **Day 5-7:** Add LICENSE, CONTRIBUTING.md, CHANGELOG.md
4. **Day 8-10:** Implement `centrevo run` command with progress bar
5. **Day 11-14:** Add CI/CD pipeline (GitHub Actions)

### Week 3-4: Essential Features
1. **Day 15-21:** Implement compression for database storage
2. **Day 22-28:** Create basic analysis module (diversity metrics)

### Week 5-8: Analysis & Visualization
1. **Weeks 5-6:** Complete analysis module implementation
2. **Week 7:** Add Python plotting capabilities
3. **Week 8:** Create example Jupyter notebooks

---

## 💡 Key Recommendations

### Immediate Priority (Do First)
1. **Fix clippy warnings** - Sets professional tone, takes <2 hours
2. **Add README.md** - Critical for GitHub visibility
3. **Implement `centrevo run`** - Makes CLI actually useful
4. **Set up CI/CD** - Prevents quality regression

### High Impact, Medium Effort
1. **Add compression** - Huge storage savings
2. **Analysis module** - Core research functionality
3. **Integration tests** - Confidence for changes
4. **Python plotting** - Publication-ready outputs

### Long-term Strategic
1. **Write methods paper** - Academic credibility
2. **Web interface** - Broader adoption
3. **GPU acceleration** - Competitive advantage

---

## 🤝 Community Building

### Attract Contributors
- Clean, well-documented code
- Good first issues tagged
- Responsive maintainers
- Clear contribution guidelines

### Build User Base
- Tutorial content
- Example use cases
- Active social media presence
- Conference presentations

### Academic Recognition
- Methods paper
- Citation guidelines
- Zenodo DOI
- Reproducible examples

---

## 📚 Resources Needed

### Development
- CI/CD credits (GitHub Actions - free tier sufficient)
- Code coverage service (free for open source)
- Documentation hosting (docs.rs for Rust, ReadTheDocs for Python)

### Research
- Benchmark datasets
- Validation against theoretical models
- Computational resources for testing

### Outreach
- Conference travel funds
- Workshop materials
- Website hosting

---

## 🎓 Learning Opportunities

This roadmap presents opportunities to learn:
- Advanced Rust patterns (async, macros, unsafe)
- Performance optimization techniques
- Database optimization
- Web development (WASM)
- Machine learning integration
- Scientific publishing

---

## ✅ Success Criteria by Phase

**Phase 1 Complete When:**
- 0 clippy warnings
- README, LICENSE, CONTRIBUTING present
- CLI can run full simulations
- CI/CD operational

**Phase 2 Complete When:**
- Compression reduces storage by >50%
- Analysis module with 10+ metrics
- 3+ example notebooks
- Test coverage >80%

**Phase 3 Complete When:**
- 5+ complex selection models
- Demographic model support
- Checkpoint/resume working
- Integration tests comprehensive

**Phase 4 Complete When:**
- Methods paper submitted
- 50+ GitHub stars
- 3+ external contributors
- Standard format support

---

## 📞 Getting Help

**Technical Questions:**
- Rust Discord server
- Stack Overflow (tag: rust, bioinformatics)

**Scientific Questions:**
- Population genetics forums
- Bioinformatics communities

**Collaboration:**
- Open GitHub issues
- Email maintainers
- Conference networking

---

**Note:** This roadmap is a living document. Priorities may shift based on user feedback, research needs, and resource availability. Update quarterly based on progress and community input.

**Version:** 1.0  
**Authors:** Development Team  
**Contact:** [Add contact information]
