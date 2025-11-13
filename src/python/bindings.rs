//! Python bindings for Centrevo using PyO3.
//!
//! This module provides Python interfaces to all major Centrevo functionality.

use pyo3::prelude::*;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::types::{PyList, PyDict};
use std::path::PathBuf;

use crate::base::{Alphabet, Nucleotide};
use crate::genome::{Chromosome, Haplotype, Individual};
use crate::simulation::{
    Population, RepeatStructure, SimulationConfig, Simulation,
    MutationConfig, RecombinationConfig, FitnessConfig, SequenceInput,
};
use crate::storage::{QueryBuilder, Recorder, RecordingStrategy, SimulationSnapshot};

/// Python module for Centrevo.
#[pymodule]
fn centrevo(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Core classes
    m.add_class::<PyNucleotide>()?;
    m.add_class::<PyAlphabet>()?;
    m.add_class::<PyChromosome>()?;
    m.add_class::<PyHaplotype>()?;
    m.add_class::<PyIndividual>()?;
    m.add_class::<PyPopulation>()?;
    m.add_class::<PyRepeatStructure>()?;
    m.add_class::<PySimulationConfig>()?;
    m.add_class::<PyRecorder>()?;
    m.add_class::<PyQueryBuilder>()?;
    m.add_class::<PyRecordingStrategy>()?;
    m.add_class::<PySimulation>()?;
    m.add_class::<PyMutationConfig>()?;
    m.add_class::<PyRecombinationConfig>()?;
    m.add_class::<PyFitnessConfig>()?;
    
    // Core functions
    m.add_function(wrap_pyfunction!(create_initial_population, m)?)?;
    
    // Analysis functions
    super::analysis::register_analysis_functions(m)?;
    
    // Export functions
    super::export::register_export_functions(m)?;
    
    Ok(())
}

/// Nucleotide base (A, C, G, T).
#[pyclass(name = "Nucleotide")]
#[derive(Clone)]
struct PyNucleotide {
    inner: Nucleotide,
}

#[pymethods]
impl PyNucleotide {
    #[new]
    fn new(base: &str) -> PyResult<Self> {
        if base.len() != 1 {
            return Err(PyValueError::new_err("Base must be a single character"));
        }

        let byte = base.as_bytes()[0];
        let inner = Nucleotide::from_ascii(byte)
            .ok_or_else(|| PyValueError::new_err(format!("Invalid nucleotide: {}", base)))?;

        Ok(Self { inner })
    }

    /// Get string representation.
    fn __str__(&self) -> String {
        self.inner.to_char().to_string()
    }

    /// Get Python representation.
    fn __repr__(&self) -> String {
        format!("Nucleotide('{}')", self.inner.to_char())
    }

    #[staticmethod]
    #[allow(non_snake_case)]
    fn A() -> Self {
        Self { inner: Nucleotide::A }
    }

    #[staticmethod]
    #[allow(non_snake_case)]
    fn C() -> Self {
        Self { inner: Nucleotide::C }
    }

    #[staticmethod]
    #[allow(non_snake_case)]
    fn G() -> Self {
        Self { inner: Nucleotide::G }
    }

    #[staticmethod]
    #[allow(non_snake_case)]
    fn T() -> Self {
        Self { inner: Nucleotide::T }
    }
}

/// Alphabet of nucleotide bases.
#[pyclass(name = "Alphabet")]
#[derive(Clone)]
struct PyAlphabet {
    inner: Alphabet,
}

#[pymethods]
impl PyAlphabet {
    #[new]
    fn new(chars: Vec<String>) -> PyResult<Self> {
        let char_vec: Result<Vec<char>, _> = chars
            .iter()
            .map(|s| {
                if s.len() == 1 {
                    Ok(s.chars().next().unwrap())
                } else {
                    Err(PyValueError::new_err("Each element must be a single character"))
                }
            })
            .collect();

        Ok(Self {
            inner: Alphabet::new(char_vec?),
        })
    }

    #[staticmethod]
    fn dna() -> Self {
        Self {
            inner: Alphabet::dna(),
        }
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __repr__(&self) -> String {
        format!("Alphabet({:?})", self.inner.chars())
    }
}

/// Chromosome with repeat structure.
#[pyclass(name = "Chromosome")]
struct PyChromosome {
    inner: Chromosome,
}

#[pymethods]
impl PyChromosome {
    #[staticmethod]
    fn uniform(
        id: String,
        base: &PyNucleotide,
        length: usize,
        ru_length: usize,
        rus_per_hor: usize,
        alphabet: &PyAlphabet,
    ) -> Self {
        Self {
            inner: Chromosome::uniform(
                id,
                base.inner,
                length,
                ru_length,
                rus_per_hor,
                alphabet.inner.clone(),
            ),
        }
    }

    #[getter]
    fn id(&self) -> String {
        self.inner.id().to_string()
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn to_string(&self) -> String {
        self.inner.to_string()
    }

    fn gc_content(&self) -> f64 {
        self.inner.gc_content()
    }

    fn __repr__(&self) -> String {
        format!("Chromosome(id='{}', length={})", self.inner.id(), self.inner.len())
    }
}

/// Haplotype (collection of chromosomes).
#[pyclass(name = "Haplotype")]
struct PyHaplotype {
    inner: Haplotype,
}

#[pymethods]
impl PyHaplotype {
    #[new]
    fn new() -> Self {
        Self {
            inner: Haplotype::new(),
        }
    }

    #[staticmethod]
    fn from_chromosomes(chromosomes: Vec<PyRef<PyChromosome>>) -> Self {
        let chrs: Vec<Chromosome> = chromosomes
            .iter()
            .map(|c| c.inner.clone())
            .collect();

        Self {
            inner: Haplotype::from_chromosomes(chrs),
        }
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __repr__(&self) -> String {
        format!("Haplotype({} chromosomes)", self.inner.len())
    }
}

/// Individual organism with diploid genome.
#[pyclass(name = "Individual")]
struct PyIndividual {
    inner: Individual,
}

#[pymethods]
impl PyIndividual {
    #[new]
    fn new(id: String, haplotype1: &PyHaplotype, haplotype2: &PyHaplotype) -> Self {
        Self {
            inner: Individual::new(id, haplotype1.inner.clone(), haplotype2.inner.clone()),
        }
    }

    #[getter]
    fn id(&self) -> String {
        self.inner.id().to_string()
    }

    #[getter]
    fn fitness(&self) -> f64 {
        self.inner.fitness()
    }

    #[setter]
    fn set_fitness(&mut self, fitness: f64) {
        self.inner.set_fitness(fitness);
    }

    fn __repr__(&self) -> String {
        format!("Individual(id='{}', fitness={})", self.inner.id(), self.inner.fitness())
    }
}

/// Population of individuals.
#[pyclass(name = "Population")]
pub(super) struct PyPopulation {
    pub(super) inner: Population,
}

#[pymethods]
impl PyPopulation {
    #[new]
    fn new(id: String, individuals: Vec<PyRef<PyIndividual>>) -> Self {
        let inds: Vec<Individual> = individuals
            .iter()
            .map(|i| i.inner.clone())
            .collect();

        Self {
            inner: Population::new(id, inds),
        }
    }

    #[getter]
    fn id(&self) -> String {
        self.inner.id().to_string()
    }

    fn size(&self) -> usize {
        self.inner.size()
    }

    fn generation(&self) -> usize {
        self.inner.generation()
    }

    fn __len__(&self) -> usize {
        self.inner.size()
    }

    fn __repr__(&self) -> String {
        format!(
            "Population(id='{}', size={}, generation={})",
            self.inner.id(),
            self.inner.size(),
            self.inner.generation()
        )
    }
}

/// Repeat structure configuration.
#[pyclass(name = "RepeatStructure")]
#[derive(Clone)]
struct PyRepeatStructure {
    inner: RepeatStructure,
}

#[pymethods]
impl PyRepeatStructure {
    #[new]
    fn new(
        alphabet: &PyAlphabet,
        init_base: &PyNucleotide,
        ru_length: usize,
        rus_per_hor: usize,
        hors_per_chr: usize,
        chrs_per_hap: usize,
    ) -> Self {
        Self {
            inner: RepeatStructure::new(
                alphabet.inner.clone(),
                init_base.inner,
                ru_length,
                rus_per_hor,
                hors_per_chr,
                chrs_per_hap,
            ),
        }
    }

    #[getter]
    fn ru_length(&self) -> usize {
        self.inner.ru_length
    }

    #[getter]
    fn rus_per_hor(&self) -> usize {
        self.inner.rus_per_hor
    }

    #[getter]
    fn hors_per_chr(&self) -> usize {
        self.inner.hors_per_chr
    }

    fn chr_length(&self) -> usize {
        self.inner.chr_length()
    }

    fn __repr__(&self) -> String {
        format!(
            "RepeatStructure(ru_length={}, rus_per_hor={}, hors_per_chr={})",
            self.inner.ru_length, self.inner.rus_per_hor, self.inner.hors_per_chr
        )
    }
}

/// Simulation configuration.
#[pyclass(name = "SimulationConfig")]
#[derive(Clone)]
struct PySimulationConfig {
    inner: SimulationConfig,
}

#[pymethods]
impl PySimulationConfig {
    #[new]
    fn new(population_size: usize, total_generations: usize, seed: Option<u64>) -> Self {
        Self {
            inner: SimulationConfig::new(population_size, total_generations, seed),
        }
    }

    #[getter]
    fn population_size(&self) -> usize {
        self.inner.population_size
    }

    #[getter]
    fn total_generations(&self) -> usize {
        self.inner.total_generations
    }

    #[getter]
    fn seed(&self) -> Option<u64> {
        self.inner.seed
    }

    fn __repr__(&self) -> String {
        format!(
            "SimulationConfig(population_size={}, total_generations={}, seed={:?})",
            self.inner.population_size, self.inner.total_generations, self.inner.seed
        )
    }
}

/// Recording strategy for simulation snapshots.
#[pyclass(name = "RecordingStrategy")]
#[derive(Clone)]
struct PyRecordingStrategy {
    inner: RecordingStrategy,
}

#[pymethods]
impl PyRecordingStrategy {
    #[staticmethod]
    fn every_n(n: usize) -> Self {
        Self {
            inner: RecordingStrategy::EveryN(n),
        }
    }

    #[staticmethod]
    fn specific(generations: Vec<usize>) -> Self {
        Self {
            inner: RecordingStrategy::Specific(generations),
        }
    }

    #[staticmethod]
    fn all() -> Self {
        Self {
            inner: RecordingStrategy::All,
        }
    }

    #[staticmethod]
    fn none() -> Self {
        Self {
            inner: RecordingStrategy::None,
        }
    }

    fn __repr__(&self) -> String {
        match &self.inner {
            RecordingStrategy::EveryN(n) => format!("RecordingStrategy.every_n({})", n),
            RecordingStrategy::Specific(gens) => format!("RecordingStrategy.specific({:?})", gens),
            RecordingStrategy::All => "RecordingStrategy.all()".to_string(),
            RecordingStrategy::None => "RecordingStrategy.none()".to_string(),
        }
    }
}

/// Simulation recorder.
#[pyclass(name = "Recorder", unsendable)]
struct PyRecorder {
    inner: Option<Recorder>,
}

#[pymethods]
impl PyRecorder {
    #[new]
    fn new(db_path: String, sim_id: String, strategy: &PyRecordingStrategy) -> PyResult<Self> {
        let recorder = Recorder::new(
            PathBuf::from(db_path),
            sim_id,
            strategy.inner.clone(),
        )
        .map_err(|e| PyRuntimeError::new_err(format!("Failed to create recorder: {}", e)))?;

        Ok(Self {
            inner: Some(recorder),
        })
    }

    fn record_metadata(&mut self, config: &PySimulationConfig) -> PyResult<()> {
        self.inner
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("Recorder has been closed"))?
            .record_metadata(&config.inner)
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to record metadata: {}", e)))?;

        Ok(())
    }

    fn record_full_config(
        &mut self,
        structure: &PyRepeatStructure,
        mutation: &PyMutationConfig,
        recombination: &PyRecombinationConfig,
        fitness: &PyFitnessConfig,
        config: &PySimulationConfig,
    ) -> PyResult<()> {
        let snapshot = SimulationSnapshot {
            structure: structure.inner.clone(),
            mutation: mutation.inner.clone(),
            recombination: recombination.inner.clone(),
            fitness: fitness.inner.clone(),
            config: config.inner.clone(),
        };

        self.inner
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("Recorder has been closed"))?
            .record_full_config(&snapshot)
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to record full config: {}", e)))?;

        Ok(())
    }

    fn record_generation(&mut self, population: &PyPopulation, generation: usize) -> PyResult<()> {
        self.inner
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("Recorder has been closed"))?
            .record_generation(&population.inner, generation)
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to record generation: {}", e)))?;

        Ok(())
    }

    fn record_checkpoint(&mut self, simulation: &PySimulation, generation: usize) -> PyResult<()> {
        let sim = simulation
            .inner
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("Simulation has been consumed"))?;

        let rng_state = sim.rng_state_bytes();

        self.inner
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("Recorder has been closed"))?
            .record_checkpoint(generation, &rng_state)
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to record checkpoint: {}", e)))?;

        Ok(())
    }

    fn finalize_metadata(&mut self) -> PyResult<()> {
        self.inner
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("Recorder has been closed"))?
            .finalize_metadata()
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to finalize metadata: {}", e)))?;

        Ok(())
    }

    fn close(&mut self) -> PyResult<()> {
        if let Some(recorder) = self.inner.take() {
            recorder
                .close()
                .map_err(|e| PyRuntimeError::new_err(format!("Failed to close recorder: {}", e)))?;
        }
        Ok(())
    }

    fn __repr__(&self) -> String {
        if self.inner.is_some() {
            "Recorder(active)".to_string()
        } else {
            "Recorder(closed)".to_string()
        }
    }
}

/// Query builder for simulation results.
#[pyclass(name = "QueryBuilder", unsendable)]
struct PyQueryBuilder {
    inner: Option<QueryBuilder>,
}

#[pymethods]
impl PyQueryBuilder {
    #[new]
    fn new(db_path: String) -> PyResult<Self> {
        let query = QueryBuilder::new(PathBuf::from(db_path))
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to open database: {}", e)))?;

        Ok(Self { inner: Some(query) })
    }

    fn list_simulations(&self, py: Python) -> PyResult<Py<PyList>> {
        let simulations = self
            .inner
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("QueryBuilder has been closed"))?
            .list_simulations()
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to list simulations: {}", e)))?;

        let list = PyList::new(py, simulations)?;
        Ok(list.into())
    }

    fn get_recorded_generations(&self, sim_id: String, py: Python) -> PyResult<Py<PyList>> {
        let generations = self
            .inner
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("QueryBuilder has been closed"))?
            .get_recorded_generations(&sim_id)
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to get generations: {}", e)))?;

        let list = PyList::new(py, generations)?;
        Ok(list.into())
    }

    fn get_simulation_info(&self, sim_id: String, py: Python) -> PyResult<Py<PyDict>> {
        let info = self
            .inner
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("QueryBuilder has been closed"))?
            .get_simulation_info(&sim_id)
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to get simulation info: {}", e)))?;

        let dict = PyDict::new(py);
        dict.set_item("sim_id", info.sim_id)?;
        dict.set_item("start_time", info.start_time)?;
        dict.set_item("end_time", info.end_time)?;
        dict.set_item("pop_size", info.pop_size)?;
        dict.set_item("num_generations", info.num_generations)?;
        dict.set_item("mutation_rate", info.mutation_rate)?;
        dict.set_item("recombination_rate", info.recombination_rate)?;
        dict.set_item("parameters_json", info.parameters_json)?;

        Ok(dict.into())
    }

    fn get_generation(&self, sim_id: String, generation: usize) -> PyResult<PyPopulation> {
        let snapshots = self
            .inner
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("QueryBuilder has been closed"))?
            .get_generation(&sim_id, generation)
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to load generation: {}", e)))?;

        if snapshots.is_empty() {
            return Err(PyRuntimeError::new_err(format!(
                "No data found for generation {}",
                generation
            )));
        }

        // Convert snapshots to individuals
        let alphabet = Alphabet::dna();
        let mut individuals = Vec::new();

        for snap in snapshots {
            let seq1 = crate::base::Sequence::from_indices(snap.haplotype1_seq, alphabet.clone());
            let seq2 = crate::base::Sequence::from_indices(snap.haplotype2_seq, alphabet.clone());

            let chr1 = Chromosome::new(snap.haplotype1_chr_id, seq1, 171, 12);
            let chr2 = Chromosome::new(snap.haplotype2_chr_id, seq2, 171, 12);

            let h1 = Haplotype::from_chromosomes(vec![chr1]);
            let h2 = Haplotype::from_chromosomes(vec![chr2]);

            let mut ind = Individual::new(snap.individual_id, h1, h2);
            ind.set_fitness(snap.fitness);
            individuals.push(ind);
        }

        Ok(PyPopulation {
            inner: Population::new(format!("{}_gen{}", sim_id, generation), individuals),
        })
    }

    fn get_fitness_history(&self, sim_id: String, py: Python) -> PyResult<Py<PyList>> {
        let history = self
            .inner
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("QueryBuilder has been closed"))?
            .get_fitness_history(&sim_id)
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to get fitness history: {}", e)))?;

        let mut records: Vec<Py<PyDict>> = Vec::new();
        for (generation, stats) in history {
            let dict = PyDict::new(py);
            dict.set_item("generation", generation)?;
            dict.set_item("mean", stats.mean)?;
            dict.set_item("std", stats.std)?;
            dict.set_item("min", stats.min)?;
            dict.set_item("max", stats.max)?;
            records.push(dict.into());
        }

        Ok(PyList::new(py, records)?.into())
    }

    fn close(&mut self) -> PyResult<()> {
        if let Some(query) = self.inner.take() {
            query
                .close()
                .map_err(|e| PyRuntimeError::new_err(format!("Failed to close database: {}", e)))?;
        }
        Ok(())
    }

    fn __repr__(&self) -> String {
        if self.inner.is_some() {
            "QueryBuilder(active)".to_string()
        } else {
            "QueryBuilder(closed)".to_string()
        }
    }
}

/// Helper function to create initial population.
#[pyfunction]
fn create_initial_population(
    size: usize,
    structure: &PyRepeatStructure,
) -> PyResult<PyPopulation> {
    let mut individuals = Vec::with_capacity(size);

    for i in 0..size {
        let chr1 = Chromosome::uniform(
            format!("ind{}_h1_chr1", i),
            structure.inner.init_base,
            structure.inner.chr_length(),
            structure.inner.ru_length,
            structure.inner.rus_per_hor,
            structure.inner.alphabet.clone(),
        );

        let chr2 = Chromosome::uniform(
            format!("ind{}_h2_chr1", i),
            structure.inner.init_base,
            structure.inner.chr_length(),
            structure.inner.ru_length,
            structure.inner.rus_per_hor,
            structure.inner.alphabet.clone(),
        );

        let h1 = Haplotype::from_chromosomes(vec![chr1]);
        let h2 = Haplotype::from_chromosomes(vec![chr2]);

        individuals.push(Individual::new(format!("ind{}", i), h1, h2));
    }

    Ok(PyPopulation {
        inner: Population::new("initial_pop", individuals),
    })
}

/// Mutation configuration for simulations.
#[pyclass(name = "MutationConfig")]
#[derive(Clone)]
struct PyMutationConfig {
    inner: MutationConfig,
}

#[pymethods]
impl PyMutationConfig {
    #[staticmethod]
    fn uniform(alphabet: &PyAlphabet, rate: f64) -> PyResult<Self> {
        let config = MutationConfig::uniform(alphabet.inner.clone(), rate)
            .map_err(|e| PyValueError::new_err(format!("Failed to create mutation config: {}", e)))?;
        Ok(Self { inner: config })
    }

    fn __repr__(&self) -> String {
        "MutationConfig(...)".to_string()
    }
}

/// Recombination configuration for simulations.
#[pyclass(name = "RecombinationConfig")]
#[derive(Clone)]
struct PyRecombinationConfig {
    inner: RecombinationConfig,
}

#[pymethods]
impl PyRecombinationConfig {
    #[staticmethod]
    fn standard(break_prob: f64, crossover_prob: f64, gc_extension_prob: f64) -> PyResult<Self> {
        let config = RecombinationConfig::standard(break_prob, crossover_prob, gc_extension_prob)
            .map_err(|e| PyValueError::new_err(format!("Failed to create recombination config: {}", e)))?;
        Ok(Self { inner: config })
    }

    fn __repr__(&self) -> String {
        "RecombinationConfig(...)".to_string()
    }
}

/// Fitness configuration for simulations.
#[pyclass(name = "FitnessConfig")]
#[derive(Clone)]
struct PyFitnessConfig {
    inner: FitnessConfig,
}

#[pymethods]
impl PyFitnessConfig {
    #[staticmethod]
    fn neutral() -> Self {
        Self {
            inner: FitnessConfig::neutral(),
        }
    }

    fn __repr__(&self) -> String {
        if self.inner.is_neutral() {
            "FitnessConfig.neutral()".to_string()
        } else {
            "FitnessConfig(...)".to_string()
        }
    }
}

/// Simulation engine for evolutionary processes.
#[pyclass(name = "Simulation", unsendable)]
struct PySimulation {
    inner: Option<Simulation>,
}

#[pymethods]
impl PySimulation {
    /// Create a new simulation with uniform initial sequences.
    #[new]
    #[pyo3(signature = (structure, mutation, recombination, fitness, config))]
    fn new(
        structure: &PyRepeatStructure,
        mutation: &PyMutationConfig,
        recombination: &PyRecombinationConfig,
        fitness: &PyFitnessConfig,
        config: &PySimulationConfig,
    ) -> PyResult<Self> {
        let sim = Simulation::new(
            structure.inner.clone(),
            mutation.inner.clone(),
            recombination.inner.clone(),
            fitness.inner.clone(),
            config.inner.clone(),
        )
        .map_err(|e| PyRuntimeError::new_err(format!("Failed to create simulation: {}", e)))?;

        Ok(Self { inner: Some(sim) })
    }

    /// Create a simulation from custom sequences.
    ///
    /// Args:
    ///     source_type: Type of input ("fasta", "json", or "database")
    ///     source_path: Path to file or JSON string
    ///     structure: Repeat structure configuration
    ///     mutation: Mutation configuration
    ///     recombination: Recombination configuration
    ///     fitness: Fitness configuration
    ///     config: Simulation configuration
    ///     sim_id: Simulation ID (required for "database" source)
    ///     generation: Generation to load (optional for "database" source)
    #[staticmethod]
    #[pyo3(signature = (source_type, source_path, structure, mutation, recombination, fitness, config, sim_id=None, generation=None))]
    fn from_sequences(
        source_type: &str,
        source_path: String,
        structure: &PyRepeatStructure,
        mutation: &PyMutationConfig,
        recombination: &PyRecombinationConfig,
        fitness: &PyFitnessConfig,
        config: &PySimulationConfig,
        sim_id: Option<String>,
        generation: Option<usize>,
    ) -> PyResult<Self> {
        let source = match source_type {
            "fasta" => SequenceInput::Fasta(source_path),
            "json" => SequenceInput::Json(source_path),
            "database" => {
                let id = sim_id.ok_or_else(|| {
                    PyValueError::new_err("sim_id required for database source")
                })?;
                SequenceInput::Database {
                    path: source_path,
                    sim_id: id,
                    generation,
                }
            }
            _ => {
                return Err(PyValueError::new_err(format!(
                    "Invalid source_type '{}'. Use 'fasta', 'json', or 'database'",
                    source_type
                )));
            }
        };

        let sim = Simulation::from_sequences(
            source,
            structure.inner.clone(),
            mutation.inner.clone(),
            recombination.inner.clone(),
            fitness.inner.clone(),
            config.inner.clone(),
        )
        .map_err(|e| PyRuntimeError::new_err(format!("Failed to create simulation: {}", e)))?;

        Ok(Self { inner: Some(sim) })
    }

    /// Resume a simulation from a checkpoint in the database.
    ///
    /// Args:
    ///     db_path: Path to the database file
    ///     sim_id: Simulation ID to resume
    #[staticmethod]
    fn from_checkpoint(db_path: String, sim_id: String) -> PyResult<Self> {
        let sim = Simulation::from_checkpoint(db_path, &sim_id)
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to resume simulation: {}", e)))?;

        Ok(Self { inner: Some(sim) })
    }

    /// Get the current population.
    fn population(&self) -> PyResult<PyPopulation> {
        let sim = self
            .inner
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("Simulation has been consumed"))?;

        Ok(PyPopulation {
            inner: sim.population().clone(),
        })
    }

    /// Get the current generation number.
    fn generation(&self) -> PyResult<usize> {
        let sim = self
            .inner
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("Simulation has been consumed"))?;

        Ok(sim.generation())
    }

    /// Advance simulation by one generation.
    fn step(&mut self) -> PyResult<()> {
        self.inner
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("Simulation has been consumed"))?
            .step()
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to step simulation: {}", e)))?;

        Ok(())
    }

    /// Run simulation for a specific number of generations.
    fn run_for(&mut self, generations: usize) -> PyResult<()> {
        self.inner
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("Simulation has been consumed"))?
            .run_for(generations)
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to run simulation: {}", e)))?;

        Ok(())
    }

    /// Run simulation for the configured number of generations.
    fn run(&mut self) -> PyResult<()> {
        self.inner
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("Simulation has been consumed"))?
            .run()
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to run simulation: {}", e)))?;

        Ok(())
    }

    fn __repr__(&self) -> String {
        if let Some(sim) = &self.inner {
            format!(
                "Simulation(generation={}, population_size={})",
                sim.generation(),
                sim.population().size()
            )
        } else {
            "Simulation(consumed)".to_string()
        }
    }
}
