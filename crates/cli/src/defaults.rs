//! Shared default values for simulation configuration.
//! These values are used by both the `init` command (via clap) and the `setup` wizard.

pub const SIMULATION_NAME: &str = "simulation";
pub const OUTPUT_DB: &str = "simulation.db";

pub const POPULATION_SIZE: usize = 100;
pub const GENERATIONS: usize = 1000;

pub const RU_LENGTH: usize = 171;
pub const RUS_PER_HOR: usize = 12;
pub const HORS_PER_CHR: usize = 100;
pub const CHRS_PER_HAP: usize = 1;

// Evolution
pub const MUTATION_RATE: f64 = 1e-5;
// Specific rates default to 1e-6 if manual entry is chosen but fields are empty,
// though typically optional fields don't have defaults in clap.
// We define them here for the setup wizard prompt defaults.
pub const SPECIFIC_RATE_DEFAULT: f64 = 1e-6;

pub const INDEL_VS_SUB_RATIO: f64 = 0.1; // Not used directly, but maybe useful context
pub const INDEL_INS_RATE: f64 = 0.0;
pub const INDEL_DEL_RATE: f64 = 0.0;
pub const INDEL_LENGTH_P: f64 = 0.5;

pub const RECOMB_RATE: f64 = 1e-6;
pub const CROSSOVER_PROB: f64 = 0.01;
pub const GC_EXTENSION_PROB: f64 = 0.95;
pub const HOMOLOGY_STRENGTH: f64 = 5.0;
pub const SEARCH_WINDOW: usize = 100;

// Fitness (defaults for prompts, not CLI args which are optional)
pub const FIT_GC_OPT: f64 = 0.4;
pub const FIT_GC_CONC: f64 = 10.0;
pub const FIT_LEN_STD: f64 = 500.0; // Optimal length derived from structure
pub const FIT_SEQ_SIM_SHAPE: f64 = 2.0;
pub const FIT_LEN_SIM_SHAPE: f64 = 2.0;

pub const RECORD_EVERY: usize = 100;
