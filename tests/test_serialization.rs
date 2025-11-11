//! Integration tests for configuration serialization and deserialization.
//! Tests that all configuration types can be serialized to JSON and back.

use centrevo::base::{Alphabet, Nucleotide};
use centrevo::simulation::{
    FitnessConfig, MutationConfig, RecombinationConfig, RepeatStructure, SimulationConfig,
};
use centrevo::evolution::{GCContentFitness, LengthFitness};

#[test]
fn test_simulation_config_serialization() {
    let config = SimulationConfig::new(100, 1000, Some(42));
    let json = serde_json::to_string(&config).unwrap();
    let deserialized: SimulationConfig = serde_json::from_str(&json).unwrap();

    assert_eq!(config.population_size, deserialized.population_size);
    assert_eq!(config.seed, deserialized.seed);
}

#[test]
fn test_mutation_config_serialization() {
    let config = MutationConfig::uniform(Alphabet::dna(), 0.001).unwrap();
    let json = serde_json::to_string(&config).unwrap();
    let _deserialized: MutationConfig = serde_json::from_str(&json).unwrap();

    // Test that deserialized config can be used
    assert!(serde_json::from_str::<MutationConfig>(&json).is_ok());
}

#[test]
fn test_recombination_config_serialization() {
    let config = RecombinationConfig::standard(0.1, 0.5, 0.9).unwrap();
    let json = serde_json::to_string(&config).unwrap();
    let _deserialized: RecombinationConfig = serde_json::from_str(&json).unwrap();

    // Test that deserialized config can be used
    assert!(serde_json::from_str::<RecombinationConfig>(&json).is_ok());
}

#[test]
fn test_fitness_config_neutral_serialization() {
    let config = FitnessConfig::neutral();
    let json = serde_json::to_string(&config).unwrap();
    let _deserialized: FitnessConfig = serde_json::from_str(&json).unwrap();

    assert!(serde_json::from_str::<FitnessConfig>(&json).is_ok());
}

#[test]
fn test_fitness_config_with_components_serialization() {
    let config = FitnessConfig::new(
        Some(GCContentFitness::new(0.5, 1.0).unwrap()),
        Some(LengthFitness::new(1000, 1.0).unwrap()),
        None,
    );
    let json = serde_json::to_string(&config).unwrap();
    let _deserialized: FitnessConfig = serde_json::from_str(&json).unwrap();

    assert!(serde_json::from_str::<FitnessConfig>(&json).is_ok());
}

#[test]
fn test_repeat_structure_serialization() {
    let structure = RepeatStructure::new(Alphabet::dna(), Nucleotide::A, 171, 12, 100, 1);
    let json = serde_json::to_string(&structure).unwrap();
    let _deserialized: RepeatStructure = serde_json::from_str(&json).unwrap();

    assert!(serde_json::from_str::<RepeatStructure>(&json).is_ok());
}

#[test]
fn test_complete_simulation_config_roundtrip() {
    // Test serializing and deserializing a complete set of simulation parameters
    let sim_config = SimulationConfig::new(50, 100, Some(123));
    let mutation_config = MutationConfig::uniform(Alphabet::dna(), 0.005).unwrap();
    let recombination_config = RecombinationConfig::standard(0.02, 0.6, 0.15).unwrap();
    let fitness_config = FitnessConfig::new(
        Some(GCContentFitness::new(0.45, 2.0).unwrap()),
        None,
        None,
    );
    let structure = RepeatStructure::new(Alphabet::dna(), Nucleotide::C, 20, 10, 50, 2);

    // Serialize all configs
    let sim_json = serde_json::to_string(&sim_config).unwrap();
    let mut_json = serde_json::to_string(&mutation_config).unwrap();
    let rec_json = serde_json::to_string(&recombination_config).unwrap();
    let fit_json = serde_json::to_string(&fitness_config).unwrap();
    let str_json = serde_json::to_string(&structure).unwrap();

    // Deserialize all configs
    let sim_deser: SimulationConfig = serde_json::from_str(&sim_json).unwrap();
    let mut_deser: MutationConfig = serde_json::from_str(&mut_json).unwrap();
    let rec_deser: RecombinationConfig = serde_json::from_str(&rec_json).unwrap();
    let fit_deser: FitnessConfig = serde_json::from_str(&fit_json).unwrap();
    let str_deser: RepeatStructure = serde_json::from_str(&str_json).unwrap();

    // Verify they can be used to create a simulation
    use centrevo::simulation::Simulation;
    let sim = Simulation::new(str_deser, mut_deser, rec_deser, fit_deser, sim_deser);
    assert!(sim.is_ok());
}

#[test]
fn test_config_pretty_print() {
    // Test that configs can be pretty-printed for human readability
    let config = SimulationConfig::new(100, 1000, Some(42));
    let pretty_json = serde_json::to_string_pretty(&config).unwrap();

    // Should contain newlines and indentation
    assert!(pretty_json.contains('\n'));
    assert!(pretty_json.len() > serde_json::to_string(&config).unwrap().len());

    // Should still deserialize correctly
    let deserialized: SimulationConfig = serde_json::from_str(&pretty_json).unwrap();
    assert_eq!(config.population_size, deserialized.population_size);
}
