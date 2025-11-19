use centrevo::base::Nucleotide;
use centrevo::evolution::{GCContentFitness, LengthFitness};
use centrevo::simulation::{
    FitnessConfig, MutationConfig, RecombinationConfig, RepeatStructure, SimulationConfig,
};

fn main() {
    println!("=== Testing Serde Serialization ===\n");

    // Test 1: SimulationConfig
    println!("1. SimulationConfig:");
    let sim_config = SimulationConfig::new(100, 1000, Some(42));
    let json = serde_json::to_string_pretty(&sim_config).unwrap();
    println!("{}", json);
    let _deserialized: SimulationConfig = serde_json::from_str(&json).unwrap();
    println!("✓ Round-trip successful\n");

    // Test 2: MutationConfig
    println!("2. MutationConfig:");
    let mutation_config = MutationConfig::uniform(0.001).unwrap();
    let json = serde_json::to_string_pretty(&mutation_config).unwrap();
    println!("{}", json);
    let _deserialized: MutationConfig = serde_json::from_str(&json).unwrap();
    println!("✓ Round-trip successful\n");

    // Test 3: RecombinationConfig
    println!("3. RecombinationConfig:");
    let recomb_config = RecombinationConfig::standard(0.1, 0.5, 0.9).unwrap();
    let json = serde_json::to_string_pretty(&recomb_config).unwrap();
    println!("{}", json);
    let _deserialized: RecombinationConfig = serde_json::from_str(&json).unwrap();
    println!("✓ Round-trip successful\n");

    // Test 4: FitnessConfig
    println!("4. FitnessConfig:");
    let fitness_config = FitnessConfig::new(
        Some(GCContentFitness::new(0.5, 1.0).unwrap()),
        Some(LengthFitness::new(1000, 1.0).unwrap()),
        None,
        None,
    );
    let json = serde_json::to_string_pretty(&fitness_config).unwrap();
    println!("{}", json);
    let _deserialized: FitnessConfig = serde_json::from_str(&json).unwrap();
    println!("✓ Round-trip successful\n");

    // Test 5: RepeatStructure
    println!("5. RepeatStructure:");
    let structure = RepeatStructure::new(Nucleotide::A, 171, 12, 100, 1);
    let json = serde_json::to_string_pretty(&structure).unwrap();
    println!("{}", json);
    let _deserialized: RepeatStructure = serde_json::from_str(&json).unwrap();
    println!("✓ Round-trip successful\n");

    println!("=== All Serialization Tests Passed! ===");
}
