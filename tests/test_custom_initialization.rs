//! Integration tests for custom sequence initialization.

use centrevo::base::{Alphabet, Nucleotide};
use centrevo::simulation::{
    Simulation, RepeatStructure, MutationConfig, RecombinationConfig,
    FitnessConfig, SimulationConfig, SequenceInput, parse_fasta,
};
use centrevo::storage::Recorder;
use std::io::Write;
use tempfile::NamedTempFile;

fn create_test_structure() -> RepeatStructure {
    RepeatStructure::new(
        Alphabet::dna(),
        Nucleotide::A,
        10,  // ru_length
        5,   // rus_per_hor
        2,   // hors_per_chr
        1,   // chrs_per_hap (1 chromosome per haplotype)
    )
}

fn create_test_configs() -> (MutationConfig, RecombinationConfig, FitnessConfig, SimulationConfig) {
    let alphabet = Alphabet::dna();
    let mutation = MutationConfig::uniform(alphabet, 0.001).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
    let fitness = FitnessConfig::neutral();
    let config = SimulationConfig::new(2, 10, Some(42)); // 2 individuals
    
    (mutation, recombination, fitness, config)
}

#[test]
fn test_initialization_from_fasta() {
    let structure = create_test_structure();
    let chr_len = structure.chr_length(); // 100 bases
    
    // Create FASTA file with sequences for 2 individuals (4 sequences total)
    let mut fasta = NamedTempFile::new().unwrap();
    writeln!(fasta, ">ind0_h1_chr0").unwrap();
    writeln!(fasta, "{}", "ACGT".repeat(chr_len / 4)).unwrap();
    writeln!(fasta, ">ind0_h2_chr0").unwrap();
    writeln!(fasta, "{}", "TGCA".repeat(chr_len / 4)).unwrap();
    writeln!(fasta, ">ind1_h1_chr0").unwrap();
    writeln!(fasta, "{}", "AAAA".repeat(chr_len / 4)).unwrap();
    writeln!(fasta, ">ind1_h2_chr0").unwrap();
    writeln!(fasta, "{}", "CCCC".repeat(chr_len / 4)).unwrap();
    fasta.flush().unwrap();
    
    // Initialize simulation from FASTA
    let (mutation, recombination, fitness, config) = create_test_configs();
    let source = SequenceInput::Fasta(fasta.path().to_string_lossy().to_string());
    
    let sim = Simulation::from_sequences(
        source,
        structure,
        mutation,
        recombination,
        fitness,
        config,
    ).unwrap();
    
    // Verify population
    assert_eq!(sim.population().size(), 2);
    assert_eq!(sim.generation(), 0);
    
    // Verify individuals have correct structure
    for ind in sim.population().individuals() {
        assert_eq!(ind.haplotype1().len(), 1);
        assert_eq!(ind.haplotype2().len(), 1);
        
        let chr1 = ind.haplotype1().get(0).unwrap();
        let chr2 = ind.haplotype2().get(0).unwrap();
        
        assert_eq!(chr1.len(), chr_len);
        assert_eq!(chr2.len(), chr_len);
    }
}

#[test]
fn test_initialization_from_json() {
    let structure = create_test_structure();
    let chr_len = structure.chr_length();
    
    // Create JSON with sequences
    let json = format!(
        r#"[
            {{"id": "ind0_h1", "seq": "{}"}},
            {{"id": "ind0_h2", "seq": "{}"}},
            {{"id": "ind1_h1", "seq": "{}"}},
            {{"id": "ind1_h2", "seq": "{}"}}
        ]"#,
        "ACGT".repeat(chr_len / 4),
        "TGCA".repeat(chr_len / 4),
        "AAAA".repeat(chr_len / 4),
        "CCCC".repeat(chr_len / 4),
    );
    
    let (mutation, recombination, fitness, config) = create_test_configs();
    let source = SequenceInput::Json(json);
    
    let sim = Simulation::from_sequences(
        source,
        structure,
        mutation,
        recombination,
        fitness,
        config,
    ).unwrap();
    
    assert_eq!(sim.population().size(), 2);
}

#[test]
fn test_initialization_from_database() {
    let structure = create_test_structure();
    let (mutation, recombination, fitness, config) = create_test_configs();
    
    // Create a simulation and run it for a few generations
    let db_path = "/tmp/test_custom_init.db";
    let _ = std::fs::remove_file(db_path);
    
    let mut sim = Simulation::new(
        structure.clone(),
        mutation.clone(),
        recombination.clone(),
        fitness.clone(),
        config.clone(),
    ).unwrap();
    
    // Record initial state
    let mut recorder = Recorder::new(db_path, "test_sim", 
        centrevo::storage::RecordingStrategy::All).unwrap();
    
    // Record full config for resumability
    let snapshot = centrevo::storage::SimulationSnapshot {
        structure: structure.clone(),
        mutation: mutation.clone(),
        recombination: recombination.clone(),
        fitness: fitness.clone(),
        config: config.clone(),
    };
    recorder.record_full_config(&snapshot).unwrap();
    recorder.record_generation(sim.population(), 0).unwrap();
    
    // Run for a few generations
    for generation in 1..=5 {
        sim.step().unwrap();
        recorder.record_generation(sim.population(), generation).unwrap();
    }
    
    recorder.close().unwrap();
    
    // Now initialize a new simulation from the database
    let source = SequenceInput::Database {
        path: db_path.to_string(),
        sim_id: "test_sim".to_string(),
        generation: Some(5), // Use generation 5
    };
    
    let sim2 = Simulation::from_sequences(
        source,
        structure,
        mutation,
        recombination,
        fitness,
        SimulationConfig::new(2, 10, Some(100)), // Different seed for new sim
    ).unwrap();
    
    // Verify new simulation has correct population
    assert_eq!(sim2.population().size(), 2);
    assert_eq!(sim2.generation(), 0); // Starts at generation 0
    
    // Clean up
    std::fs::remove_file(db_path).ok();
}

#[test]
fn test_validation_wrong_sequence_count() {
    let structure = create_test_structure();
    let chr_len = structure.chr_length();
    
    // Only 2 sequences instead of 4
    let json = format!(
        r#"[
            {{"id": "ind0_h1", "seq": "{}"}},
            {{"id": "ind0_h2", "seq": "{}"}}
        ]"#,
        "ACGT".repeat(chr_len / 4),
        "TGCA".repeat(chr_len / 4),
    );
    
    let (mutation, recombination, fitness, config) = create_test_configs();
    let source = SequenceInput::Json(json);
    
    let result = Simulation::from_sequences(
        source,
        structure,
        mutation,
        recombination,
        fitness,
        config,
    );
    
    assert!(result.is_err());
    let err_msg = result.unwrap_err();
    assert!(err_msg.contains("Expected 4 sequences"));
}

#[test]
fn test_validation_wrong_sequence_length() {
    let structure = create_test_structure();
    
    // Wrong length sequences
    let json = r#"[
        {"id": "ind0_h1", "seq": "ACGT"},
        {"id": "ind0_h2", "seq": "TGCA"},
        {"id": "ind1_h1", "seq": "AAAA"},
        {"id": "ind1_h2", "seq": "CCCC"}
    ]"#;
    
    let (mutation, recombination, fitness, config) = create_test_configs();
    let source = SequenceInput::Json(json.to_string());
    
    let result = Simulation::from_sequences(
        source,
        structure,
        mutation,
        recombination,
        fitness,
        config,
    );
    
    assert!(result.is_err());
    let err_msg = result.unwrap_err();
    assert!(err_msg.contains("has length 4, expected 100"));
}

#[test]
fn test_validation_invalid_base() {
    let structure = create_test_structure();
    let chr_len = structure.chr_length(); // 100 bases
    
    // Create a sequence with 'N' - make sure sequence is exactly chr_len
    // chr_len = 100, so we need 100 bases total
    // "ACGT" = 4 bases, repeated 24 times = 96 bases, then add "ACGN" = 100 bases
    let mut seq = "ACGT".repeat(24);
    seq.push_str("ACGN"); // Add 4 more chars to make it exactly 100, with 'N' being invalid
    
    assert_eq!(seq.len(), chr_len, "Test sequence should be exactly {} bases", chr_len);
    
    let json = format!(
        r#"[
            {{"id": "ind0_h1", "seq": "{}"}},
            {{"id": "ind0_h2", "seq": "{}"}},
            {{"id": "ind1_h1", "seq": "{}"}},
            {{"id": "ind1_h2", "seq": "{}"}}
        ]"#,
        seq,
        "TGCA".repeat(chr_len / 4),
        "AAAA".repeat(chr_len / 4),
        "CCCC".repeat(chr_len / 4),
    );
    
    let (mutation, recombination, fitness, config) = create_test_configs();
    let source = SequenceInput::Json(json);
    
    let result = Simulation::from_sequences(
        source,
        structure,
        mutation,
        recombination,
        fitness,
        config,
    );
    
    assert!(result.is_err());
    let err_msg = result.unwrap_err();
    // The error should mention either "invalid" or "Invalid"
    assert!(err_msg.to_lowercase().contains("invalid"), "Expected 'invalid' in error message, got: {}", err_msg);
}

#[test]
fn test_fasta_parse_multiline_sequences() {
    let structure = create_test_structure();
    let chr_len = structure.chr_length(); // 100 bases
    
    // Create FASTA with multiline sequences
    let mut fasta = NamedTempFile::new().unwrap();
    // Split first sequence into 4 lines of 25 bases each
    let seq_line = "ACGT".repeat(25 / 4); // 24 bases - repeats "ACGT" 6 times
    
    writeln!(fasta, ">ind0_h1_chr0").unwrap();
    // Write 4 lines to get close to 100 (96 bases)
    for _ in 0..4 {
        writeln!(fasta, "{}", seq_line).unwrap();
    }
    // Add remaining 4 bases to get exactly 100
    writeln!(fasta, "ACGT").unwrap();
    
    writeln!(fasta, ">ind0_h2_chr0").unwrap();
    writeln!(fasta, "{}", "TGCA".repeat(chr_len / 4)).unwrap();
    writeln!(fasta, ">ind1_h1_chr0").unwrap();
    writeln!(fasta, "{}", "AAAA".repeat(chr_len / 4)).unwrap();
    writeln!(fasta, ">ind1_h2_chr0").unwrap();
    writeln!(fasta, "{}", "CCCC".repeat(chr_len / 4)).unwrap();
    fasta.flush().unwrap();
    
    let entries = parse_fasta(fasta.path()).unwrap();
    assert_eq!(entries.len(), 4);
    // First sequence should combine all lines
    assert_eq!(entries[0].seq.len(), chr_len, "First sequence should be {} bases", chr_len);
}

#[test]
fn test_simulation_runs_with_custom_sequences() {
    let structure = create_test_structure();
    let chr_len = structure.chr_length();
    
    // Create JSON with sequences
    let json = format!(
        r#"[
            {{"id": "ind0_h1", "seq": "{}"}},
            {{"id": "ind0_h2", "seq": "{}"}},
            {{"id": "ind1_h1", "seq": "{}"}},
            {{"id": "ind1_h2", "seq": "{}"}}
        ]"#,
        "ACGT".repeat(chr_len / 4),
        "TGCA".repeat(chr_len / 4),
        "AAAA".repeat(chr_len / 4),
        "CCCC".repeat(chr_len / 4),
    );
    
    let (mutation, recombination, fitness, config) = create_test_configs();
    let source = SequenceInput::Json(json);
    
    let mut sim = Simulation::from_sequences(
        source,
        structure,
        mutation,
        recombination,
        fitness,
        config,
    ).unwrap();
    
    // Run simulation for a few generations
    sim.run_for(5).unwrap();
    
    // Verify it ran successfully
    assert_eq!(sim.generation(), 5);
    assert_eq!(sim.population().size(), 2);
}
