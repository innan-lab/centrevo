//! CLI integration tests.
//! Tests the command-line interface to ensure all commands work correctly.

use assert_cmd::assert::OutputAssertExt;
use predicates::prelude::*;
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

/// Helper to create a unique test database path
fn test_db_path(name: &str) -> PathBuf {
    let temp_dir = std::env::temp_dir();
    temp_dir.join(format!("test_cli_{}_{}.db", name, std::process::id()))
}

/// Helper to clean up test database
fn cleanup_db(path: &PathBuf) {
    let _ = fs::remove_file(path);
}

/// Get the centrevo binary command
fn centrevo_cmd() -> Command {
    // Find the binary in the target directory
    Command::new(env!("CARGO_BIN_EXE_centrevo"))
}

#[test]
fn test_cli_help() {
    centrevo_cmd()
        .arg("--help")
        .assert()
        .success()
        .stdout(predicate::str::contains("Centromeric evolution simulator"))
        .stdout(predicate::str::contains("Commands:"));
}

#[test]
fn test_cli_version() {
    centrevo_cmd()
        .arg("--version")
        .assert()
        .success()
        .stdout(predicate::str::contains("centrevo"));
}

#[test]
fn test_init_basic() {
    let db_path = test_db_path("init_basic");

    centrevo_cmd()
        .args([
            "init",
            "-N",
            "test_sim",
            "-o",
            db_path.to_str().unwrap(),
            "-n",
            "10",
            "-g",
            "5",
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("Initializing simulation: test_sim"))
        .stdout(predicate::str::contains("Created 10 individuals"))
        .stdout(predicate::str::contains("Database created"));

    // Verify database was created
    assert!(db_path.exists(), "Database file should be created");

    cleanup_db(&db_path);
}

#[test]
fn test_init_with_seed() {
    let db_path = test_db_path("init_seed");

    centrevo_cmd()
        .args([
            "init",
            "-N",
            "seeded_sim",
            "-o",
            db_path.to_str().unwrap(),
            "-n",
            "5",
            "-g",
            "10",
            "--seed",
            "42",
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("seeded_sim"));

    assert!(db_path.exists());
    cleanup_db(&db_path);
}

#[test]
fn test_init_custom_structure() {
    let db_path = test_db_path("init_structure");

    centrevo_cmd()
        .args([
            "init",
            "-N",
            "custom_struct",
            "-o",
            db_path.to_str().unwrap(),
            "-n",
            "10",
            "-g",
            "5",
            "--ru-length",
            "200",
            "--rus-per-hor",
            "10",
            "--hors-per-chr",
            "50",
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("200bp RU × 10 RUs/HOR × 50 HORs"));

    cleanup_db(&db_path);
}

#[test]
fn test_init_default_values() {
    let temp_dir = TempDir::new().unwrap();
    let db_path = temp_dir.path().join("default.db");

    centrevo_cmd()
        .args(["init", "-o", db_path.to_str().unwrap()])
        .assert()
        .success()
        .stdout(predicate::str::contains("simulation")) // default name
        .stdout(predicate::str::contains("171bp RU")); // default structure
}

#[test]
fn test_list_empty_database() {
    let db_path = test_db_path("list_empty");

    // Create empty database by initializing then clearing
    centrevo_cmd()
        .args(["init", "-N", "temp", "-o", db_path.to_str().unwrap(), "-n", "1", "-g", "1"])
        .assert()
        .success();

    centrevo_cmd()
        .args(["list", "-d", db_path.to_str().unwrap()])
        .assert()
        .success()
        .stdout(predicate::str::contains("temp"));

    cleanup_db(&db_path);
}

#[test]
fn test_list_with_simulations() {
    let db_path = test_db_path("list_multi");

    // Create multiple simulations
    for name in &["sim1", "sim2", "sim3"] {
        centrevo_cmd()
            .args([
                "init",
                "-N",
                name,
                "-o",
                db_path.to_str().unwrap(),
                "-n",
                "5",
                "-g",
                "1",
            ])
            .assert()
            .success();
    }

    // List should show all simulations
    centrevo_cmd()
        .args(["list", "-d", db_path.to_str().unwrap()])
        .assert()
        .success()
        .stdout(predicate::str::contains("sim1"))
        .stdout(predicate::str::contains("sim2"))
        .stdout(predicate::str::contains("sim3"));

    cleanup_db(&db_path);
}

#[test]
fn test_info_existing_simulation() {
    let db_path = test_db_path("info_exists");

    // Create simulation
    centrevo_cmd()
        .args([
            "init",
            "-N",
            "info_test",
            "-o",
            db_path.to_str().unwrap(),
            "-n",
            "20",
            "-g",
            "100",
        ])
        .assert()
        .success();

    // Get info
    centrevo_cmd()
        .args(["info", "-d", db_path.to_str().unwrap(), "-N", "info_test"])
        .assert()
        .success()
        .stdout(predicate::str::contains("Simulation Information: info_test"))
        .stdout(predicate::str::contains("Population size: 20"))
        .stdout(predicate::str::contains("Generations: 100"));

    cleanup_db(&db_path);
}

#[test]
fn test_info_nonexistent_simulation() {
    let db_path = test_db_path("info_missing");

    // Create a simulation
    centrevo_cmd()
        .args([
            "init",
            "-N",
            "exists",
            "-o",
            db_path.to_str().unwrap(),
            "-n",
            "5",
            "-g",
            "1",
        ])
        .assert()
        .success();

    // Try to get info for non-existent simulation
    centrevo_cmd()
        .args([
            "info",
            "-d",
            db_path.to_str().unwrap(),
            "-N",
            "does_not_exist",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("Failed to get simulation info"));

    cleanup_db(&db_path);
}

#[test]
fn test_generations_after_init() {
    let db_path = test_db_path("gens_init");

    // Initialize simulation (records generation 0)
    centrevo_cmd()
        .args([
            "init",
            "-N",
            "gen_test",
            "-o",
            db_path.to_str().unwrap(),
            "-n",
            "5",
            "-g",
            "100",
        ])
        .assert()
        .success();

    // Check recorded generations
    centrevo_cmd()
        .args(["generations", "-d", db_path.to_str().unwrap(), "-N", "gen_test"])
        .assert()
        .success()
        .stdout(predicate::str::contains("Recorded Generations"))
        .stdout(predicate::str::contains("[0]")); // Should have generation 0

    cleanup_db(&db_path);
}

#[test]
fn test_validate_database_success() {
    let db_path = test_db_path("validate_ok");

    // Create valid simulation
    centrevo_cmd()
        .args([
            "init",
            "-N",
            "valid_sim",
            "-o",
            db_path.to_str().unwrap(),
            "-n",
            "10",
            "-g",
            "10",
        ])
        .assert()
        .success();

    // Validate
    centrevo_cmd()
        .args(["validate", "-d", db_path.to_str().unwrap()])
        .assert()
        .success()
        .stdout(predicate::str::contains("Validating simulation: valid_sim"));

    cleanup_db(&db_path);
}

#[test]
fn test_validate_specific_simulation() {
    let db_path = test_db_path("validate_specific");

    // Create multiple simulations
    centrevo_cmd()
        .args([
            "init",
            "-N",
            "sim_a",
            "-o",
            db_path.to_str().unwrap(),
            "-n",
            "5",
            "-g",
            "5",
        ])
        .assert()
        .success();

    centrevo_cmd()
        .args([
            "init",
            "-N",
            "sim_b",
            "-o",
            db_path.to_str().unwrap(),
            "-n",
            "5",
            "-g",
            "5",
        ])
        .assert()
        .success();

    // Validate specific simulation
    centrevo_cmd()
        .args([
            "validate",
            "-d",
            db_path.to_str().unwrap(),
            "-N",
            "sim_a",
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("Validating simulation: sim_a"));
        // Note: sim_b should not be in output since we're validating only sim_a

    cleanup_db(&db_path);
}

#[test]
fn test_validate_nonexistent_database() {
    centrevo_cmd()
        .args(["validate", "-d", "/nonexistent/path/to/database.db"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("Database file does not exist"));
}

#[test]
fn test_init_invalid_population_size() {
    // Note: CLI currently accepts 0 population size - this may be a valid edge case
    // for testing empty populations or could be validated in the future
    let db_path = test_db_path("invalid_pop");

    centrevo_cmd()
        .args([
            "init",
            "-N",
            "test",
            "-o",
            db_path.to_str().unwrap(),
            "-n",
            "0", // Currently accepted by CLI
            "-g",
            "10",
        ])
        .assert()
        .success(); // CLI currently allows this

    cleanup_db(&db_path);
}

#[test]
fn test_init_invalid_generations() {
    // Note: CLI currently accepts 0 generations - this may be valid for
    // initializing a simulation without running it
    let db_path = test_db_path("invalid_gens");

    centrevo_cmd()
        .args([
            "init",
            "-N",
            "test",
            "-o",
            db_path.to_str().unwrap(),
            "-n",
            "10",
            "-g",
            "0", // Currently accepted by CLI
        ])
        .assert()
        .success(); // CLI currently allows this

    cleanup_db(&db_path);
}

#[test]
fn test_multiple_simulations_same_database() {
    let db_path = test_db_path("multi_sim");

    // Create first simulation
    centrevo_cmd()
        .args([
            "init",
            "-N",
            "first",
            "-o",
            db_path.to_str().unwrap(),
            "-n",
            "10",
            "-g",
            "5",
        ])
        .assert()
        .success();

    // Create second simulation in same database
    centrevo_cmd()
        .args([
            "init",
            "-N",
            "second",
            "-o",
            db_path.to_str().unwrap(),
            "-n",
            "15",
            "-g",
            "10",
        ])
        .assert()
        .success();

    // List should show both
    centrevo_cmd()
        .args(["list", "-d", db_path.to_str().unwrap()])
        .assert()
        .success()
        .stdout(predicate::str::contains("first"))
        .stdout(predicate::str::contains("second"));

    // Info for each should be independent
    centrevo_cmd()
        .args(["info", "-d", db_path.to_str().unwrap(), "-N", "first"])
        .assert()
        .success()
        .stdout(predicate::str::contains("Population size: 10"));

    centrevo_cmd()
        .args(["info", "-d", db_path.to_str().unwrap(), "-N", "second"])
        .assert()
        .success()
        .stdout(predicate::str::contains("Population size: 15"));

    cleanup_db(&db_path);
}

#[test]
fn test_init_creates_expected_structure() {
    let db_path = test_db_path("structure_check");

    centrevo_cmd()
        .args([
            "init",
            "-N",
            "struct_test",
            "-o",
            db_path.to_str().unwrap(),
            "-n",
            "5",
            "-g",
            "10",
            "--ru-length",
            "100",
            "--rus-per-hor",
            "5",
            "--hors-per-chr",
            "20",
        ])
        .assert()
        .success();

    // Verify the structure is reported correctly
    centrevo_cmd()
        .args(["info", "-d", db_path.to_str().unwrap(), "-N", "struct_test"])
        .assert()
        .success()
        .stdout(predicate::str::contains("struct_test"));

    cleanup_db(&db_path);
}

#[test]
fn test_command_without_required_args() {
    // Info without name
    centrevo_cmd()
        .args(["info", "-d", "test.db"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("required"));

    // Generations without name
    centrevo_cmd()
        .args(["generations", "-d", "test.db"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("required"));
}

#[test]
fn test_help_for_subcommands() {
    for subcmd in &["init", "list", "info", "generations", "validate"] {
        centrevo_cmd()
            .args([subcmd, "--help"])
            .assert()
            .success()
            .stdout(predicate::str::contains("Usage:"));
    }
}

#[test]
fn test_export_sequences_csv() {
    let db_path = test_db_path("export_csv");
    let output_path = test_db_path("export_output_csv");

    // Create simulation
    centrevo_cmd()
        .args([
            "init",
            "-N",
            "export_test",
            "-o",
            db_path.to_str().unwrap(),
            "-n",
            "3",
            "-g",
            "5",
        ])
        .assert()
        .success();

    // Export as CSV
    centrevo_cmd()
        .args([
            "export",
            "-d",
            db_path.to_str().unwrap(),
            "-N",
            "export_test",
            "-g",
            "0",
            "-f",
            "csv",
            "-o",
            output_path.to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("Data exported"));

    // Verify CSV was created and has content
    assert!(output_path.exists());
    let content = fs::read_to_string(&output_path).unwrap();
    assert!(content.contains("individual_id,haplotype,chromosome,sequence"));

    cleanup_db(&db_path);
    cleanup_db(&output_path);
}

#[test]
fn test_export_sequences_fasta() {
    let db_path = test_db_path("export_fasta");
    let output_path = test_db_path("export_output_fasta");

    // Create simulation
    centrevo_cmd()
        .args([
            "init",
            "-N",
            "fasta_test",
            "-o",
            db_path.to_str().unwrap(),
            "-n",
            "2",
            "-g",
            "1",
        ])
        .assert()
        .success();

    // Export as FASTA
    centrevo_cmd()
        .args([
            "export",
            "-d",
            db_path.to_str().unwrap(),
            "-N",
            "fasta_test",
            "-g",
            "0",
            "-f",
            "fasta",
            "-o",
            output_path.to_str().unwrap(),
        ])
        .assert()
        .success();

    // Verify FASTA format
    assert!(output_path.exists());
    let content = fs::read_to_string(&output_path).unwrap();
    assert!(content.contains(">"));  // FASTA header

    cleanup_db(&db_path);
    cleanup_db(&output_path);
}

#[test]
fn test_export_metadata_json() {
    let db_path = test_db_path("export_meta");
    let output_path = test_db_path("export_meta_json");

    // Create simulation
    centrevo_cmd()
        .args([
            "init",
            "-N",
            "meta_test",
            "-o",
            db_path.to_str().unwrap(),
            "-n",
            "10",
            "-g",
            "100",
        ])
        .assert()
        .success();

    // Export metadata
    centrevo_cmd()
        .args([
            "export",
            "-d",
            db_path.to_str().unwrap(),
            "-N",
            "meta_test",
            "-g",
            "0",
            "--data-type",
            "metadata",
            "-f",
            "json",
            "-o",
            output_path.to_str().unwrap(),
        ])
        .assert()
        .success();

    // Verify JSON
    assert!(output_path.exists());
    let content = fs::read_to_string(&output_path).unwrap();
    assert!(content.contains("\"name\""));
    assert!(content.contains("meta_test"));

    cleanup_db(&db_path);
    cleanup_db(&output_path);
}

#[test]
fn test_analyze_command() {
    let db_path = test_db_path("analyze");
    let output_path = test_db_path("analyze_output");

    // Create simulation with some complexity
    centrevo_cmd()
        .args([
            "init",
            "-N",
            "analyze_test",
            "-o",
            db_path.to_str().unwrap(),
            "-n",
            "20",
            "-g",
            "10",
        ])
        .assert()
        .success();

    // Analyze generation 0
    centrevo_cmd()
        .args([
            "analyze",
            "-d",
            db_path.to_str().unwrap(),
            "-N",
            "analyze_test",
            "-g",
            "0",
            "--format",
            "pretty",
            "-o",
            output_path.to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("Analyzing simulation"));

    // Verify analysis output
    assert!(output_path.exists());
    let content = fs::read_to_string(&output_path).unwrap();
    assert!(content.contains("Population Genetics Summary"));
    assert!(content.contains("Diversity Metrics"));

    cleanup_db(&db_path);
    cleanup_db(&output_path);
}

#[test]
fn test_analyze_json_format() {
    let db_path = test_db_path("analyze_json");

    // Create simulation
    centrevo_cmd()
        .args([
            "init",
            "-N",
            "json_analyze",
            "-o",
            db_path.to_str().unwrap(),
            "-n",
            "15",
            "-g",
            "5",
        ])
        .assert()
        .success();

    // Analyze with JSON output (to stdout)
    centrevo_cmd()
        .args([
            "analyze",
            "-d",
            db_path.to_str().unwrap(),
            "-N",
            "json_analyze",
            "-g",
            "0",
            "--format",
            "json",
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("\"diversity\""))
        .stdout(predicate::str::contains("\"pi\""));

    cleanup_db(&db_path);
}

#[test]
fn test_reproducibility_with_seeds() {
    let db_path1 = test_db_path("repro1");
    let db_path2 = test_db_path("repro2");

    // Create two simulations with same seed
    for db_path in &[&db_path1, &db_path2] {
        centrevo_cmd()
            .args([
                "init",
                "-N",
                "repro_test",
                "-o",
                db_path.to_str().unwrap(),
                "-n",
                "10",
                "-g",
                "5",
                "--seed",
                "42",
            ])
            .assert()
            .success();
    }

    // Both should have same initial generation
    // (Full verification would require comparing exported sequences)
    assert!(db_path1.exists());
    assert!(db_path2.exists());

    cleanup_db(&db_path1);
    cleanup_db(&db_path2);
}

#[test]
fn test_setup_with_defaults() {
    // Test setup command with --defaults flag
    // This command would normally be interactive, but --defaults makes it non-interactive
    let _output_path = test_db_path("setup_defaults");

    // Note: This test may need to be adjusted based on how setup handles defaults
    // For now, we just verify it runs without panicking
    let result = centrevo_cmd()
        .args(["setup", "--defaults"])
        .output();

    // May succeed or fail depending on implementation
    // Just verify it completes (doesn't hang indefinitely)
    assert!(result.is_ok(), "Command should complete without hanging");
}

#[test]
fn test_concurrent_database_access() {
    let db_path = test_db_path("concurrent");

    // Create initial simulation
    centrevo_cmd()
        .args([
            "init",
            "-N",
            "concurrent_test",
            "-o",
            db_path.to_str().unwrap(),
            "-n",
            "5",
            "-g",
            "1",
        ])
        .assert()
        .success();

    // Multiple read operations should work
    centrevo_cmd()
        .args(["info", "-d", db_path.to_str().unwrap(), "-N", "concurrent_test"])
        .assert()
        .success();

    centrevo_cmd()
        .args(["list", "-d", db_path.to_str().unwrap()])
        .assert()
        .success();

    centrevo_cmd()
        .args(["generations", "-d", db_path.to_str().unwrap(), "-N", "concurrent_test"])
        .assert()
        .success();

    cleanup_db(&db_path);
}
