use assert_cmd::Command;
use predicates::prelude::*;
use tempfile::tempdir;

#[test]
fn test_init_creates_database() {
    let temp = tempdir().unwrap();
    let db_path = temp.path().join("test_sim.db");

    let mut cmd = Command::cargo_bin("centrevo").unwrap();
    cmd.arg("init")
        .arg("--name")
        .arg("test_sim_default")
        .arg("--output")
        .arg(&db_path)
        .assert()
        .success()
        .stdout(predicate::str::contains(
            "Simulation initialized successfully!",
        ));

    assert!(db_path.exists());
}

#[test]
fn test_init_population_param() {
    let temp = tempdir().unwrap();
    let db_path = temp.path().join("test_pop.db");

    let mut cmd = Command::cargo_bin("centrevo").unwrap();
    cmd.arg("init")
        .arg("--name")
        .arg("test_pop_size")
        .arg("--output")
        .arg(&db_path)
        .arg("--population-size")
        .arg("10")
        .assert()
        .success()
        .stdout(predicate::str::contains("Population size: 10"));
}

#[test]
fn test_init_generations_param() {
    let temp = tempdir().unwrap();
    let db_path = temp.path().join("test_gen.db");

    let mut cmd = Command::cargo_bin("centrevo").unwrap();
    cmd.arg("init")
        .arg("--name")
        .arg("test_gen_count")
        .arg("--output")
        .arg(&db_path)
        .arg("--generations")
        .arg("50")
        .assert()
        .success()
        .stdout(predicate::str::contains("Generations: 50"));
}

#[test]
fn test_init_mutation_rate_param() {
    let temp = tempdir().unwrap();
    let db_path = temp.path().join("test_mut.db");

    let mut cmd = Command::cargo_bin("centrevo").unwrap();
    cmd.arg("init")
        .arg("--name")
        .arg("test_mut")
        .arg("--output")
        .arg(&db_path)
        .arg("--mutation-rate")
        .arg("1e-4")
        .assert()
        .success();
}

#[test]
fn test_init_advanced_mutation_params() {
    let temp = tempdir().unwrap();
    let db_path = temp.path().join("test_adv_mut.db");

    let mut cmd = Command::cargo_bin("centrevo").unwrap();
    cmd.arg("init")
        .arg("--name")
        .arg("test_adv_mut")
        .arg("--output")
        .arg(&db_path)
        .arg("--rate-ac")
        .arg("1e-6")
        .arg("--rate-ag")
        .arg("1e-6")
        .arg("--rate-at")
        .arg("1e-6")
        .arg("--rate-cg")
        .arg("1e-6")
        .arg("--rate-ct")
        .arg("1e-6")
        .arg("--rate-gt")
        .arg("1e-6")
        .assert()
        .success();
}

#[test]
fn test_init_error_partial_mutation_rates() {
    let temp = tempdir().unwrap();
    let db_path = temp.path().join("test_fail_mut.db");

    let mut cmd = Command::cargo_bin("centrevo").unwrap();
    cmd.arg("init")
        .arg("--name")
        .arg("test_fail_mut")
        .arg("--output")
        .arg(&db_path)
        .arg("--rate-ac")
        .arg("1e-6")
        // Missing other rates
        .assert()
        .failure()
        .stderr(predicate::str::contains("must be provided"));
}

#[test]
fn test_init_indel_params() {
    let temp = tempdir().unwrap();
    let db_path = temp.path().join("test_indel.db");

    let mut cmd = Command::cargo_bin("centrevo").unwrap();
    cmd.arg("init")
        .arg("--name")
        .arg("test_indel")
        .arg("--output")
        .arg(&db_path)
        .arg("--indel-ins-rate")
        .arg("0.01")
        .arg("--indel-del-rate")
        .arg("0.01")
        .assert()
        .success();
}

#[test]
fn test_init_recombination_params() {
    let temp = tempdir().unwrap();
    let db_path = temp.path().join("test_recomb.db");

    let mut cmd = Command::cargo_bin("centrevo").unwrap();
    cmd.arg("init")
        .arg("--name")
        .arg("test_recomb")
        .arg("--output")
        .arg(&db_path)
        .arg("--recomb-rate")
        .arg("0.05")
        .assert()
        .success();
}

#[test]
fn test_init_fitness_params() {
    let temp = tempdir().unwrap();
    let db_path = temp.path().join("test_fit.db");

    let mut cmd = Command::cargo_bin("centrevo").unwrap();
    cmd.arg("init")
        .arg("--name")
        .arg("test_fit")
        .arg("--output")
        .arg(&db_path)
        .arg("--fit-gc-opt")
        .arg("0.4")
        .arg("--fit-gc-conc")
        .arg("10.0")
        .assert()
        .success();
}

#[test]
fn test_init_defaults() {
    let temp = tempdir().unwrap();
    // Move to temp dir so default output "simulation.db" is written there
    // But Command::cargo_bin executes the binary, so we need to set current_dir
    let mut cmd = Command::cargo_bin("centrevo").unwrap();
    cmd.current_dir(&temp)
        .arg("init")
        .assert()
        .success()
        .stdout(predicate::str::contains("Name: simulation"));

    assert!(temp.path().join("simulation.db").exists());
}

#[test]
fn test_run_execution() {
    let temp = tempdir().unwrap();
    let db_path = temp.path().join("test_run.db");

    // Init first
    let mut cmd_init = Command::cargo_bin("centrevo").unwrap();
    cmd_init
        .arg("init")
        .arg("--name")
        .arg("test_run_exec")
        .arg("--output")
        .arg(&db_path)
        .arg("--generations")
        .arg("2")
        .arg("--population-size")
        .arg("5")
        // Use small values
        .arg("--ru-length")
        .arg("10")
        .arg("--rus-per-hor")
        .arg("2")
        .arg("--hors-per-chr")
        .arg("2")
        .assert()
        .success();

    // Then run
    let mut cmd_run = Command::cargo_bin("centrevo").unwrap();
    cmd_run
        .arg("run")
        .arg("--name")
        .arg("test_run_exec")
        .arg("--database")
        .arg(&db_path)
        .assert()
        .success()
        .stdout(predicate::str::contains("Simulation complete!"));
}

#[test]
fn test_run_error_missing_db() {
    let temp = tempdir().unwrap();
    let db_path = temp.path().join("non_existent.db");

    let mut cmd = Command::cargo_bin("centrevo").unwrap();
    cmd.arg("run")
        .arg("--name")
        .arg("test_run_missing")
        .arg("--database")
        .arg(&db_path)
        .assert()
        .failure()
        .stderr(predicate::str::contains("Failed to load configuration"));
}

#[test]
fn test_setup_flag_defaults() {
    // This requires interaction unless defaults is true
    // But setup creates a file in current dir if not specified, or asks validation.
    // Setup mostly wraps init but interactively.
    // verify setup --defaults runs init non-interactively.

    let temp = tempdir().unwrap();
    let _db_path = temp.path().join("setup_defaults.db");

    // We can't easy change current dir for the arbitrary output of setup unless setup supports --output arg which init does but setup wizard asks for it.
    // setup --defaults uses default name 'simulation' and output 'simulation.db' in CWD.
    // To avoid cluttering repo, we should run this in temp dir?
    // assert_cmd current_dir can set CWD.

    let mut cmd = Command::cargo_bin("centrevo").unwrap();
    cmd.current_dir(&temp)
        .arg("setup")
        .arg("--defaults")
        .assert()
        .success();

    assert!(temp.path().join("simulation.db").exists());
}

#[test]
fn test_list_displays_simulations() {
    let temp = tempdir().unwrap();
    let db_path = temp.path().join("test_list.db");

    let mut cmd_init = Command::cargo_bin("centrevo").unwrap();
    cmd_init
        .arg("init")
        .arg("--name")
        .arg("my_simulation_name")
        .arg("--output")
        .arg(&db_path)
        .assert()
        .success();

    let mut cmd_list = Command::cargo_bin("centrevo").unwrap();
    cmd_list
        .arg("list")
        .arg("--database")
        .arg(&db_path)
        .assert()
        .success()
        .stdout(predicate::str::contains("my_simulation_name"));
}
