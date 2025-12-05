use centrevo_sim::base::Sequence;
use centrevo_sim::evolution::SubstitutionModel;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::str::FromStr;
use std::time::Instant;

fn main() {
    let len = 10_000_000;
    let rate = 1e-5;

    println!("Generating sequence of length {}...", len);
    // Create a random sequence to be more realistic (though 'A's would show the same issue)
    // Actually, let's just use 'A's for speed of generation in this script
    let seq_str = "A".repeat(len);
    let sequence = Sequence::from_str(&seq_str).unwrap();

    let model = SubstitutionModel::uniform(rate).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

    println!("Benchmarking Direct Mutation...");
    let start = Instant::now();
    let mut seq_direct = sequence.clone();
    let count_direct = model.mutate_sequence(&mut seq_direct, &mut rng);
    let duration_direct = start.elapsed();
    println!("Direct: {:?} ({} mutations)", duration_direct, count_direct);

    println!("Benchmarking Sparse Mutation...");
    let start = Instant::now();
    let mut seq_sparse = sequence.clone();
    let count_sparse = model.mutate_sequence_sparse(&mut seq_sparse, &mut rng);
    let duration_sparse = start.elapsed();
    println!("Sparse: {:?} ({} mutations)", duration_sparse, count_sparse);

    let speedup = duration_direct.as_secs_f64() / duration_sparse.as_secs_f64();
    println!("Speedup: {:.2}x", speedup);
}
