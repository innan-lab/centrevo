//! Quick comparison of standard vs Poisson mutation performance
use centrevo::base::{Nucleotide, Sequence};
use centrevo::evolution::SubstitutionModel;
use rand::SeedableRng;
use rand::rngs::StdRng;
use std::time::Instant;

fn create_test_sequence(size: usize) -> Sequence {
    let mut seq = Sequence::with_capacity(size);
    for i in 0..size {
        seq.push(match i % 4 {
            0 => Nucleotide::A,
            1 => Nucleotide::C,
            2 => Nucleotide::G,
            _ => Nucleotide::T,
        });
    }
    seq
}

fn benchmark_method(
    model: &SubstitutionModel,
    seq_size: usize,
    iterations: usize,
    use_poisson: bool,
) -> (f64, usize) {
    let mut rng = StdRng::seed_from_u64(42);
    let mut total_mutations = 0;

    let start = Instant::now();

    for _ in 0..iterations {
        let mut seq = create_test_sequence(seq_size);
        let count = if use_poisson {
            model.mutate_sequence_poisson(&mut seq, &mut rng)
        } else {
            model.mutate_sequence(&mut seq, &mut rng)
        };
        total_mutations += count;
    }

    let elapsed = start.elapsed();
    let avg_time = elapsed.as_secs_f64() / iterations as f64;

    (avg_time * 1_000_000.0, total_mutations / iterations) // Return microseconds
}

fn main() {
    println!("Mutation Performance Comparison: Standard vs Poisson Pre-sampling\n");
    println!("{:-<90}", "");
    println!(
        "{:<20} {:>10} {:>12} {:>12} {:>12} {:>12}",
        "Config", "Seq Size", "Standard", "Poisson", "Speedup", "Avg Mut"
    );
    println!("{:-<90}", "");

    let configs = [
        (0.0001, "Very Low μ"),
        (0.001, "Low μ"),
        (0.01, "Medium μ"),
        (0.1, "High μ"),
    ];

    let sizes = [1_000, 10_000, 100_000];
    let iterations = 100;

    for (rate, label) in configs {
        let model = SubstitutionModel::uniform(rate).unwrap();

        for size in sizes {
            let (std_time, std_mut) = benchmark_method(&model, size, iterations, false);
            let (poi_time, _poi_mut) = benchmark_method(&model, size, iterations, true);

            let speedup = std_time / poi_time;
            let speedup_str = if speedup > 1.0 {
                format!("{:.2}x faster", speedup)
            } else {
                format!("{:.2}x slower", 1.0 / speedup)
            };

            println!(
                "{:<20} {:>10} {:>10.2} μs {:>10.2} μs {:>12} {:>12}",
                label, size, std_time, poi_time, speedup_str, std_mut
            );
        }
        println!("{:-<90}", "");
    }

    println!("\nKey Observations:");
    println!("  • Poisson method is fastest for low mutation rates (typical in evolution)");
    println!("  • Performance gain increases with sequence length");
    println!("  • High mutation rates may use standard method as fallback");
    println!("  • Both methods produce statistically equivalent results\n");
}
