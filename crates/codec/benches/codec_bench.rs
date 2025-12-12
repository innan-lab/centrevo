use centrevo_codec::{BitPackedRS, Codec, CodecError, ParallelBitPackedRS, UnpackedRS};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rand::Rng;
use std::hint::black_box;

// Baseline: Unpacked No RS (Simulates raw Vec<u8> storage/access)
struct UnpackedNoRS;
impl Codec for UnpackedNoRS {
    fn encode(&self, seq: &[u8]) -> Result<Vec<u8>, CodecError> {
        Ok(seq.to_vec())
    }
    fn decode(&self, data: &[u8]) -> Result<Vec<u8>, CodecError> {
        Ok(data.to_vec())
    }
}

fn bench_codecs(c: &mut Criterion) {
    let mut rng = rand::thread_rng();

    // Strategy Definitions
    let strategies: Vec<(&str, Box<dyn Codec>)> = vec![
        ("UnpackedNoRS", Box::new(UnpackedNoRS)),
        ("UnpackedRS", Box::new(UnpackedRS)),
        ("BitPackedRS", Box::new(BitPackedRS)),
        ("ParallelBitPackedRS", Box::new(ParallelBitPackedRS)),
    ];

    // Test Cases (Name, Size)
    // Block alignment = 14,272 bytes (64 * 223)
    let block_size = 14_272;

    let sizes = vec![
        ("Small_Worst", 1),                   // 99.9% padding
        ("Small_Best", block_size),           // 0% padding
        ("Medium_Best", 7 * block_size),      // Exact fit
        ("Medium_Worst", 7 * block_size + 1), // Just over, forces 2x shard size padding
        ("Large_Best", 70 * block_size),
        ("Large_Worst", 70 * block_size + 1),
        ("Huge_Best", 700 * block_size),   // ~10MB
        ("Giant_Best", 3500 * block_size), // ~50MB (Limited to keep bench time reasonable)
    ];

    // Group by Operation then Size
    for (size_name, size) in sizes {
        // Generate random input data
        let input: Vec<u8> = (0..size).map(|_| rng.gen_range(0..4)).collect();

        // Benchmark Encode (Write)
        let mut group_encode = c.benchmark_group(format!("Encode_{size_name}"));
        group_encode.throughput(Throughput::Bytes(size as u64));

        for (strategy_name, strategy) in &strategies {
            group_encode.bench_with_input(
                BenchmarkId::new(*strategy_name, size),
                &input,
                |b, i| b.iter(|| strategy.encode(black_box(i)).unwrap()),
            );
        }
        group_encode.finish();

        // Benchmark Decode (Read)
        // Pre-encode data for each strategy to benchmark decode
        let mut group_decode = c.benchmark_group(format!("Decode_{size_name}"));
        group_decode.throughput(Throughput::Bytes(size as u64));

        for (strategy_name, strategy) in &strategies {
            let encoded = strategy.encode(&input).unwrap();
            group_decode.bench_with_input(
                BenchmarkId::new(*strategy_name, size),
                &encoded,
                |b, e| b.iter(|| strategy.decode(black_box(e)).unwrap()),
            );
        }
        group_decode.finish();
    }
}

criterion_group!(benches, bench_codecs);
criterion_main!(benches);
