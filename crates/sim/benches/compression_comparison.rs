use centrevo_codec::CodecStrategy;
use criterion::{Criterion, Throughput, black_box, criterion_group, criterion_main};
use rand::Rng;

fn generate_dna(len: usize) -> Vec<u8> {
    let mut rng = rand::rng();
    (0..len).map(|_| rng.random_range(0..4)).collect()
}

fn generate_repetitive_dna(len: usize) -> Vec<u8> {
    let mut rng = rand::rng();
    // Simulate 171bp monomer repeats with slight variation
    let monomer: Vec<u8> = (0..171).map(|_| rng.random_range(0..4)).collect();
    let mut seq = Vec::with_capacity(len);
    while seq.len() < len {
        // Copy monomer with 1% mutation rate to be realistic
        let mut variance = monomer.clone();
        for base in variance.iter_mut() {
            if rng.random_bool(0.01) {
                *base = rng.random_range(0..4);
            }
        }
        seq.extend(variance);
    }
    seq.truncate(len);
    seq
}

fn bench_compression(c: &mut Criterion) {
    let mut group = c.benchmark_group("compression_comparison");
    let len = 1_000_000; // 1MB sequence
    let random_data = generate_dna(len);
    let repetitive_data = generate_repetitive_dna(len);

    group.throughput(Throughput::Bytes(len as u64));

    // Print sizes for comparison
    let packed_random = CodecStrategy::BitPackedRS.encode(&random_data).unwrap();
    let unpacked_z_random = CodecStrategy::UnpackedZ.encode(&random_data).unwrap();
    let unpacked_rsz_random = CodecStrategy::UnpackedZRS.encode(&random_data).unwrap();

    println!("Original: {len} bytes");
    println!(
        "BitPackedRS: {} bytes ({:.2}%)",
        packed_random.len(),
        (packed_random.len() as f64 / len as f64) * 100.0
    );
    println!(
        "UnpackedZ:   {} bytes ({:.2}%)",
        unpacked_z_random.len(),
        (unpacked_z_random.len() as f64 / len as f64) * 100.0
    );
    println!(
        "UnpackedZRS: {} bytes ({:.2}%)",
        unpacked_rsz_random.len(),
        (unpacked_rsz_random.len() as f64 / len as f64) * 100.0
    );

    // Benchmarks for Random Data
    group.bench_function("bitpacked_rs_encode", |b| {
        b.iter(|| {
            CodecStrategy::BitPackedRS
                .encode(black_box(&random_data))
                .unwrap()
        })
    });
    group.bench_function("unpacked_z_encode", |b| {
        b.iter(|| {
            CodecStrategy::UnpackedZ
                .encode(black_box(&random_data))
                .unwrap()
        })
    });
    group.bench_function("unpacked_rsz_encode", |b| {
        b.iter(|| {
            CodecStrategy::UnpackedZRS
                .encode(black_box(&random_data))
                .unwrap()
        })
    });

    // Decoding benchmarks
    group.bench_function("bitpacked_rs_decode", |b| {
        b.iter(|| {
            CodecStrategy::BitPackedRS
                .decode(black_box(&packed_random))
                .unwrap()
        })
    });
    group.bench_function("unpacked_z_decode", |b| {
        b.iter(|| {
            CodecStrategy::UnpackedZ
                .decode(black_box(&unpacked_z_random))
                .unwrap()
        })
    });

    // --- REPETITIVE DATA ---
    println!("\n--- REPETITIVE DATA (1MB) ---");
    let packed_rep = CodecStrategy::BitPackedRS.encode(&repetitive_data).unwrap();
    let unpacked_z_rep = CodecStrategy::UnpackedZ.encode(&repetitive_data).unwrap();
    let unpacked_rsz_rep = CodecStrategy::UnpackedZRS.encode(&repetitive_data).unwrap();

    println!("Original: {len} bytes");
    println!(
        "BitPackedRS: {} bytes ({:.2}%)",
        packed_rep.len(),
        (packed_rep.len() as f64 / len as f64) * 100.0
    );
    println!(
        "UnpackedZ:   {} bytes ({:.2}%)",
        unpacked_z_rep.len(),
        (unpacked_z_rep.len() as f64 / len as f64) * 100.0
    );
    println!(
        "UnpackedZRS: {} bytes ({:.2}%)",
        unpacked_rsz_rep.len(),
        (unpacked_rsz_rep.len() as f64 / len as f64) * 100.0
    );

    // 1. BitPackedRS Encode
    group.bench_function("rep_bitpacked_rs_encode", |b| {
        b.iter(|| {
            CodecStrategy::BitPackedRS
                .encode(black_box(&repetitive_data))
                .unwrap()
        })
    });

    // 2. UnpackedZ Encode
    group.bench_function("rep_unpacked_z_encode", |b| {
        b.iter(|| {
            CodecStrategy::UnpackedZ
                .encode(black_box(&repetitive_data))
                .unwrap()
        })
    });

    // 3. UnpackedZRS Encode
    group.bench_function("rep_unpacked_rsz_encode", |b| {
        b.iter(|| {
            CodecStrategy::UnpackedZRS
                .encode(black_box(&repetitive_data))
                .unwrap()
        })
    });

    // Decoding benchmarks

    // 1. BitPackedRS Decode
    group.bench_function("rep_bitpacked_rs_decode", |b| {
        b.iter(|| {
            CodecStrategy::BitPackedRS
                .decode(black_box(&packed_rep))
                .unwrap()
        })
    });

    // 2. UnpackedZ Decode
    group.bench_function("rep_unpacked_z_decode", |b| {
        b.iter(|| {
            CodecStrategy::UnpackedZ
                .decode(black_box(&unpacked_z_rep))
                .unwrap()
        })
    });

    group.finish();
}

criterion_group!(benches, bench_compression);
criterion_main!(benches);
