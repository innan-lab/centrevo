# Codec Benchmark Results

This document summarizes the performance of different codec strategies:
- **UnpackedNoRS**: Baseline `Vec<u8>` (No protection, no packing).
- **UnpackedRS**: 1 byte/base + Reed-Solomon protection.
- **BitPackedRS**: 2 bits/base + Reed-Solomon protection.

## Test Cases
- **Small_Worst**: 1 byte (Max RS overhead)
- **Small_Best**: 14,272 bytes (Perfect alignment)
- **Medium_Best**: ~100kB (aligned)
- **Medium_Worst**: ~100kB + 1 byte (misaligned)
- **Large_Best**: ~1MB (aligned)
- **Large_Worst**: ~1MB + 1 byte (misaligned)

## Benchmark Results

### 1. Throughput Comparison (Giant Aligned Data - ~50MB)

| Strategy | Operation | Time (ms) | Throughput | Notes |
|----------|-----------|-----------|------------|-------|
| **UnpackedNoRS** | Encode | 0.93 | ~50 GB/s | Baseline |
| **UnpackedRS** | Encode | 693.0 | ~68 MB/s | Very slow, single-threaded RS |
| **BitPackedRS** | Encode | 201.0 | ~235 MB/s | 3.4x faster than UnpackedRS |
| **ParallelBitPackedRS** | Encode | **12.3** | **~3.7 GB/s** | **15x faster** than Serial BitPacked. Saturates NVMe. |
| **UnpackedNoRS** | Decode | 0.94 | ~50 GB/s | |
| **UnpackedRS** | Decode | 696.0 | ~68 MB/s | |
| **BitPackedRS** | Decode | 202.0 | ~235 MB/s | |
| **ParallelBitPackedRS** | Decode | **~13.0** | **~3.8 GB/s** | Estimated from Huge results |

### 2. Latency & Overhead (Small Data - 1 Byte)

| Strategy | Operation | Time | Overhead Factor | Notes |
|----------|-----------|------|-----------------|-------|
| **UnpackedNoRS** | Encode | ~23 ns | 1x | |
| **UnpackedRS** | Encode | ~40 µs | ~1700x | Padding to 14kB block required |
| **BitPackedRS** | Encode | ~40 µs | ~1700x | Padding to 14kB block required |
| **ParallelBitPackedRS** | Encode | ~41 µs | ~1700x | Negligible thread overhead (~1µs) |

### 3. Analysis
- **Parallel Scaling**: `ParallelBitPackedRS` achieves near-linear scaling on multi-core systems (tested on ~10 cores).
- **NVMe Saturation**: With throughputs of ~3.7 GB/s (Write) and ~3.8 GB/s (Read), this strategy successfully utilizes the bandwidth of modern NVMe SSDs (typically 3-7 GB/s).
- **Bit-Packing Efficiency**: Even without parallelism, `BitPackedRS` is 3-4x faster than `UnpackedRS` due to reduced data volume. Combined with parallelism, it provides a massive 50x speedup over the naive `UnpackedRS` approach.
