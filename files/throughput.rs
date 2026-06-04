// =============================================================================
//  HelixArchive — Throughput benchmarks
//  benches/throughput.rs
// =============================================================================

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use helix_archive::*;
use std::io::Cursor;

fn bench_encode_bytes(c: &mut Criterion) {
    let mut group = c.benchmark_group("encode_bytes");
    for size in [64usize, 1024, 65_536, 1_048_576] {
        let data: Vec<u8> = (0..size).map(|i| (i * 7 + 13) as u8).collect();
        group.throughput(Throughput::Bytes(size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), &data, |b, d| {
            b.iter(|| encode_bytes(black_box(d)))
        });
    }
    group.finish();
}

fn bench_decode_bytes(c: &mut Criterion) {
    let mut group = c.benchmark_group("decode_bytes");
    for size in [64usize, 1024, 65_536, 1_048_576] {
        let data: Vec<u8> = (0..size).map(|i| (i * 7 + 13) as u8).collect();
        let encoded = encode_bytes(&data);
        group.throughput(Throughput::Bytes(size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), &encoded, |b, enc| {
            b.iter(|| decode_bytes(black_box(enc)))
        });
    }
    group.finish();
}

fn bench_encode_stream(c: &mut Criterion) {
    let mut group = c.benchmark_group("encode_stream");
    for size in [65_536usize, 1_048_576] {
        let data: Vec<u8> = (0..size).map(|i| (i * 7 + 13) as u8).collect();
        group.throughput(Throughput::Bytes(size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), &data, |b, d| {
            b.iter(|| {
                let mut out = Vec::with_capacity(d.len() * 4 + 32);
                encode_stream(Cursor::new(black_box(d)), Cursor::new(&mut out), CHUNK_SIZE).unwrap();
                out
            })
        });
    }
    group.finish();
}

fn bench_decode_stream(c: &mut Criterion) {
    let mut group = c.benchmark_group("decode_stream");
    for size in [65_536usize, 1_048_576] {
        let data: Vec<u8> = (0..size).map(|i| (i * 7 + 13) as u8).collect();
        let mut enc = Vec::new();
        encode_stream(Cursor::new(&data), Cursor::new(&mut enc), CHUNK_SIZE).unwrap();
        group.throughput(Throughput::Bytes(size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), &enc, |b, e| {
            b.iter(|| {
                let mut out = Vec::with_capacity(size);
                decode_stream(Cursor::new(black_box(e)), Cursor::new(&mut out), CHUNK_SIZE).unwrap();
                out
            })
        });
    }
    group.finish();
}

fn bench_gc_content(c: &mut Criterion) {
    let strand: Vec<u8> = encode_bytes(&(0u8..=255).cycle().take(65_536).collect::<Vec<_>>());
    c.bench_function("gc_content_64k", |b| {
        b.iter(|| gc_content(black_box(&strand)))
    });
}

fn bench_homopolymer(c: &mut Criterion) {
    let strand: Vec<u8> = encode_bytes(&(0u8..=255).cycle().take(65_536).collect::<Vec<_>>());
    c.bench_function("homopolymer_64k", |b| {
        b.iter(|| find_homopolymer_runs(black_box(&strand), 5))
    });
}

criterion_group!(
    benches,
    bench_encode_bytes,
    bench_decode_bytes,
    bench_encode_stream,
    bench_decode_stream,
    bench_gc_content,
    bench_homopolymer,
);
criterion_main!(benches);
