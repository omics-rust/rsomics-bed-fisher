use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_bed_fisher::fisher;
use std::fmt::Write as FmtWrite;
use std::io::Cursor;

fn bench_fisher(c: &mut Criterion) {
    // 500 A intervals on chr1, sorted by start
    let a_data: String = (0..500u64).fold(String::new(), |mut s, i| {
        let start = i * 2000;
        let _ = writeln!(s, "chr1\t{start}\t{}", start + 500);
        s
    });

    // 500 B intervals on chr1, offset so ~half overlap, sorted by start
    let b_data: String = (0..500u64).fold(String::new(), |mut s, i| {
        let start = i * 2000 + 300;
        let _ = writeln!(s, "chr1\t{start}\t{}", start + 500);
        s
    });

    let genome = "chr1\t1000000\n";

    c.bench_function("fisher_500x500", |bencher| {
        bencher.iter(|| {
            fisher(
                Cursor::new(a_data.as_str()),
                Cursor::new(b_data.as_str()),
                Cursor::new(genome),
                false,
            )
            .unwrap()
        });
    });
}

criterion_group!(benches, bench_fisher);
criterion_main!(benches);
