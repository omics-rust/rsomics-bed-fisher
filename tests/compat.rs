//! Compatibility tests: compare rsomics-bed-fisher output against bedtools fisher.
//!
//! Requires `bedtools` in PATH. Tests are skipped if bedtools is absent.

use std::process::Command;

fn bedtools_available() -> bool {
    Command::new("bedtools")
        .arg("--version")
        .output()
        .is_ok_and(|o| o.status.success())
}

fn run_our(a: &str, b: &str, genome: &str) -> String {
    let bin = env!("CARGO_BIN_EXE_rsomics-bed-fisher");
    let out = Command::new(bin)
        .args(["-a", a, "-b", b, "-g", genome])
        .output()
        .expect("failed to run rsomics-bed-fisher");
    assert!(
        out.status.success(),
        "rsomics-bed-fisher failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn run_bedtools(a: &str, b: &str, genome: &str) -> String {
    let out = Command::new("bedtools")
        .args(["fisher", "-a", a, "-b", b, "-g", genome])
        .output()
        .expect("failed to run bedtools fisher");
    assert!(
        out.status.success(),
        "bedtools fisher failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn extract_pvalue_line(output: &str) -> &str {
    output
        .lines()
        .find(|l| !l.starts_with('#') && !l.starts_with("left\t"))
        .unwrap_or("")
}

fn golden(name: &str) -> String {
    format!("{}/tests/golden/{}", env!("CARGO_MANIFEST_DIR"), name)
}

#[test]
fn header_line_counts_match() {
    if !bedtools_available() {
        return;
    }
    let a = golden("a.bed");
    let b = golden("b.bed");
    let g = golden("genome.txt");

    let ours = run_our(&a, &b, &g);
    let theirs = run_bedtools(&a, &b, &g);

    // Compare the "Number of overlaps" line
    let our_overlap = ours
        .lines()
        .find(|l| l.contains("Number of overlaps"))
        .unwrap_or("");
    let their_overlap = theirs
        .lines()
        .find(|l| l.contains("Number of overlaps"))
        .unwrap_or("");
    assert_eq!(our_overlap, their_overlap, "overlap count mismatch");
}

#[test]
fn n_possible_matches() {
    if !bedtools_available() {
        return;
    }
    let a = golden("a.bed");
    let b = golden("b.bed");
    let g = golden("genome.txt");

    let ours = run_our(&a, &b, &g);
    let theirs = run_bedtools(&a, &b, &g);

    let our_possible = ours
        .lines()
        .find(|l| l.contains("Number of possible"))
        .unwrap_or("");
    let their_possible = theirs
        .lines()
        .find(|l| l.contains("Number of possible"))
        .unwrap_or("");
    assert_eq!(our_possible, their_possible, "n_possible mismatch");
}

#[test]
fn p_values_match() {
    if !bedtools_available() {
        return;
    }
    let a = golden("a.bed");
    let b = golden("b.bed");
    let g = golden("genome.txt");

    let ours = run_our(&a, &b, &g);
    let theirs = run_bedtools(&a, &b, &g);

    // Parse the numeric output line (left/right/two-tail/ratio).
    let our_line = extract_pvalue_line(&ours);
    let their_line = extract_pvalue_line(&theirs);

    // Compare left, right, two-tail p-values (first 3 tab-separated fields)
    // with tolerance for floating-point differences.
    let parse_vals = |line: &str| -> Vec<f64> {
        line.split('\t')
            .take(3)
            .filter_map(|s| s.parse::<f64>().ok())
            .collect()
    };
    let our_vals = parse_vals(our_line);
    let their_vals = parse_vals(their_line);

    assert_eq!(
        our_vals.len(),
        3,
        "expected 3 p-values in our output, got: {our_line:?}"
    );
    assert_eq!(
        their_vals.len(),
        3,
        "expected 3 p-values in bedtools output, got: {their_line:?}"
    );

    for (i, (o, t)) in our_vals.iter().zip(their_vals.iter()).enumerate() {
        let diff = (o - t).abs();
        assert!(
            diff < 1e-4,
            "p-value[{i}] mismatch: ours={o:.6} bedtools={t:.6}"
        );
    }
}

#[test]
fn no_overlap_produces_p_right_near_one() {
    if !bedtools_available() {
        return;
    }
    let a_content = "chr1\t100\t200\n";
    let b_content = "chr1\t500\t600\n";
    let g_content = "chr1\t1000000\n";

    let tmp_a = "/tmp/fisher_compat_a_nooverlap.bed";
    let tmp_b = "/tmp/fisher_compat_b_nooverlap.bed";
    let tmp_g = "/tmp/fisher_compat_g.txt";

    std::fs::write(tmp_a, a_content).unwrap();
    std::fs::write(tmp_b, b_content).unwrap();
    std::fs::write(tmp_g, g_content).unwrap();

    let ours = run_our(tmp_a, tmp_b, tmp_g);
    let pvalue_line = extract_pvalue_line(&ours);
    let parts: Vec<f64> = pvalue_line
        .split('\t')
        .take(3)
        .filter_map(|s| s.parse().ok())
        .collect();
    assert_eq!(parts.len(), 3, "expected 3 p-values");
    // P(X >= 0) = 1 — right p-value is 1 when n11=0 (every outcome qualifies).
    assert!(
        (parts[1] - 1.0).abs() < 1e-9,
        "p_right={} (expected 1.0 when n11=0)",
        parts[1]
    );
}

#[test]
fn contingency_table_n11_equals_overlap_count() {
    if !bedtools_available() {
        return;
    }
    let a = golden("a.bed");
    let b = golden("b.bed");
    let g = golden("genome.txt");

    let ours = run_our(&a, &b, &g);

    let n_overlaps: i64 = ours
        .lines()
        .find(|l| l.contains("Number of overlaps"))
        .and_then(|l| l.rsplit(':').next())
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(-1);

    let n11: i64 = ours
        .lines()
        .find(|l| l.contains("in -a |") && !l.contains("not in"))
        .and_then(|l| {
            // "  #     in -a | 2            | 1            |"
            let after = l.split('|').nth(1)?;
            after.trim().parse().ok()
        })
        .unwrap_or(-2);

    assert_eq!(
        n_overlaps, n11,
        "n11 ({n11}) must equal overlap count ({n_overlaps})"
    );
}
