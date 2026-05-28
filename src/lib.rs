//! Fisher's exact test for interval overlap significance between two BED sets.
//!
//! ## Algorithm
//!
//! Counts, for each query interval in A, how many database intervals in B it
//! overlaps (counting each A-B pair once). Builds a 2×2 contingency table:
//!
//! ```text
//!             | in B  | not in B |
//!     in A    |  n11  |   n12    |
//! not in A    |  n21  |   n22    |
//! ```
//!
//! - `n11` = total overlap pairs (A intervals × B intervals that overlap them)
//! - `n12` = max(0, |A| − n11)
//! - `n21` = max(0, |B| − n11)
//! - `n22_full` (displayed as "Number of possible intervals") =
//!   max(n11+n12+n21, `genome_size` / `bMean`), where
//!   `bMean` = (1 + `mean_len_A`) + (1 + `mean_len_B`)
//! - `n22` = max(0, `n22_full` − n11 − n12 − n21)
//!
//! Applies Fisher's exact test (hypergeometric, three alternatives) to the
//! resulting table.
//!
//! Both files must be sorted by chromosome then start position.
//!
//! ## Reference
//!
//! `BEDTools` fisher — Quinlan & Hall (2010). `BEDTools`: a flexible suite of
//! utilities for comparing genomic features. Bioinformatics 26(6): 841–842.
//! DOI: 10.1093/bioinformatics/btq033

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read, Write};

use rsomics_common::{Result, RsomicsError};
use rsomics_stats::{Alternative, StatsError, fisher_exact_2x2};

/// Parsed BED3 record.
#[derive(Debug, Clone)]
struct Bed3 {
    chrom: String,
    start: i64,
    end: i64,
}

fn parse_bed3(line: &str) -> Option<Bed3> {
    let t = line.trim_end_matches(['\n', '\r']);
    if t.is_empty() || t.starts_with('#') || t.starts_with("track") || t.starts_with("browser") {
        return None;
    }
    let mut cols = t.splitn(4, '\t');
    let chrom = cols.next()?.to_owned();
    let start: i64 = cols.next()?.parse().ok()?;
    let end: i64 = cols.next()?.parse().ok()?;
    Some(Bed3 { chrom, start, end })
}

fn read_bed<R: Read>(r: R) -> Result<Vec<Bed3>> {
    let mut records = Vec::new();
    for line in BufReader::new(r).lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if let Some(rec) = parse_bed3(&line) {
            records.push(rec);
        }
    }
    Ok(records)
}

/// Parse a genome-sizes file (chrom TAB length per line).
fn read_genome<R: Read>(r: R) -> Result<i64> {
    let mut total: i64 = 0;
    for line in BufReader::new(r).lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let mut cols = t.splitn(2, '\t');
        cols.next(); // chrom name
        if let Some(len) = cols.next().and_then(|s| s.trim().parse::<i64>().ok()) {
            total += len;
        }
    }
    Ok(total)
}

/// Result of the Fisher test.
#[derive(Debug, Clone)]
pub struct FisherResult {
    pub query_count: u64,
    pub db_count: u64,
    pub overlap_count: u64,
    pub n_possible: i64,
    /// Contingency table: [[n11, n12], [n21, n22]]
    pub table: [[i64; 2]; 2],
    pub p_left: f64,
    pub p_right: f64,
    pub p_two_tail: f64,
    pub ratio: f64,
}

/// Compute Fisher's exact test for overlap between sorted BED files `a` and `b`
/// against the given genome size.
pub fn fisher<RA: Read, RB: Read, RG: Read>(
    a: RA,
    b: RB,
    genome: RG,
    merge_before: bool,
) -> Result<FisherResult> {
    let mut a_recs = read_bed(a)?;
    let mut b_recs = read_bed(b)?;
    let genome_size = read_genome(genome)?;

    if merge_before {
        a_recs = merge_intervals(a_recs);
        b_recs = merge_intervals(b_recs);
    }

    if a_recs.is_empty() || b_recs.is_empty() {
        #[allow(clippy::cast_possible_wrap)]
        let n12 = a_recs.len() as i64;
        #[allow(clippy::cast_possible_wrap)]
        let n21 = b_recs.len() as i64;
        let n11 = 0i64;
        let n22_full = n11 + n12 + n21;
        let n22 = 0i64;
        let (p_left, p_right, p_two_tail) = compute_fisher(n11, n12, n21, n22);
        let ratio = fisher_ratio(n11, n12, n21, n22);
        return Ok(FisherResult {
            #[allow(clippy::cast_possible_wrap)]
            query_count: a_recs.len() as u64,
            #[allow(clippy::cast_possible_wrap)]
            db_count: b_recs.len() as u64,
            overlap_count: 0,
            n_possible: n22_full,
            table: [[n11, n12], [n21, n22]],
            p_left,
            p_right,
            p_two_tail,
            ratio,
        });
    }

    let query_union: i64 = a_recs.iter().map(|r| r.end - r.start).sum();
    let db_union: i64 = b_recs.iter().map(|r| r.end - r.start).sum();
    #[allow(clippy::cast_possible_wrap)]
    let query_counts = a_recs.len() as i64;
    #[allow(clippy::cast_possible_wrap)]
    let db_counts = b_recs.len() as i64;

    // Group B by chromosome for fast lookup.
    let mut b_by_chrom: HashMap<&str, (usize, usize)> = HashMap::new();
    {
        let mut i = 0;
        while i < b_recs.len() {
            let chrom = b_recs[i].chrom.as_str();
            let start = i;
            while i < b_recs.len() && b_recs[i].chrom == chrom {
                i += 1;
            }
            b_by_chrom.insert(chrom, (start, i));
        }
    }

    // Count overlap pairs: for each A interval, count how many B intervals it overlaps.
    let mut overlap_counts: i64 = 0;
    for a_rec in &a_recs {
        let Some(&(b_start, b_end)) = b_by_chrom.get(a_rec.chrom.as_str()) else {
            continue;
        };
        let b_slice = &b_recs[b_start..b_end];

        // First B with start < a_rec.end may overlap.
        let upper = b_slice.partition_point(|b| b.start < a_rec.end);
        for b_rec in &b_slice[..upper] {
            if b_rec.end > a_rec.start {
                overlap_counts += 1;
            }
        }
    }

    // Contingency table following bedtools' formula.
    let n11 = overlap_counts;
    let n12 = (query_counts - overlap_counts).max(0);
    let n21 = (db_counts - overlap_counts).max(0);

    // bMean = (1 + mean_len_A) + (1 + mean_len_B)
    #[allow(clippy::cast_precision_loss)]
    let q_mean = 1.0 + query_union as f64 / query_counts as f64;
    #[allow(clippy::cast_precision_loss)]
    let d_mean = 1.0 + db_union as f64 / db_counts as f64;
    let b_mean = q_mean + d_mean;

    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let n22_full = (n11 + n12 + n21).max((genome_size as f64 / b_mean) as i64);
    let n22 = (n22_full - n12 - n21 - n11).max(0);

    let (p_left, p_right, p_two_tail) = compute_fisher(n11, n12, n21, n22);
    let ratio = fisher_ratio(n11, n12, n21, n22);

    Ok(FisherResult {
        #[allow(clippy::cast_sign_loss)]
        query_count: query_counts as u64,
        #[allow(clippy::cast_sign_loss)]
        db_count: db_counts as u64,
        #[allow(clippy::cast_sign_loss)]
        overlap_count: overlap_counts as u64,
        n_possible: n22_full,
        table: [[n11, n12], [n21, n22]],
        p_left,
        p_right,
        p_two_tail,
        ratio,
    })
}

fn p_or_fallback(
    result: std::result::Result<rsomics_stats::hypothesis::TestResult, StatsError>,
) -> f64 {
    result.map_or_else(
        |e| {
            if matches!(e, StatsError::Empty) {
                1.0
            } else {
                f64::NAN
            }
        },
        |r| r.p_value,
    )
}

fn compute_fisher(n11: i64, n12: i64, n21: i64, n22: i64) -> (f64, f64, f64) {
    let a = n11.max(0) as u64;
    let b = n12.max(0) as u64;
    let c = n21.max(0) as u64;
    let d = n22.max(0) as u64;
    let left = p_or_fallback(fisher_exact_2x2(a, b, c, d, Alternative::Less));
    let right = p_or_fallback(fisher_exact_2x2(a, b, c, d, Alternative::Greater));
    let two = p_or_fallback(fisher_exact_2x2(a, b, c, d, Alternative::TwoSided));
    (left, right, two)
}

fn fisher_ratio(n11: i64, n12: i64, n21: i64, n22: i64) -> f64 {
    if n12 == 0 || n21 == 0 {
        f64::INFINITY
    } else {
        #[allow(clippy::cast_precision_loss)]
        let r = (n11 as f64 * n22 as f64) / (n12 as f64 * n21 as f64);
        r
    }
}

/// Merge overlapping/adjacent intervals within a sorted BED record list.
fn merge_intervals(mut recs: Vec<Bed3>) -> Vec<Bed3> {
    if recs.is_empty() {
        return recs;
    }
    recs.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));
    let mut merged: Vec<Bed3> = Vec::with_capacity(recs.len());
    for rec in recs {
        if let Some(last) = merged.last_mut()
            && last.chrom == rec.chrom
            && last.end >= rec.start
        {
            last.end = last.end.max(rec.end);
            continue;
        }
        merged.push(rec);
    }
    merged
}

/// Format a float using C `%.5g` semantics (5 significant figures, scientific
/// notation when the exponent is < −4 or >= 5). The exponent always uses at
/// least 2 digits (e.g. `e-07` not `e-7`), matching C `printf` behaviour.
fn fmt_g5(v: f64) -> String {
    if v == 0.0 {
        return "0".to_owned();
    }
    if v == 1.0 {
        return "1".to_owned();
    }
    let exp = v.abs().log10().floor() as i32;
    if (-4..5).contains(&exp) {
        // Fixed notation with enough sig-figs.
        let decimals = (4 - exp).max(0) as usize;
        format!("{:.prec$}", v, prec = decimals)
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_owned()
    } else {
        // Scientific notation: 5 sig figs → 4 decimal places in mantissa.
        // Rust emits single-digit exponents (e-7); C uses at least 2 (e-07).
        let s = format!("{:.4e}", v);
        normalise_sci_exp(&s)
    }
}

/// Pad single-digit exponents to match C printf behaviour: `e-7` → `e-07`.
fn normalise_sci_exp(s: &str) -> String {
    if let Some(e_pos) = s.find('e') {
        let (mantissa, exp_part) = s.split_at(e_pos + 1); // includes 'e'
        let sign = if exp_part.starts_with('-') { "-" } else { "+" };
        let digits = exp_part.trim_start_matches(['-', '+']);
        if digits.len() == 1 {
            return format!("{mantissa}{sign}0{digits}");
        }
    }
    s.to_owned()
}

/// Write Fisher result to a writer, matching bedtools fisher output format.
pub fn write_fisher<W: Write>(result: &FisherResult, w: &mut W) -> Result<()> {
    let [[n11, n12], [n21, n22]] = result.table;
    writeln!(w, "# Number of query intervals: {}", result.query_count).map_err(RsomicsError::Io)?;
    writeln!(w, "# Number of db intervals: {}", result.db_count).map_err(RsomicsError::Io)?;
    writeln!(w, "# Number of overlaps: {}", result.overlap_count).map_err(RsomicsError::Io)?;
    writeln!(
        w,
        "# Number of possible intervals (estimated): {}",
        result.n_possible,
    )
    .map_err(RsomicsError::Io)?;
    writeln!(
        w,
        "# phyper({} - 1, {}, {} - {}, {}, lower.tail=F)",
        n11, result.query_count, result.n_possible, result.query_count, result.db_count
    )
    .map_err(RsomicsError::Io)?;
    writeln!(w, "# Contingency Table Of Counts").map_err(RsomicsError::Io)?;
    writeln!(w, "#_________________________________________").map_err(RsomicsError::Io)?;
    writeln!(w, "#           | {:<12} | {:<12} |", " in -b", "not in -b")
        .map_err(RsomicsError::Io)?;
    writeln!(w, "#     in -a | {:<12} | {:<12} |", n11, n12).map_err(RsomicsError::Io)?;
    writeln!(w, "# not in -a | {:<12} | {:<12} |", n21, n22).map_err(RsomicsError::Io)?;
    writeln!(w, "#_________________________________________").map_err(RsomicsError::Io)?;
    writeln!(w, "# p-values for fisher's exact test").map_err(RsomicsError::Io)?;
    writeln!(w, "left\tright\ttwo-tail\tratio").map_err(RsomicsError::Io)?;
    let pl = fmt_g5(result.p_left);
    let pr = fmt_g5(result.p_right);
    let pt = fmt_g5(result.p_two_tail);
    if result.ratio.is_nan() {
        writeln!(w, "{pl}\t{pr}\t{pt}\t-nan").map_err(RsomicsError::Io)?;
    } else if result.ratio.is_infinite() {
        writeln!(w, "{pl}\t{pr}\t{pt}\tinf").map_err(RsomicsError::Io)?;
    } else {
        writeln!(w, "{pl}\t{pr}\t{pt}\t{:.3}", result.ratio).map_err(RsomicsError::Io)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn genome_1m() -> Cursor<&'static str> {
        Cursor::new("chr1\t1000000\n")
    }

    fn run(a: &str, b: &str) -> FisherResult {
        fisher(Cursor::new(a), Cursor::new(b), genome_1m(), false).unwrap()
    }

    #[test]
    fn no_overlaps_p_values_in_range() {
        // With 0 observed overlaps, p_right = P(X >= 0) = 1 (trivially satisfied).
        let r = run("chr1\t100\t200\n", "chr1\t300\t400\n");
        assert_eq!(r.overlap_count, 0);
        assert_eq!(r.p_right, 1.0, "p_right must be 1 when n11=0 (P(X>=0)=1)");
        assert!((0.0..=1.0).contains(&r.p_left), "p_left={}", r.p_left);
    }

    #[test]
    fn full_overlap_p_right_is_low() {
        let r = run(
            "chr1\t100\t200\nchr1\t400\t500\n",
            "chr1\t150\t250\nchr1\t450\t550\n",
        );
        assert_eq!(r.overlap_count, 2);
        assert!(r.p_right < 0.01, "p_right={}", r.p_right);
    }

    #[test]
    fn different_chroms_no_overlap() {
        let r = run("chr1\t100\t200\n", "chr2\t100\t200\n");
        assert_eq!(r.overlap_count, 0);
    }

    #[test]
    fn empty_a_returns_zero_overlaps() {
        let r = run("", "chr1\t100\t200\n");
        assert_eq!(r.overlap_count, 0);
    }

    #[test]
    fn contingency_table_matches_bedtools() {
        let a = "chr1\t100\t200\nchr1\t400\t500\nchr1\t700\t800\n";
        let b = "chr1\t150\t250\n";
        let r = run(a, b);
        assert_eq!(r.table[0][0], 1, "n11"); // 1 overlap
        assert_eq!(r.table[0][1], 2, "n12"); // 2 A intervals not overlapping
        assert_eq!(r.table[1][0], 0, "n21"); // 0 B intervals not overlapping
    }

    #[test]
    fn ratio_is_inf_when_n12_zero() {
        let r = run("chr1\t100\t200\n", "chr1\t150\t250\n");
        assert!(r.ratio.is_infinite(), "ratio={}", r.ratio);
    }

    #[test]
    fn p_values_in_range() {
        let a = "chr1\t100\t200\nchr1\t300\t400\nchr1\t600\t700\n";
        let b = "chr1\t150\t250\nchr1\t350\t450\n";
        let r = run(a, b);
        assert!((0.0..=1.0).contains(&r.p_left), "p_left={}", r.p_left);
        assert!((0.0..=1.0).contains(&r.p_right), "p_right={}", r.p_right);
        assert!(
            (0.0..=1.0).contains(&r.p_two_tail),
            "p_two={}",
            r.p_two_tail
        );
    }
}
