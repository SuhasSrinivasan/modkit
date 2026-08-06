use crate::common::{
    check_against_expected_text_file, check_legal_csv, run_modkit,
};
use std::fs;
use std::path::Path;
use std::process::{Command, Output};

mod common;

#[test]
fn test_dmr_helps() {
    let _ = run_modkit(&["dmr", "pair", "--help"])
        .expect("failed to run modkit dmr pair help");
    let _ = run_modkit(&["dmr", "multi", "--help"])
        .expect("failed to run modkit dmr multi help");
}

#[test]
fn test_dmr_regression() {
    let out_bed = std::env::temp_dir().join("test_dmr_regression.bed");
    let _ = run_modkit(&[
        "dmr",
        "pair",
        "-a",
        "../tests/resources/\
         lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-b",
        "../tests/resources/\
         lung_00733-m_primary-tumour_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-o",
        out_bed.to_str().unwrap(),
        "-r",
        "../tests/resources/cpg_chr20_with_orig_names_selection.bed",
        "--ref",
        "../tests/resources/GRCh38_chr20.fa",
        "--header",
        "-f",
        "--base",
        "C",
    ])
    .expect("failed to run modkit dmr");

    check_legal_csv::<{ '\t' as u8 }>(&out_bed).expect("should be a legal CSV");
    check_against_expected_text_file(
        out_bed.to_str().unwrap(),
        "../tests/resources/test_output_chr20-2.bed",
    );

    let out_bed =
        std::env::temp_dir().join("foo").join("test_dmr_regression_2.bed");

    let _ = run_modkit(&[
        "dmr",
        "pair",
        "-a",
        "../tests/resources/\
         lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-b",
        "../tests/resources/\
         lung_00733-m_primary-tumour_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-o",
        out_bed.to_str().unwrap(),
        "-r",
        "../tests/resources/cpg_chr20_with_orig_names_selection.bed",
        "--ref",
        "../tests/resources/GRCh38_chr20.fa",
        "-f",
        "--header",
        "--base",
        "C",
    ])
    .expect("failed to run modkit dmr");

    check_legal_csv::<{ '\t' as u8 }>(&out_bed).expect("should be a legal CSV");
    check_against_expected_text_file(
        out_bed.to_str().unwrap(),
        "../tests/resources/test_output_chr20-2.bed",
    );
}

fn run_dmr_with_prior(output: &Path, alpha: &str, beta: &str) -> Output {
    Command::new(env!("CARGO_BIN_EXE_modkit"))
        .args([
            "dmr",
            "pair",
            "-a",
            "../tests/resources/lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz",
            "-b",
            "../tests/resources/lung_00733-m_primary-tumour_5mc-5hmc_chr20_cpg_pileup.bed.gz",
            "-o",
            output.to_str().unwrap(),
            "--ref",
            "../tests/resources/GRCh38_chr20.fa",
            "--base",
            "C",
            "--prior",
            alpha,
            beta,
            "--delta",
            "1",
            "--max-coverages",
            "100",
            "100",
            "--threads",
            "1",
            "--io-threads",
            "1",
            "--suppress-progress",
            "--force",
        ])
        .output()
        .unwrap()
}

#[test]
fn dmr_prior_cli_accepts_boundary_and_interior_but_rejects_invalid_inputs() {
    let temp_dir = tempfile::tempdir().unwrap();

    for (label, alpha, beta) in
        [("boundary", "0.5", "0.5"), ("interior", "0.55", "0.55")]
    {
        let output_path = temp_dir.path().join(format!("{label}.bed"));
        let output = run_dmr_with_prior(&output_path, alpha, beta);
        assert!(
            output.status.success(),
            "{label} prior ({alpha}, {beta}) was rejected: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    let below_boundary = run_dmr_with_prior(
        &temp_dir.path().join("below-boundary.bed"),
        "0.4",
        "0.5",
    );
    assert!(!below_boundary.status.success());
    assert!(String::from_utf8_lossy(&below_boundary.stderr)
        .contains("alpha + beta must be >= 1.0 for numerical stability"));

    let non_positive =
        run_dmr_with_prior(&temp_dir.path().join("non-positive.bed"), "0", "1");
    assert!(!non_positive.status.success());
    assert!(String::from_utf8_lossy(&non_positive.stderr)
        .contains("invalid beta parameters 0, 1"));
}

fn run_segmented_dmr(
    output: &Path,
    segments: &Path,
    threads: &str,
    io_threads: &str,
) {
    run_modkit(&[
        "dmr",
        "pair",
        "-a",
        "../tests/resources/\
         lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-b",
        "../tests/resources/\
         lung_00733-m_primary-tumour_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-o",
        output.to_str().unwrap(),
        "--segment",
        segments.to_str().unwrap(),
        "--ref",
        "../tests/resources/GRCh38_chr20.fa",
        "--header",
        "--base",
        "C",
        "--max-coverages",
        "100",
        "100",
        "--threads",
        threads,
        "--io-threads",
        io_threads,
        "--suppress-progress",
        "--force",
    ])
    .expect("segmented DMR run should succeed");
}

#[test]
fn segmentation_includes_last_site_and_is_thread_deterministic() {
    let temp_dir = tempfile::tempdir().unwrap();
    let sites_one = temp_dir.path().join("sites-threads-1.bed");
    let segments_one = temp_dir.path().join("segments-threads-1.bed");
    let sites_four = temp_dir.path().join("sites-threads-4.bed");
    let segments_four = temp_dir.path().join("segments-threads-4.bed");

    run_segmented_dmr(&sites_one, &segments_one, "1", "1");
    run_segmented_dmr(&sites_four, &segments_four, "4", "2");

    let sites_one = fs::read_to_string(sites_one).unwrap();
    let sites_four = fs::read_to_string(sites_four).unwrap();
    let segments_one = fs::read_to_string(segments_one).unwrap();
    let segments_four = fs::read_to_string(segments_four).unwrap();
    assert_eq!(sites_one.as_bytes(), sites_four.as_bytes());
    assert_eq!(segments_one.as_bytes(), segments_four.as_bytes());

    let site_rows = sites_one.lines().skip(1).collect::<Vec<_>>();
    let segment_rows = segments_one.lines().skip(1).collect::<Vec<_>>();
    assert_eq!(site_rows.len(), 17_271);

    let segment_site_total = segment_rows
        .iter()
        .map(|row| row.split('\t').nth(5).unwrap().parse::<usize>().unwrap())
        .sum::<usize>();
    assert_eq!(segment_site_total, site_rows.len());

    let last_site_end = site_rows.last().unwrap().split('\t').nth(2).unwrap();
    let last_segment_end =
        segment_rows.last().unwrap().split('\t').nth(2).unwrap();
    assert_eq!(last_site_end, "10804378");
    assert_eq!(last_segment_end, last_site_end);
}

// todo
//  test pair with explicit index
//  test multi
