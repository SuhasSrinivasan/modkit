use crate::common::{check_against_expected_text_file, run_modkit};
use rust_htslib::bam::{self, Read};
use std::fs::{self, File};
use std::io::{Read as IoRead, Seek, SeekFrom};
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

mod common;

#[test]
fn test_pileup_hemi_help() {
    let pileup_help_args = ["pileup-hemi", "--help"];
    let _out = run_modkit(&pileup_help_args).unwrap();
}

fn make_truncated_duplex_bam() -> (TempDir, PathBuf) {
    const SOURCE: &str = "../tests/resources/duplex_modcalls_sort.bam";
    let mut reader = bam::IndexedReader::from_path(SOURCE).unwrap();
    let tid = reader.header().tid(b"chr20").unwrap();
    reader.fetch((tid, 22_613_835, 22_640_468)).unwrap();
    let mut first_record = bam::Record::new();
    reader
        .read(&mut first_record)
        .expect("region should contain a record")
        .unwrap();
    let first_data_block = (reader.tell() as u64) >> 16;
    assert!(first_data_block > 0);

    let temp_dir = tempfile::tempdir().unwrap();
    let bam_fp = temp_dir.path().join("truncated.bam");
    fs::copy(SOURCE, &bam_fp).unwrap();
    fs::copy(format!("{SOURCE}.bai"), bam_fp.with_extension("bam.bai"))
        .unwrap();

    let mut bam = File::options().read(true).write(true).open(&bam_fp).unwrap();
    let len = bam.metadata().unwrap().len();
    bam.seek(SeekFrom::Start(first_data_block + 16)).unwrap();
    let mut block_size = [0u8; 2];
    IoRead::read_exact(&mut bam, &mut block_size).unwrap();
    let block_size = u16::from_le_bytes(block_size) as u64 + 1;
    let truncate_at = first_data_block + (block_size / 2);
    assert!(truncate_at < len);
    bam.set_len(truncate_at).unwrap();
    (temp_dir, bam_fp)
}

#[test]
fn test_pileup_hemi_rejects_truncated_bam() {
    let (_temp_dir, bam_fp) = make_truncated_duplex_bam();
    let out_fp = bam_fp.with_extension("bed");
    let output = Command::new(env!("CARGO_BIN_EXE_modkit"))
        .arg("pileup-hemi")
        .arg(&bam_fp)
        .args(["-o", out_fp.to_str().unwrap()])
        .args(["-r", "../tests/resources/GRCh38_chr20.fa"])
        .args(["--cpg", "--no-filtering", "--threads", "1"])
        .args(["--chunk-size", "1"])
        .args(["--region", "chr20:22,613,835-22,640,468"])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(
        !output.status.success(),
        "truncated BAM unexpectedly succeeded; stderr: {stderr}"
    );
    assert!(
        stderr.contains("failed to read BAM pileup during hemi pileup"),
        "missing BAM pileup context in stderr: {stderr}"
    );
}

#[test]
fn test_pileup_hemi_hm() {
    let temp_file = std::env::temp_dir().join("test_pileup_hemi_hm.bed");
    let args = [
        "pileup-hemi",
        "../tests/resources/duplex_modcalls_sort.bam",
        "-o",
        temp_file.to_str().unwrap(),
        "-r",
        "../tests/resources/GRCh38_chr20.fa",
        "--motif",
        "CG",
        "0",
        "--region",
        "chr20:22,613,835-22,640,468",
        "--no-filtering",
        "--mixed-delim",
    ];

    run_modkit(&args).unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/duplex_hemi_nofilt.bed",
    );
}

#[test]
fn test_pileup_hemi_preset() {
    let temp_file = std::env::temp_dir().join("test_pileup_hemi_preset.bed");
    let args = [
        "pileup-hemi",
        "../tests/resources/duplex_modcalls_sort.bam",
        "-o",
        temp_file.to_str().unwrap(),
        "-r",
        "../tests/resources/GRCh38_chr20.fa",
        "--cpg",
        "--region",
        "chr20:22,613,835-22,640,468",
        "--mixed-delim",
    ];

    run_modkit(&args).unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/duplex_hemi.bed",
    );
}

// todo test with combine mods
