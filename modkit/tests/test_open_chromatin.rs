use crate::common::{check_against_expected_text_file, run_modkit};
use mod_kit::util::GenomeRegion;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

mod common;

#[test]
fn test_open_chromatin_help() {
    run_modkit(&["open-chromatin", "--help"]).expect("should run help");
    run_modkit(&["oc", "--help"]).expect("should run (abbrev) help");
}

struct BedGraphRecord {
    chrom: String,
    start: u32,
    end: u32,
    score: f32,
}

#[test]
fn test_open_chromatin_regression() {
    fn parse_bedgraph_line(l: &str) -> BedGraphRecord {
        let parts = l.split_ascii_whitespace().collect::<Vec<&str>>();
        let chrom = parts[0].to_string();
        let start = parts[1].parse::<u32>().unwrap();
        let end = parts[2].parse::<u32>().unwrap();
        let score = parts[3].parse::<f32>().unwrap();
        BedGraphRecord { chrom, start, end, score }
    }

    fn parse_bedgraph_file<T: AsRef<Path>>(fp: T) -> Vec<BedGraphRecord> {
        BufReader::new(File::open(fp).unwrap())
            .lines()
            .skip_while(|r| {
                r.as_ref().map(|l| l.starts_with('#')).unwrap_or(false)
            })
            .map(|res| res.unwrap())
            .map(|line| parse_bedgraph_line(&line))
            .collect::<Vec<BedGraphRecord>>()
    }
    let temp_file =
        std::env::temp_dir().join("test_open_chromatin_regression.bg");
    let args = [
        "oc",
        "predict",
        "../tests/resources/cs_test_hac52.bam",
        "-o",
        temp_file.to_str().unwrap(),
        "--super-batch-size",
        "1",
        "--batch-size",
        "16",
        "--force",
        "--model",
        "../ochm/models/r1041_e82_400bps_hac_v5.2.0@v0.1.0",
        "--region",
        "chr19:12,936,464-12,940,246",
    ];

    run_modkit(&args).unwrap();

    let expected_records = parse_bedgraph_file(
        "../tests/resources/test_open_chromatin_regression_expected.bg",
    );
    let observed_records = parse_bedgraph_file(&temp_file);
    assert_eq!(
        expected_records.len(),
        observed_records.len(),
        "number of records is not the same"
    );
    let diff = expected_records
        .iter()
        .zip(observed_records)
        .map(|(y, y_hat)| (y.score - y_hat.score).abs())
        .sum::<f32>();
    assert!(diff < 0.0001)
}
