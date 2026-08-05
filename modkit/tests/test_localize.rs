use crate::common::run_modkit;
use std::fs::{read_to_string, write};

mod common;

#[test]
fn test_localise_helps() {
    let _ = run_modkit(&["localize", "--help"])
        .expect("failed to run modkit localize help");
    let _ = run_modkit(&["localise", "--help"])
        .expect("failed to run modkit localise help");
}

#[test]
fn test_localize_empty_chart_does_not_create_or_truncate_file() {
    const EMPTY_TABLE: &str =
        "mod_code\toffset\tn_valid\tn_mod\tpercent_modified\n";

    let temp_dir = tempfile::tempdir().unwrap();
    let regions = temp_dir.path().join("regions.bed");
    let genome_sizes = temp_dir.path().join("genome-sizes.tsv");
    let table = temp_dir.path().join("localize.tsv");
    let new_chart = temp_dir.path().join("new-localize.html");
    let existing_chart = temp_dir.path().join("existing-localize.html");

    write(&regions, "chr20\t1\t2\n").unwrap();
    write(&genome_sizes, "chr20\t100000000\n").unwrap();

    let run = |chart: &std::path::Path, force: bool| {
        let mut args = vec![
            "localize",
            "../tests/resources/lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz",
            "--regions",
            regions.to_str().unwrap(),
            "--genome-sizes",
            genome_sizes.to_str().unwrap(),
            "--window",
            "1",
            "--out-file",
            table.to_str().unwrap(),
            "--chart",
            chart.to_str().unwrap(),
        ];
        if force {
            args.push("--force");
        }
        run_modkit(&args)
    };

    let create_result = run(&new_chart, false);
    assert!(create_result.is_err());
    assert_eq!(read_to_string(&table).unwrap(), EMPTY_TABLE);

    let sentinel = "keep existing chart content\n";
    write(&existing_chart, sentinel).unwrap();
    let truncate_result = run(&existing_chart, true);

    assert!(truncate_result.is_err());
    assert_eq!(read_to_string(&table).unwrap(), EMPTY_TABLE);
    let existing_content = read_to_string(existing_chart).unwrap();
    assert!(
        !new_chart.exists() && existing_content == sentinel,
        "failed chart rendering mutated chart paths: new_exists={}, \
         existing_content={existing_content:?}",
        new_chart.exists()
    );
}
