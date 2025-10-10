use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{anyhow, Context, Result};
use assert_cmd::prelude::*;
use csv::ReaderBuilder;
use parquet::file::reader::{FileReader, SerializedFileReader};
use parquet::record::RowAccessor;
use serde::Deserialize;
use tempfile::TempDir;

#[derive(Debug, Deserialize)]
struct GoldenObservation {
    scannr: u32,
    charge: i32,
    #[allow(dead_code)]
    mass: f64,
    intensity: f64,
    #[allow(dead_code)]
    isotope: i32,
    #[allow(dead_code)]
    rt: f64,
}

struct LfqAggregate {
    run_names: Vec<String>,
    intensities: HashMap<(usize, i32), f64>,
    presence: HashSet<(usize, i32)>,
}

#[derive(Debug, Deserialize)]
struct GoldenData {
    run_names: Vec<String>,
    observations: Vec<GoldenObservation>,
}

fn read_lfq_tsv(lfq_path: &Path) -> Result<LfqAggregate> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(lfq_path)
        .with_context(|| format!("failed to read {}", lfq_path.display()))?;

    let headers = reader
        .headers()
        .with_context(|| format!("failed to read headers from {}", lfq_path.display()))?
        .clone();

    let run_names: Vec<String> = headers.iter().skip(6).map(|s| s.to_string()).collect();

    let mut intensities: HashMap<(usize, i32), f64> = HashMap::new();
    let mut presence: HashSet<(usize, i32)> = HashSet::new();

    for record in reader.records() {
        let record =
            record.with_context(|| format!("failed to read row from {}", lfq_path.display()))?;

        let charge: i32 = record
            .get(1)
            .ok_or_else(|| anyhow!("missing charge column in {}", lfq_path.display()))?
            .parse()
            .with_context(|| "failed to parse charge value")?;

        for (idx, _) in run_names.iter().enumerate() {
            let value = record
                .get(6 + idx)
                .ok_or_else(|| anyhow!("missing intensity column in {}", lfq_path.display()))?;
            let intensity = if value.trim().is_empty() {
                0.0
            } else {
                value.parse::<f64>().with_context(|| {
                    format!(
                        "failed to parse intensity for run {} (charge {})",
                        idx, charge
                    )
                })?
            };
            if intensity > 0.0 {
                presence.insert((idx, charge));
            }
            *intensities.entry((idx, charge)).or_insert(0.0) += intensity;
        }
    }

    Ok(LfqAggregate {
        run_names,
        intensities,
        presence,
    })
}

fn read_lfq_parquet(
    parquet_path: &Path,
    run_names: &[String],
) -> Result<HashMap<(usize, i32), f64>> {
    let run_lookup: HashMap<String, usize> = run_names
        .iter()
        .enumerate()
        .map(|(idx, name)| (name.clone(), idx))
        .collect();

    let file = fs::File::open(parquet_path)
        .with_context(|| format!("failed to open {}", parquet_path.display()))?;
    let reader = SerializedFileReader::new(file)
        .with_context(|| format!("failed to read parquet from {}", parquet_path.display()))?;
    let mut rows = reader
        .get_row_iter(None)
        .with_context(|| format!("failed to iterate rows in {}", parquet_path.display()))?;

    let mut intensities: HashMap<(usize, i32), f64> = HashMap::new();

    while let Some(row) = rows.next() {
        let row = row.with_context(|| {
            format!("failed to read parquet row from {}", parquet_path.display())
        })?;

        let filename = row
            .get_string(6)
            .with_context(|| format!("missing filename column in {}", parquet_path.display()))?
            .to_string();

        let charge = match row.get_int(2) {
            Ok(value) => value,
            Err(_) => continue,
        };

        let intensity = row
            .get_double(7)
            .or_else(|_| row.get_float(7).map(f64::from))
            .with_context(|| format!("missing intensity column in {}", parquet_path.display()))?;

        let run_idx = run_lookup
            .get(&filename)
            .ok_or_else(|| anyhow!("unexpected run '{filename}' in parquet output"))?;

        if intensity > 0.0 {
            *intensities.entry((*run_idx, charge)).or_insert(0.0) += intensity;
        }
    }

    Ok(intensities)
}

fn approx_equal(actual: f64, expected: f64) -> bool {
    let scale = actual.abs().max(expected.abs()).max(1.0);
    (actual - expected).abs() <= scale * 1.0e-5
}

#[test]
fn lfq_outputs_match_golden() -> Result<()> {
    let crate_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let repo_root = crate_dir
        .parent()
        .and_then(Path::parent)
        .ok_or_else(|| anyhow!("unable to resolve repository root"))?;
    let fixtures = repo_root.join("tests");

    let config_template = fixtures.join("config.json");
    let config_contents = fs::read_to_string(&config_template)
        .with_context(|| format!("failed to read {}", config_template.display()))?;

    let mut config: serde_json::Value = serde_json::from_str(&config_contents)
        .with_context(|| format!("failed to parse {}", config_template.display()))?;

    let mzml_paths = [
        fixtures.join("lfq_run1.mzML"),
        fixtures.join("lfq_run2.mzML"),
    ];
    let fasta_path = fixtures.join("Q99536.fasta");

    config["database"]["fasta"] = serde_json::Value::String(fasta_path.to_string_lossy().into());
    config["mzml_paths"] = serde_json::Value::Array(
        mzml_paths
            .iter()
            .map(|p| serde_json::Value::String(p.to_string_lossy().into()))
            .collect(),
    );
    config["quant"] = serde_json::json!({
        "lfq": true,
        "lfq_settings": {
            "combine_charge_states": false,
            "min_peptide_q": 1.0,
            "max_precursor_q": 1.0
        }
    });

    let temp_dir = TempDir::new().context("failed to create temporary directory")?;
    let config_path = temp_dir.path().join("config.json");
    fs::write(&config_path, serde_json::to_vec_pretty(&config)?).with_context(|| {
        format!(
            "failed to write temporary configuration to {}",
            config_path.display()
        )
    })?;

    let tsv_dir = temp_dir.path().join("tsv");
    fs::create_dir_all(&tsv_dir)
        .with_context(|| format!("failed to create output directory {}", tsv_dir.display()))?;

    let mut cmd = Command::cargo_bin("sage")?;
    cmd.current_dir(&repo_root)
        .arg(&config_path)
        .arg("--output_directory")
        .arg(&tsv_dir)
        .arg("--disable-telemetry-i-dont-want-to-improve-sage");
    cmd.assert().success();

    let lfq_path = tsv_dir.join("lfq.tsv");
    let lfq = read_lfq_tsv(&lfq_path)?;

    let golden_path = fixtures.join("lfq_test.json");
    let golden_data: GoldenData = serde_json::from_reader(
        fs::File::open(&golden_path)
            .with_context(|| format!("failed to open {}", golden_path.display()))?,
    )
    .with_context(|| format!("failed to parse {}", golden_path.display()))?;

    assert_eq!(
        lfq.run_names, golden_data.run_names,
        "LFQ TSV run headers did not match the golden run names"
    );

    if golden_data.observations.is_empty() {
        assert!(
            lfq.presence.is_empty(),
            "LFQ TSV unexpectedly contained quantified entries"
        );
        assert!(
            lfq.intensities.values().all(|v| approx_equal(*v, 0.0)),
            "LFQ TSV contained non-zero intensities without matching golden data"
        );
    } else {
        let mut unique_scans: Vec<u32> =
            golden_data.observations.iter().map(|g| g.scannr).collect();
        unique_scans.sort_unstable();
        unique_scans.dedup();

        let scan_to_index: HashMap<u32, usize> = unique_scans
            .iter()
            .enumerate()
            .map(|(idx, scan)| (*scan, idx))
            .collect();

        let mut expected_intensity: HashMap<(usize, i32), f64> = HashMap::new();
        let mut expected_presence: HashSet<(usize, i32)> = HashSet::new();

        for entry in &golden_data.observations {
            let run_idx = scan_to_index[&entry.scannr];
            *expected_intensity
                .entry((run_idx, entry.charge))
                .or_insert(0.0) += entry.intensity;
            expected_presence.insert((run_idx, entry.charge));
        }

        assert_eq!(
            lfq.presence, expected_presence,
            "LFQ TSV peptides did not match the golden peptide identifiers"
        );

        for (key, expected) in &expected_intensity {
            let actual = lfq.intensities.get(key).copied().unwrap_or_default();
            assert!(
                approx_equal(actual, *expected),
                "intensity mismatch for run {} charge {}: expected {}, observed {}",
                key.0,
                key.1,
                expected,
                actual
            );
        }
    }

    let parquet_dir = temp_dir.path().join("parquet");
    fs::create_dir_all(&parquet_dir).with_context(|| {
        format!(
            "failed to create output directory {}",
            parquet_dir.display()
        )
    })?;

    let mut parquet_cmd = Command::cargo_bin("sage")?;
    parquet_cmd
        .current_dir(&repo_root)
        .arg(&config_path)
        .arg("--output_directory")
        .arg(&parquet_dir)
        .arg("--parquet")
        .arg("--disable-telemetry-i-dont-want-to-improve-sage");
    parquet_cmd.assert().success();

    let parquet_path = parquet_dir.join("lfq.parquet");
    let parquet_map = read_lfq_parquet(&parquet_path, &lfq.run_names)?;

    for key in &lfq.presence {
        let tsv_intensity = lfq.intensities.get(key).copied().unwrap_or_default();
        let parquet_intensity = parquet_map.get(key).copied().unwrap_or_default();
        assert!(
            approx_equal(tsv_intensity, parquet_intensity),
            "mismatch between TSV and parquet for run {} charge {}: {} vs {}",
            key.0,
            key.1,
            tsv_intensity,
            parquet_intensity
        );
    }

    Ok(())
}
