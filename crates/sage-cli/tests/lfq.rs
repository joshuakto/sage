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
    run_name: String,
    #[allow(dead_code)]
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

#[derive(Debug, Deserialize)]
struct GoldenProteinData {
    run_names: Vec<String>,
    proteins: Vec<GoldenProtein>,
}

#[derive(Debug, Deserialize)]
struct GoldenProtein {
    accessions: Vec<String>,
    min_q_value: f32,
    has_intensities: bool,
    min_coverage_count: usize,
}

#[derive(Debug, Clone)]
struct ProteinRecord {
    accessions: Vec<String>,
    q_value: f32,
    total_precursors: usize,
    unique_peptides: usize,
    lfq_peptide_count: usize,
    intensities: Vec<f64>,
    coverage: Vec<bool>,
    contributors: Vec<usize>,
}

struct ProteinLfqData {
    run_names: Vec<String>,
    proteins: Vec<ProteinRecord>,
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

fn map_dynamic_columns_to_runs(
    headers: &csv::StringRecord,
    run_names: &[String],
    prefix: &str,
    path: &Path,
) -> Result<Vec<usize>> {
    let run_lookup: HashMap<&str, usize> = run_names
        .iter()
        .enumerate()
        .map(|(idx, name)| (name.as_str(), idx))
        .collect();

    let mut ordered: Vec<Option<usize>> = vec![None; run_names.len()];
    for (column_idx, header) in headers.iter().enumerate() {
        if let Some(name) = header.strip_prefix(prefix) {
            let run_idx = *run_lookup.get(name).ok_or_else(|| {
                anyhow!(
                    "unexpected {} run '{}' in {}",
                    prefix.trim_end_matches(':'),
                    name,
                    path.display()
                )
            })?;
            ordered[run_idx] = Some(column_idx);
        }
    }

    ordered
        .into_iter()
        .enumerate()
        .map(|(run_idx, column)| {
            column.ok_or_else(|| {
                anyhow!(
                    "missing {} column for run '{}' in {}",
                    prefix.trim_end_matches(':'),
                    run_names[run_idx],
                    path.display()
                )
            })
        })
        .collect()
}

fn sort_proteins_by_accession(records: &mut [ProteinRecord]) {
    records.sort_by(|a, b| a.accessions.cmp(&b.accessions));
}

fn read_lfq_proteins_tsv(path: &Path) -> Result<ProteinLfqData> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("failed to read {}", path.display()))?;

    let headers = reader
        .headers()
        .with_context(|| format!("failed to read headers from {}", path.display()))?
        .clone();

    let mut run_names: Vec<String> = Vec::new();
    let mut intensity_columns: Vec<usize> = Vec::new();
    let mut coverage_columns: Vec<usize> = Vec::new();
    let mut contributor_columns: Vec<usize> = Vec::new();
    for (idx, header) in headers.iter().enumerate() {
        if let Some(name) = header.strip_prefix("intensity:") {
            run_names.push(name.to_string());
            intensity_columns.push(idx);
        } else if header.starts_with("coverage:") {
            coverage_columns.push(idx);
        } else if header.starts_with("contributors:") {
            contributor_columns.push(idx);
        }
    }

    if run_names.is_empty() {
        return Err(anyhow!("no run names found in {}", path.display()));
    }

    if coverage_columns.len() != run_names.len() || contributor_columns.len() != run_names.len() {
        return Err(anyhow!(
            "unexpected LFQ protein column layout in {}",
            path.display()
        ));
    }

    let ordered_coverage = map_dynamic_columns_to_runs(&headers, &run_names, "coverage:", path)?;
    let ordered_contributors =
        map_dynamic_columns_to_runs(&headers, &run_names, "contributors:", path)?;

    let mut proteins: Vec<ProteinRecord> = Vec::new();
    for record in reader.records() {
        let record =
            record.with_context(|| format!("failed to read row from {}", path.display()))?;

        let accessions = record
            .get(0)
            .ok_or_else(|| anyhow!("missing proteins column in {}", path.display()))?
            .split(';')
            .filter(|value| !value.trim().is_empty())
            .map(|value| value.to_string())
            .collect::<Vec<_>>();

        let parse_usize = |value: Option<&str>, column: &str| -> Result<usize> {
            let text =
                value.ok_or_else(|| anyhow!("missing {column} column in {}", path.display()))?;
            if text.trim().is_empty() {
                Ok(0)
            } else {
                text.parse::<usize>()
                    .with_context(|| format!("failed to parse {column}"))
            }
        };

        let q_value: f32 = record
            .get(1)
            .ok_or_else(|| anyhow!("missing q_value column in {}", path.display()))?
            .parse()
            .with_context(|| "failed to parse q_value")?;
        let total_precursors = parse_usize(record.get(2), "total_precursors")?;
        let unique_peptides = parse_usize(record.get(3), "unique_peptides")?;
        let lfq_peptide_count = parse_usize(record.get(4), "lfq_peptide_count")?;

        let mut intensities = Vec::with_capacity(run_names.len());
        for &column_idx in &intensity_columns {
            let value = record
                .get(column_idx)
                .ok_or_else(|| anyhow!("missing intensity column in {}", path.display()))?;
            let intensity = if value.trim().is_empty() {
                0.0
            } else {
                value
                    .parse::<f64>()
                    .with_context(|| format!("failed to parse intensity value '{}'", value))?
            };
            intensities.push(intensity);
        }

        let mut coverage = Vec::with_capacity(run_names.len());
        for (run_idx, &column_idx) in ordered_coverage.iter().enumerate() {
            let value = record
                .get(column_idx)
                .ok_or_else(|| anyhow!("missing coverage column in {}", path.display()))?;
            let trimmed = value.trim();
            let covered = match trimmed {
                "1" | "true" | "TRUE" => true,
                "0" | "false" | "FALSE" | "" => false,
                _ => {
                    return Err(anyhow!(
                        "unexpected coverage value '{}' for run '{}' (column {}) in {}",
                        trimmed,
                        run_names[run_idx],
                        column_idx,
                        path.display()
                    ))
                }
            };
            coverage.push(covered);
        }

        let mut contributors = Vec::with_capacity(run_names.len());
        for &column_idx in &ordered_contributors {
            let value = record
                .get(column_idx)
                .ok_or_else(|| anyhow!("missing contributors column in {}", path.display()))?;
            let parsed = if value.trim().is_empty() {
                0
            } else {
                value
                    .parse::<usize>()
                    .with_context(|| format!("failed to parse contributors value '{}'", value))?
            };
            contributors.push(parsed);
        }

        proteins.push(ProteinRecord {
            accessions,
            q_value,
            total_precursors,
            unique_peptides,
            lfq_peptide_count,
            intensities,
            coverage,
            contributors,
        });
    }

    sort_proteins_by_accession(&mut proteins);

    Ok(ProteinLfqData {
        run_names,
        proteins,
    })
}

fn read_lfq_proteins_parquet(path: &Path, run_names: &[String]) -> Result<Vec<ProteinRecord>> {
    let run_lookup: HashMap<String, usize> = run_names
        .iter()
        .enumerate()
        .map(|(idx, name)| (name.clone(), idx))
        .collect();

    let file =
        fs::File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader = SerializedFileReader::new(file)
        .with_context(|| format!("failed to read parquet from {}", path.display()))?;
    let mut rows = reader
        .get_row_iter(None)
        .with_context(|| format!("failed to iterate rows in {}", path.display()))?;

    let mut proteins: HashMap<String, ProteinRecord> = HashMap::new();

    while let Some(row) = rows.next() {
        let row =
            row.with_context(|| format!("failed to read parquet row from {}", path.display()))?;

        let proteins_field = row
            .get_string(0)
            .with_context(|| format!("missing proteins column in {}", path.display()))?
            .to_string();
        let q_value_raw = row
            .get_double(1)
            .or_else(|_| row.get_float(1).map(f64::from))
            .with_context(|| format!("missing q_value column in {}", path.display()))?;
        if !(f32::MIN as f64..=f32::MAX as f64).contains(&q_value_raw) {
            return Err(anyhow!(
                "protein '{}' has out-of-range q_value {} in {}",
                proteins_field,
                q_value_raw,
                path.display()
            ));
        }
        let q_value = q_value_raw as f32;
        let total_precursors =
            usize::try_from(row.get_int(2).with_context(|| {
                format!("missing total_precursors column in {}", path.display())
            })?)
            .map_err(|_| {
                anyhow!(
                    "protein '{}' has negative total_precursors in {}",
                    proteins_field,
                    path.display()
                )
            })?;
        let unique_peptides =
            usize::try_from(row.get_int(3).with_context(|| {
                format!("missing unique_peptides column in {}", path.display())
            })?)
            .map_err(|_| {
                anyhow!(
                    "protein '{}' has negative unique_peptides in {}",
                    proteins_field,
                    path.display()
                )
            })?;
        let lfq_peptide_count =
            usize::try_from(row.get_int(4).with_context(|| {
                format!("missing lfq_peptide_count column in {}", path.display())
            })?)
            .map_err(|_| {
                anyhow!(
                    "protein '{}' has negative lfq_peptide_count in {}",
                    proteins_field,
                    path.display()
                )
            })?;
        let filename = row
            .get_string(5)
            .with_context(|| format!("missing filename column in {}", path.display()))?
            .to_string();
        let run_idx = run_lookup
            .get(&filename)
            .copied()
            .ok_or_else(|| anyhow!("unexpected run '{filename}' in {}", path.display()))?;
        let intensity = row
            .get_double(6)
            .or_else(|_| row.get_float(6).map(f64::from))
            .with_context(|| format!("missing intensity column in {}", path.display()))?;
        let covered = row
            .get_bool(7)
            .with_context(|| format!("missing covered column in {}", path.display()))?;
        let contributors = usize::try_from(
            row.get_int(8)
                .with_context(|| format!("missing contributors column in {}", path.display()))?,
        )
        .map_err(|_| {
            anyhow!(
                "protein '{}' has negative contributors count in {}",
                proteins_field,
                path.display()
            )
        })?;

        let entry = proteins
            .entry(proteins_field.clone())
            .or_insert_with(|| ProteinRecord {
                accessions: proteins_field
                    .split(';')
                    .filter(|value| !value.trim().is_empty())
                    .map(|value| value.to_string())
                    .collect(),
                q_value,
                total_precursors,
                unique_peptides,
                lfq_peptide_count,
                intensities: vec![0.0; run_names.len()],
                coverage: vec![false; run_names.len()],
                contributors: vec![0; run_names.len()],
            });

        // Sanity check: parquet rows should contain consistent metadata for each protein group.
        if !approx_equal_f32(entry.q_value, q_value) {
            return Err(anyhow!(
                "protein '{}' has inconsistent q_value: {} vs {} in {}",
                proteins_field,
                entry.q_value,
                q_value,
                path.display()
            ));
        }
        if entry.total_precursors != total_precursors {
            return Err(anyhow!(
                "protein '{}' has inconsistent total_precursors: {} vs {} in {}",
                proteins_field,
                entry.total_precursors,
                total_precursors,
                path.display()
            ));
        }
        if entry.unique_peptides != unique_peptides {
            return Err(anyhow!(
                "protein '{}' has inconsistent unique_peptides: {} vs {} in {}",
                proteins_field,
                entry.unique_peptides,
                unique_peptides,
                path.display()
            ));
        }
        if entry.lfq_peptide_count != lfq_peptide_count {
            return Err(anyhow!(
                "protein '{}' has inconsistent lfq_peptide_count: {} vs {} in {}",
                proteins_field,
                entry.lfq_peptide_count,
                lfq_peptide_count,
                path.display()
            ));
        }

        entry.intensities[run_idx] = intensity;
        entry.coverage[run_idx] = covered;
        entry.contributors[run_idx] = contributors;
    }

    let mut proteins: Vec<ProteinRecord> = proteins.into_values().collect();
    sort_proteins_by_accession(&mut proteins);

    Ok(proteins)
}

fn approx_equal(actual: f64, expected: f64) -> bool {
    let scale = actual.abs().max(expected.abs()).max(1.0);
    (actual - expected).abs() <= scale * 1.0e-5
}

fn approx_equal_f32(actual: f32, expected: f32) -> bool {
    approx_equal(f64::from(actual), f64::from(expected))
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
        },
        "lfq_proteins": {
            "enabled": true,
            "min_peptides": 1,
            "min_samples": 1,
            "normalize": false,
            "reference_sample": null
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
        let run_lookup: HashMap<&str, usize> = lfq
            .run_names
            .iter()
            .enumerate()
            .map(|(idx, name)| (name.as_str(), idx))
            .collect();

        let mut expected_intensity: HashMap<(usize, i32), f64> = HashMap::new();
        let mut expected_presence: HashSet<(usize, i32)> = HashSet::new();

        for entry in &golden_data.observations {
            let run_idx = *run_lookup
                .get(entry.run_name.as_str())
                .ok_or_else(|| anyhow!("unknown run '{}' in golden data", entry.run_name))?;
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

    let proteins_tsv = tsv_dir.join("lfq_proteins.tsv");
    assert!(proteins_tsv.exists(), "lfq_proteins.tsv not created");

    let ProteinLfqData {
        run_names: protein_run_names,
        proteins,
    } = read_lfq_proteins_tsv(&proteins_tsv)?;
    assert_eq!(
        protein_run_names, lfq.run_names,
        "LFQ protein TSV run headers did not match the LFQ run names"
    );
    assert!(!proteins.is_empty(), "no proteins quantified");

    for (idx, protein) in proteins.iter().enumerate() {
        assert!(
            protein.q_value <= 1.0,
            "protein {} q-value out of range",
            idx
        );
        assert_eq!(protein.intensities.len(), lfq.run_names.len());
        assert_eq!(protein.coverage.len(), lfq.run_names.len());
        assert_eq!(protein.contributors.len(), lfq.run_names.len());
        assert!(
            protein.total_precursors >= protein.lfq_peptide_count,
            "protein {} has more quantified peptides than total precursors",
            idx
        );
        assert!(
            protein.unique_peptides >= 1,
            "protein {} missing unique peptides",
            idx
        );

        for (run_idx, &covered) in protein.coverage.iter().enumerate() {
            if covered {
                assert!(
                    protein.contributors[run_idx] > 0,
                    "covered run {} has zero contributors",
                    run_idx
                );
            }
        }
    }

    let proteins_golden_path = fixtures.join("lfq_proteins_test.json");
    let golden_proteins: GoldenProteinData = serde_json::from_reader(
        fs::File::open(&proteins_golden_path)
            .with_context(|| format!("failed to open {}", proteins_golden_path.display()))?,
    )
    .with_context(|| format!("failed to parse {}", proteins_golden_path.display()))?;

    assert_eq!(
        golden_proteins.run_names, lfq.run_names,
        "LFQ protein golden run names did not match observed runs"
    );
    assert_eq!(
        proteins.len(),
        golden_proteins.proteins.len(),
        "LFQ protein count mismatch with golden data"
    );

    for golden in &golden_proteins.proteins {
        let key = golden.accessions.clone();
        let actual = proteins
            .iter()
            .find(|record| record.accessions == key)
            .ok_or_else(|| anyhow!("missing protein {:?} in LFQ output", key))?;

        assert!(
            actual.q_value >= golden.min_q_value,
            "protein {:?} q-value {} below golden threshold {}",
            golden.accessions,
            actual.q_value,
            golden.min_q_value
        );
        if golden.has_intensities {
            assert!(
                actual.intensities.iter().any(|value| *value > 0.0),
                "protein {:?} missing non-zero intensity",
                golden.accessions
            );
        }
        let covered_runs = actual.coverage.iter().filter(|covered| **covered).count();
        assert!(
            covered_runs >= golden.min_coverage_count,
            "protein {:?} covered runs {} below expected {}",
            golden.accessions,
            covered_runs,
            golden.min_coverage_count
        );
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

    let proteins_parquet = parquet_dir.join("lfq_proteins.parquet");
    assert!(
        proteins_parquet.exists(),
        "lfq_proteins.parquet not created"
    );
    let parquet_proteins = read_lfq_proteins_parquet(&proteins_parquet, &lfq.run_names)?;
    assert_eq!(
        proteins.len(),
        parquet_proteins.len(),
        "protein parquet row count mismatch"
    );

    let mut parquet_sorted = parquet_proteins;
    sort_proteins_by_accession(&mut parquet_sorted);

    for (tsv_protein, parquet_protein) in proteins.iter().zip(parquet_sorted.iter()) {
        assert_eq!(tsv_protein.accessions, parquet_protein.accessions);
        assert!(
            approx_equal_f32(tsv_protein.q_value, parquet_protein.q_value),
            "protein {:?} q-value mismatch between TSV and parquet",
            tsv_protein.accessions
        );
        assert_eq!(
            tsv_protein.total_precursors, parquet_protein.total_precursors,
            "protein {:?} total_precursors mismatch",
            tsv_protein.accessions
        );
        assert_eq!(
            tsv_protein.unique_peptides, parquet_protein.unique_peptides,
            "protein {:?} unique_peptides mismatch",
            tsv_protein.accessions
        );
        assert_eq!(
            tsv_protein.lfq_peptide_count, parquet_protein.lfq_peptide_count,
            "protein {:?} peptide count mismatch",
            tsv_protein.accessions
        );
        for (run_idx, (tsv_intensity, parquet_intensity)) in tsv_protein
            .intensities
            .iter()
            .zip(parquet_protein.intensities.iter())
            .enumerate()
        {
            assert!(
                approx_equal(*tsv_intensity, *parquet_intensity),
                "protein {:?} intensity mismatch for run {}",
                tsv_protein.accessions,
                run_idx
            );
        }
        assert_eq!(tsv_protein.coverage, parquet_protein.coverage);
        assert_eq!(tsv_protein.contributors, parquet_protein.contributors);
    }

    Ok(())
}
