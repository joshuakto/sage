//! Use the low-level `parquet` file writer API to serialize Sage results
//!
//! Modifying the file formats here requires some digging into documentation
//! about Dremel definition and repetition levels and the Parquet file format
//! https://akshays-blog.medium.com/wrapping-head-around-repetition-and-definition-levels-in-dremel-powering-bigquery-c1a33c9695da
//! https://blog.twitter.com/engineering/en_us/a/2013/dremel-made-simple-with-parquet
//! https://github.com/apache/parquet-format/blob/master/LogicalTypes.md

#![cfg(feature = "parquet")]

use std::collections::HashMap;
use std::hash::BuildHasher;

use parquet::data_type::{BoolType, ByteArray, FloatType, Int64Type};
use parquet::file::writer::SerializedColumnWriter;
use parquet::{
    basic::ZstdLevel,
    data_type::{ByteArrayType, DataType, Int32Type},
    file::{properties::WriterProperties, writer::SerializedFileWriter},
    schema::types::Type,
};
use sage_core::database::IndexedDatabase;
use sage_core::ion_series::Kind;
use sage_core::lfq::{PeptideQuantTrace, PrecursorId, ProteinRollupResult};
use sage_core::scoring::Feature;
use sage_core::tmt::TmtQuant;

pub fn build_schema() -> Result<Type, parquet::errors::ParquetError> {
    let msg = r#"
        message schema {
            required int64 psm_id;
            required byte_array filename (utf8);
            required byte_array scannr (utf8);
            required byte_array peptide (utf8);
            required byte_array stripped_peptide (utf8);
            required byte_array proteins (utf8);
            required int32 num_proteins;
            required int32 rank;
            required boolean is_decoy;
            required float expmass;
            required float calcmass;
            required int32 charge;
            required int32 peptide_len;
            required int32 missed_cleavages;
            required boolean semi_enzymatic;
            required float ms2_intensity;
            required float isotope_error;
            required float precursor_ppm;
            required float fragment_ppm;
            required float hyperscore;
            required float delta_next;
            required float delta_best;
            required float rt;
            required float aligned_rt;
            required float predicted_rt;
            required float delta_rt_model;
            required float ion_mobility;
            required float predicted_mobility;
            required float delta_mobility;
            required int32 matched_peaks;
            required int32 longest_b;
            required int32 longest_y;
            required float longest_y_pct;
            required float matched_intensity_pct;
            required int32 scored_candidates;
            required float poisson;
            required float sage_discriminant_score;
            required float posterior_error;
            required float spectrum_q;
            required float peptide_q;
            required float protein_q;
            optional group reporter_ion_intensity (LIST) {
                repeated group list {
                    optional float element;
                }
            }
        }
    "#;
    parquet::schema::parser::parse_message_type(msg)
}

/// Caller must guarantee that `reporter_ions` is not an empty slice
fn write_reporter_ions(
    mut column: SerializedColumnWriter,
    features: &[Feature],
    reporter_ions: &[TmtQuant],
) -> parquet::errors::Result<()> {
    let mut scan_map = HashMap::new();

    for r in reporter_ions {
        scan_map.entry((r.file_id, &r.spec_id)).or_insert(r);
    }

    // Caller guarantees `reporter_ions` is not empty
    let channels = reporter_ions[0].peaks.len();

    // https://docs.rs/parquet/44.0.0/parquet/column/index.html
    // Using the low level API here is not very pleasant...
    let def_levels = vec![3; channels];
    let mut rep_levels = vec![1; channels];
    rep_levels[0] = 0;

    let col = column.typed::<FloatType>();
    for feature in features {
        if let Some(rs) = scan_map.get(&(feature.file_id, &feature.spec_id)) {
            col.write_batch(&rs.peaks, Some(&def_levels), Some(&rep_levels))?;
        } else {
            col.write_batch(&[], Some(&[0]), Some(&[0]))?;
        }
    }

    column.close()?;
    Ok(())
}

fn write_null_column(
    mut column: SerializedColumnWriter,
    length: usize,
) -> Result<usize, parquet::errors::ParquetError> {
    let levels = vec![0i16; length];
    let wrote = column
        .typed::<FloatType>()
        .write_batch(&[], Some(&levels), Some(&levels))?;
    column.close().map(|_| wrote)
}

pub fn serialize_features(
    features: &[Feature],
    reporter_ions: &[TmtQuant],
    filenames: &[String],
    database: &IndexedDatabase,
) -> Result<Vec<u8>, parquet::errors::ParquetError> {
    let schema = build_schema()?;

    let options = WriterProperties::builder()
        .set_compression(parquet::basic::Compression::ZSTD(ZstdLevel::try_new(3)?))
        .build();

    let buf = Vec::new();
    let mut writer = SerializedFileWriter::new(buf, schema.into(), options.into())?;

    for features in features.chunks(65536) {
        let mut rg = writer.next_row_group()?;
        macro_rules! write_col {
            ($field:ident, $ty:ident) => {
                if let Some(mut col) = rg.next_column()? {
                    col.typed::<$ty>().write_batch(
                        &features
                            .iter()
                            .map(|f| f.$field as <$ty as DataType>::T)
                            .collect::<Vec<_>>(),
                        None,
                        None,
                    )?;
                    col.close()?;
                }
            };
            ($lambda:expr, $ty:ident) => {
                if let Some(mut col) = rg.next_column()? {
                    col.typed::<$ty>().write_batch(
                        &features.iter().map($lambda).collect::<Vec<_>>(),
                        None,
                        None,
                    )?;
                    col.close()?;
                }
            };
        }

        write_col!(|f: &Feature| f.psm_id as i64, Int64Type);
        write_col!(
            |f: &Feature| filenames[f.file_id].as_str().into(),
            ByteArrayType
        );
        write_col!(|f: &Feature| f.spec_id.as_str().into(), ByteArrayType);
        write_col!(
            |f: &Feature| database[f.peptide_idx].to_string().as_bytes().into(),
            ByteArrayType
        );
        write_col!(
            |f: &Feature| database[f.peptide_idx].sequence.as_ref().into(),
            ByteArrayType
        );
        write_col!(
            |f: &Feature| database[f.peptide_idx]
                .proteins(&database.decoy_tag, database.generate_decoys)
                .as_str()
                .into(),
            ByteArrayType
        );
        write_col!(
            |f: &Feature| database[f.peptide_idx].proteins.len() as i32,
            Int32Type
        );
        write_col!(rank, Int32Type);
        write_col!(|f: &Feature| f.label == -1, BoolType);
        write_col!(expmass, FloatType);
        write_col!(calcmass, FloatType);
        write_col!(charge, Int32Type);
        write_col!(peptide_len, Int32Type);
        write_col!(missed_cleavages, Int32Type);
        write_col!(
            |f: &Feature| database[f.peptide_idx].semi_enzymatic,
            BoolType
        );
        write_col!(ms2_intensity, FloatType);
        write_col!(isotope_error, FloatType);
        write_col!(delta_mass, FloatType);
        write_col!(average_ppm, FloatType);
        write_col!(hyperscore, FloatType);
        write_col!(delta_next, FloatType);
        write_col!(delta_best, FloatType);
        write_col!(rt, FloatType);
        write_col!(aligned_rt, FloatType);
        write_col!(predicted_rt, FloatType);
        write_col!(delta_rt_model, FloatType);
        write_col!(ims, FloatType);
        write_col!(predicted_ims, FloatType);
        write_col!(delta_ims_model, FloatType);
        write_col!(matched_peaks, Int32Type);
        write_col!(longest_b, Int32Type);
        write_col!(longest_y, Int32Type);
        write_col!(longest_y_pct, FloatType);
        write_col!(matched_intensity_pct, FloatType);
        write_col!(scored_candidates, Int32Type);
        write_col!(poisson, FloatType);
        write_col!(discriminant_score, FloatType);
        write_col!(posterior_error, FloatType);
        write_col!(spectrum_q, FloatType);
        write_col!(peptide_q, FloatType);
        write_col!(protein_q, FloatType);

        if let Some(col) = rg.next_column()? {
            if reporter_ions.is_empty() {
                write_null_column(col, features.len())?;
            } else {
                write_reporter_ions(col, features, reporter_ions)?;
            }
        }

        rg.close()?;
    }
    writer.into_inner()
}

pub fn build_matched_fragment_schema() -> parquet::errors::Result<Type> {
    let msg = r#"
        message schema {
            required int64 psm_id;
            required byte_array fragment_type (utf8);
            required int32 fragment_ordinals;
            required int32 fragment_charge;
            required float fragment_mz_experimental;
            required float fragment_mz_calculated;
            required float fragment_intensity;
        }
    "#;

    parquet::schema::parser::parse_message_type(msg)
}

pub fn serialize_matched_fragments(
    features: &[Feature],
) -> Result<Vec<u8>, parquet::errors::ParquetError> {
    let schema = build_matched_fragment_schema()?;

    let options = WriterProperties::builder()
        .set_compression(parquet::basic::Compression::ZSTD(ZstdLevel::try_new(3)?))
        .build();

    let buf = Vec::new();

    let mut writer = SerializedFileWriter::new(buf, schema.into(), options.into())?;

    for features in features.chunks(65536) {
        let mut rg = writer.next_row_group()?;

        if let Some(mut col) = rg.next_column()? {
            let psm_ids = features
                .iter()
                .flat_map(|f| {
                    std::iter::repeat(f.psm_id as i64).take(
                        f.fragments
                            .as_ref()
                            .map(|fragments| fragments.fragment_ordinals.len())
                            .unwrap_or_default(),
                    )
                })
                .collect::<Vec<_>>();

            col.typed::<Int64Type>().write_batch(&psm_ids, None, None)?;
            col.close()?;
        }

        if let Some(mut col) = rg.next_column()? {
            let fragment_types = features
                .iter()
                .flat_map(|f| {
                    f.fragments
                        .as_ref()
                        .map(|fragments| fragments.kinds.iter().copied())
                })
                .flatten()
                .map(|kind| match kind {
                    Kind::A => "a".as_bytes().into(),
                    Kind::B => "b".as_bytes().into(),
                    Kind::C => "c".as_bytes().into(),
                    Kind::X => "x".as_bytes().into(),
                    Kind::Y => "y".as_bytes().into(),
                    Kind::Z => "z".as_bytes().into(),
                })
                .collect::<Vec<ByteArray>>();

            col.typed::<ByteArrayType>()
                .write_batch(&fragment_types, None, None)?;
            col.close()?;
        }

        if let Some(mut col) = rg.next_column()? {
            let fragment_ordinals = features
                .iter()
                .flat_map(|f| {
                    f.fragments
                        .as_ref()
                        .map(|fragments| fragments.fragment_ordinals.iter().copied())
                })
                .flatten()
                .collect::<Vec<_>>();

            col.typed::<Int32Type>()
                .write_batch(&fragment_ordinals, None, None)?;
            col.close()?;
        }

        if let Some(mut col) = rg.next_column()? {
            let fragment_charge = features
                .iter()
                .flat_map(|f| {
                    f.fragments
                        .as_ref()
                        .map(|fragments| fragments.charges.iter().copied())
                })
                .flatten()
                .collect::<Vec<i32>>();

            col.typed::<Int32Type>()
                .write_batch(&fragment_charge, None, None)?;
            col.close()?;
        }

        if let Some(mut col) = rg.next_column()? {
            let fragment_mz_experimental = features
                .iter()
                .flat_map(|f| {
                    f.fragments
                        .as_ref()
                        .map(|fragments| fragments.mz_experimental.iter().copied())
                })
                .flatten()
                .collect::<Vec<_>>();

            col.typed::<FloatType>()
                .write_batch(&fragment_mz_experimental, None, None)?;
            col.close()?;
        }

        if let Some(mut col) = rg.next_column()? {
            let fragment_mz_calculated = features
                .iter()
                .flat_map(|f| {
                    f.fragments
                        .as_ref()
                        .map(|fragments| fragments.mz_calculated.iter().copied())
                })
                .flatten()
                .collect::<Vec<_>>();

            col.typed::<FloatType>()
                .write_batch(&fragment_mz_calculated, None, None)?;
            col.close()?;
        }

        if let Some(mut col) = rg.next_column()? {
            let fragment_intensity = features
                .iter()
                .flat_map(|f| {
                    f.fragments
                        .as_ref()
                        .map(|fragments| fragments.intensities.iter().copied())
                })
                .flatten()
                .collect::<Vec<_>>();

            col.typed::<FloatType>()
                .write_batch(&fragment_intensity, None, None)?;
            col.close()?;
        }

        rg.close()?;
    }

    writer.into_inner()
}

pub fn build_lfq_schema() -> parquet::errors::Result<Type> {
    let msg = r#"
        message schema {
            required byte_array peptide (utf8);
            required byte_array stripped_peptide (utf8);
            optional int32 charge;
            required byte_array proteins (utf8);
            required boolean is_decoy;
            required float q_value;
            required byte_array filename (utf8);
            required float intensity;
        }
    "#;
    parquet::schema::parser::parse_message_type(msg)
}

pub fn build_lfq_protein_schema() -> parquet::errors::Result<Type> {
    let msg = r#"
        message schema {
            required byte_array proteins (utf8);
            required float q_value;
            required int32 total_precursors;
            required int32 unique_peptides;
            required int32 lfq_peptide_count;
            required byte_array filename (utf8);
            required float intensity;
            required boolean covered;
            required int32 contributors;
        }
    "#;
    parquet::schema::parser::parse_message_type(msg)
}

pub fn serialize_lfq<H: BuildHasher>(
    areas: &HashMap<(PrecursorId, bool), PeptideQuantTrace, H>,
    filenames: &[String],
    database: &IndexedDatabase,
) -> parquet::errors::Result<Vec<u8>> {
    let schema = build_lfq_schema()?;

    let options = WriterProperties::builder()
        .set_compression(parquet::basic::Compression::ZSTD(ZstdLevel::try_new(3)?))
        .build();

    let buf = Vec::new();
    let mut writer = SerializedFileWriter::new(buf, schema.into(), options.into())?;
    let mut rg = writer.next_row_group()?;

    if let Some(mut col) = rg.next_column()? {
        let values = areas
            .iter()
            .flat_map(|(_, trace)| {
                let val = database[trace.peptide].to_string().as_bytes().into();
                std::iter::repeat(val).take(filenames.len())
            })
            .collect::<Vec<_>>();

        col.typed::<ByteArrayType>()
            .write_batch(&values, None, None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = areas
            .iter()
            .flat_map(|(_, trace)| {
                let val = database[trace.peptide].sequence.as_ref().into();
                std::iter::repeat(val).take(filenames.len())
            })
            .collect::<Vec<_>>();

        col.typed::<ByteArrayType>()
            .write_batch(&values, None, None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let mut values = Vec::with_capacity(areas.len() * filenames.len());
        let mut def_levels = Vec::with_capacity(areas.len() * filenames.len());

        for (_, trace) in areas.iter() {
            match trace.precursor {
                PrecursorId::Combined(_) => {
                    def_levels.extend(std::iter::repeat(0).take(filenames.len()));
                }
                PrecursorId::Charged((_, charge)) => {
                    values.extend(std::iter::repeat(charge as i32).take(filenames.len()));
                    def_levels.extend(std::iter::repeat(1).take(filenames.len()));
                }
            }
        }

        col.typed::<Int32Type>()
            .write_batch(&values, Some(&def_levels), None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = areas
            .iter()
            .flat_map(|(_, trace)| {
                let val = database[trace.peptide]
                    .proteins(&database.decoy_tag, database.generate_decoys)
                    .as_str()
                    .into();
                std::iter::repeat(val).take(filenames.len())
            })
            .collect::<Vec<_>>();

        col.typed::<ByteArrayType>()
            .write_batch(&values, None, None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = areas
            .iter()
            .flat_map(|(_, trace)| std::iter::repeat(trace.decoy).take(filenames.len()))
            .collect::<Vec<_>>();

        col.typed::<BoolType>().write_batch(&values, None, None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = areas
            .iter()
            .flat_map(|(_, trace)| std::iter::repeat(trace.peak.q_value).take(filenames.len()))
            .collect::<Vec<_>>();

        col.typed::<FloatType>().write_batch(&values, None, None)?;
        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = areas
            .iter()
            .flat_map(|(_, trace)| {
                (0..trace.intensities.len()).map(|idx| filenames[idx].as_bytes().into())
            })
            .collect::<Vec<_>>();

        col.typed::<ByteArrayType>()
            .write_batch(&values, None, None)?;

        col.close()?;
    }

    if let Some(mut col) = rg.next_column()? {
        let values = areas
            .iter()
            .flat_map(|(_, trace)| trace.intensities.iter().copied().map(|v| v as f32))
            .collect::<Vec<_>>();

        col.typed::<FloatType>().write_batch(&values, None, None)?;
        col.close()?;
    }

    rg.close()?;
    writer.into_inner()
}

pub fn serialize_lfq_proteins(
    proteins: &[ProteinRollupResult],
    filenames: &[String],
) -> parquet::errors::Result<Vec<u8>> {
    let schema = build_lfq_protein_schema()?;

    let options = WriterProperties::builder()
        .set_compression(parquet::basic::Compression::ZSTD(ZstdLevel::try_new(3)?))
        .build();

    let buf = Vec::new();
    let mut writer = SerializedFileWriter::new(buf, schema.into(), options.into())?;
    let mut rg = writer.next_row_group()?;

    // Column 1: proteins
    if let Some(mut col) = rg.next_column()? {
        let mut values = Vec::with_capacity(proteins.len() * filenames.len());
        for protein in proteins {
            let joined = protein
                .accessions
                .iter()
                .map(|s| s.as_ref())
                .collect::<Vec<_>>()
                .join(";");
            for _ in filenames {
                values.push(ByteArray::from(joined.as_str()));
            }
        }
        col.typed::<ByteArrayType>()
            .write_batch(&values, None, None)?;
        col.close()?;
    }

    // Column 2: q_value
    if let Some(mut col) = rg.next_column()? {
        let mut values = Vec::with_capacity(proteins.len() * filenames.len());
        for protein in proteins {
            for _ in filenames {
                values.push(protein.q_value);
            }
        }
        col.typed::<FloatType>().write_batch(&values, None, None)?;
        col.close()?;
    }

    // Column 3: total_precursors
    if let Some(mut col) = rg.next_column()? {
        let mut values = Vec::with_capacity(proteins.len() * filenames.len());
        for protein in proteins {
            for _ in filenames {
                values.push(protein.total_precursors as i32);
            }
        }
        col.typed::<Int32Type>().write_batch(&values, None, None)?;
        col.close()?;
    }

    // Column 4: unique_peptides
    if let Some(mut col) = rg.next_column()? {
        let mut values = Vec::with_capacity(proteins.len() * filenames.len());
        for protein in proteins {
            for _ in filenames {
                values.push(protein.unique_peptides as i32);
            }
        }
        col.typed::<Int32Type>().write_batch(&values, None, None)?;
        col.close()?;
    }

    // Column 5: lfq_peptide_count
    if let Some(mut col) = rg.next_column()? {
        let mut values = Vec::with_capacity(proteins.len() * filenames.len());
        for protein in proteins {
            for _ in filenames {
                values.push(protein.quant.peptide_count as i32);
            }
        }
        col.typed::<Int32Type>().write_batch(&values, None, None)?;
        col.close()?;
    }

    // Column 6: filename
    if let Some(mut col) = rg.next_column()? {
        let mut values = Vec::with_capacity(proteins.len() * filenames.len());
        for _ in proteins {
            for filename in filenames {
                values.push(filename.as_bytes().into());
            }
        }
        col.typed::<ByteArrayType>()
            .write_batch(&values, None, None)?;
        col.close()?;
    }

    // Column 7: intensity
    if let Some(mut col) = rg.next_column()? {
        let mut values = Vec::with_capacity(proteins.len() * filenames.len());
        for protein in proteins {
            for sample_idx in 0..filenames.len() {
                let intensity = protein
                    .quant
                    .lfq_intensities
                    .get(sample_idx)
                    .copied()
                    .unwrap_or_default();
                values.push(intensity);
            }
        }
        col.typed::<FloatType>().write_batch(&values, None, None)?;
        col.close()?;
    }

    // Column 8: covered
    if let Some(mut col) = rg.next_column()? {
        let mut values = Vec::with_capacity(proteins.len() * filenames.len());
        for protein in proteins {
            for sample_idx in 0..filenames.len() {
                let covered = protein
                    .quant
                    .sample_coverage
                    .get(sample_idx)
                    .copied()
                    .unwrap_or(false);
                values.push(covered);
            }
        }
        col.typed::<BoolType>().write_batch(&values, None, None)?;
        col.close()?;
    }

    // Column 9: contributors
    if let Some(mut col) = rg.next_column()? {
        let mut values = Vec::with_capacity(proteins.len() * filenames.len());
        for protein in proteins {
            for sample_idx in 0..filenames.len() {
                let contributors = protein
                    .run_coverage
                    .get(sample_idx)
                    .copied()
                    .unwrap_or_default();
                values.push(contributors as i32);
            }
        }
        col.typed::<Int32Type>().write_batch(&values, None, None)?;
        col.close()?;
    }

    rg.close()?;
    writer.into_inner()
}

#[cfg(test)]
mod tests {
    use super::*;
    use bytes::Bytes;
    use fnv::FnvHasher;
    use parquet::column::reader::ColumnReader;
    use parquet::file::reader::{FileReader, SerializedFileReader};
    use parquet::record::{Field, RowAccessor};
    use sage_core::database::PeptideIx;
    use sage_core::enzyme::Position;
    use sage_core::lfq::{Peak, ProteinQuantResult, ProteinRollupResult};
    use sage_core::ml::matrix::Matrix;
    use sage_core::peptide::Peptide;
    use std::hash::BuildHasherDefault;
    use std::sync::Arc;

    fn make_peptide(sequence: &str, proteins: &[&str], decoy: bool) -> Peptide {
        Peptide {
            decoy,
            sequence: Arc::from(sequence.as_bytes()),
            modifications: vec![0.0; sequence.len()],
            nterm: None,
            cterm: None,
            monoisotopic: 0.0,
            missed_cleavages: 0,
            semi_enzymatic: false,
            position: Position::Full,
            proteins: proteins
                .iter()
                .map(|protein| Arc::<str>::from(*protein))
                .collect(),
        }
    }

    fn make_trace(
        precursor: PrecursorId,
        peptide: PeptideIx,
        decoy: bool,
        q_value: f32,
        intensities: Vec<f64>,
    ) -> PeptideQuantTrace {
        PeptideQuantTrace {
            precursor,
            peptide,
            decoy,
            peak: Peak {
                rt: 10,
                spectral_angle: 0.0,
                score: 0.0,
                q_value,
            },
            intensities,
            reference_file_id: 0,
            dot_product: Matrix::zeros(0, 0),
            spectral_angle: Matrix::zeros(0, 0),
            isotope_traces: Matrix::zeros(0, 0),
            raw_isotope_traces: Matrix::zeros(0, 0),
            isotopic_distribution: [0.0; 3],
            time_warps: Vec::new(),
        }
    }

    #[derive(Debug, Clone, PartialEq)]
    struct ExpectedRow {
        peptide: String,
        stripped: String,
        precursor: PrecursorId,
        filename: String,
        intensity: f32,
    }

    #[test]
    fn serialize_lfq_optional_charge_levels() -> parquet::errors::Result<()> {
        let filenames = vec!["run_a".to_string(), "run_b".to_string()];

        let peptides = vec![
            make_peptide("PEPTIDEK", &["PROT1"], false),
            make_peptide("PEPTIDER", &["PROT2"], false),
        ];

        let database = IndexedDatabase {
            peptides,
            fragments: Vec::new(),
            ion_kinds: Vec::new(),
            min_value: Vec::new(),
            potential_mods: Vec::new(),
            bucket_size: 1,
            generate_decoys: false,
            decoy_tag: "rev_".into(),
        };

        let combined_precursor = PrecursorId::Combined(PeptideIx(0));
        let charged_precursor = PrecursorId::Charged((PeptideIx(1), 3));

        let mut areas: HashMap<_, _, BuildHasherDefault<FnvHasher>> = HashMap::default();
        areas.insert(
            (combined_precursor, false),
            make_trace(
                combined_precursor,
                PeptideIx(0),
                false,
                0.01,
                vec![111.25, 222.5],
            ),
        );
        areas.insert(
            (charged_precursor, false),
            make_trace(
                charged_precursor,
                PeptideIx(1),
                false,
                0.02,
                vec![333.75, 444.5],
            ),
        );

        let expected_rows: Vec<ExpectedRow> = areas
            .iter()
            .flat_map(|((precursor, _), trace)| {
                let peptide = database[trace.peptide].to_string();
                let stripped = std::str::from_utf8(&database[trace.peptide].sequence)
                    .unwrap()
                    .to_owned();
                filenames
                    .iter()
                    .zip(trace.intensities.iter())
                    .map(move |(file, intensity)| ExpectedRow {
                        peptide: peptide.clone(),
                        stripped: stripped.clone(),
                        precursor: *precursor,
                        filename: file.clone(),
                        intensity: *intensity as f32,
                    })
            })
            .collect();

        let parquet_buffer = serialize_lfq(&areas, &filenames, &database)?;

        let reader = SerializedFileReader::new(Bytes::from(parquet_buffer))?;
        let actual_rows = reader
            .get_row_iter(None)?
            .collect::<parquet::errors::Result<Vec<_>>>()?;

        assert_eq!(actual_rows.len(), expected_rows.len());

        for (expected, row) in expected_rows.iter().zip(actual_rows.iter()) {
            assert_eq!(row.get_string(0)?, &expected.peptide);
            assert_eq!(row.get_string(1)?, &expected.stripped);

            let charge_field = row
                .get_column_iter()
                .nth(2)
                .expect("charge column present")
                .1;
            match expected.precursor {
                PrecursorId::Combined(_) => {
                    assert!(matches!(charge_field, Field::Null));
                }
                PrecursorId::Charged((_, charge)) => {
                    assert!(matches!(charge_field, Field::Int(value) if *value == charge as i32));
                }
            }

            assert_eq!(row.get_string(6)?, &expected.filename);
            let intensity = row.get_float(7)?;
            assert!(intensity > 0.0, "intensity should retain magnitude");
            assert!((intensity - expected.intensity).abs() < 1e-4);
        }

        let row_group = reader.get_row_group(0)?;
        let num_rows = row_group.metadata().num_rows() as usize;
        assert_eq!(num_rows, expected_rows.len());

        if let ColumnReader::Int32ColumnReader(mut charge_reader) =
            row_group.get_column_reader(2)?
        {
            let mut def_levels = Vec::with_capacity(num_rows);
            let mut charge_values = Vec::with_capacity(num_rows);
            let mut def_batch = Vec::with_capacity(num_rows);
            let mut value_batch = Vec::with_capacity(num_rows);

            while def_levels.len() < num_rows {
                def_batch.clear();
                value_batch.clear();
                let (batch_values, batch_levels) = charge_reader.read_batch(
                    num_rows - def_levels.len(),
                    Some(&mut def_batch),
                    None,
                    &mut value_batch,
                )?;

                assert!(batch_levels > 0, "no definition levels read");
                def_levels.extend_from_slice(&def_batch[..batch_levels]);
                charge_values.extend_from_slice(&value_batch[..batch_values]);
            }

            assert_eq!(def_levels.len(), num_rows);

            let mut extracted = Vec::with_capacity(num_rows);
            let mut value_idx = 0;
            for level in def_levels.into_iter() {
                if level == 0 {
                    extracted.push(None);
                } else {
                    extracted.push(Some(charge_values[value_idx]));
                    value_idx += 1;
                }
            }

            assert_eq!(value_idx, charge_values.len());

            let expected_charge: Vec<Option<i32>> = expected_rows
                .iter()
                .map(|row| match row.precursor {
                    PrecursorId::Combined(_) => None,
                    PrecursorId::Charged((_, charge)) => Some(charge as i32),
                })
                .collect();

            assert_eq!(extracted, expected_charge);
        } else {
            panic!("charge column reader not available");
        }

        if let ColumnReader::ByteArrayColumnReader(mut filename_reader) =
            row_group.get_column_reader(6)?
        {
            let mut filename_values = Vec::with_capacity(num_rows);
            let mut filename_batch = Vec::with_capacity(num_rows);
            let mut levels_read = 0usize;

            while levels_read < num_rows {
                filename_batch.clear();
                let (batch_values, batch_levels) = filename_reader.read_batch(
                    num_rows - levels_read,
                    None,
                    None,
                    &mut filename_batch,
                )?;

                assert!(batch_levels > 0, "no filenames read");
                levels_read += batch_levels;
                filename_values.extend(filename_batch.drain(0..batch_values));
            }

            assert_eq!(levels_read, num_rows);

            let filenames_from_column: Vec<String> = filename_values
                .into_iter()
                .map(|ba| std::str::from_utf8(ba.data()).unwrap().to_owned())
                .collect();

            for (expected, actual) in expected_rows.iter().zip(filenames_from_column.iter()) {
                assert_eq!(actual, &expected.filename);
            }
        } else {
            panic!("filename column reader not available");
        }

        if let ColumnReader::FloatColumnReader(mut intensity_reader) =
            row_group.get_column_reader(7)?
        {
            let mut intensity_values = Vec::with_capacity(num_rows);
            let mut intensity_batch = Vec::with_capacity(num_rows);
            let mut levels_read = 0usize;

            while levels_read < num_rows {
                intensity_batch.clear();
                let (batch_values, batch_levels) = intensity_reader.read_batch(
                    num_rows - levels_read,
                    None,
                    None,
                    &mut intensity_batch,
                )?;

                assert!(batch_levels > 0, "no intensities read");
                levels_read += batch_levels;
                intensity_values.extend_from_slice(&intensity_batch[..batch_values]);
            }

            assert_eq!(levels_read, num_rows);

            for (expected, actual) in expected_rows.iter().zip(intensity_values.iter()) {
                assert!((*actual - expected.intensity).abs() < 1e-4);
                assert!(*actual > 0.0);
            }
        } else {
            panic!("intensity column reader not available");
        }

        Ok(())
    }

    #[test]
    fn serialize_lfq_proteins_round_trip() -> parquet::errors::Result<()> {
        let filenames = vec!["run_a".to_string(), "run_b".to_string()];
        
        // Build ProteinRollupResult with new schema fields
        let proteins = vec![ProteinRollupResult {
            accessions: vec![Arc::<str>::from("P1"), Arc::<str>::from("P2")],
            quant: ProteinQuantResult {
                protein_ids: vec!["P1".to_string(), "P2".to_string()],
                lfq_intensities: vec![123.0, 456.0],
                peptide_count: 3,
                sample_coverage: vec![true, false],  // run_b not covered by MaxLFQ
            },
            q_value: 0.123,
            total_precursors: 5,
            unique_peptides: 4,
            run_coverage: vec![2, 0],  // run_b has 0 contributors (not covered)
        }];

        let parquet_buffer = serialize_lfq_proteins(&proteins, &filenames)?;
        let reader = SerializedFileReader::new(Bytes::from(parquet_buffer))?;
        let mut rows = reader.get_row_iter(None)?;

        // Verify first row (run_a)
        // Schema: proteins, q_value, total_precursors, unique_peptides, lfq_peptide_count,
        //         filename, intensity, covered, contributors
        let first = rows.next().expect("protein row")?;
        assert_eq!(first.get_string(0)?, "P1;P2");                           // proteins
        assert!((first.get_float(1)? - 0.123).abs() < f32::EPSILON);         // q_value
        assert_eq!(first.get_int(2)?, 5);                                    // total_precursors
        assert_eq!(first.get_int(3)?, 4);                                    // unique_peptides
        assert_eq!(first.get_int(4)?, 3);                                    // lfq_peptide_count
        assert_eq!(first.get_string(5)?, "run_a");                           // filename
        assert!((first.get_float(6)? - 123.0).abs() < f32::EPSILON);         // intensity
        assert!(first.get_bool(7)?);                                         // covered
        assert_eq!(first.get_int(8)?, 2);                                    // contributors

        // Verify second row (run_b)
        let second = rows.next().expect("second protein row")?;
        assert_eq!(second.get_string(0)?, "P1;P2");
        assert!((second.get_float(1)? - 0.123).abs() < f32::EPSILON);
        assert_eq!(second.get_int(2)?, 5);
        assert_eq!(second.get_int(3)?, 4);
        assert_eq!(second.get_int(4)?, 3);
        assert_eq!(second.get_string(5)?, "run_b");
        assert!((second.get_float(6)? - 456.0).abs() < f32::EPSILON);
        assert!(!second.get_bool(7)?);                                       // not covered
        assert_eq!(second.get_int(8)?, 0);                                   // contributors (0 because not covered)

        assert!(rows.next().is_none());

        Ok(())
    }
}
