use err::UmiVarCalError;
use noodles::{bam, fasta};
use std::fs::File;
use std::io::BufReader;

pub mod err;

pub mod commons;
use commons::StrandBiasMethod;

pub mod treat_reads;

pub mod pileup;
use pileup::Pileup;

#[allow(unused)]
pub fn call(
    input: &std::path::Path,
    bed: &std::path::Path,
    fasta: &std::path::Path,
    min_base_quality: u16,
    min_mapping_quality: u16,
    min_read_quality: u16,
    min_variant_umi: u16,
    min_variant_depth: u16,
    strand_bias_method: StrandBiasMethod,
    max_strand_bias: f64,
    pileup: Option<&std::path::Path>,
    output: Option<&std::path::Path>,
    cores: usize,
    alpha: f64,
    max_hp_length: u16,
    gvcf: bool,
    keep_pileup: bool,
    black_list: Option<&std::path::Path>,
    white_list: Option<&std::path::Path>,
    min_phase_umi: u16,
    min_phase_vaf_ratio: f64,
    max_phase_distance: u16,
    compute_coverage_stats: bool,
) -> Result<(), UmiVarCalError> {
    let mut cores = cores;

    // Sample name is the file name without the extension. Used to name the output files.
    let sample_name = input
        .file_stem()
        .expect("No file name!")
        .to_str()
        .expect("Path cannot be converted to string!");

    // Output is a directory. If not provided, use the current directory.
    let output = output.unwrap_or(std::path::Path::new("."));
    // Create the output directory if it does not exist.
    std::fs::create_dir_all(output).expect("Could not create output directory!");

    let f = fasta::Reader::new(BufReader::new(File::open(fasta)?));
    let mut reader = bam::io::reader::Builder::default().build_from_path(input)?;
    reader.read_header()?;
    let read_count = reader.records().count();
    // Exhausted iterator, drop it to avoid accidental use.
    drop(reader);

    let mut valid_reads = 0;

    let rebuild = pileup.is_none();

    let mut pileup = match pileup {
        None => Pileup::pileup_from_bed(bed),
        Some(pileup) => Pileup::load_pileup(pileup),
    };

    let pileup_len = pileup.len();
    if cores > 1 {
        if (pileup_len > 1000000 && pileup_len <= 2000000 && read_count < 5000000)
            || (pileup_len > 2000000 && pileup_len <= 5000000 && read_count < 20000000)
            || (pileup_len > 5000000 && pileup_len <= 10000000 && read_count < 50000000)
            || (pileup_len > 10000000)
        {
            eprintln!("Warning: Using more cores to analyze the provided data will not result in any significant performance gains!\nLaunching UMI-VarCal on one core only...\n");
            cores = 1;
        }
    }
    if rebuild {
        valid_reads = treat_reads::treat_reads(
            &input,
            &mut pileup,
            read_count,
            &output,
            min_base_quality,
            min_read_quality,
            min_mapping_quality,
        )?;
        if keep_pileup {
            pileup.save_pileup_bed(output.join(sample_name).to_str().unwrap(), bed);
        }
    }

    Ok(())
}
