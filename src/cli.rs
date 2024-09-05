use rusty_umivarcal::commons::{self, UmiSource};

#[derive(clap::ValueEnum, std::fmt::Debug, Clone)]
pub enum StrandBiasMethod {
    Default,
    Torrent,
}

#[derive(clap::Parser, std::fmt::Debug)]
#[group(multiple = false)]
struct Umi {
    /// Get UMI from the read name
    #[clap(long = "umi-from-read-name")]
    from_read_name: bool,

    /// Get UMI from a BAM tag
    #[clap(long = "umi-from-tag")]
    tag: Option<String>,

    /// Get UMI from the read sequence, with a fixed length
    #[clap(long = "umi-from-length")]
    from_length: Option<usize>,
}

#[derive(clap::Parser, std::fmt::Debug)]
#[command(
    name = "umivarcal",
    author = "Charles Monod-Broca",
    version,
    about = "Rust implementation of UMIVarCal (see https://gitlab.com/vincent-sater/umi-varcal/-/tree/master for the original implementation)"
)]
pub struct Command {
    /// Fasta file
    #[clap(short = 'f', long = "fasta")]
    fasta: std::path::PathBuf,
    /// Bed file
    #[clap(short = 'b', long = "bed")]
    bed: std::path::PathBuf,
    /// Input path
    #[clap(short = 'i', long = "input")]
    input: std::path::PathBuf,

    /// Pileup file
    #[clap(short = 'p', long = "pileup")]
    pileup: Option<std::path::PathBuf>,

    /// Output path (directory)
    #[clap(short = 'o', long = "output")]
    output: Option<std::path::PathBuf>,

    /// Number of threads
    #[clap(short = 'c', long = "cores", default_value = "1")]
    cores: usize,

    /// Minimum base quality
    #[clap(long = "min_base_quality", default_value = "10")]
    min_base_quality: u16,

    /// Minimum read quality
    #[clap(long = "min_read_quality", default_value = "20")]
    min_read_quality: u16,

    /// Minimum mapping quality
    #[clap(long = "min_mapping_quality", default_value = "20")]
    min_mapping_quality: u8,

    /// Minimum variant UMI
    #[clap(long = "min_variant_umi", default_value = "5")]
    min_variant_umi: u16,

    /// Minimum variant depth
    #[clap(long = "min_variant_depth", default_value = "5")]
    min_variant_depth: u16,

    /// Alpha
    #[clap(long = "alpha", default_value = "0.05")]
    alpha: f64,

    /// Strand bias method
    #[clap(long = "strand_bias_method")]
    strand_bias_method: Option<StrandBiasMethod>,

    /// Max strand bias
    #[clap(long = "max_strand_bias")]
    max_strand_bias: Option<f64>,

    /// Max homopolymer length
    #[clap(long = "max_hp_length", default_value = "10")]
    max_hp_length: u16,

    /// Whether to output GVCF
    #[clap(long = "gvcf")]
    gvcf: bool,

    /// Whether to keep pileup
    #[clap(short = 'k', long = "keep_pileup")]
    keep_pileup: bool,

    /// Whether to compute coverage stats
    #[clap(long = "compute_coverage_stats")]
    compute_coverage_stats: bool,

    /// Blacklist file
    #[clap(long = "black_list")]
    black_list: Option<std::path::PathBuf>,

    /// White list file
    #[clap(long = "white_list")]
    white_list: Option<std::path::PathBuf>,

    /// Minimum phase UMI
    #[clap(long = "min_phase_umi", default_value = "3")]
    min_phase_umi: u16,

    /// Minimum phase vaf ratio
    #[clap(long = "min_phase_vaf_ratio", default_value = "0.8")]
    min_phase_vaf_ratio: f64,

    /// Maximum phase distance
    #[clap(long = "max_phase_distance", default_value = "100")]
    max_phase_distance: u16,

    /// Whether to save messagepack objects with or without names. Please use for debugging purposes only (increases file size a lot).
    #[clap(long = "msgpack-with-names")]
    msgpack_with_names: bool,

    /// Where to get UMI barcodes from. Three optional mutually exclusive flags:
    /// 1. `--umi-from-read-name` to get UMI from the read name.
    /// 2. `--umi-from-tag` to get UMI from a BAM tag (accepts tag name, but defaults to `RX`).
    /// 3. `--umi-from-length` to get UMI from the read sequence, with a fixed length.
    #[command(flatten)]
    umi: Umi,
}

impl Command {
    pub fn parse() -> Self {
        <Self as clap::Parser>::parse()
    }

    pub fn fasta(&self) -> &std::path::Path {
        &self.fasta
    }

    pub fn bed(&self) -> &std::path::Path {
        &self.bed
    }

    pub fn input(&self) -> &std::path::Path {
        &self.input
    }

    pub fn pileup(&self) -> Option<&std::path::Path> {
        self.pileup.as_deref()
    }

    pub fn output(&self) -> Option<&std::path::Path> {
        self.output.as_deref()
    }

    pub fn cores(&self) -> usize {
        self.cores
    }

    pub fn min_base_quality(&self) -> u16 {
        self.min_base_quality
    }

    pub fn min_read_quality(&self) -> u16 {
        self.min_read_quality
    }

    pub fn min_mapping_quality(&self) -> u8 {
        self.min_mapping_quality
    }

    pub fn min_variant_umi(&self) -> u16 {
        self.min_variant_umi
    }

    pub fn min_variant_depth(&self) -> u16 {
        self.min_variant_depth
    }

    pub fn alpha(&self) -> f64 {
        self.alpha
    }

    pub fn strand_bias_method(&self) -> commons::StrandBiasMethod {
        if let Some(value) = &self.strand_bias_method {
            match value {
                StrandBiasMethod::Default => commons::StrandBiasMethod::Default,
                StrandBiasMethod::Torrent => commons::StrandBiasMethod::Torrent,
            }
        } else {
            commons::StrandBiasMethod::Default
        }
    }

    pub fn max_strand_bias(&self) -> f64 {
        match self.max_strand_bias {
            Some(value) => value,
            None => match self.strand_bias_method() {
                commons::StrandBiasMethod::Default => 1.0,
                commons::StrandBiasMethod::Torrent => 0.743,
            },
        }
    }

    pub fn max_hp_length(&self) -> u16 {
        self.max_hp_length
    }

    pub fn gvcf(&self) -> bool {
        self.gvcf
    }

    pub fn keep_pileup(&self) -> bool {
        self.keep_pileup
    }

    pub fn compute_coverage_stats(&self) -> bool {
        self.compute_coverage_stats
    }

    pub fn black_list(&self) -> Option<&std::path::Path> {
        self.black_list.as_deref()
    }

    pub fn white_list(&self) -> Option<&std::path::Path> {
        self.white_list.as_deref()
    }

    pub fn min_phase_umi(&self) -> u16 {
        self.min_phase_umi
    }

    pub fn min_phase_vaf_ratio(&self) -> f64 {
        self.min_phase_vaf_ratio
    }

    pub fn max_phase_distance(&self) -> u16 {
        self.max_phase_distance
    }

    pub fn msgpack_with_names(&self) -> bool {
        self.msgpack_with_names
    }

    pub fn umi_source(&self) -> UmiSource {
        match &self.umi {
            Umi {
                from_read_name: true,
                tag: None,
                from_length: None,
            } => UmiSource::ReadName,
            Umi {
                from_read_name: false,
                tag: Some(tag),
                from_length: None,
            } => UmiSource::Tag(tag.clone()),
            Umi {
                from_read_name: false,
                tag: None,
                from_length: Some(length),
            } => UmiSource::Length(*length),
            Umi {
                from_read_name: false,
                tag: None,
                from_length: None,
            } => UmiSource::Tag("RX".to_string()),
            _ => panic!("Arguments are mutually exclusive!"),
        }
    }
}
