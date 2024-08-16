mod cli;

use rusty_umivarcal as lib;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli_options = cli::Command::parse();
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli_options.cores())
        .build_global()
        .unwrap();
    lib::call(
        cli_options.input(),
        cli_options.bed(),
        cli_options.fasta(),
        cli_options.min_base_quality(),
        cli_options.min_mapping_quality(),
        cli_options.min_read_quality(),
        cli_options.min_variant_umi(),
        cli_options.min_variant_depth(),
        cli_options.strand_bias_method(),
        cli_options.max_strand_bias(),
        cli_options.pileup(),
        cli_options.output(),
        cli_options.cores(),
        cli_options.alpha(),
        cli_options.max_hp_length(),
        cli_options.gvcf(),
        cli_options.keep_pileup(),
        cli_options.black_list(),
        cli_options.white_list(),
        cli_options.min_phase_umi(),
        cli_options.min_phase_vaf_ratio(),
        cli_options.max_phase_distance(),
        cli_options.compute_coverage_stats(),
    )?;
    Ok(())
}
