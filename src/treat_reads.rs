use std::sync::{Arc, Mutex};

use noodles::sam::alignment::Record;

use crate::{err::UmiVarCalError, pileup::Pileup};
use noodles::bam;
use rayon::iter::{ParallelBridge, ParallelIterator};

pub fn treat_reads(
    bam_file: &std::path::Path,
    pileup: &mut Pileup,
    n_reads: usize,
    output: &std::path::Path,
    min_base_quality: u16,
    min_read_quality: u16,
    min_mapping_quality: u16,
) -> Result<usize, UmiVarCalError> {
    let mut valid_reads = 0;

    let mut current_line: u64 = 1;

    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file)?;

    let header = reader.read_header()?;

    let pileup = Arc::new(Mutex::new(pileup));

    //Read bam_file in parallel, and update pileup
    reader.records().par_bridge().for_each(move |record| {
        let record = record.expect("Could not read record!");
        let strand = record.flags().is_reverse_complemented() as u8;
        let first_in_pair = record.flags().is_first_segment();
        let chromosome: String = record
            .reference_sequence(&header)
            .transpose()
            .expect("Could not get reference sequence!")
            .map(|(name, _)| name.to_string())
            .unwrap();
        let position: usize = record
            .alignment_start()
            .expect("Could not get alignment start!")
            .unwrap()
            .get();
        let mapq: u16 = record
            .mapping_quality()
            .expect("Could not get mapping quality!")
            .get() as u16;
        let cigar: &[u8] = &record.cigar().as_ref();
        let sequence: &[u8] = record.sequence().as_ref();
        let base_quals: &[u8] = record.quality_scores().as_ref();

        let mate_pos: usize = match record.mate_alignment_start() {
            Some(Ok(pos)) => pos.get() as usize,
            //Early return if mate is not mapped. TODO: Handle this case for single-end reads.
            _ => return,
        };
        let umi = match record.data().get(b"RX") {
            Some(Ok(umi)) => format!("{:?}", umi),
            // No UMI, early return. Maybe it should panic?
            _ => return,
        };
        eprintln!("UMI: {}", umi);

        //umi_pos: last 4 digits of the position
        let umi_pos = if strand == 0 {
            position % 10000
        } else {
            mate_pos % 10000
        };

        let ini_umi = format!("{}-{}", umi, umi_pos);
        let mut found = false;

        for i in umi_pos - 3..umi_pos + 3 {
            let test_umi = format!("{}-{}", umi, i);
        }
    });

    Ok(valid_reads)
}
