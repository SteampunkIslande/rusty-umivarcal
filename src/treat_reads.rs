use std::{
    ops::Deref,
    sync::{Arc, Mutex},
};

use indexmap::IndexMap;

use noodles::sam::alignment::Record;

use crate::{commons::read_is_valid, err::UmiVarCalError, pileup::Pileup};
use noodles::bam;
use rayon::iter::{ParallelBridge, ParallelIterator};

pub fn treat_reads(
    bam_file: &std::path::Path,
    pileup: &mut Pileup,
    n_reads: usize,
    output: &std::path::Path,
    min_base_quality: u16,
    min_read_quality: u16,
    min_mapping_quality: u8,
) -> Result<usize, UmiVarCalError> {
    let valid_reads = Arc::new(Mutex::new(0 as usize));

    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file)?;

    let header = reader.read_header()?;

    let pileup = Arc::new(Mutex::new(pileup));

    //Create a set with all UMIs (contained in a Arc of Mutex)
    let all_umis: Arc<Mutex<IndexMap<String, u32>>> = Arc::new(Mutex::new(IndexMap::new()));

    //Read bam_file in parallel, and update pileup
    reader.records().par_bridge().for_each(|record| {
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
        let mapq: u8 = record
            .mapping_quality()
            .expect("Could not get mapping quality!")
            .get();
        let cigar = record.cigar();
        let sequence = record.sequence();
        let base_quals = record.quality_scores();
        let flag = record.flags().bits();

        let mate_pos: usize = match record.mate_alignment_start() {
            Some(Ok(pos)) => pos.get() as usize,
            //Early return if mate is not mapped. TODO: Handle this case for single-end reads.
            _ => return,
        };
        let mut umi = match record.data().get(b"RX") {
            Some(Ok(umi)) => format!("{:?}", umi),
            // No UMI, early return. Maybe it should panic? This read won't be used for variant calling.
            _ => return,
        };
        eprintln!("UMI: {}", umi);

        //umi_pos: last 4 digits of the position
        let umi_pos = if strand == 0 {
            position % 10000
        } else {
            mate_pos % 10000
        };

        let mut found = false;

        for i in umi_pos - 3..umi_pos + 3 {
            let test_umi = format!("{}-{}", umi, i);
            if all_umis.lock().unwrap().contains_key(&test_umi) {
                umi = test_umi;
                found = true;
                break;
            }
        }
        if !found {
            umi = format!("{}-{}", umi, umi_pos);
        }
        all_umis
            .lock()
            .unwrap()
            .entry(umi.clone())
            .and_modify(|e| *e += 1)
            .or_insert(1);
        if read_is_valid(
            flag,
            mapq,
            min_mapping_quality,
            base_quals.as_ref(),
            min_read_quality,
        ) {
            let mut pileup = pileup.lock().unwrap();
            pileup.add_read(
                &umi,
                strand,
                &chromosome,
                position as u32,
                cigar.as_ref(),
                sequence.as_ref(),
                base_quals.as_ref(),
                min_base_quality,
            );
            *valid_reads.lock().unwrap() += 1;
        }
    });
    let valid_reads = valid_reads.lock().unwrap().deref().clone();

    Ok(valid_reads)
}
