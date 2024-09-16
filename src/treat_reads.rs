use std::collections::HashMap;

use indexmap::IndexMap;

use noodles::sam::alignment::Record;

use crate::{
    commons::{read_is_valid, UmiSource},
    err::{PileupError, UmiVarCalError},
    pileup::Pileup,
};
use noodles::bam;

pub fn treat_reads(
    bam_file: &std::path::Path,
    pileup: &mut Pileup,
    n_reads: usize,
    output: &std::path::Path,
    min_base_quality: u16,
    min_read_quality: u16,
    min_mapping_quality: u8,
    umi_source: UmiSource,
) -> Result<(usize, HashMap<String, i32>), UmiVarCalError> {
    let mut valid_reads = 0;

    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file)?;

    let header = reader.read_header()?;

    //Create a set with all UMIs (contained in a Arc of Mutex)
    let mut all_umis: IndexMap<String, u32> = IndexMap::new();

    let mut encountered_errors = HashMap::new();

    //Read bam_file in parallel, and update pileup
    for record in reader.records() {
        let record = record.expect("Could not read record!");

        let strand = record.flags().is_reverse_complemented() as u8;
        let first_in_pair = record.flags().is_first_segment();
        let chromosome: Option<String> = record
            .reference_sequence(&header)
            .transpose()
            .ok()
            .flatten()
            .map(|(name, _)| name.to_string());
        let position: Option<usize> = record
            .alignment_start()
            .transpose()
            .ok()
            .flatten()
            .map(|pos| pos.get());
        let mapq = record.mapping_quality().map(|q| q.get()).unwrap_or(0);
        let cigar = record.cigar();
        let sequence = record.sequence();
        let base_quals = record.quality_scores();
        let flag = record.flags();
        if position.is_none() {
            continue;
        }
        if chromosome.is_none() {
            continue;
        }
        let position = position.unwrap();
        let chromosome = chromosome.unwrap();

        let mate_pos: usize = match record.mate_alignment_start() {
            Some(Ok(pos)) => pos.get() as usize,
            //Early return if mate is not mapped. TODO: Handle this case for single-end reads.
            _ => {
                eprintln!("TODO: Handle single-end reads!");
                return Err(UmiVarCalError::SingleEndReadsNotSupported);
            }
        };
        let mut umi = match &umi_source {
            UmiSource::ReadName => {
                let read_name = record.name().expect("Could not get read name!");
                let umi = read_name.rsplitn(2, |s| *s == b'_').next().unwrap();
                String::from_utf8(umi.to_vec()).unwrap()
            }
            UmiSource::Tag(tag) => {
                match record.data().get::<[u8; 2]>(
                    tag.as_bytes()[..2]
                        .try_into()
                        .expect(&format!("Cannot read tag {:?}!", tag)),
                ) {
                    Some(Ok(umi)) => format!("{:?}", umi),
                    // No UMI, early return. Maybe it should panic? This read won't be used for variant calling.
                    _ => {
                        return Err(UmiVarCalError::NoUmi(
                            record
                                .name()
                                .map_or("Unnamed".to_string(), |n| n.to_string()),
                        ))
                    }
                }
            }
            UmiSource::Length(len) => {
                let umi = sequence.iter().take(*len).collect::<Vec<u8>>();
                String::from_utf8(umi).unwrap()
            }
        };

        //umi_pos: last 4 digits of the position
        let umi_pos = if strand == 0 {
            position % 10000
        } else {
            mate_pos % 10000
        };

        let mut found = false;

        for i in umi_pos - 3..umi_pos + 3 {
            let test_umi = format!("{}-{}", umi, i);
            if all_umis.contains_key(&test_umi) {
                umi = test_umi;
                found = true;
                break;
            }
        }
        if !found {
            umi = format!("{}-{}", umi, umi_pos);
        }
        all_umis
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
            valid_reads += 1;
            match pileup.add_read(
                &umi,
                strand,
                &chromosome,
                position as u32,
                &cigar,
                &sequence,
                base_quals.as_ref(),
                min_base_quality,
            ) {
                Ok(_) => {}
                Err(e) => match e {
                    UmiVarCalError::NoReferenceSequence => {
                        encountered_errors
                            .entry("No reference sequence found".to_string())
                            .and_modify(|e| *e += 1)
                            .or_insert(1);
                    }
                    UmiVarCalError::PileupError(PileupError::ChromosomeNotFound(c)) => {
                        encountered_errors
                            .entry(format!("Chromosome not found: {}", c))
                            .and_modify(|e| *e += 1)
                            .or_insert(1);
                    }
                    UmiVarCalError::PileupError(PileupError::PositionNotFound(_)) => {
                        encountered_errors
                            .entry(format!("Position not found"))
                            .and_modify(|e| *e += 1)
                            .or_insert(1);
                    }
                    _ => {}
                },
            }
        }
    }

    Ok((valid_reads, encountered_errors))
}
