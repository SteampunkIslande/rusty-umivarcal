use crate::pileup::Pileup;
use std::path::Path;

use std::io::BufRead;

use noodles::core::Region;

use noodles::fasta;

use crate::err::UmiVarCalError;

fn get_ref_from_fasta_region<R>(
    reader: &mut fasta::IndexedReader<R>,
    chromosome: &str,
    position: usize,
    length: Option<usize>,
) -> Result<String, UmiVarCalError>
where
    R: BufRead + std::io::Seek,
{
    let length = length.unwrap_or(1);
    let query_start = reader.index().query(
        &format!("{}:{}-{}", chromosome, position, position + length)
            .parse::<Region>()
            .unwrap(),
    )?;
    let mut sequence: Vec<u8> = vec![0; length];
    reader
        .get_mut()
        .seek(std::io::SeekFrom::Start(query_start))?;
    reader.get_mut().read_exact(&mut sequence)?;
    let result = String::from_utf8(sequence)?;
    Ok(result)
}

pub fn add_depth_noise_ref_hp(pileup: &mut Pileup, fasta: &Path) -> Result<(), UmiVarCalError> {
    let mut reader = fasta::io::indexed_reader::Builder::default().build_from_path(fasta)?;

    let total_lines: usize = pileup
        .pileup()
        .iter()
        .map(|(_chromosome, infos)| infos.len())
        .sum();

    let mut current_line = 0;

    for (chromosome, infos) in pileup.pileup_mut().iter_mut() {
        for (position, composition) in infos.iter_mut() {
            // Add reference base to pileup composition
            let reference_base =
                get_ref_from_fasta_region(&mut reader, chromosome, *position as usize, None)?
                    .to_uppercase();
            composition.set_reference(&reference_base.bytes().next().unwrap());

            // Add homopolymer length to pileup composition
            let mut hp = 1;
            let seq_around_pos = get_ref_from_fasta_region(
                &mut reader,
                chromosome,
                (*position as usize) - 20,
                Some(41),
            )?;
            hp += seq_around_pos
                .chars()
                .rev()
                .skip(21)
                .take_while(|&c| c == reference_base.chars().next().unwrap())
                .count();
            hp += seq_around_pos
                .chars()
                .skip(21)
                .take_while(|&c| c == reference_base.chars().next().unwrap())
                .count();
            composition.set_homopolymer(hp as u16);

            composition.compute_depth();

            composition.compute_mean_qscore();

            composition.compute_mean_base_qscore();

            current_line += 1;
            eprint!("\rComputing noise...{}%", current_line * 100 / total_lines);
        }
    }

    Ok(())
}
