use crate::pileup::Pileup;
use std::path::Path;

use std::io::{BufReader, Read};

use noodles::core::Region;

use noodles::fasta;

use crate::err::UmiVarCalError;

pub fn add_depth_noise_ref_hp(pileup: &mut Pileup, fasta: &Path) -> Result<(), UmiVarCalError> {
    let mut reader = fasta::io::indexed_reader::Builder::default().build_from_path(fasta)?;

    for (chromosome, infos) in pileup.pileup().iter_mut() {
        for (position, composition) in infos.iter_mut() {
            let reference = reader.index().query(
                &format!("{}:{}-{}", chromosome, position, position + 1)
                    .parse::<Region>()
                    .unwrap(),
            )?;
        }
    }

    Ok(())
}
