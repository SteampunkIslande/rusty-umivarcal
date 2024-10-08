use core::str;
use indexmap::IndexMap;
use noodles::bam::record::{self, Sequence};
use noodles::sam::alignment::record::cigar;
use std::collections::HashSet;
use std::fs::File;

use crate::commons::BedRecord;
use crate::err::{PileupError, UmiVarCalError};

use std::cmp;

#[derive(Debug, Default, serde::Deserialize, serde::Serialize)]
pub struct InsertionDict {
    insertions: IndexMap<String, NucleotideCounter>,

    mean_qscore: Option<f32>,
    qscore_stdev: Option<f32>,
}

#[derive(Debug, Default, serde::Deserialize, serde::Serialize)]
pub struct DeletionDict {
    deletions: IndexMap<u32, NucleotideCounter>,
}

#[derive(Debug, Default, serde::Deserialize, serde::Serialize)]
pub struct NucleotideCounter {
    forward: u32,
    reverse: u32,
    umis: HashSet<String>,
    qscore: Option<u8>,
}

#[derive(Debug, Default, serde::Deserialize, serde::Serialize)]
pub struct PileupCounter {
    // Reference base at this position
    reference: char,

    //Keep track of the number of each nucleotide at this position (plus strand), the set of UMIs, and the sum of the quality scores.
    a: NucleotideCounter,
    c: NucleotideCounter,
    g: NucleotideCounter,
    t: NucleotideCounter,
    insertions: InsertionDict,
    deletions: DeletionDict,
    base_error_probability: f32,

    // Homopolymer length around this position
    hp: u16,

    // Number of reads at this position
    depth: u32,

    // Quality scores at this position: do not serialize
    #[serde(skip)]
    qscore: Vec<f32>,

    // Mean quality score at this position.
    // If None, this means that the quality scores have not been computed yet.
    mean_qscore: Option<f32>,

    // Standard deviation of quality scores at this position.
    // If None, this means that the quality scores have not been computed yet.
    qscore_stdev: Option<f32>,
}

impl PileupCounter {
    pub fn add_nucleotide(&mut self, nuc: &u8, umi: &str, qscore: u8, strand: u8) {
        match (nuc, strand) {
            (b'A', 0) => {
                self.a.add_forward(umi, qscore);
            }
            (b'C', 0) => {
                self.c.add_forward(umi, qscore);
            }
            (b'G', 0) => {
                self.g.add_forward(umi, qscore);
            }
            (b'T', 0) => {
                self.t.add_forward(umi, qscore);
            }
            (b'A', 1) => {
                self.a.add_reverse(umi, qscore);
            }
            (b'C', 1) => {
                self.c.add_reverse(umi, qscore);
            }
            (b'G', 1) => {
                self.g.add_reverse(umi, qscore);
            }
            (b'T', 1) => {
                self.t.add_reverse(umi, qscore);
            }
            _ => todo!("Found {:?} in pileup!", *nuc),
        }
        self.qscore.push(qscore as f32);
    }

    pub fn add_insertion(&mut self, seq: &str, umi: &str, qscore: u8, strand: u8) {
        self.insertions.insert(seq, umi, qscore, strand);
    }

    pub fn add_deletion(&mut self, position: u32, umi: &str, strand: u8) {
        self.deletions.add_deletion(position, umi, strand);
    }

    pub fn set_reference(&mut self, nuc: &u8) {
        self.reference = *(nuc) as char;
    }

    pub fn set_homopolymer(&mut self, hp: u16) {
        self.hp = hp;
    }

    pub fn compute_depth(&mut self) {
        self.depth = 0
            + self.get_a_depth()
            + self.get_c_depth()
            + self.get_g_depth()
            + self.get_t_depth()
            + self.get_deletions_count()
            + self.get_insertions_count();
    }

    pub fn compute_mean_base_qscore(&mut self) {
        let mut qscore = 0f32;
        let mut n = 0;
        if self.a.qscore.is_some() {
            qscore += self.a.qscore.unwrap() as f32;
            n += self.a.forward + self.a.reverse;
        }
        if self.c.qscore.is_some() {
            qscore += self.c.qscore.unwrap() as f32;
            n += self.c.forward + self.c.reverse;
        }
        if self.g.qscore.is_some() {
            qscore += self.g.qscore.unwrap() as f32;
            n += self.g.forward + self.g.reverse;
        }
        if self.t.qscore.is_some() {
            qscore += self.t.qscore.unwrap() as f32;
            n += self.t.forward + self.t.reverse;
        }
        if n > 0 {
            self.mean_qscore = Some(qscore / n as f32);
        } else {
            self.mean_qscore = None;
        }
    }

    pub fn compute_mean_qscore(&mut self) {
        self.mean_qscore = Some(self.qscore.iter().sum::<f32>() / self.qscore.len() as f32);
        self.qscore_stdev = Some(
            ((self
                .qscore
                .iter()
                .map(|x| (x - self.mean_qscore.unwrap()).powi(2))
                .sum::<f32>()
                / self.qscore.len() as f32)
                .sqrt()
                * 1e6)
                .round()
                / 1e6,
        );
        self.base_error_probability =
            (10f32.powf(-self.mean_qscore.unwrap().floor() / 10f32) * 1e6).round() / 1e6;
    }

    fn get_a_depth(&self) -> u32 {
        self.a.forward + self.a.reverse
    }

    fn get_c_depth(&self) -> u32 {
        self.c.forward + self.c.reverse
    }

    fn get_g_depth(&self) -> u32 {
        self.g.forward + self.g.reverse
    }

    fn get_t_depth(&self) -> u32 {
        self.t.forward + self.t.reverse
    }

    fn get_insertions_count(&self) -> u32 {
        self.insertions.insertions.len() as u32
    }

    fn get_deletions_count(&self) -> u32 {
        self.deletions.deletions.len() as u32
    }
}

impl NucleotideCounter {
    pub fn add_forward(&mut self, umi: &str, qscore: u8) {
        self.forward += 1;
        self.umis.insert(umi.to_string());
        self.qscore = Some(qscore + self.qscore.unwrap_or(0));
    }
    pub fn add_reverse(&mut self, umi: &str, qscore: u8) {
        self.reverse += 1;
        self.umis.insert(umi.to_string());
        self.qscore = Some(qscore + self.qscore.unwrap_or(0));
    }
}

impl InsertionDict {
    pub fn insert(&mut self, nuc: &str, umi: &str, qscore: u8, strand: u8) {
        match strand {
            0 => self
                .insertions
                .entry(nuc.to_string())
                .or_insert(NucleotideCounter::default())
                .add_forward(umi, qscore),
            1 => self
                .insertions
                .entry(nuc.to_string())
                .or_insert(NucleotideCounter::default())
                .add_reverse(umi, qscore),
            _ => (),
        }
    }
}

impl DeletionDict {
    pub fn add_deletion(&mut self, position: u32, umi: &str, strand: u8) {
        match strand {
            0 => self
                .deletions
                .entry(position)
                .or_insert(NucleotideCounter::default())
                .add_forward(umi, 0),
            1 => self
                .deletions
                .entry(position)
                .or_insert(NucleotideCounter::default())
                .add_reverse(umi, 0),
            _ => (),
        }
    }
}

#[derive(Debug, Default, serde::Deserialize, serde::Serialize)]
pub struct Pileup {
    pileup: IndexMap<String, IndexMap<u32, PileupCounter>>,
}

fn cmp_chromosome(chr1: &str, chr2: &str) -> cmp::Ordering {
    let chr1: u32 = chr1
        .chars()
        .filter(|c| c.is_digit(10))
        .collect::<String>()
        .parse()
        .unwrap_or(0);
    let chr2: u32 = chr2
        .chars()
        .filter(|c| c.is_digit(10))
        .collect::<String>()
        .parse()
        .unwrap_or(0);
    chr1.cmp(&chr2)
}

impl Pileup {
    pub fn save_pileup(
        &self,
        output_prefix: &str,
        chrom: Option<&str>,
        msgpack_save_with_names: bool,
    ) -> Result<(), UmiVarCalError> {
        let file = File::create(match chrom {
            Some(chrom) => format!("{}.{}.pileup", output_prefix, chrom),
            None => format!("{}.pileup", output_prefix),
        })
        .unwrap();
        let mut wr = std::io::BufWriter::new(&file);

        let mut obj = IndexMap::new();
        for (pileup_chrom, infos) in self.pileup.iter() {
            if chrom.is_some() && chrom.unwrap() != pileup_chrom {
                continue;
            }
            obj.insert(pileup_chrom.to_owned(), infos);
        }
        match msgpack_save_with_names {
            true => rmp_serde::encode::write_named(&mut wr, &obj)?,
            false => rmp_serde::encode::write(&mut wr, &obj)?,
        }
        Ok(())
    }

    /// Create new, empty pileup from bed file.
    /// Prepares all the necessary data structures from the bed file.
    pub fn new_pileup_from_bed(bed: &std::path::Path) -> Self {
        let mut pileup = IndexMap::new();
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .comment(Some(b'#'))
            .from_path(bed)
            .expect("Could not read bed file!");

        for record in reader.records() {
            let record = record.expect("Could not read record!");
            let bed_record: BedRecord = record
                .deserialize(None)
                .expect("Could not deserialize record!");
            let chrom_entry = pileup.entry(bed_record.chrom).or_insert(IndexMap::new());
            for i in bed_record.start..=bed_record.end {
                chrom_entry.insert(i, PileupCounter::default());
            }
        }

        let mut result = Self { pileup };
        result.sort();
        result
    }

    /// Load pileup from serialized messagepack file.
    pub fn load_pileup(pileup: &std::path::Path) -> Self {
        let file = File::open(pileup).expect("Could not open pileup file!");
        let mut rd = std::io::BufReader::new(&file);
        let pileup: IndexMap<String, IndexMap<u32, PileupCounter>> =
            rmp_serde::decode::from_read(&mut rd).expect("Could not decode pileup file!");
        Self { pileup }
    }

    pub fn pileup_mut(&mut self) -> &mut IndexMap<String, IndexMap<u32, PileupCounter>> {
        &mut self.pileup
    }

    pub fn pileup(&self) -> &IndexMap<String, IndexMap<u32, PileupCounter>> {
        &self.pileup
    }

    pub fn len(&self) -> usize {
        self.pileup.iter().map(|(_, v)| v.len()).sum()
    }

    pub fn sort(&mut self) {
        self.pileup
            .par_sort_by(|k1, _v1, k2, _v2| cmp_chromosome(k1, k2));
        for (_, pileup) in self.pileup.iter_mut() {
            pileup.par_sort_keys();
        }
    }

    fn add_matches(
        &mut self,
        umi: &str,
        chromosome: &str,
        start: u32,
        sequence: &Sequence,
        strand: u8,
        cursor_pos: u32,
        cursor_seq: usize,
        op_length: u32,
        base_qualities: &[u8],
        min_base_quality: u16,
    ) -> Result<(u32, usize, u32), UmiVarCalError> {
        for position in start + cursor_pos..start + cursor_pos + op_length {
            if cursor_seq >= sequence.len() {
                break;
            }
            let base = sequence.get(cursor_seq).unwrap();
            let base_quality = base_qualities[cursor_seq] - 33;
            if base_quality as u16 >= min_base_quality {
                //Access pileup[chromosome][position]
                self.pileup
                    .get_mut(chromosome)
                    .ok_or(PileupError::ChromosomeNotFound(chromosome.to_string()))?
                    .get_mut(&position)
                    .ok_or(PileupError::PositionNotFound(position))?
                    .add_nucleotide(&base, umi, base_quality, strand);
            }
        }
        Ok((start + cursor_pos + op_length - 1, cursor_seq, cursor_pos))
    }

    fn add_insertions(
        &mut self,
        umi: &str,
        chromosome: &str,
        position: u32,
        sequence: &Sequence,
        strand: u8,
        cursor_pos: u32,
        cursor_seq: usize,
        op_length: u32,
        base_qualities: &[u8],
        min_base_quality: u16,
    ) -> Result<(u32, usize, u32), UmiVarCalError> {
        //Check if chromosome and position+1 exist in pileup
        let mut cursor_pos = cursor_pos;
        let mut cursor_seq = cursor_seq;
        if self
            .pileup
            .get(chromosome)
            .map(|x| x.get(&(position + 1)))
            .flatten()
            .is_some()
        {
            let mut inserted_sequence: Vec<u8> = Vec::new();
            let mut inserted_qscore = 0f32;

            for _ in 0..op_length {
                let base = sequence.get(cursor_seq).unwrap();
                let base_quality = base_qualities[cursor_seq] - 33;
                inserted_sequence.push(base);
                inserted_qscore += base_quality as f32;

                cursor_seq += 1;
            }
            inserted_qscore /= op_length as f32;

            if inserted_qscore as u16 >= min_base_quality {
                self.pileup
                    .get_mut(chromosome)
                    .ok_or(PileupError::ChromosomeNotFound(chromosome.to_string()))?
                    .get_mut(&(position + 1))
                    .ok_or(PileupError::PositionNotFound(position + 1))?
                    .add_insertion(
                        &String::from_utf8(inserted_sequence)?,
                        umi,
                        inserted_qscore as u8,
                        strand,
                    );
                cursor_pos += 1;
                cursor_seq += 1;
            }
        } else {
            cursor_seq += op_length as usize;
        }
        Ok((position, cursor_seq, cursor_pos))
    }

    fn add_deletions(
        &mut self,
        umi: &str,
        chromosome: &str,
        start: u32,
        strand: u8,
        cursor_pos: u32,
        cursor_seq: usize,
        op_length: u32,
    ) -> Result<(u32, usize, u32), UmiVarCalError> {
        let mut del_cursor = 0;
        let mut cursor_pos = cursor_pos;
        for position in start + cursor_pos..start + cursor_pos + op_length {
            if self
                .pileup
                .get(chromosome)
                .map(|x| x.get(&position))
                .flatten()
                .is_some()
            {
                self.pileup
                    .get_mut(chromosome)
                    .ok_or(PileupError::ChromosomeNotFound(chromosome.to_string()))?
                    .get_mut(&position)
                    .ok_or(PileupError::PositionNotFound(position))?
                    .add_deletion(op_length - del_cursor, umi, strand);
                del_cursor += 1;
            }
        }
        cursor_pos += op_length;
        Ok((start + cursor_pos + op_length - 1, cursor_seq, cursor_pos))
    }

    pub fn add_read(
        &mut self,
        umi: &str,
        strand: u8,
        chromosome: &str,
        start: u32,
        cigar: &record::Cigar,
        sequence: &Sequence,
        base_qualities: &[u8],
        min_base_quality: u16,
    ) -> Result<(), UmiVarCalError> {
        let mut cursor_pos = 0;
        let mut cursor_seq = 0;
        let mut position = start;
        let mut start = start;

        for (op_type, op_length) in cigar
            .iter()
            .filter_map(|op| op.ok())
            .map(|op| (op.kind(), op.len() as u32))
        {
            match op_type {
                cigar::op::Kind::Match => {
                    (position, cursor_seq, cursor_pos) = self.add_matches(
                        umi,
                        chromosome,
                        start,
                        sequence,
                        strand,
                        cursor_pos,
                        cursor_seq,
                        op_length,
                        base_qualities,
                        min_base_quality,
                    )?;
                }
                cigar::op::Kind::Insertion => {
                    (position, cursor_seq, cursor_pos) = self.add_insertions(
                        umi,
                        chromosome,
                        position,
                        sequence,
                        strand,
                        cursor_pos,
                        cursor_seq,
                        op_length,
                        base_qualities,
                        min_base_quality,
                    )?;
                }
                cigar::op::Kind::Deletion => {
                    (position, cursor_seq, cursor_pos) = self.add_deletions(
                        umi, chromosome, start, strand, cursor_pos, cursor_seq, op_length,
                    )?;
                }
                cigar::op::Kind::SoftClip => {
                    if cursor_pos == 0 {
                        start -= op_length;
                        cursor_pos += op_length;
                        cursor_seq += op_length as usize;
                    }
                }
                _ => {}
            }
        }
        Ok(())
    }
}
