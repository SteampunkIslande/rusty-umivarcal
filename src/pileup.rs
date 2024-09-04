use core::str;
use indexmap::IndexMap;
use std::collections::HashSet;
use std::fs::File;

use crate::commons::BedRecord;
use crate::err::{PileupError, UmiVarCalError};

use regex::Regex;

use std::cmp;

struct CIGAROperation {
    length: u32,
    operation: char,
}

#[derive(Debug, Default, serde::Deserialize, serde::Serialize)]
pub struct InsertionDict {
    insertions: IndexMap<String, NucleotideCounter>,
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
    reference: u8,

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

    // Mean quality score at this position
    qscore: f32,

    // Standard deviation of quality scores at this position
    qscore_stdev: f32,
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
            _ => (),
        }
    }
    pub fn add_insertion(&mut self, seq: &str, umi: &str, qscore: u8, strand: u8) {
        self.insertions.insert(seq, umi, qscore, strand);
    }
    pub fn add_deletion(&mut self, position: u32, umi: &str, strand: u8) {
        self.deletions.add_deletion(position, umi, strand);
    }
    pub fn set_reference(&mut self, nuc: &u8) {
        self.reference = *nuc;
    }

    pub fn set_homopolymer(&mut self, hp: u16) {
        self.hp = hp;
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
    pub fn save_pileup_chromosome(
        &self,
        output_prefix: &str,
        chrom: &str,
    ) -> Result<(), UmiVarCalError> {
        let pileup = self.pileup.get(chrom).unwrap();
        let file = File::create(format!("{}.{}.pileup", output_prefix, chrom)).unwrap();
        let mut wr = std::io::BufWriter::new(&file);

        let mut obj = IndexMap::new();
        obj.insert(&chrom, pileup);
        rmp_serde::encode::write(&mut wr, &obj)?;
        Ok(())
    }

    /// Save pileup for a range of positions in a chromosome. If no start and end are provided, the entire chromosome is saved.
    /// The output file is named as follows: `output_prefix.chrom.start.end.pileup`.
    pub fn save_pileup_range(
        &self,
        output_prefix: &str,
        chrom: &str,
        start: Option<u32>,
        end: Option<u32>,
    ) {
        //Get pileup for the range, as a IndexMap of chromosome => position => PileupCounter.

        let pileup_range: IndexMap<u32, &PileupCounter> = self
            .pileup
            .get(chrom)
            .map(|x| {
                x.iter()
                    .filter(|(k, _)| match (start, end) {
                        (Some(start), Some(end)) => k >= &&start && k < &&end,
                        (Some(_), None) => {
                            todo!("End position not provided! Either provide both or none!")
                        }
                        (None, Some(_)) => {
                            todo!("Starts position not provided! Either provide both or none!")
                        }
                        (None, None) => true,
                    })
                    .map(|(k, v)| (*k, v))
                    .collect()
            })
            .unwrap();

        let start = start.unwrap_or(
            *pileup_range
                .keys()
                .min()
                .expect("Could not get minimum position!"),
        );
        let end = end.unwrap_or(
            *pileup_range
                .keys()
                .max()
                .expect("Could not get maximum position!"),
        );

        let mut obj = IndexMap::new();
        obj.insert(chrom, &pileup_range);

        let file = File::create(format!(
            "{}.{}.{}.{}.pileup",
            output_prefix, chrom, start, end
        ))
        .unwrap();

        let mut wr = std::io::BufWriter::new(&file);
        rmp_serde::encode::write(&mut wr, &pileup_range).unwrap();
    }

    /// Save pileup for a range of positions in a bed file.
    pub fn save_pileup_bed(
        &self,
        output_prefix: &str,
        bed: &std::path::Path,
    ) -> Result<(), UmiVarCalError> {
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
            let chrom = bed_record.chrom;
            let start = bed_record.start;
            let end = bed_record.end;
            let file_name = format!("{}.{}.{}.{}.pileup", output_prefix, chrom, start, end);
            let file = File::create(&file_name).expect("Could not create pileup file!");
            let mut wr = std::io::BufWriter::new(&file);
            let pileup_range: Option<IndexMap<&u32, &PileupCounter>> =
                self.pileup.get(&chrom).map(|x| {
                    x.iter()
                        .filter(|(k, _)| k >= &&start && *k < &&end)
                        .collect()
                });
            match pileup_range {
                None => {
                    //remove file
                    std::fs::remove_file(&file_name).expect("Could not remove pileup file!");
                    continue;
                }
                Some(pileup_range) => {
                    let mut obj = IndexMap::new();
                    obj.insert(&chrom, pileup_range);
                    rmp_serde::encode::write(&mut wr, &obj)?
                }
            }
        }
        Ok(())
    }

    /// Create new, empty pileup from bed file.
    /// Prepares all the necessary data structures from the bed file.
    pub fn pileup_from_bed(bed: &std::path::Path) -> Self {
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

        Self { pileup }
    }

    /// Load pileup from serialized messagepack file.
    pub fn load_pileup(pileup: &std::path::Path) -> Self {
        let file = File::open(pileup).expect("Could not open pileup file!");
        let mut rd = std::io::BufReader::new(&file);
        let pileup: IndexMap<String, IndexMap<u32, PileupCounter>> =
            rmp_serde::decode::from_read(&mut rd).expect("Could not decode pileup file!");
        Self { pileup }
    }

    pub fn pileup(&mut self) -> &mut IndexMap<String, IndexMap<u32, PileupCounter>> {
        &mut self.pileup
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
        sequence: &[u8],
        strand: u8,
        cursor_pos: u32,
        cursor_seq: usize,
        op_length: u32,
        base_qualities: &[u8],
        min_base_quality: u16,
    ) -> Result<(u32, usize, u32), UmiVarCalError> {
        for position in start + cursor_pos..start + cursor_pos + op_length {
            let base = sequence[cursor_seq];
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
        sequence: &[u8],
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
                let base = sequence[cursor_seq];
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
        cigar: &[u8],
        sequence: &[u8],
        base_qualities: &[u8],
        min_base_quality: u16,
    ) -> Result<(), UmiVarCalError> {
        let re = Regex::new(r"(\d+)([MSIDHN])").unwrap();
        let mut operations = Vec::new();
        for cap in
            re.captures_iter(str::from_utf8(cigar).expect("Could not convert cigar to string"))
        {
            let length = cap[1].parse::<u32>().unwrap(); //Safe to unwrap because the regex guarantees that the capture is a number.
            let operation = cap[2].chars().next().unwrap(); //Safe to unwrap because the regex guarantees that the capture is a single character.
            operations.push(CIGAROperation { length, operation });
        }
        // Check if operations and sequence length match
        if operations
            .iter()
            .map(|op| op.length as usize)
            .sum::<usize>()
            != sequence.len()
        {
            panic!("Invalid CIGAR string: CIGAR and sequence length do not match!");
        }
        //Remove N operations
        operations.retain(|op| op.operation != 'N');

        let mut cursor_pos = 0;
        let mut cursor_seq = 0;
        let mut position = start;
        let mut start = start;

        for (op_length, op_type) in operations.iter().map(|op| (op.length, op.operation)) {
            match op_type {
                'M' => {
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
                'I' => {
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
                'D' => {
                    (position, cursor_seq, cursor_pos) = self.add_deletions(
                        umi, chromosome, start, strand, cursor_pos, cursor_seq, op_length,
                    )?;
                }
                'S' => todo!(),
                _ => {}
            }
        }
        Ok(())
    }
}
