use std::collections::{HashMap, HashSet};
use std::fs::File;

use crate::commons::BedRecord;
use crate::err::UmiVarCalError;

#[derive(Debug, Default, serde::Deserialize, serde::Serialize)]
pub struct InsertionDict {
    insertions: HashMap<String, NucleotideCounter>,
}

#[derive(Debug, Default, serde::Deserialize, serde::Serialize)]
pub struct DeletionDict {
    deletions: HashMap<u32, NucleotideCounter>,
}

#[derive(Debug, Default, serde::Deserialize, serde::Serialize)]
pub struct NucleotideCounter {
    forward: u32,
    reverse: u32,
    umis: HashSet<String>,
    qscore: Option<u16>,
}

#[derive(Debug, Default, serde::Deserialize, serde::Serialize)]
pub struct PileupCounter {
    a: NucleotideCounter,
    c: NucleotideCounter,
    g: NucleotideCounter,
    t: NucleotideCounter,
    insertions: InsertionDict,
    deletions: DeletionDict,
    base_error_probability: f32,
}

#[derive(Debug, Default, serde::Deserialize, serde::Serialize)]
pub struct Pileup {
    pileup: HashMap<String, HashMap<u32, PileupCounter>>,
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

        let mut obj = HashMap::new();
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
        //Get pileup for the range, as a HashMap of chromosome => position => PileupCounter.

        let pileup_range: HashMap<u32, &PileupCounter> = self
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

        let mut obj = HashMap::new();
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
            let pileup_range: Option<HashMap<&u32, &PileupCounter>> =
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
                    let mut obj = HashMap::new();
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
        let mut pileup = HashMap::new();
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
            let chrom_entry = pileup.entry(bed_record.chrom).or_insert(HashMap::new());
            for i in bed_record.start..=bed_record.end {
                chrom_entry.insert(i, PileupCounter::default());
            }
        }

        Self { pileup }
    }

    /// Load pileup from serialized RMP file.
    pub fn load_pileup(pileup: &std::path::Path) -> Self {
        let file = File::open(pileup).expect("Could not open pileup file!");
        let mut rd = std::io::BufReader::new(&file);
        let pileup: HashMap<String, HashMap<u32, PileupCounter>> =
            rmp_serde::decode::from_read(&mut rd).expect("Could not decode pileup file!");
        Self { pileup }
    }

    pub fn len(&self) -> usize {
        self.pileup.len()
    }

    pub fn add_read(
        &mut self,
        umi: &str,
        strand: u8,
        chromosome: &str,
        position: u32,
        cigar: &[u8],
        sequence: &[u8],
        base_qualities: &[u8],
        min_base_quality: u16,
    ) {
    }
}
