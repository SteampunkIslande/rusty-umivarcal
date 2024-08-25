use indexmap::IndexMap;
use noodles::bam::record::Cigar;
use std::collections::HashSet;
use std::fs::File;

use crate::commons::BedRecord;
use crate::err::UmiVarCalError;

use std::cmp;

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
    qscore: Option<u16>,
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
}

impl PileupCounter {
    // pub fn add_nucleotide(&mut self, nuc: &u8, umi: &str, qscore: u16) {
    //     match nuc {
    //         b'A' => {
    //             self.a.add_forward(umi, qscore);
    //         }
    //         b'C' => {
    //             self.c.add_forward(umi, qscore);
    //         }
    //         b'G' => {
    //             self.g.add_forward(umi, qscore);
    //         }
    //         b'T' => {
    //             self.t.add_forward(umi, qscore);
    //         }
    //         _ => (),
    //     }
    // }
    pub fn set_reference(&mut self, nuc: &u8) {
        self.reference = *nuc;
    }

    pub fn set_homopolymer(&mut self, hp: u16) {
        self.hp = hp;
    }
}

impl NucleotideCounter {
    // pub fn add_forward(&mut self, umi: &str, qscore: u16) {
    //     self.forward += 1;
    //     self.umis.insert(umi.to_string());
    //     self.qscore = Some(qscore + self.qscore.unwrap_or(0));
    // }
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

    /// Load pileup from serialized RMP file.
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
        self.pileup.len()
    }

    pub fn sort(&mut self) {
        self.pileup
            .par_sort_by(|k1, _v1, k2, _v2| cmp_chromosome(k1, k2));
        for (_, pileup) in self.pileup.iter_mut() {
            pileup.par_sort_keys();
        }
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
