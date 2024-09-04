pub enum StrandBiasMethod {
    Default,
    Torrent,
}

pub enum UmiSource {
    ReadName,
    Tag(String),
    Length(usize),
}

#[derive(Debug, Default, serde::Deserialize, serde::Serialize)]
pub struct BedRecord {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
}

/// Check if a read is valid based on its base qualities. If the mean base quality is below the threshold, the read is discarded.
pub fn qual_is_valid(read_base_qualities: &[u8], min_read_quality: u16) -> bool {
    let threshold = min_read_quality * read_base_qualities.len() as u16;
    read_base_qualities
        .iter()
        .map(|q| (*q - 33) as u16)
        .sum::<u16>()
        >= threshold
}

pub fn read_is_valid(
    flag: u16,
    mapq: u8,
    min_mapping_quality: u8,
    read_base_qualities: &[u8],
    min_read_quality: u16,
) -> bool {
    let orphan = flag & 0x2 == 0x0;
    if orphan {
        eprintln!("Orphan read found!");
        return false;
    }
    if !qual_is_valid(read_base_qualities, min_read_quality) {
        eprintln!(
            "Read quality is below threshold! {:?} < {:?}",
            read_base_qualities, min_read_quality
        );
        return false;
    }
    if mapq < min_mapping_quality || mapq == 255 {
        eprintln!(
            "Mapping quality is below threshold! {:?} < {:?}",
            mapq, min_mapping_quality
        );
        return false;
    }
    true
}
