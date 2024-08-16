pub enum StrandBiasMethod {
    Default,
    Torrent,
}

#[derive(Debug, Default, serde::Deserialize, serde::Serialize)]
pub struct BedRecord {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
}
