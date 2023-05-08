use bedanno;
use bio::io::gff;
use flate2::read::GzDecoder;
use std::env;
use std::fs;
use std::io;

fn main() {
    let mut args = env::args();
    args.next();

    let query: Box<dyn io::Read> = match args.next() {
        Some(path) => {
            let reader = fs::File::open(&path).unwrap();
            if path.ends_with(".gz") {
                Box::new(GzDecoder::new(reader))
            } else {
                Box::new(reader)
            }
        }
        None => Box::new(io::stdin()),
    };

    let mut gff_path: String = args.next().unwrap_or("hg38".to_string());
    if gff_path == "hg38" {
        gff_path = "data/hg38/gencode.v43.basic.annotation.gtf.gz".to_string();
    }
    let gff_type = {
        let p = gff_path.trim_end_matches(".gz").to_owned();
        if p.ends_with(".gff") || p.ends_with(".gff3") {
            gff::GffType::GFF3
        } else if p.ends_with(".gff2") {
            gff::GffType::GFF2
        } else if p.ends_with(".gtf") {
            gff::GffType::GTF2
        } else {
            panic!("Reference must be a GFF or GTF file, or genome name (hg38 is supported)")
        }
    };

    let reader = fs::File::open(&gff_path).expect("Cannot open GTF/GFF file");
    let target: Box<dyn io::Read> = if gff_path.ends_with(".gz") {
        Box::new(GzDecoder::new(reader))
    } else {
        Box::new(reader)
    };

    let output = io::stdout();

    bedanno::run(query, target, output, gff_type).expect("Should be able to annotate BED file")
}
