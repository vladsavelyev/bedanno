use anyhow;
use bio::io::bed;
use bio::io::gff;
use sorted_list::SortedList;
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::io;

#[derive(Clone, Eq, PartialEq, PartialOrd, Ord, Debug)]
struct Interval {
    start: u64,
    end: u64,
}

impl Interval {
    fn new(start: u64, end: u64) -> Self {
        Self { start, end }
    }

    fn overlapping(&self, other: &Self) -> bool {
        return other.start <= self.start && self.start < other.end
            || self.start <= other.start && other.start <= self.end;
    }
}

impl fmt::Display for Interval {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}, {})", self.start, self.end)
    }
}

fn find_annotations<'a>(
    queries: &SortedList<Interval, ()>,
    targets: &'a SortedList<Interval, gff::Record>,
) -> Vec<Vec<&'a gff::Record>> {
    //
    let mut result: Vec<Vec<&gff::Record>> = Vec::new();

    let mut targets = targets.values();
    let mut prev_overlaps: Vec<&gff::Record> = vec![];

    let mut trg: Option<&'a gff::Record> = None;
    for q in queries.keys() {
        // skip targets that are before the query
        while let Some(rec) = trg.or_else(|| targets.next()) {
            let t = Interval::new(rec.start() - 1, *rec.end());
            if t.end >= q.start {
                trg = Some(rec);
                break;
            }
            trg = None;
        }
        // take overlapping targets
        let mut new_overlaps = vec![];
        while let Some(rec) = trg.or_else(|| targets.next()) {
            trg = None;
            let t = Interval::new(rec.start() - 1, *rec.end());
            if t.overlapping(&q) {
                new_overlaps.push(rec);
            } else {
                trg = Some(rec);
                break;
            }
        }

        result.push(
            prev_overlaps
                .iter()
                .skip_while(|rec| {
                    let t = Interval::new(rec.start() - 1, *rec.end());
                    !q.overlapping(&t)
                })
                .chain(new_overlaps.iter())
                .map(|rec| *rec)
                .collect(),
        );
        prev_overlaps = new_overlaps;
    }
    result
}

fn resolve_all_overlaps<'a>(
    queries: &'a SortedList<Interval, ()>,
    annotations: &'a mut Vec<Vec<&'a gff::Record>>,
) -> Vec<Option<&'a gff::Record>> {
    let feature_type_map = [
        "CDS",
        "stop_codon",
        "start_codon",
        "UTR",
        "exon",
        "transcript",
        "gene",
    ]
    .iter()
    .enumerate()
    .map(|(i, x)| (x.to_owned().to_owned(), i))
    .collect::<HashMap<String, usize>>();

    let tsl_map = ["1", "2", "3", "4", "5", "NA"]
        .iter()
        .enumerate()
        .map(|(i, x)| (x.to_owned().to_owned(), i))
        .collect::<HashMap<String, usize>>();

    queries
        .keys()
        .enumerate()
        .map(|(i, q)| {
            annotations[i].sort_by_key(|t| {
                // chr1    HAVANA  exon    12010   12057   .       +       .       gene_id "ENSG00000223972.6"; transcript_id "ENST00000450305.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-201"; exon_number 1; exon_id "ENSE00001948541.1"; level 2; transcript_support_level "NA"; hgnc_id "HGNC:37102"; ont "PGO:0000005"; ont "PGO:0000019"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000002844.2";
                let feature_type = feature_type_map
                    .get(t.feature_type())
                    .unwrap_or(&255);
                let mane: u8 = t
                    .attributes()
                    .get("tag")
                    .map(|x| if x == "MANE_Select" { 0 } else { 1 })
                    .unwrap_or(255);
                let tsl = t
                    .attributes()
                    .get("transcript_support_level")
                    .map(|x| tsl_map.get(x).unwrap())
                    .unwrap_or(&255);
                let level = t
                    .attributes()
                    .get("level")
                    .map(|x| x.parse::<u8>().unwrap_or(255))
                    .unwrap_or(255);
                let coding: u8 = t
                    .attributes()
                    .get("transcript_type")
                    .map(|x| if x == "protein_coding" { 0 } else { 1 })
                    .unwrap_or(255);
                let overlap = if q.start < *t.start() - 1 {
                    q.end - *t.start() - 1
                } else {
                    *t.end() - q.start
                };
                (feature_type, mane, tsl, level, coding, overlap)
            });
            annotations[i].first().map(|x| *x)
        })
        .collect()
}

pub fn run(
    qreader: impl io::Read,
    treader: impl io::Read,
    writer: impl io::Write,
    gff_type: gff::GffType,
) -> anyhow::Result<()> {
    let mut qreader = bed::Reader::new(qreader);
    let mut treader = gff::Reader::new(treader, gff_type);
    let mut writer = bed::Writer::new(writer);

    let mut qreader = qreader
        .records()
        .enumerate()
        .map(|(i, x)| (i, x.expect(&format!("Parsing BED line {}", i))));

    let mut treader = treader
        .records()
        .enumerate()
        .map(|(i, x)| (i, x.expect(&format!("Parsing GTF line {}", i))));

    let mut seen_qry_contigs: HashSet<String> = HashSet::new();
    let mut next_qry_contig = "".to_string();
    // Buffering into this HashMap until we meet the contig we are looking for:
    let mut targets_buffer: HashMap<String, SortedList<Interval, gff::Record>> = HashMap::new();
    // Buffering one last BED and GTF record so we can do 2 parallel loops over BED and GTF:
    let mut buf_qry = None;
    let mut buf_trg = None;
    loop {
        // Contigs loop
        let mut cur_contig = next_qry_contig.clone();
        let mut cur_queries: SortedList<Interval, ()> = SortedList::new();

        while let Some((i, rec)) = buf_qry.clone().or_else(|| qreader.next()) {
            buf_qry = None;
            if cur_contig == "" {
                cur_contig = rec.chrom().to_owned();
            }
            if rec.chrom() == &cur_contig {
                let interval = Interval::new(rec.start(), rec.end());
                cur_queries.insert(interval, ());
            } else {
                // New contig
                if seen_qry_contigs.contains(rec.chrom()) {
                    return Err(anyhow::anyhow!(
                        "Parsing BED record {i}: contig {} is out of order",
                        rec.chrom()
                    ));
                }
                seen_qry_contigs.insert(cur_contig.clone());
                next_qry_contig = rec.chrom().to_owned();
                buf_qry = Some((i, rec.clone()));
                break;
            }
        }
        if cur_contig == "" {
            return Err(anyhow::anyhow!("No records found in BED file"));
        }
        if cur_queries.len() == 0 {
            return Ok(());
        }

        // Now that we collected queries for qry_contig, looking for corresponding targets
        let cur_targets = if targets_buffer.contains_key(&cur_contig) {
            // If we already recorded current contig before, taking it from the buffer:
            targets_buffer.remove(&cur_contig).unwrap()
        } else {
            // Else searching in the file:
            let mut res: SortedList<Interval, gff::Record> = SortedList::new();
            while let Some((i, rec)) = buf_trg.clone().or_else(|| treader.next()) {
                buf_trg = None;
                if rec.seqname() != cur_contig && res.len() > 0 {
                    buf_trg = Some((i, rec));
                    break; // finished collecting contig data
                }
                let interval = Interval::new(*rec.start() - 1, *rec.end());
                if rec.seqname() != &cur_contig {
                    // Buffering
                    targets_buffer
                        .entry(rec.seqname().to_owned())
                        .or_insert(SortedList::new())
                } else {
                    &mut res
                }
                .insert(interval, rec.clone());
            }
            res
        };

        let mut annotations = find_annotations(&mut cur_queries, &cur_targets);
        let annotations = resolve_all_overlaps(&cur_queries, &mut annotations);
        for (q, anno) in cur_queries.keys().zip(annotations) {
            let mut rec = bed::Record::new();
            rec.set_chrom(&cur_contig);
            rec.set_start(q.start);
            rec.set_end(q.end);
            let gene = anno
                .and_then(|x| x.attributes().get("gene_name"))
                .unwrap_or(&".".to_owned())
                .to_owned();
            rec.set_name(&gene);
            writer.write(&rec)?;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn to_str(lines: &[&str]) -> String {
        lines
            .into_iter()
            .map(|x| {
                x.split_whitespace()
                    .collect::<Vec<&str>>()
                    .join("\t")
            })
            .collect::<Vec<String>>()
            .join("\n")
    }

    #[test]
    fn test_annotate() {
        let queries = to_str(&[
            "chr1  10  50",
            "chr1  100 150",
            "chr1  400 500",
            "chr1  600 700",
            "chr10 10  50",
            "chr10 55  60",
            "chr2  5   55",
            "chrX  100 200",
            "chrM  0   200",
        ]);
        let targets = to_str(&[
            "chr1  havana gene 1    5    . + . gene_name=SKIP1;",
            "chr1  havana gene 21   60   . + . gene_name=LOW1;",
            "chr1  havana CDS  31   40   . + . gene_name=KEEP1;",
            "chr1  havana CDS  91   170  . + . gene_name=KEEP2;",
            "chr1  havana gene 201  300  . + . gene_name=SKIP2;",
            "chr1  havana gene 601  700  . + . gene_name=LOW4;level=2;",
            "chr1  havana gene 551  650  . + . gene_name=KEEP4;level=1;",
            "chr1  havana gene 801  900  . + . gene_name=SKIP5;",
            "chr2  havana gene 1    500  . + . gene_name=KEEPchr2;",
            "chr10 havana gene 1    500  . + . gene_name=LOWchr10;transcript_support_level=2;",
            "chr10 havana gene 1    500  . + . gene_name=KEEPchr10;transcript_support_level=1;",
            "chrX  havana gene 1    500  . + . gene_name=KEEPchrX;",
        ]);
        let expected = to_str(&[
            "chr1  10  50  KEEP1",
            "chr1  100 150 KEEP2",
            "chr1  400 500 .",
            "chr1  600 700 KEEP4",
            "chr10 10  50  KEEPchr10",
            "chr10 55  60  KEEPchr10",
            "chr2  5   55  KEEPchr2",
            "chrX  100 200 KEEPchrX",
            "chrM  0   200 .",
        ]);
        let mut output = vec![];

        run(
            Box::new(queries.as_bytes()),
            Box::new(targets.as_bytes()),
            Box::new(&mut output),
            gff::GffType::GFF3,
        )
        .expect("Cannot annotate BED file");

        let output = String::from_utf8(output).unwrap();
        assert_eq!(&expected.trim(), &output.trim())
    }
}
