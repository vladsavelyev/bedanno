use anyhow;
use bio::io::bed;
use sorted_list::SortedList;
use std::collections::{HashMap, HashSet};
use std::io::{self, BufRead};

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

#[derive(Clone, PartialEq)]
struct GffLine {
    contig: String,
    interval: Interval,
    feature_type: String,
    annotation: String,
}

#[derive(Clone)]
struct Annotation {
    interval: Interval,
    gene_name: Option<String>,
    feature_type: String,
    mane: bool,
    tsl: String,
    level: u8,
    transcript_type: String,
}

impl Annotation {
    fn from_gff_line(line: &GffLine) -> Self {
        let mut attributes: HashMap<String, String> = HashMap::new();
        for attr in line
            .annotation
            .split(';')
            .filter(|x| !x.is_empty())
        {
            let mut kv = attr
                .trim_matches(' ')
                .split(|x| x == '=' || x == ' ');
            let key = kv.next().expect(attr).to_string();
            let value = kv.next().expect(attr).trim_matches('"').to_string();
            attributes.insert(key, value);
        }
        let mane = attributes
            .get("tag")
            .map_or(false, |x| x == "MANE_Select");
        let tsl = attributes
            .get("transcript_support_level")
            .map_or("NA".to_string(), |x| x.to_string());
        let level = attributes
            .get("level")
            .map_or(255, |x| x.parse::<u8>().unwrap_or(255));
        let transcript_type = attributes
            .get("transcript_type")
            .map_or("".to_string(), |x| x.to_string());
        Self {
            interval: line.interval.clone(),
            gene_name: attributes.get("gene_name").cloned(),
            feature_type: line.feature_type.clone(),
            mane,
            tsl,
            level,
            transcript_type,
        }
    }
}

impl GffLine {
    fn from_line(line: &str) -> Self {
        let tokens: Vec<&str> = line.split("\t").collect::<Vec<&str>>();
        GffLine {
            contig: tokens[0].to_string(),
            interval: Interval::new(
                tokens[3].parse::<u64>().unwrap() - 1,
                tokens[4].parse::<u64>().unwrap(),
            ),
            feature_type: tokens[2].to_string(),
            annotation: tokens[8].to_string(),
        }
    }
}

fn find_overlaps<'a>(
    queries: &SortedList<Interval, ()>,
    targets: &'a SortedList<Interval, GffLine>,
) -> Vec<Vec<&'a GffLine>> {
    let mut result: Vec<Vec<&GffLine>> = Vec::new();

    let mut targets = targets.values();
    let mut prev_overlaps: Vec<&GffLine> = vec![];

    let mut trg: Option<&'a GffLine> = None;
    for q in queries.keys() {
        // skip targets that are before the query
        while let Some(rec) = trg.or_else(|| targets.next()) {
            if rec.interval.end >= q.start {
                trg = Some(rec);
                break;
            }
            trg = None;
        }
        // take overlapping targets
        let mut new_overlaps = vec![];
        while let Some(rec) = trg.or_else(|| targets.next()) {
            trg = None;
            if q.overlapping(&rec.interval) {
                new_overlaps.push(rec);
            } else {
                trg = Some(rec);
                break;
            }
        }

        result.push(
            prev_overlaps
                .iter()
                .skip_while(|rec| !q.overlapping(&rec.interval))
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
    gfflines: &'a mut Vec<Vec<&'a GffLine>>,
) -> Vec<Option<Annotation>> {
    let feature_type_rank = [
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

    let tsl_rank = ["1", "2", "3", "4", "5", "NA"]
        .iter()
        .enumerate()
        .map(|(i, x)| (x.to_owned().to_owned(), i))
        .collect::<HashMap<String, usize>>();

    queries
        .keys()
        .enumerate()
        .map(|(i, q)| {
            let mut annotations: Vec<Annotation> = gfflines[i]
                .iter()
                .map(|&t| Annotation::from_gff_line(t))
                .collect();
            annotations.sort_by_key(|anno| {
                // chr1    HAVANA  exon    12010   12057   .       +       .       gene_id "ENSG00000223972.6"; transcript_id "ENST00000450305.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-201"; exon_number 1; exon_id "ENSE00001948541.1"; level 2; transcript_support_level "NA"; hgnc_id "HGNC:37102"; ont "PGO:0000005"; ont "PGO:0000019"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000002844.2";
                let feature_type = feature_type_rank
                    .get(&anno.feature_type)
                    .unwrap_or(&255);
                let mane: u8 = if anno.mane { 0 } else { 255 };
                let tsl = tsl_rank.get(&anno.tsl).unwrap_or(&255);
                let level = anno.level;
                let coding: u8 = (anno.transcript_type != "protein_coding") as u8;
                assert!(dbg!(&q).overlapping(dbg!(&anno.interval)));
                let overlap = if q.start < anno.interval.start {
                    q.end - anno.interval.start
                } else {
                    anno.interval.end - q.start
                };
                (feature_type, mane, tsl, level, coding, overlap)
            });
            annotations.first().cloned()
        })
        .collect::<Vec<Option<Annotation>>>()
}

pub fn run(
    qreader: impl io::Read,
    treader: impl io::Read,
    writer: impl io::Write,
) -> anyhow::Result<()> {
    // let mut treader = gff::Reader::new(treader, gff_type);

    let mut qreader = bed::Reader::new(qreader);
    let mut qreader = qreader
        .records()
        .enumerate()
        .map(|(i, x)| (i, x.expect(&format!("Parsing BED line {}", i))));

    let mut treader = io::BufReader::new(treader)
        .lines()
        .map(|x| x.expect("Reading GTF line"))
        .skip_while(|x| x.starts_with('#'))
        .enumerate();

    let mut writer = bed::Writer::new(writer);

    let mut seen_qry_contigs: HashSet<String> = HashSet::new();
    let mut next_qry_contig = "".to_string();
    // Buffering into this HashMap until we meet the contig we are looking for:
    let mut targets_buffer: HashMap<String, SortedList<Interval, GffLine>> = HashMap::new();
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
        if cur_queries.len() == 0 {
            return Ok(());
        }
        eprintln!("Processing contig, {}", cur_contig);

        // Now that we collected queries for qry_contig, looking for corresponding targets
        let cur_targets = if targets_buffer.contains_key(&cur_contig) {
            // If we already recorded current contig before, taking it from the buffer:
            targets_buffer.remove(&cur_contig).unwrap()
        } else {
            // Else searching in the file:
            let mut res: SortedList<Interval, GffLine> = SortedList::new();
            while let Some((i, line)) = buf_trg.clone().or_else(|| treader.next()) {
                if i % 100_000 == 0 {
                    eprintln!("GTF line number {i}");
                }
                buf_trg = None;

                let rec = GffLine::from_line(&line);

                if rec.contig != cur_contig && res.len() > 0 {
                    buf_trg = Some((i, line));
                    break; // finished collecting contig data
                }
                if rec.contig != cur_contig {
                    // Buffering
                    targets_buffer
                        .entry(rec.contig.to_owned())
                        .or_insert(SortedList::new())
                } else {
                    &mut res
                }
                .insert(rec.interval.clone(), rec);
            }
            res
        };

        let mut annotations = find_overlaps(&mut cur_queries, &cur_targets);
        let annotations = resolve_all_overlaps(&cur_queries, &mut annotations);
        for (q, anno) in cur_queries.keys().zip(annotations) {
            let mut rec = bed::Record::new();
            rec.set_chrom(&cur_contig);
            rec.set_start(q.start);
            rec.set_end(q.end);
            let gene = anno
                .map_or(".".to_string(), |x| x.gene_name.unwrap_or(".".to_string()))
                .to_string();
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

    // #[test]
    fn test1() {
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
        )
        .expect("Cannot annotate BED file");

        let output = String::from_utf8(output).unwrap();
        assert_eq!(&expected.trim(), &output.trim())
    }

    #[test]
    fn test2() {
        let queries = to_str(&["chr1	11144641	11144761", "chr1	11144947	11145067"]);
        let targets = &[
            "chr1\tHAVANA\texon\t11144667\t11144886\t.\t+\t.\tgene_id \"ENSG00000225602.5\"; transcript_id \"ENST00000445982.5\"; gene_type \"lncRNA\"; gene_name \"MTOR-AS1\"; transcript_type \"lncRNA\"; transcript_name \"MTOR-AS1-202\"; exon_number 2; exon_id \"ENSE00001736483.1\"; level 2; transcript_support_level \"5\"; hgnc_id \"HGNC:40242\"; tag \"dotter_confirmed\"; tag \"basic\"; tag \"Ensembl_canonical\"; havana_gene \"OTTHUMG00000002003.1\"; havana_transcript \"OTTHUMT00000005565.1\";"
        ].join("\n");
        let expected = to_str(&[
            "chr1	11144641	11144761	MTOR-AS1",
            "chr1	11144947	11145067	.",
        ]);
        let mut output = vec![];

        run(
            Box::new(queries.as_bytes()),
            Box::new(targets.as_bytes()),
            Box::new(&mut output),
        )
        .expect("Cannot annotate BED file");

        let output = String::from_utf8(output).unwrap();
        assert_eq!(&expected.trim(), &output.trim())
    }
}
