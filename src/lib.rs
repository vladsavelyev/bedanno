use std::io;
use bio::io::bed;
use bio::io::gff;

pub fn overlap(qry: &bed::Record, trg: &gff::Record) -> bool {
    let qs = qry.start();
    let qe = qry.end();
    let ts = *trg.start();
    let te = *trg.end();
    return ts <= qs && qs < te || qs <= ts && ts <= qe;
}

pub fn run(
    input: impl io::Read,
    reference: impl io::Read,
    output: impl io::Write,
) -> io::Result<()> {
    // let mut reader = BufReader::new(reader);
    let mut reader = bed::Reader::new(input);
    let mut writer = bed::Writer::new(output);
    let mut reference = gff::Reader::new(
        reference, gff::GffType::GTF2,
    );

    let mut trg = reference.records().next();
    'bed: while let Some(q) = reader.records().next() {
        let mut q = q.unwrap();
        while let Some(Ok(t)) = trg {
            if *t.end() >= q.start() {
                if overlap(&q, &t) {
                    if let Some(gene) = t.attributes().get("gene_name") {
                        q.set_name(gene);
                    }
                } else {
                    q.set_name(".");
                }
                trg = Some(Ok(t));
                break;
            }
            trg = reference.records().next();
        }
        writer.write(&q)?;
        continue 'bed;  // continue to next region
    }
    Ok(())
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_annotate() {
        let input = vec![
            "chr1\t100\t500\tname1\t0.5\n",
            "chr1\t1000\t1500\tname2\t1.0\n",
            "chr1\t4000\t5000\tname3\t2.0\n",
            "chr1\t6000\t7000\tname4\t1.0\n",
        ];
        let reference = vec![
            "chr1\thavana\tgene\t1\t1200\t.\t+\t.\tgene_id \"ENSG1\"; gene_name \"GENE1\";\n",
            "chr1\thavana\tgene\t2000\t3000\t.\t+\t.\tgene_id \"ENSG2\"; gene_name \"GENE2\";\n",
            "chr1\thavana\tgene\t6500\t6600\t.\t+\t.\tgene_id \"ENSG3\"; gene_name \"GENE3\";\n",
        ];
        let expected = vec![
            "chr1\t100\t500\tGENE1\t0.5\n",
            "chr1\t1000\t1500\tGENE1\t1.0\n",
            "chr1\t4000\t5000\t.\t2.0\n",
            "chr1\t6000\t7000\tGENE3\t1.0\n",
        ];
        let input = input.into_iter().collect::<String>();
        let reference = reference.into_iter().collect::<String>();
        let expected = expected.into_iter().collect::<String>();
        let mut output = vec![];

        run(
            Box::new(input.as_bytes()),
            Box::new(reference.as_bytes()),
            Box::new(&mut output),
        ).expect("Cannot annotate BED file");

        let output = String::from_utf8(output).unwrap();
        assert_eq!(&expected.trim(), &output.trim())
    }
}
