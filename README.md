# BedAnno

Assigns gene names to regions in a BED file.

```sh
cat regions.bed
chr1   10   50
chr1   100  150
chr1   400  500
chr1   600  700
chr10  10   50
chr10  55   60
chr2   5    55
chrX   100  200
chrM   0    200
```

```sh
bedanno regions.bed [gencode.gff.gz] > regions.anno.bed
```

```sh
cat regions.anno.bed
chr1   10   50   GENEA
chr1   100  150  GENEB
chr1   400  500  .
chr1   600  700  GENED
chr10  10   50   GENE10
chr10  55   60   GENE10
chr2   5    55   GENE2
chrX   100  200  KEEPX
chrM   0    200  .
```

When there are multiple genes overlapping a region, the gene with the highest priority will be used. The priority is
defined as follows:

1. Feature type (`CDS` > `stop_codon` > `start_codon` > `UTR` > `exon` > `transcript` > `gene` > everything else)
2. Whether the annotation is selected in `MANE` as the primary transcript.
3. Transcript support level (`1` (all splice junctions are supported by at least one non-suspect mRNA) > `2` (best
   supporting mRNA is flagged as suspect or the support is from multiple ESTs) > `3` (only support is from a single
   EST) > `4` (best supporting EST is flagged as suspect) > `5` (no single transcript supports the model
   structure) > `NA` (transcript that was not analysed for TSL).
4. Confidence level (`1` (feature is validated) > `2` (manual annotation) > `3` (automated annotation)
5. Whether transcript type is `protein_coding`.
6. Overlap % between GTF and BED region.

The tool works fastest when chromosomes are sorted and in the same order between BED and GTF, however doesn't require
it.

Regions within chromosomes also don't need to be sorted, however the tool would benefit if they are.

Note, however, that chromosomes can't be intermixed, e.g. if `chr2` is met after `chr1`, meeting `chr1` again will cause
an error.

If GTF/GFF not provided, assumes that regions are hg38, and uses a built-in GTF
file `gencode.v43.basic.annotation.gtf.gz`.
