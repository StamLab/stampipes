# Description

This process implements our Alt-seq processing pipeline

## Usage:

TODO

## Developing

TODO

## Notes

We use STARsolo, which is like CellRanger. Designed for single-cell analysis, we use it by telling the software that the barcode internal to the pool is a "cell barcode." All references by STAR and output files to "cell barcodes" actually refer to individual library barcodes. That is, each "cell barcode" corresponds to a DSnumber and LNnumber.

### STARsolo output directories

We get 4 output directories for STARSolo, each with their own way of counting results.

`Gene/`
: counts reads that are fully concordant with gene transcript.

`GeneFull/`
: counts all reads that overlap gene loci, including exonic and intronic reads.

`GeneFull_Ex50pAS/`
: "Gene Full, Except 50% anti-sense." Excludes reads that map >50% in the antisense direction

`GeneFull_ExonOverIntron`
: "GeneFull, prefer exon over intron." (This is important, for instance, for a read that overlap gene A exons and gene B introns (i.e. gene A exons are located within gene B introns). With GeneFull option, such a read will be ambiguous and not counted. With GeneFull_ExonOverIntron it will be counted towards gene A.)

### CellReads.stats interpretation

The meaning of `CellReads.stats` columns is not obvious, but it is documented in a github issue [here](https://github.com/alexdobin/STAR/issues/1501).
For convenience, I've collected the answers below:

Column descriptions:

`CB`
: cell barcode

`cbMatch`
: number of reads that matched the cell barcode

`cbPerfect` 
: number of perfect match on cell barcode

`exonic` 
: number of reads mapping on exonic (only for `GeneFull_Ex50pAS` and `GeneFull_ExonOverIntron`)

`intronic` 
: number of reads mapping on intronic (only for `GeneFull_Ex50pAS` and `GeneFull_ExonOverIntron`)

`mito` 
: number of reads mapping on mitochondrial genome

`genomeU`
: number of reads mapping to one locus in the genome

`genomeM`
: number of reads mapping to multiple loci in the genome

`featureU`
: number of reads mapping to one feature (Gene, GeneFull, etc)

`featureM`
: number of reads mapping to multiple features

`cbMMunique`
: number of reads with cell barcodes that map with mismatches to one barcode in the passlist

`cbMMmultiple`
: number of reads with cell barcodes that map with mismatches to multiple barcodes in the passlist

`exonicAS`
: number of reads mapping antisense to annotated exons (only for `GeneFull_Ex50pAS`)

`intronicAS`
: number of reads mapping antisense to annotated introns (only for `GeneFull_Ex50pAS`)

`countedU`
: number of unique-gene reads that were used in counting UMIs (!= number of UMIs), i.e. reads with valid CB/UMI/gene

`countedM`
: number of multi-gene reads that were used in counting UMIs (!= number of UMIs), i.e. reads with valid CB/UMI/gene

`nUMIunique` 
: total number of counted UMI

`nGenesUnique` 
: number of genes having non 0 counts

`nUMImulti`
: number of UMI for multi-gene reads, if requested

`nGenesMulti`
: number of genes supported by just multi-gene reads, if requested

NB:
> All columns are read counts, except for CB and the last 4: nUMIunique nGenesUnique nUMImulti nGenesMulti
