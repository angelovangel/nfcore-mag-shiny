# nfcore-mag-shiny
A Shiny GUI for the [nf-core/mag](https://nf-co.re/mag) pipeline. The frontend collects the required parameters and
executes the nf-core/mag nextflow pipeline in the background.

The pipeline is used for assembly, binning and annotation of metagenomes. It supports **both short and long** reads, quality trims the reads and adapters with `fastp` and `Porechop` and performs basic QC with `FastQC`.

The pipeline then:

- assigns taxonomy to reads using `kraken2`
- performs assembly using `megahit` and `spades`
- checks their quality using `quast`
- performs metagenome binning using `metabat`
- checks the quality of the genome bins with `busco`

For usage and how to interpret the results
please see the [nf-core/mag documentation](https://nf-co.re/mag/docs/output)