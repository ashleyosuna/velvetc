# BabyVelvet

This project introduces BabyVelvet, a lightweight re-implementation of the Velvet algorithm for short read de Bruijn graph assembly in the Codon programming language. Codon is a domain-extensible compiler and DSL framework with Python’s syntax and semantics. Codon provides the benefits of readability, while achieving performance similar to low level languages such as C or C++. The goal of this project was to implement the core ideas of the original Velvet program [1].

BabyVelvet constructs a deBruijn graph, performs concatenation of nodes, clipping short tips, and bubble resolution (inspired by Velvet's Tour Bus bubble-resolution method), and outputs the assembled contigs.

## Usage

BabyVelvet can be run using the following command:
codon run [-release] main.py <output_prefix> <input_file> filetype=<FASTA|FASTQ> k=<hash_length>

For example, running "codon run main.py test_output ecoli_3x_75k_region.fq filetype=FASTQ k=21" will generate contigs using 21-mers from the provided FAST1 dataset and write the assembled contigs to the file test_output.txt

### Notes

- BabyVelvet currently supports only FASTA and FASTQ files.

[1] Zerbino DR, Birney E. Velvet: algorithms for de novo short read assembly using de Bruijn graphs. Genome Res. 2008 May;18(5):821-9. doi: 10.1101/gr.074492.107. Epub 2008 Mar 18. PMID: 18349386; PMCID: PMC2336801.
