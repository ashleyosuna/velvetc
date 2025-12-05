# BabyVelvet: A Reimplementation of Velvet in the Codon Programming Language

## ABSTRACT
Velvet is a set of algorithms to create and manipulate de Bruijn graphs for the purpose of genomic sequence assembly. Velvet is able to produce significantly long and useful contiguous sequences even with very short reads. The original implementation of Velvet is written in C. Here we attempt to re-implement Velvet in Codon, a high performance Python implementation that compiles to native machine code, to evaluate whether Codon can match or improve on the efficiency achieved through Velvet’s C implementation.

## INTRODUCTION
The field of genomics depends on accurate and efficient genome sequence assembly, for everything from sequencing RNA molecules to determining the genome of a new species. Genome sequencing is done through the piecing together of DNA fragments randomly extracted from the sample, and assembling them into a set of contiguous sequences, known as contigs. Longer contigs provide better knowledge of the DNA sequence, and are more informative for analysis.

At the time of Velvet’s publishing in 2007, recent sequencing technologies like Solexa were producing shorter reads (35bp) than previous technologies (Bentley, 2006). These shorter reads were not well suited for the widely-used overlap-layout-consensus approach for sequence assembly, which represents each read as a node (Batzoglou 2005). This is because shorter reads correspond to more nodes, increasing the time and space complexity of the resulting graph.

Alternatively, Velvet relies on a de Bruijn graph, with nodes organized around k-mers. In this data representation, reads are mapped as paths through the graph. In this way, high redundancy is handled by the graph without an increase in the number of nodes. These de Bruijn graphs require efficient and robust methods for eliminating errors and resolving repeats, which is what Velvet’s algorithms do.

Codon is a domain-extensible compiler and DSL framework with Python’s syntax and semantics (Shajii, 2023). High-level languages like Python are growing more common in many domains, including scientific computing, due to their ease and flexibility. Like a high-level programming language, Codon provides the benefits of readability and an easier learning curve, while achieving performance equal to or better than low level languages such as C or C++. For more information, see: https://docs.exaloop.io/.

## MATERIALS AND METHODS

Graph construction begins with identification of k-mers. Each k-mer is extracted from the reads and stored in a hash table in its canonical form, the lexicographically smaller of itself and its reverse complement. A node is created for each of the unique canonical k-mers, along with a twin referencing the k-mer’s reverse complement.

Edges are then created by treading once more through each read sequence and making cuts whenever a k-mer is found to overlap (i.e., is found to appear in multiple reads). A node is created for the uninterrupted sequence of k-mers before the overlap, which achieves a reduction in nodes, making graph traversals more computationally efficient. If the overlapping k-mer is the first occurrence of this k-mer, a new pre-node is created. If one already exists, the older one is used. Every time a node is created or reused, an edge is created from the previous node to the current node. For every node, a twin-node is created for its reverse complement. Nodes are identified by their id, with twin nodes identified by the negative of that id.

The graph is then simplified by concatenating pre-nodes according to the number of edges. If two nodes, Node A and Node B, are connected by an edge, a check is done: if Node A has only one outgoing edge, and if Node B has only one incoming edge, the nodes (and their twins) are merged.

The graph is cleaned up by removing “tips,” dead-end paths that are shorter than 2k. They are considered to be the result of sequencing errors, and are iteratively removed. The concatenation process is run again on this cleaned up graph. 

Paths that start and end at the same node, known as bubbles, are identified through the TourBus algorithm, a Dijkstra-like breadth first search. The algorithm starts from a select few start nodes and finds the closest common ancestor of two paths that reach the same node. If the paths are determined to be similar enough (over 50% of bases match), they are merged according to a global alignment. 

Once more, errors are removed through tip clipping, and the graph is simplified once again through concatenation.

Contigs are identified as long nodes in the graph (of length >= k). From the set of contigs, we calculate the N50, a common statistical measure that describes the quality of an assembly by its continuity. In our implementation, we consider the length of a contig to be the number of k-mers collapsed into it, rather than its nucleotide sequence length. We used this definition to produce a metric as comparable as possible to Velvet’s, as we believe this is how Velvet measures contig length. However, we cannot guarantee that the two metrics are exactly identical, so small differences might remain.

Following this, Velvet resolves repeats using read-pair information through a module it calls ‘Breadcrumb.’ We chose not to implement this due to time constraints. We believe that BabyVelvet encompasses the main idea of Velvet by assembling short reads into long contigs; Breadcrumb is an added, but unnecessary in our case, bonus.

## RESULTS

In order to compare our implementation to Velvet, we ran several tests. Initially we tested on a small dataset of three randomly generated reads with length of 50bp. Testing on such a small dataset allowed us to debug during development, and ensure our mental model of the implementation was similar to Velvet’s. We ran this test with a hash length of 5. From this initial testing we obtained encouraging results: Velvet resulted in a final graph with 30 nodes, outputting 8 final contigs with n50 of 8 and maximum contig length of 17, while BabyVelvet resulted in a final graph with 35 nodes, outputting 8 final contigs with n50 of 9 and maximum contig length of 17. Furthermore, upon closer examination of the final contigs generated, BabyVelvet’s final contigs corresponded to either a contig reported by Velvet, or to its reverse complement.

We then decided to run both Velvet and BabyVelvet on larger datasets. To do this, we simulated 35bp-long reads from a 50-Kb region within the Escherichia coli genome with a coverage value of 3. We then ran both Velvet and BabyVelvet using a hash length of 21. We selected these values for the length of the reads and the hash length to test our implementation under similar conditions to those under which Velvet was tested in the original article [cite?]. Table 1 summarizes the results of this experiment. 

Prior to graph simplification, both assemblers produced graphs with a very large number of k-mers (with Velvet producing ~5k initial prenodes and our assembler ~37k initial nodes). The difference in the number of initial nodes is due to the fact that our assembler initially creates a node for each canonical k-mer before proceeding to concatenate, while Velvet creates pre-nodes for sequences of uninterrupted k-mers. However, they both demonstrated a significant reduction in this number of initial nodes after concatenating and pruning tips off the graph. For example, the graph produced by our assembler went from ~37k initial nodes down to 1416 nodes in this test.


Although both assemblers exhibit a similar behavior of a drastic graph compaction through concatenation and pruning, our assembler consistently results in a much larger number of nodes in the final graph. Figure 1 shows the final number of nodes in the graphs generated by Velvet and our assembler for different datasets with varying numbers of reads. This discrepancy is likely due to differences in graph-pruning heuristics and bubble-handling. In particular, Velvet uses much more sophisticated methods to aggressively merge redundant paths, while our assembler uses a much simpler heuristic of comparing nodes in similar paths locally to decide whether to merge them, resulting in a more conservative merging of paths. Furthermore, due to time constraints, we did not implement coverage-pruning, which Velvet uses to further simplify the graph. 

![Figure 1]([http://url/to/img.png](https://github.com/ashleyosuna/velvetc/main/images/figure1.png))
Figure 1. Number of nodes in the final graph for different numbers of reads in the input data.


Figure 2. N50 metric for the final contigs found in the final graphs for different numbers of reads in the input data.

Despite the difference in graph sizes, both Velvet and BabyVelvet exhibit similar trends in the reported N50 metrics for datasets of varying sizes. Figure 2 shows how the N50 metric changed as the numbers of reads contained in the datasets increased. For both, the N50 remained relatively stable even as the size of the dataset increased. Velvet consistently achieves an N50 hovering around 60, while BabyVelvet’s N50 values hover around 40. This is consistent with the fact that our assembler is a lot more conservative when merging paths. This result is still encouraging, as it shows that the continuity reported by our implementation scales similarly to that of Velvet’s. 

While the results of our preliminary experiments indicate that BabyVelvet behaves similarly to Velvet, it is significantly slower in practice. For example, when running both programs on the 100-Kb dataset, Velvet completed in 1 second while BabyVelvet finished in 21 seconds. This difference in performance may be partially explained by the use of higher-level data structures and our more conservative approach to merging paths. Additionally, C is known for its speed, and so some of this performance discrepancy may also be attributed to the use of a higher-level programming language. We are hopeful that with some performance optimizations, BabyVelvet could achieve improved runtimes and approximate Velvet’s performance.

## CHALLENGES

The main challenge while attempting to re-implement Velvet in Codon was due to the robustness and complexity of the source code. Since Velvet was originally written in C, which handles more complex structures, such as arrays and strings, in a much different way than languages such as Python and Codon do, the source code was hard to read, trace, and interpret. We had initially planned on porting the code, but the difficulty of reading the code made it essentially impossible under the time constraints to do so. While we did manage to understand the main flow of Velvet and frequently referenced the source code throughout the development of BabyVelvet, we eventually decided that it was more realistic to try to re-implement the core ideas of the algorithm described in the original paper as opposed to rewriting the code line-by-line.

Another issue we encountered had to do with testing. Since our implementation ended up being slower than Velvet, we could not run larger tests locally, and were limited to smaller test cases. Even with smaller datasets (those generated from the Escherichia coli genome), we were unable to run the tests locally, and so we were limited to testing using Github actions.

## CONCLUSION
This project introduces BabyVelvet, a re-implementation of the Velvet algorithms for short read de Bruijn graph assembly. Though we were able to replicate some of Velvet’s behaviour, we were not able to match its efficiency. However, there is significant potential to improve the efficiency of our algorithms and bring us closer to Velvet’s results. Our work demonstrates that Velvet’s algorithms are effective in producing long contigs from short reads, and that implementing them in Codon provides an opportunity to improve readability while maintaining functionality.

## REFERENCES
Batzoglou, S. 2005. Algorithmic challenges in mammalian genome sequence assembly. In Encyclopedia of genomics, proteomics and bioinformatics (eds. M. Dunn et al.), Part 4. John Wiley and Sons, New York.

Bentley, D.R. 2006. Whole-genome re-sequencing. Curr. Opin. Genet. Dev. 16: 545–552

Shajii, A., Ramirez, G., Smajlović, H., Ray, J., Berger, B., Amarasinghe, S., and Numanagić, I. Codon: A Compiler for High-Performance Pythonic Applications and DSLs. In Proceedings of the 32nd ACM SIGPLAN International Conference on Compiler Construction (CC ’23), February 25–26, 2023, Montréal, QC, Canada. ACM, New York, NY, USA, 12 pages. https://doi.org/10.1145/3578360.3580275

Zerbino DR, Birney E. Velvet: algorithms for de novo short read assembly using de Bruijn graphs. Genome Res. 2008 May;18(5):821-9. doi: 10.1101/gr.074492.107. Epub 2008 Mar 18. PMID: 18349386; PMCID: PMC2336801.

