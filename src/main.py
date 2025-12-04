import sys
from readSet import parse_and_read_file, output_contigs
from utils import settings
from kmerOccurrences import kmer_occurrences
from graph import Graph
from tourbus import tourbus
from utils import n50

# getting command-line arguments
argv = sys.argv
output_file, filepath, hash_length, file_format, read_type = settings(argv[1:])

# Read input file (focus on FASTA/FASQ files)
reads = parse_and_read_file(filepath, file_format)

kmers = kmer_occurrences(reads, hash_length)

print(len(kmers), "kmers found")

graph = Graph(hash_length)

graph.create_init_nodes(kmers)

print("Initial nodes", graph.node_count)

graph.map_through_reads(reads)

graph.concatenate_nodes()

print("After initial concatenation", graph.node_count)

graph.clip_tips()

print("After clipping tips", graph.node_count)

tourbus(graph)

print("After running tourbus", graph.node_count)

graph.concatenate_nodes()

print("After concatenating", graph.node_count)

contigs = graph.get_contigs()

print(contigs)

print(n50(contigs, hash_length))

# output nodes to file
output_contigs(output_file, contigs)