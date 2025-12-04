import sys
from readSet import parse_and_read_file
from utils import settings
from kmerOccurrences import kmer_occurrences
from graph import Graph
from tourbus import tourbus
from utils import n50

# getting command-line arguments
argv = sys.argv
output_dir, filepath, hash_length, file_format, read_type = settings(argv[1:])

# Read input file (focus on FASTA/FASQ files)
reads = parse_and_read_file(filepath, file_format)

kmers = kmer_occurrences(reads, hash_length)

graph = Graph(hash_length)

graph.create_init_nodes(kmers)

graph.map_through_reads(reads)

graph.concatenate_nodes()

graph.clip_tips()

tourbus(graph)

graph.concatenate_nodes()

contigs = graph.get_contigs()

print(contigs)

print(n50(contigs, hash_length))

# output nodes to file