import sys
from utils import settings, canonical_form
from readSet import parse_and_read_file
from kmerOccurrenceTable import kmer_occurrences
import graph
from tourbus import tourbus

# getting command-line arguments
argv = sys.argv
output_dir, filepath, hash_length, file_format, read_type = settings(argv[1:])

# 1. Read input file (focus on FASTA/FASQ files)
reads = parse_and_read_file(filepath, file_format)

# 2. Build kmer hash table
    # For each k-mer observed in the set of reads, the hash table records the ID of the first read encountered containing that k-mer and the position of its occurrence within that read. 
    # Each k-mer is recorded simultaneously to its reverse complement. 
hash_length = 5
kmer_table = kmer_occurrences(reads, hash_length)

# 3. Build graph
pre_graph = graph.Graph()

graph.create_pre_nodes(reads, kmer_table, hash_length, pre_graph)

graph.concatenate_nodes(pre_graph)

# print(pre_graph.nodes, pre_graph.starts)

for n in pre_graph.nodes:
    print(n)
    # print(n[0].out_edges, n[0].in_edges) # does the source node not have any in edges or out edges??
#     print(n[1].out_edges, n[1].in_edges, '\n\n')

# 4. Remove bubbles using Tour Bus
unbubbled = tourbus(pre_graph, pre_graph.nodes[1])
print(unbubbled)

# 7. Write contigs and graph stats to output file specified?