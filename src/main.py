import sys
from utils import settings, canonical_form
from readSet import parse_and_read_file
from kmerOccurrenceTable import kmer_occurrences
import graph
from tourbus import tourbus

# getting command-line arguments
argv = sys.argv
output_dir, filepath, hash_length, file_format, read_type = settings(argv[1:])

pre_graph = graph.Graph(hash_length=hash_length)

# 1. Read input file (focus on FASTA/FASQ files)
reads = parse_and_read_file(filepath, file_format)

# 2. Build kmer hash table
    # For each k-mer observed in the set of reads, the hash table records the ID of the first read encountered containing that k-mer and the position of its occurrence within that read. 
    # Each k-mer is recorded simultaneously to its reverse complement. 
hash_length = 3
# print(pre_graph.enable_twin) # for debugging
kmer_table = kmer_occurrences(reads, hash_length, pre_graph.enable_twin)

# 3. Build graph
# pre_graph = graph.Graph(hash_length=hash_length)
graph.create_pre_nodes(reads, kmer_table, hash_length, pre_graph)
graph.concatenate_nodes(pre_graph)

print(pre_graph.nodes, pre_graph.starts)

# for n in pre_graph.nodes:
#     print(pre_graph.nodes[n].in_edges, pre_graph.nodes[n], pre_graph.nodes[n].out_edges)

# 4. Remove bubbles using Tour Bus
unbubbled = tourbus(pre_graph)
print(unbubbled)
# for n in pre_graph.nodes:
#     print(n[0].out_edges, n[0].in_edges)
#     print(n[1].out_edges, n[1].in_edges, '\n\n')

# graph.clip_tips(pre_graph)

# print(pre_graph.nodes, pre_graph.starts)

# 7. Write contigs and graph stats to output file specified?