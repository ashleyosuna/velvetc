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
hash_length = 3
kmer_table = kmer_occurrences(reads, hash_length)

# 3. Build graph
pre_graph = graph.Graph(hash_length=hash_length)

graph.create_pre_nodes(reads, kmer_table, hash_length, pre_graph)

graph.concatenate_nodes(pre_graph)

print(pre_graph.nodes, pre_graph.starts)

# no in edges or out edges?
for n in pre_graph.nodes:
    print(pre_graph.nodes[n].in_edges, pre_graph.nodes[n], pre_graph.nodes[n].out_edges)

# hardcoded example
# node_1 = graph.Node(
#     seq = "ATC",
#     id = 1,
#     in_edges =
#     out_edges = 
# )
# node_2 = graph.Node(
#     seq="TCA",
#     id=2,
#     in_edges=
#     out_edges=
# )
# node_3 = graph.Node(
#     seq="CAA",
#     id=3,
#     in_edges=
#     out_edges=
# )
# node_4 = graph.Node(
#     seq="AAT",
#     id=4,
#     in_edges=
#     out_edges=
# )
# node_5 = graph.Node(
#     seq="TCC",
#     id=5,
#     in_edges=
#     out_edges=
# )
# node_6 = graph.Node(
#     seq="CCA",
#     id=6,
#     in_edges=
#     out_edges=
# )
# node_7 = graph.Node(
#     seq="CAA",
#     id=7,
#     in_edges=
#     out_edges=
# )
# hardcoded_graph = graph.Graph()


# 4. Remove bubbles using Tour Bus
# for node_id in pre_graph.nodes:
#     node = pre_graph.nodes[node_id]
#     print(node_id, node.out_edges)
# print(pre_graph.nodes)

print('\n\n')
unbubbled = tourbus(pre_graph)
print(unbubbled)
# for n in pre_graph.nodes:
#     print(n[0].out_edges, n[0].in_edges)
#     print(n[1].out_edges, n[1].in_edges, '\n\n')

# graph.clip_tips(pre_graph)

# print(pre_graph.nodes, pre_graph.starts)

# 7. Write contigs and graph stats to output file specified?