import sys

argv = sys.argv
k = argv[1]

# Rough outline
# 1. Read input file (focus on FASTA/FASQ files)
# 2. Build kmer hash table
    # For each k-mer observed in the set of reads, the hash table records the ID of the first read encountered containing that k-mer and the position of its occurrence within that read. 
    # Each k-mer is recorded simultaneously to its reverse complement. 
# 3. Build roadmaps
    # rewrite each read as a set of original k-mers combined with overlaps with previously hashed reads
# 4. Build second database
    # A second database is created with the opposite information. It records, for each read, which of its original k-mers are overlapped by subsequent reads.
# 5. Build graph
# 6. Simplify the graph
# 7. Write contigs and graph stats to output file specified?