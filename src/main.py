import sys
from utils import settings
from readSet import parse_and_read_file
from kmerOccurrenceTable import kmer_occurrences

# getting command-line arguments
argv = sys.argv
output_dir, filepath, hash_length, file_format, read_type = settings(argv[1:])

# 1. Read input file (focus on FASTA/FASQ files)
reads = parse_and_read_file(filepath, file_format)

# 2. Build kmer hash table
    # For each k-mer observed in the set of reads, the hash table records the ID of the first read encountered containing that k-mer and the position of its occurrence within that read. 
    # Each k-mer is recorded simultaneously to its reverse complement. 
kmers, reverse_kmers = kmer_occurrences(reads, hash_length)

print(kmers, reverse_kmers) 
# 3. Build roadmaps
    # rewrite each read as a set of original k-mers combined with overlaps with previously hashed reads
# 4. Build second database
    # A second database is created with the opposite information. It records, for each read, which of its original k-mers are overlapped by subsequent reads.
# 5. Build graph
# 6. Simplify the graph
# 7. Write contigs and graph stats to output file specified?