from typing import Dict, Tuple
from utils import canonical_form, reverse_complement

# REVERSE_COMPLEMENTS = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

# class KmerOccurrenceTable:
#     """
#     Keeps a record of kmers observed within reads.
#     For each kmer, we keep track of the id of the read where it was first encountered and its
#     position within the read.
#     """
#     _kmerOccurrences: Dict[str, Tuple[str, int]]
    
#     def __init__(self):
#         self._kmerOccurrences = {}
    
#     def _add_kmer(self, kmer, readId, coordinate):
#         self._kmerOccurrences.setdefault(kmer, (readId, coordinate))
#         self._add_reverse(kmer, readId, coordinate)
    
#     def _add_reverse(self, kmer, readId, coordinate):
#         reverse = "".join([REVERSE_COMPLEMENTS[base] for base in kmer])
#         self._kmerOccurrences.setdefault(reverse, (readId, coordinate))

def kmer_occurrences(reads, k):
    """Records all k-mer occurrences {'kmer': [(read_id, position)]}"""
    kmers = {}
    reverse_kmers = {}

    for read_name, seq in reads.items():
        for i in range(len(seq) - k + 1):
            canonical_kmer = canonical_form(seq[i:i+k])
            rev_kmer = reverse_complement(canonical_kmer) # what is rev_kmer for?

            # kmers.setdefault(canonical_kmer, (read_name, i))
            # reverse_kmers.setdefault(rev_kmer, (read_name, i))

            if canonical_kmer in kmers:
                kmers[canonical_kmer].append((read_name, i))
            else:
                kmers[canonical_kmer] = [(read_name, i)]
    
    return kmers, reverse_kmers

def kmer_to_id(kmers):
    """Returns khash: {id, 'kmer'}"""
    khash = {}
    count = 0
    for kmer in kmers:
        khash[kmer] = count
        count += 1

    return khash