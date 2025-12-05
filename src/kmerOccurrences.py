from typing import Dict, Tuple
from utils import canonical_form, reverse_complement
from collections import defaultdict

def kmer_occurrences(reads, k):
    # kmers = set()

    # for i in range(len(reads)):
    #     seq = reads[i]

    #     # initial window
    #     kmer = seq[:k]
    #     for j in range(k - 1, len(seq)):
    #         # slide window to include new nucleotide if not initial window
    #         if j >= k: kmer = kmer[1:] + seq[j]

    #         canonical_kmer, dir = canonical_form(kmer)
    #         # kmer_table[canonical_kmer].append((i, dir, j - k + 1)) 
    #         kmers.add(canonical_kmer)

    # return list(kmers)
    kmer_table = defaultdict(list)

    for i in range(len(reads)):
        seq = reads[i]

        # initial window
        kmer = seq[:k]
        for j in range(k - 1, len(seq)):
            # slide window to include new nucleotide if not initial window
            if j >= k: kmer = kmer[1:] + seq[j]

            canonical_kmer, dir = canonical_form(kmer)
            kmer_table[canonical_kmer].append((i, dir, j - k + 1)) 

    return kmer_table