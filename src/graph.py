from utils import canonical_form

# TODO: convert this into a graph class

def create_pre_nodes(reads, kmer_table, hash_length):
    pre_nodes = []

    for i in range(len(reads)):
        seq = reads[i]
        start = 0
        
        # get initial kmer
        new_kmer = seq[start:hash_length]

        # initialize sequence of uninterrupted kmers
        consecutive_seq = seq[start:hash_length - 1]

        for end in range(hash_length - 1, len(seq)):
            # if not in initial kmer, slide kmer window
            if end >= hash_length: new_kmer = new_kmer[1:] + seq[end]
            
            can_kmer, dir = canonical_form(new_kmer)
            first_occurrence = kmer_table[can_kmer][0]

            # if newly added kmer to the window overlaps with other reads
            if len(kmer_table[can_kmer]) > 1:
                # create a new node for the previously uninterrupted sequence of kmers
                if (end - 1) - start + 1 >= hash_length:
                    pre_nodes.append(consecutive_seq)
                
                # create a node for the overlapping kmer
                if (i, dir, end - hash_length + 1) == first_occurrence:
                    pre_nodes.append(new_kmer)
                
                # slide window
                start = end - hash_length + 2
                consecutive_seq = consecutive_seq[1:]
            
            consecutive_seq += seq[end]
            
            # if we have reached the end of the read and there is no overlap create a node
            # for this rightmost run of uninterrupted kmers
            if end == len(seq) - 1 and len(kmer_table[can_kmer]) == 1:
                pre_nodes.append(consecutive_seq)
    return pre_nodes