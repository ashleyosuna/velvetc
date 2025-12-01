import sys
import numpy as np

def reverse_complement(seq: str):
    copy = seq[::-1]
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join([complements[base] for base in copy])

def canonical_form(seq: str):
    rev_seq = reverse_complement(seq)
    # return min(seq, rev_seq)
    return ((seq, 1) if seq < rev_seq else (rev_seq, -1))

def compare_sequences(u, v):
    MAX_GAPS = 3
    MAX_DIVERGENCE = 0.2
    INDEL = 0
    MATCH = 1

    u = u.upper()
    v = v.upper()
    
    # CONSTRUCT MATRIX
    n, m = len(u) + 1, len(v) + 1
    network = np.full((n, m), 0, dtype=float)

    # fill first row
    for i in range(n): network[i][0] = INDEL * i

    # fill first column
    for j in range(m): network[0][j] = INDEL * j

    # fill matrix row by row
    for i in range(1, n):
        for j in range(1, m):
            network[i][j] = max(
                network[i-1][j] + INDEL,
                network[i][j-1] + INDEL,
                network[i-1][j-1] + (MATCH if u[i-1] == v[j-1] else -MATCH)
            )

    # BACKTRACK TO CONSTRUCT ALIGNMENT
    i, j = n - 1, m - 1

    while i > 0 or j > 0:
        # vertical gap, 'move' along u
        if network[i][j] == network[i-1][j] + INDEL:
            i -= 1
        # horizontal gap, 'move' along v
        elif network[i][j] == network[i][j-1] + INDEL:
            j -= 1
        # if match or mismatch, 'move' along both sequences
        else:
            i -= 1
            j -= 1
        
    max_score = network[n-1][m-1]
    max_len = max(u, v)

    if max_score < max_len - MAX_GAPS: return False
    if (1 - max_score / max_len) > MAX_DIVERGENCE: return False

    return True

def settings(args):
    if len(args) < 2:
        print("The following arguments must be provided: output_dir and filepath.")
        sys.exit(0)
    
    output_dir = args[0]
    filepath = args[1]

    # setting default options
    k = 31
    filetype = "FASTA"
    read_type = "short"

    for option in args[2:]:
        key, val = option.split("=")
        if key == "k": k = int(val)
        elif key == "filetype": filetype = val
        elif key == "read_type": read_type = val
        else:
            print(f"Option {option} is not recognized.")
            sys.exit(0)
    
    # value for k should be odd, if even it is decremented
    if k % 2 == 0: k -= 1

    return output_dir, filepath, k, filetype, read_type