import sys

def reverse_complement(seq: str):
    copy = seq[::-1]
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join([complements[base] for base in copy])

def canonical_form(seq: str):
    rev_seq = reverse_complement(seq)
    # return min(seq, rev_seq)
    return ((seq, 1) if seq < rev_seq else (rev_seq, -1))

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