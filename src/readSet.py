import gzip

"""
Functions associated with parsing input files.
"""
def parse_and_read_file(filename, filetype="FASTA", double_strand = False, no_hash = False):
    if filetype in ["FASTA_GZ", "FASTQ_GZ"]:
        with gzip.open(filename, "rt") as f:
            content = f.readlines()
        filetype = filetype.rstrip("_GZ")
    elif filetype in ["FASTA", "FASTQ"]:
        with open(filename) as f:
            content = f.readlines()
    else:
        print(f"Filetype {filetype} is not yet supported.")
        exit(0)

    # verifying file is in correct format
    if filetype in ["FASTA_GZ", "FASTA"]:
        if len(content) and len(content[0]) and content[0][0] != '>':
            print(f"{filename} does not seem to be in FASTA format.")
            exit(0)
    elif filetype in ["FASTQ_GZ", "FASTQ"]:
        if len(content) and len(content[0]) and content[0][0] != '@':
            print(f"{filename} does not seem to be in FASTQ format.")
            exit(0)
    
    # parsing read ids and their corresponding sequence
    names, seqs = [], []
    for line in content:
        line = line.strip()

        if len(line) and line[0] in ['>', '@']:
            names.append(line[1:])
        elif len(line):
            seqs.append(line)

    return dict(zip(names, seqs))
