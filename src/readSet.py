import gzip
import sys

"""
Functions associated with parsing input files.
"""

def read_fastq(lines):
    i = 0
    seqs = []
    while i < len(lines) - 3:
        header = lines[0].strip()
        if not header: break

        seq = lines[i+1].strip()
        seqs.append(seq)
        i += 4

    return seqs

def parse_and_read_file(filename, filetype="FASTA"):
    if filetype in ["FASTA_GZ", "FASTQ_GZ"]:
        with gzip.open(filename, "rt") as f:
            content = f.readlines()
        filetype = filetype.rstrip("_GZ")
    elif filetype in ["FASTA", "FASTQ"]:
        with open(filename) as f:
            content = f.readlines()
    else:
        print(f"Filetype {filetype} is not yet supported.")
        sys.exit(0)

    # verifying file is in correct format
    if filetype in ["FASTA_GZ", "FASTA"]:
        if len(content) and len(content[0]) and content[0][0] != '>':
            print(f"{filename} does not seem to be in FASTA format.")
            sys.exit(0)
    elif filetype in ["FASTQ_GZ", "FASTQ"]:
        if len(content) and len(content[0]) and content[0][0] != '@':
            print(f"{filename} does not seem to be in FASTQ format.")
            sys.exit(0)
        
        return read_fastq(content)
    
    # parsing read ids and their corresponding sequence
    names, seqs = [], []
    for line in content:
        line = line.strip()

        if len(line) and line[0] in ['>', '@']:
            names.append(line[1:])
        elif len(line):
            seqs.append(line)
    
    # return dict(zip(names, seqs))
    return seqs

def output_contigs(filename, contigs):
    with open(f"{filename}.txt", "w") as f:
        lines = """"""

        for contig in contigs:
            line = f"{contig}\n"
            lines += line

        f.write(lines)
    
    return