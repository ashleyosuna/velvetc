REF=ecoli.fasta
REGION_LEN=5000000
COVERAGE=50
READLEN=100
OUT_PREFIX=ecoli_50x

samtools faidx $REF

art_illumina \
  -ss HS25 \
  -sam \
  -i $REF \
  -l $READLEN \
  -f $COVERAGE \
  -o ${OUT_PREFIX} \
  -s 12345   # seed for reproducibility