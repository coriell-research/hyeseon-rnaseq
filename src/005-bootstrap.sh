#!/bin/bash
#
# Create boostrap replicates of Hyseon fastq files
#
# -----------------------------------------------------------------------------
FASTQ_DIR=/home/gcalendo/data/projects/hyeSeon-rnaseq/data/fastq/uncompressed
SAMPLE_NAMES=/home/gcalendo/data/projects/hyeSeon-rnaseq/data/fastq/uncompressed/sample-names.txt
BOOT_EXE=/home/gcalendo/data/projects/hyeSeon-rnaseq/src/004-bootstrap-fastq.py
OUT_DIR=/home/gcalendo/data/projects/hyeSeon-rnaseq/data/fastq/uncompressed/boot
REPLICATES=5
SEED=123


for samp in $(cat $SAMPLE_NAMES); 
    do python $BOOT_EXE $FASTQ_DIR/${samp}_R1.fastq $FASTQ_DIR/${samp}_R2.fastq \
                        --replicates $REPLICATES \
                        --outDir $OUT_DIR \
                        --seed $SEED; 
done
