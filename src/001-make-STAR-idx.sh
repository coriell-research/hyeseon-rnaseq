#!/bin/bash
#
# create STAR index for GRCm38 genome
#
# Old STAR index did not have splice junctions
# To run newest version of STAR:
#  conda activate bioinfo
# -----------------------------------------------------------------------------
GENOME_FASTA=/mnt/data/gdata/mouse/mm38_mm10/Mus_musculus.GRCm38.dna.primary_assembly.fa
ANNOTATION_GTF=/mnt/data/gdata/mouse/mm38_mm10/Mus_musculus.GRCm38.100.gtf
THREADS=12

STAR --runThreadN $THREADS \
     --runMode genomeGenerate \
     --genomeDir /home/gcalendo/data/projects/hyeSeon-rnaseq/data/GRCm38_STAR_idx \
     --genomeFastaFiles $GENOME_FASTA \
     --sjdbGTFfile $ANNOTATION_GTF
