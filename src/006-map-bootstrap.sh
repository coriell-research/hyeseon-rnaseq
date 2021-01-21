#!/bin/bash
#
# Map boostrapped reads to genome
#
# conda activate bioinfo for STAR
# -----------------------------------------------------------------------------
GENOME_DIR=/mnt/data/gdata/mouse/mm38_mm10/GRCm38_STAR_idx
RESULTS_DIR=/home/gcalendo/data/projects/hyeSeon-rnaseq/results/boot
SAMPLE_LIST=/home/gcalendo/data/projects/hyeSeon-rnaseq/data/fastq/uncompressed/sample-list.txt
FASTQ_DIR=/home/gcalendo/data/projects/hyeSeon-rnaseq/data/fastq/uncompressed/boot
THREADS=12
REPLICATES=5


for samp in $(cat $SAMPLE_LIST);
do
    for rep in $(seq 1 $REPLICATES);
    do
        STAR --runThreadN $THREADS \
             --genomeDir $GENOME_DIR \
             --readFilesIn $FASTQ_DIR/${samp}_R1.boot.${rep}.fastq $FASTQ_DIR/${samp}_R2.boot.${rep}.fastq \
             --outFilterType BySJout \
             --outFilterMultimapNmax 20 \
             --alignSJoverhangMin 8 \
             --alignSJDBoverhangMin 1 \
             --outFilterMismatchNmax 999 \
             --outFilterMismatchNoverReadLmax 0.04 \
             --alignIntronMin 20 \
             --alignIntronMax 1000000 \
             --alignMatesGapMax 1000000 \
             --outFileNamePrefix $RESULTS_DIR/${samp}.boot.${rep}. \
             --outMultimapperOrder Random \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode GeneCounts;
    done
done

