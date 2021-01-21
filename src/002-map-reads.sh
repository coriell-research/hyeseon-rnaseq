#!/bin/bash
#
# Map reads to reference genome using star index created in previous step
#
# conda activate bioinfo for STAR
#
# -----------------------------------------------------------------------------
GENOME_DIR=/mnt/data/gdata/mouse/mm38_mm10/GRCm38_STAR_idx
RESULTS_DIR=/home/gcalendo/data/projects/hyeSeon-rnaseq/results
FILE_LIST=/home/gcalendo/data/projects/hyeSeon-rnaseq/data/fastq/file-list.txt
THREADS=12


while read -r samp fqFiles; do
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_DIR \
     --readFilesIn $fqFiles \
     --outFilterType BySJout \
     --outFilterMultimapNmax 20 \
     --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 \
     --outFilterMismatchNoverReadLmax 0.04 \
     --alignIntronMin 20 \
     --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000 \
     --outFileNamePrefix $RESULTS_DIR/${samp}_star/${samp}. \
     --outMultimapperOrder Random \
     --outSAMtype BAM SortedByCoordinate \
     --readFilesCommand zcat \
     --quantMode GeneCounts;
done < $FILE_LIST
