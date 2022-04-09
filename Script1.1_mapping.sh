#!/bin/bash

############## 1.1. MAPPING
##Align reads and Add read groups – 5h (+ sam-->bam; index bam)

PATH_ASSEMBLY_SORTED="path/to/assembly/GRCh37.dna.primary_assembly.fa"
flowcell="1";
lane="1";
library="1";
population="Human";

sample=$1

PATH_FASTQ="/path/to/fastq_files"
PATH_OUT_BAM="/path/to/output_folder_bam_files"

bwa mem -M -R @RG\tID:${flowcell}_${lane}\tLB:${library}\tPL:ILLUMINA\tSM:${sample}\tPU:1234 $PATH_ASSEMBLY_SORTED <(zcat ${PATH_FASTQ}/${sample}_1_paired.fastq.gz) <(zcat ${PATH_FASTQ}/${sample}_2_paired.fastq.gz) | samtools view -Sbh - | java -Xmx8g -jar picard.jar SortSam I=/dev/stdin O=${PATH_OUT_BAM}/${sample}/${sample}.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=STRICT

samtools index ${PATH_OUT_BAM}/${sample}/${sample}.bam

##PICARD OPTIONS
# - QUIET: Whether to suppress job-summary info on System.err. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false}
# - VALIDATION_STRINGENCY: Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded. Default value: STRICT. This option can be set to 'null' to clear the default value. Possible values: {STRICT, LENIENT, SILENT} 
# - COMPRESSION_LEVEL:	Compression level for all compressed files created (e.g. BAM and GELI). Default value: 5. This option can be set to 'null' to clear the default value.
# - SORT_ORDER: Optional sort order to output in. If not supplied OUTPUT is in the same order as INPUT. Default value: null. Possible values: {unsorted, queryname, coordinate, duplicate, unknown}