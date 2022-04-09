#!/bin/bash

############## 0.1. FASTQ QC
PATH_FASTQ="/path/to/fastq_files"

all_samples=$(ls $PATH_FASTQ | sed 's/_.*$//g' | sort | uniq) 
PATH_OUT_QC="0.1.FastQC" 

mkdir PATH_OUT_QC
for sample in $all_samples;
do 
	fastqc ${PATH_FASTQ}/${sample}_1.fastq.gz  ${PATH_FASTQ}/${sample}_2.fastq.gz -o $PATH_OUT_QC
done

############## 0.1a CONCATENATE FASTQ (when same sample in multiple lanes)
cd $PATH_FASTQ 
samples=$(ls *R1* | grep -v "JOINED" | sed 's/_L.*$//g' | uniq) 
for sample in $samples; do
	cat ${sample}_*R1* > ${sample}_R1_JOINED.fq.gz
	cat ${sample}_*R2* > ${sample}_R2_JOINED.fq.gz
done

##Check the cat output
mkdir checking;
for sample in $samples; do
	wc -l ${sample}_*R1_001* | tail -n 1 > checking/${sample}1.txt; wc -l ${sample}_*R1_JOINED* >> checking/${sample}1.txt
	wc -l ${sample}_*R2_001* | tail -n 1 > checking/${sample}2.txt; wc -l ${sample}_*R2_JOINED* >> checking/${sample}2.txt
done

############## 0.2. TRIMMING 
for sample in $all_samples; do
  java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -trimlog outs/${sample}.log -phred33 ${PATH_FASTQ}/${sample}_1.fastq.gz ${PATH_FASTQ}/${sample}_2.fastq.gz ${sample}_1_paired.fastq.gz ${sample}_1_upaired.fastq.gz ${sample}_2_paired.fastq.gz ${sample}_2_unpaired.fastq.gz ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:25 MINLEN:30
done

#--> Then, repeat fastQC

########################################################################################################
############## 1. FROM FASTQ to BAM (Mapping)
########################################################################################################

##############
############## 1.1. MAPPING
##############

PATH_ASSEMBLY_SORTED="path/to/assembly/GRCh37.dna.primary_assembly.fa"

#### Index
bwa index -a bwtsw $PATH_ASSEMBLY_SORTED
samtools faidx $PATH_ASSEMBLY_SORTED
java -jar picard.jar CreateSequenceDictionary REFERENCE=$PATH_ASSEMBLY_SORTED OUTPUT=${PATH_ASSEMBLY_SORTED}.dict

####Align reads and Add read groups – 5h (+ sam-->bam; sort bam; index bam)
PATH_FASTQ="/path/to/fastq_files"
PATH_OUT_BAM="/path/to/output_folder_bam_files"

all_samples=$(ls $PATH_FASTQ | sed 's/_.*$//g' | sort | uniq) 
for sample in $all_samples; do 
	mkdir -p $PATH_OUT_BAM/${sample}
	sh /scripts/mapping.sh $sample
done 

##############
############## 1.2 REORDER BAM (reorders reads in BAM file to match the contig ordering in a provided reference file)
##############
PATH_ASSEMBLY_SORTED="path/to/assembly/GRCh37.dna.primary_assembly.fa"
PATH_BAM=$PATH_OUT_BAM

all_samples=$(ls $PATH_BAM | sed 's/outs//g')
for sample in $all_samples;
do 
	tmpFolder="./"
	java -Xmx15g -XX:+UseSerialGC -Djava.io.tmpdir=$tmpFolder -jar ${PICARDPATH}/picard.jar ReorderSam I=${PATH_BAM}/${sample}/${sample}.bam O=${PATH_BAM}/${sample}/${sample}.reorder.bam REFERENCE=$PATH_ASSEMBLY_SORTED VALIDATION_STRINGENCY=STRICT CREATE_INDEX=true
done

##############
############## 1.3. MARK AND REMOVE DUPLICATES
##############
PATH_ASSEMBLY_SORTED="path/to/assembly/GRCh37.dna.primary_assembly.fa"
PATH_BAM=$PATH_OUT_BAM

all_samples=$(ls $PATH_BAM | sed 's/outs//g')
for sample in $all_samples;
do 
 tmpFolder="./"
 java -Xmx8g -XX:+UseSerialGC -Djava.io.tmpdir=$tmpFolder -jar ${PICARDPATH}/picard.jar MarkDuplicates I=${PATH_BAM}/${sample}/${sample}.reorder.bam O=${PATH_BAM}/${sample}/${sample}.rmDup.bam METRICS_FILE=${PATH_BAM}/${sample}/${sample}.rmDup.stats REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=STRICT CREATE_INDEX=true
done

#Summary 1.3: MarkDuplicates picard
# It locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g. library construction using PCR. Duplicate reads can also result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument. These duplication artifacts are referred to as optical duplicates. The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in a SAM/BAM file. After duplicate reads are collected, the tool differentiates the primary and duplicate reads using an algorithm that ranks reads by the sums of their base-quality scores (default method). MarkDuplicates also produces a metrics file indicating the numbers of duplicates for both single-and paired-end reads.
# During the sequencing process, the same DNA molecules can be sequenced several times. The resulting duplicate reads are not informative and should not be counted as additional evidence for or against a putative variant. 
# - REMOVE_DUPLICATES: If true do not write duplicates to the output file instead of writing them with appropriate flags set. Default value: false. 
# - CREATE_INDEX (Boolean)	Whether to create a BAM index when writing a coordinate-sorted BAM file. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false}

##############
############## 1.4.COVERAGE AND MAPPING STATISTICS
##############

############## 1.4.1. COVERAGE (INPUT: .rmDup.bam)
PATH_ASSEMBLY_SORTED="path/to/assembly/GRCh37.dna.primary_assembly.fa"
PATH_BAM=$PATH_OUT_BAM
PATH_OUT_COVERAGE="/path/to/Stats_output/Coverage"
PATH_INTERVALS="/path/to/intervals_exome_capture/intervals_covered.bed" #Intervals file used in the exome capture

all_samples=$(ls $PATH_BAM | sed 's/outs//g')
for sample in $all_samples;
do 
 # Coverage without a minBaseQuality and minMappingQuality
java -Xmx8g -XX:+UseSerialGC -jar ${EBROOTGATK}/GenomeAnalysisTK.jar -T DepthOfCoverage -R $PATH_ASSEMBLY_SORTED -o ${PATH_OUT_COVERAGE}/${sample}_coverage_V6 -I ${PATH_BAM}/${sample}/${sample}.rmDup.bam -L ${PATH_INTERVALS}
 # Coverage with a minBaseQuality and minMappingQuality
 java -Xmx8g -XX:+UseSerialGC -jar ${EBROOTGATK}/GenomeAnalysisTK.jar -T DepthOfCoverage -R $PATH_ASSEMBLY_SORTED -o ${PATH_OUT_COVERAGE}/${sample}_coverage_minQual -I ${PATH_BAM}/${sample}/${sample}.rmDup.bam --minBaseQuality 20 --minMappingQuality 20 -L ${PATH_INTERVALS}
 # DiagnoseTargets
java -Xmx8g -XX:+UseSerialGC -jar ${EBROOTGATK}/GenomeAnalysisTK.jar -T DiagnoseTargets -R $PATH_ASSEMBLY_SORTED -I ${PATH_BAM}/${sample}/${sample}.rmDup.bam -L ${PATH_INTERVALS} -o ${PATH_OUT_COVERAGE}/${sample}.DiagnoseTargets.vcf
done

# Merge each sample COV summary results in one file and nº of PASS intervals (DiagnoseTargets)
echo "sample_id	total mean granular_third_quartile granular_median granular_first_quartile %_bases_above_15" > COVERAGE_allsamples.txt
echo "sample_id	total mean granular_third_quartile granular_median granular_first_quartile %_bases_above_15" > COVERAGE_minQual_allsamples.txt
echo "sample_id PASS_intervals" > DiagnoseTargets_PASS_intervals.txt
for sample in $all_samples; do 
 awk -v var=$sample '{$1=""; if (NR==2) print var $0}' ${PATH_OUT_COVERAGE}/${sample}_coverage_V6.sample_summary >> COVERAGE_allsamples.txt
 awk -v var=$sample '{$1=""; if (NR==2) print var $0}' ${PATH_OUT_COVERAGE}/${sample}_coverage_minQual.sample_summary >> COVERAGE_minQual_allsamples.txt
 echo $sample $(grep -v "#" ${PATH_OUT_COVERAGE}/${sample}.DiagnoseTargets.vcf | grep "PASS" | wc -l)  >> DiagnoseTargets_PASS_intervals.txt
done

############## 1.4.2 MAPPING STATISTICS
############## Get nº of mapped and unmapped reads (efficiency of mapping)
PATH_ASSEMBLY_SORTED="path/to/assembly/GRCh37.dna.primary_assembly.fa"
PATH_BAM=$PATH_OUT_BAM
PATH_OUT_STATS="/path/to/Stats_output/"

all_samples=$(ls $PATH_BAM | sed 's/outs//g') 
echo "Sample Mapped Unmapped Mapped/Total" > ${PATH_OUT_STATS}/mapped_unmapped.txt
for sample in $all_samples; do
	samtools idxstats ${PATH_BAM}/${sample}/${sample}.rmDup.bam | awk -v sample="$sample" '{m+=$3; u+=$4} END {print sample" " m" "u" "m/(m+u) }' >> ${PATH_OUT_STATS}/mapped_unmapped.txt
done

echo "Sample Mapped" > ${PATH_OUT_STATS}/mapped_reads_dup.txt
for sample in $all_samples; do
	samtools idxstats ${PATH_BAM}/${sample}/${sample}.reorder.bam | awk -v sample="$sample" '{m+=$3; u+=$4} END {print sample" " m}' >> ${PATH_OUT_STATS}/mapped_reads_dup.txt
done

echo "Sample Mapped" > ${PATH_OUT_STATS}/mapped_reads_nodup.txt
for sample in $all_samples; do
	samtools idxstats ${PATH_BAM}/${sample}/${sample}.rmDup.bam | awk -v sample="$sample" '{m+=$3; u+=$4} END {print sample" " m}' >> ${PATH_OUT_STATS}/mapped_reads_nodup.txt
done

############## More statistics
PATH_ASSEMBLY_SORTED="path/to/assembly/GRCh37.dna.primary_assembly.fa"
PATH_BAM=$PATH_OUT_BAM
PATH_OUT_STATS="/path/to/Stats_output/"
PATH_INTERVALS="/path/to/intervals_exome_capture/intervals_covered.bed"

all_samples=$(ls $PATH_BAM | sed 's/outs//g')
for sample in $all_samples;
do 
 # GC content
	java -Xmx4g -jar picard.jar CollectGcBiasMetrics INPUT=${PATH_BAM}/${sample}/${sample}.rmDup.bam OUTPUT=${PATH_OUT_STATS}/GC/${sample}.gccov.stats CHART=${PATH_OUT_STATS}/GC/${sample}.gc_bias.pdf R=${PATH_ASSEMBLY_SORTED} S=${PATH_OUT_STATS}/GC/${sample}.summary_metrics.txt VALIDATION_STRINGENCY=STRICT ASSUME_SORTED=false
 # Insert size
	java -Xmx30g -jar picard.jar CollectInsertSizeMetrics INPUT=${PATH_BAM}/${sample}/${sample}.rmDup.bam OUTPUT=${PATH_OUT_STATS}/InsertSize/${sample}.insert_size_metrics H=${PATH_OUT_STATS}/InsertSize/${sample}.insert_size_histogram.pdf VALIDATION_STRINGENCY=STRICT
 # Alignment summary metrics
	java -Xmx25g -jar picard.jar CollectAlignmentSummaryMetrics INPUT=${PATH_BAM}/${sample}/${sample}.rmDup.bam REFERENCE_SEQUENCE=${PATH_ASSEMBLY_SORTED} OUTPUT=${PATH_OUT_STATS}/AlignmentSummary/${sample}.summary.stats VALIDATION_STRINGENCY=STRICT
done

# Get strand bias for each sample from Alignment summary metrics
all_samples=$(ls $PATH_BAM | sed 's/outs//g')
echo "Sample Strand_Balance" > ${PATH_OUT_STATS}/strand_balance.txt
for sample in $all_samples; do
 awk -v var=$sample '{if (NR==10) print var" "$22}' ${PATH_OUT_STATS}/AlignmentSummary/${sample}.summary.stats >> ${PATH_OUT_STATS}/strand_balance.txt
done

# Hs metrics
java -jar picard.jar BedToIntervalList I=${PATH_INTERVALS} O=${PATH_OUT_STATS}/HsMetrics/intervals SD=$PATH_ASSEMBLY_SORTED

for sample in $all_samples;
do 
	java -Xmx25g -jar picard.jar CollectHsMetrics INPUT=${PATH_BAM}/${sample}/${sample}.rmDup.bam REFERENCE_SEQUENCE=${PATH_ASSEMBLY_SORTED} OUTPUT=${PATH_OUT_STATS}/HsMetrics/${sample}.summary.stats BAIT_INTERVALS=${PATH_OUT_STATS}/HsMetrics/intervals TARGET_INTERVALS=${PATH_OUT_STATS}/HsMetrics/intervals VALIDATION_STRINGENCY=STRICT
done

# Get some Hs metrics for each sample
all_samples=$(ls $PATH_BAM | sed 's/outs//g')
echo "Sample %SELECTED_BASES %OFF_BAIT %TARGET_BASES_1X %TARGET_BASES_2X %TARGET_BASES_10X %TARGET_BASES_20X %TARGET_BASES_30X %TARGET_BASES_40X %TARGET_BASES_50X %TARGET_BASES_100X" > ${PATH_OUT_STATS}/hsMetrics.txt
for sample in $all_samples; do
 awk -v var=$sample '{if (NR==8) print var" "$19" "$20" "$36" "$37" "$38" "$39" "$40" "$41" "$42" "$43}' ${PATH_OUT_STATS}/HsMetrics/${sample}.summary.stats >> ${PATH_OUT_STATS}/hsMetrics.txt
done


##############
############## 1.5. INDEL REALIGNMENT
##############
# Reads that align on the edges of indels often get mapped with mismatching bases that might look like evidence for SNPs, but are actually mapping artifacts. The realignment process identifies the most consistent placement of the reads with respect to the indel in order to clean up these artifacts. First the program identifies intervals that need to be realigned, then, in the second step, it determines the optimal consensus sequence and performs the actual realignment of reads
PATH_ASSEMBLY="path/to/assembly"
PATH_ASSEMBLY_SORTED="path/to/assembly/GRCh37.dna.primary_assembly.fa"
PATH_BAM=$PATH_OUT_BAM

PATH_KNOWN_INDELS="/path/to/1000G_phase1.indels.b37.vcf"
PATH_KNOWN_INDELS_UPDATE="/path/to/1000G_phase1.indels.b37.update.vcf"
PATH_KNOWN_INDELS_SORTED="/path/to/1000G_phase1.indels.b37.sorted.vcf"

#Update and sort known.vcf
java -Xmx25g -jar picard.jar UpdateVcfSequenceDictionary I=${PATH_KNOWN_INDELS} O=${PATH_KNOWN_INDELS_UPDATE} SEQUENCE_DICTIONARY=${PATH_ASSEMBLY_SORTED}
java -Xmx25g -jar picard.jar SortVcf I=${PATH_KNOWN_INDELS_UPDATE} O=${PATH_KNOWN_INDELS_SORTED} SEQUENCE_DICTIONARY=${PATH_ASSEMBLY}/GRCh37.dna.primary_assembly_sorted.dict

all_samples=$(ls $PATH_BAM | sed 's/outs//g')
for sample in $all_samples;
do 
 # 1.5.1 RealignerTargetCreator 
 java -Xmx8g -XX:+UseSerialGC -jar ${EBROOTGATK}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${PATH_ASSEMBLY_SORTED} -I ${PATH_BAM}/$sample/$sample.rmDup.bam -known ${PATH_KNOWN_INDELS_SORTED} -o ${PATH_BAM}/$sample/$sample.realigner.intervals
 echo $jid1
 # 1.5.2 IndelRealigner 
java -Xmx8g -XX:+UseSerialGC -jar ${EBROOTGATK}/GenomeAnalysisTK.jar -T IndelRealigner -R ${PATH_ASSEMBLY_SORTED} -I ${PATH_BAM}/$sample/$sample.rmDup.bam -known ${PATH_KNOWN_INDELS_SORTED} -targetIntervals ${PATH_BAM}/$sample/$sample.realigner.intervals -o ${PATH_BAM}/$sample/$sample.realigned.bam
done


##############
############## 1.6. BQSR
##############
# The per-base estimate of error known as the base quality score is the foundation upon which all statistical calling algorithms are based. We have found that the estimates provided by the sequencing machines are often inaccurate and/or biased. The recalibration process applies an empirically accurate error model to the bases, producing a BAM file that is suitable for analysis.

##!!! First!!! Sort and update known sites (as in 1.5 IndelRealignment step)

PATH_KNOWN_INDELS_SORTED="/path/to/1000G_phase1.indels.b37.sorted.vcf"
PATH_KNOWN_1000G_SORTED="/path/to/1000G_omni2.5.b37.sorted.vcf" 
PATH_KNOWN_dbSNP_SORTED="/path/to/dbsnp_138.b37.sorted.vcf"
PATH_KNOWN_HapMap_SORTED="/path/to/hapmap_3.3.b37.sorted.vcf"
PATH_KNOWN_Mills_SORTED="/path/to/Mills_and_1000G_gold_standard.indels.b37.sorted.vcf"

PATH_ASSEMBLY_SORTED="path/to/assembly/GRCh37.dna.primary_assembly.fa"
PATH_BAM=$PATH_OUT_BAM

all_samples=$(ls $PATH_BAM | sed 's/outs//g')
for sample in $all_samples;
do  
 #1.6.1 Analyze patterns of covariation
 java -Xmx8g -XX:+UseSerialGC -jar ${EBROOTGATK}/GenomeAnalysisTK.jar  -T BaseRecalibrator -R ${PATH_ASSEMBLY_SORTED} -I ${PATH_BAM}/$sample/$sample.realigned.bam -knownSites ${PATH_KNOWN_1000G_SORTED} -knownSites ${PATH_KNOWN_INDELS_SORTED} -knownSites ${PATH_KNOWN_dbSNP_SORTED} -knownSites ${PATH_KNOWN_HapMap_SORTED} -knownSites ${PATH_KNOWN_Mills_SORTED} -o ${PATH_BAM}/$sample/$sample.RecalibrationFile.grp
  echo $jid2
 #1.6.2 Apply the recalibration
java -Xmx8g -XX:+UseSerialGC -jar ${EBROOTGATK}/GenomeAnalysisTK.jar -T PrintReads -R ${PATH_ASSEMBLY_SORTED} -I ${PATH_BAM}/$sample/$sample.realigned.bam -BQSR ${PATH_BAM}/$sample/$sample.RecalibrationFile.grp -o ${PATH_BAM}/$sample.recalibrated.bam
done

##If the sample in BAM ReadGroups was WRONG (always ${sample}); do:
flowcell="1"; lane="1"; library="1";
for sample in $all_samples;do java -Xmx8g -jar picard.jar AddOrReplaceReadGroups I=$PATH_BAM}/$sample.recalibrated.bam O=${PATH_BAM}/$sample.recalibrated.readgroups.bam RGLB=$library RGPL=Illumina RGPU=${flowcell}_${lane} RGSM=${sample} COMPRESSION_LEVEL=9; done
for sample in $all_samples; do samtools index ${PATH_BAM}/$sample.recalibrated.readgroups.bam; done

########################################################################################################
############## 2. FROM BAM TO VCF (Variant Calling)
########################################################################################################

PATH_ASSEMBLY_SORTED="path/to/assembly/GRCh37.dna.primary_assembly.fa"
PATH_KNOWN_dbSNP_SORTED="/path/to/dbsnp_138.b37.sorted.vcf"
PATH_INTERVALS="/path/to/intervals_exome_capture/intervals_covered.bed"
PATH_BAM=$PATH_OUT_BAM
PATH_VCF="/path/to/output_folder_VCF_files"

##############
############## 2.1. HAPLOTYPE CALLER
##############
all_samples=$(ls $PATH_BAM | sed 's/outs//g')
for sample in $all_samples;
do 
	java -Xmx8g -XX:+UseSerialGC -jar ${EBROOTGATK}/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${PATH_ASSEMBLY_SORTED} -I ${PATH_BAM}/$sample.recalibrated.readgroups.bam --dbsnp $PATH_KNOWN_dbSNP_SORTED --genotyping_mode DISCOVERY  -L $PATH_INTERVALS --emitRefConfidence GVCF -o $PATH_VCF/${sample}.HC.raw_variants.g.vcf
done

##############
############## 2.2. GENOTYPE GVCFs
##############

samples=$(ls $PATH_VCF/*HC*g.vcf | sed 's/^\/homes/--variant \/homes/g' ) # get all files ending with g.vcf and add --variant before it

 java -Xmx8g -XX:+UseSerialGC -jar ${EBROOTGATK}/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ${PATH_ASSEMBLY_SORTED}  --dbsnp $PATH_KNOWN_dbSNP_SORTED -stand_call_conf 20.0 -L $PATH_INTERVALS -o $PATH_VCF/all_samples.WES.GC.vcf $(echo $samples)

##############
############## 2.3. VQSR
##############
PATH_KNOWN_1000G_SNPS_SORTED="/path/to/1000G_phase1.snps.high_confidence.b37.sorted.vcf"
PATH_KNOWN_1000G_SORTED="/path/to/1000G_omni2.5.b37.sorted.vcf" 
PATH_KNOWN_dbSNP_SORTED="/path/to/dbsnp_138.b37.sorted.vcf"
PATH_KNOWN_HapMap_SORTED="/path/to/hapmap_3.3.b37.sorted.vcf"
PATH_KNOWN_Mills_SORTED="/path/to/Mills_and_1000G_gold_standard.indels.b37.sorted.vcf"

PATH_ASSEMBLY_SORTED="path/to/assembly/GRCh37.dna.primary_assembly.fa"
PATH_VCF="/path/to/output_folder_VCF_files"

############# 2.3.1 FOR SNPS
#############2.3.1.1 VARIANT RECALIBRATOR 

java -Xmx8g -XX:+UseSerialGC -jar ${EBROOTGATK}/GenomeAnalysisTK.jar -T VariantRecalibrator -R $PATH_ASSEMBLY_SORTED -input ${PATH_VCF}/all_samples.WES.GC.vcf -recalFile ${PATH_VCF}/VQSR/all_samples.WES.joint_variants_raw_SNPs.recal -tranchesFile ${PATH_VCF}/VQSR/joint_variants_Raw_SNPs.tranches -nt 4  -mode SNP -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $PATH_KNOWN_HapMap_SORTED -resource:omni,known=false,training=true,truth=true,prior=12.0 $PATH_KNOWN_1000G_SORTED -resource:1000G,known=false,training=true,truth=false,prior=10.0 $PATH_KNOWN_1000G_SNPS_SORTED -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $PATH_KNOWN_dbSNP_SORTED -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an InbreedingCoeff -rscriptFile ${PATH_VCF}/VQSR/recal_snp.plots.R

#############2.3.1.2 APPLY RECALIBRATION 
java -Xmx8g -XX:+UseSerialGC -jar ${EBROOTGATK}/GenomeAnalysisTK.jar  -T ApplyRecalibration -R $PATH_ASSEMBLY_SORTED -input ${PATH_VCF}/all_samples.WES.GC.vcf -recalFile ${PATH_VCF}/VQSR/all_samples.WES.joint_variants_raw_SNPs.recal -tranchesFile ${PATH_VCF}/VQSR/joint_variants_Raw_SNPs.tranches -mode SNP --ts_filter_level 99.5 -o ${PATH_VCF}/all_samples.WES.snp_filtered.vcf

############# 2.3.2 FOR INDELS
#############2.3.2.1 VARIANT RECALIBRATOR
java -Xmx8g -XX:+UseSerialGC -jar ${EBROOTGATK}/GenomeAnalysisTK.jar -T VariantRecalibrator -R $PATH_ASSEMBLY_SORTED -input ${PATH_VCF}/all_samples.WES.snp_filtered.vcf -recalFile ${PATH_VCF}/VQSR/all_samples.WES.joint_variants_raw_INDELS.recal -tranchesFile ${PATH_VCF}/VQSR/joint_variants_Raw_INDELS.tranches -nt 4 --maxGaussians 4 -mode INDEL -resource:mills,known=false,training=true,truth=true,prior=12.0 $PATH_KNOWN_Mills_SORTED -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $PATH_KNOWN_dbSNP_SORTED -an QD -an SOR -an MQRankSum -an ReadPosRankSum -an FS -an InbreedingCoeff -rscriptFile ${PATH_VCF}/VQSR/recal_indels.plots.R
#############2.3.2.2 APPLY RECALIBRATION 
java -Xmx8g -XX:+UseSerialGC -jar ${EBROOTGATK}/GenomeAnalysisTK.jar -T ApplyRecalibration -R $PATH_ASSEMBLY_SORTED -input ${PATH_VCF}/all_samples.WES.snp_filtered.vcf -recalFile ${PATH_VCF}/VQSR/all_samples.WES.joint_variants_raw_INDELS.recal -tranchesFile ${PATH_VCF}/VQSR/joint_variants_Raw_INDELS.tranches -mode INDEL --ts_filter_level 99.0 -o ${PATH_VCF}/all_samples.WES.snp_indel_filtered.vcf

