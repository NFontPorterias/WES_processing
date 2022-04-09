#############################################################################################################
##########################################     VCF FILTERING-1     ##########################################
#############################################################################################################
PATH_VCF="/path/to/output_folder_VCF_files"
VCF_filename="filename"
PATH_VCF_OUT="path/to/out_VCF_filtered_and_QC"

### 1. We kept only PASS variants. 
vcftools --vcf ${PATH_VCF}/${VCF_filename}.WES.snp_indel_filtered.vcf --remove-filtered-all --recode --out ${PATH_VCF_OUT}/${VCF_filename}.WES.PASS

### 2. We exclude all indels. 
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}.WES.PASS.recode.vcf --remove-indels --recode --recode-INFO-all --out ${PATH_VCF_OUT}/${VCF_filename}.WES.PASS_noIndels

### 3. We exclude X,Y chromosomes 
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}.WES.PASS_noIndels.recode.vcf --not-chr X --not-chr Y --recode --out ${PATH_VCF_OUT}/${VCF_filename}.WES.PASS_noIndels_noXY

### 4. We keep only biallelic sites. 
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}.WES.PASS_noIndels_noXY.recode.vcf --min-alleles 2 --max-alleles 2 --recode --out ${PATH_VCF_OUT}/${VCF_filename}.WES.PASS_noIndels_noXY_biallelic

### 5. We keep only not fixed sites In order not to bias the measures of depth of coverage and missingness with SNPs specific to the excluded populations, we applied a first round  of SNPs filtering and excluded all sites that were fixed (from the --freq2 file)
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}.WES.PASS_noIndels_noXY_biallelic.recode.vcf --freq2 --out ${PATH_VCF_OUT}/FREQ_${VCF_filename}.WES.PASS_noIndels_noXY_biallelic

awk '{if ($5==0) print}' ${PATH_VCF_OUT}/FREQ_${VCF_filename}.WES.PASS_noIndels_noXY_biallelic.frq > ${PATH_VCF_OUT}/FIXED_${VCF_filename}.WES.PASS_noIndels_noXY_biallelic.txt
awk '{if ($6==0) print}' ${PATH_VCF_OUT}/FREQ_${VCF_filename}.WES.PASS_noIndels_noXY_biallelic.frq >> ${PATH_VCF_OUT}/FIXED_${VCF_filename}.WES.PASS_noIndels_noXY_biallelic.txt

vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}.WES.PASS_noIndels_noXY_biallelic.recode.vcf --exclude-positions ${PATH_VCF_OUT}/FIXED_${VCF_filename}.WES.PASS_noIndels_noXY_biallelic.txt --recode --out ${PATH_VCF_OUT}/${VCF_filename}.WES.PASS_noIndels_noXY_biallelic_noFixed

#############################################################################################################
##########################################   VCF QUALITY CONTROL   ##########################################
#############################################################################################################

### Depth of Coverage per site
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}.WES.PASS_noIndels_noXY_biallelic_noFixed.recode.vcf --site-mean-depth --out ${PATH_VCF_OUT}/DEPTH_${VCF_filename}.WES
### Depth of Coverage per individual
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}.WES.PASS_noIndels_noXY_biallelic_noFixed.recode.vcf --depth --out ${PATH_VCF_OUT}/DEPTH_${VCF_filename}.WES
### Missingness per individual
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}.WES.PASS_noIndels_noXY_biallelic_noFixed.recode.vcf --missing-indv --out ${PATH_VCF_OUT}/MISSING_${VCF_filename}.WES
### Missingness per site
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}.WES.PASS_noIndels_noXY_biallelic_noFixed.recode.vcf --missing-site --out ${PATH_VCF_OUT}/MISSING_${VCF_filename}.WES
### Genotype Quality
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}.WES.PASS_noIndels_noXY_biallelic_noFixed.recode.vcf --extract-FORMAT-info GQ --out ${PATH_VCF_OUT}/GQ_${VCF_filename}.WES
### Ts/Tv ratio
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}.WES.PASS_noIndels_noXY_biallelic_noFixed.recode.vcf --TsTv-summary --out ${PATH_VCF_OUT}/TsTv_${VCF_filename}.WES

###########################################################
#####   REPEAT QC AFTER REMOVING SITES MISS > 0.10   ##### 
###########################################################
# Remove sites with site missingness > 10%. 
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}.WES.PASS_noIndels_noXY_biallelic_noFixed.recode.vcf --max-missing 0.90 --recode --out ${PATH_VCF_OUT}/${VCF_filename}_mis10

### Depth of Coverage per site
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10.recode.vcf --site-mean-depth  --out ${PATH_VCF_OUT}/DEPTH_${VCF_filename}.WES_mis10
### Depth of Coverage per individual
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10.recode.vcf --depth --out ${PATH_VCF_OUT}/DEPTH_${VCF_filename}.WES_mis10
### Missingness per individual
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10.recode.vcf --missing-indv --out ${PATH_VCF_OUT}/MISSING_${VCF_filename}.WES_mis10
### Missingness per site
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10.recode.vcf --missing-site --out ${PATH_VCF_OUT}/MISSING_${VCF_filename}.WES_mis10
### Genotype Quality
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10.recode.vcf --extract-FORMAT-info GQ --out ${PATH_VCF_OUT}/GQ_${VCF_filename}.WES_mis10
### Ts/Tv ratio
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10.recode.vcf --TsTv-summary --out ${PATH_VCF_OUT}/TsTv_${VCF_filename}.WES_mis10


#############################################################################################################
##########################################     VCF FILTERING-2     ##########################################
#############################################################################################################

#6. Remove sites with DP < 5X, GQ < 20 (these values will be treated as missing) and site missingness > 5% 
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10.recode.vcf --minDP 5 --recode --recode-INFO-all --out ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5.recode.vcf --minGQ 20 --recode --recode-INFO-all --out ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20.recode.vcf  --max-missing 0.95 --recode --out ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5

#7. Remove bad quality individuals. 
##	7.1 Missigness > 5% 
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5.recode.vcf --missing-indv --out ${PATH_VCF_OUT}/MISSING_${VCF_filename}_mis10_DP5_GQ20_mis5
awk '{if ($5>0.05) print }' ${PATH_VCF_OUT}/MISSING_${VCF_filename}_mis10_DP5_GQ20_mis5.imiss > ${PATH_VCF_OUT}/ind_mis_higher5.txt

##	7.2 Mean depth of coverage <40
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5.recode.vcf --depth --out ${PATH_VCF_OUT}/DEPTH_${VCF_filename}_mis10_DP5_GQ20_mis5
awk '{if ($3<40) print $1}' ${PATH_VCF_OUT}/DEPTH_${VCF_filename}_mis10_DP5_GQ20_mis5.idepth > ${PATH_VCF_OUT}/ind_cov_below40.txt

## 7.3 Remove bad quality inds (miss > 5% + cov < 40)
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5.recode.vcf --remove ${PATH_VCF_OUT}/ind_cov_below40.txt --recode --out ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5_inds

#8. Calculate heterozigosity
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5_inds.recode.vcf --het --out ${PATH_VCF_OUT}/HET_${VCF_filename}_mis10_DP5_GQ20_mis5_inds
##--> Remove those inds with het < 4 SD together with related inds (step 12).

##	8.1 Inbreeding coefficient taking pop structure into account (using the individual allele frequencies predicted by the PCA)
plink --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5_inds.recode.vcf --make-bed --out ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5_inds
PCANGSD='python pcangsd/pcangsd.py'
$PCANGSD -plink ${VCF_filename}_mis10_DP5_GQ20_mis5_inds -inbreed 2 -o ${VCF_filename}_mis10_DP5_GQ20_mis5_inds.nonhom

#9. Calculate Ts/Tv ratio
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5_inds.recode.vcf --TsTv-summary --out ${PATH_VCF_OUT}/TsTv_${VCF_filename}_mis10_DP5_GQ20_mis5_inds

#10. hwe pvalue 10-3 
##	10.1 VCFTools 
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5_inds.recode.vcf --hwe 0.001 --recode --out ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5_inds_hweVCFtools

##	10.2 PCAngsd: taking pop structure into account 
PCANGSD='python pcangsd/pcangsd.py'
$PCANGSD -plink ${VCF_filename}_mis10_DP5_GQ20_mis5_inds -inbreedSites -minMaf 0 -o ${VCF_filename}_mis10_DP5_GQ20_mis5_inds -sites_save
$PCANGSD -plink ${VCF_filename}_mis10_DP5_GQ20_mis5_inds -HWE_filter ${VCF_filename}_mis10_DP5_GQ20_mis5_inds.lrt.sites.gz -o ${VCF_filename}_mis10_DP5_GQ20_mis5_inds -minMaf 0 -HWE_tole 0.001 -sites_save
tr '-' '\t' < ${VCF_filename}_mis10_DP5_GQ20_mis5_inds.sites > ${VCF_filename}_mis10_DP5_GQ20_mis5_inds_exclude.sites

## 10.3 VCFTools per POP 
pops="list of population in vcf"
for pop in $pops; do
	vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5_inds.recode.vcf --keep ${PATH_VCF_OUT}/inds_${pop}.txt --hardy --out ${PATH_VCF_OUT}/$pop
done
echo "chr pos pop" > ${PATH_VCF_OUT}/excluded_snps.txt
for pop in $pops; do awk -v pop=$pop '{if ($6 < 0.001) print $1" "$2" "pop}' ${PATH_VCF_OUT}/$pop.hwe >> ${PATH_VCF_OUT}/excluded_snps.txt; done 
cut -d' ' -f 1,2 ${PATH_VCF_OUT}/excluded_snps.txt | sort -u > ${PATH_VCF_OUT}/excluded_snps_uniq.txt 

vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5_inds.recode.vcf --exclude-positions ${PATH_VCF_OUT}/excluded_snps_uniq.txt --recode --out ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5_inds_hwe


#############################################################################################################
##########################################  RELATEDNESS ESTIMATION  #########################################
#############################################################################################################

#11. MAF 1%
vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5_inds_hwe.recode.vcf --maf 0.01 --recode --out ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5_inds_hwe_maf

#12. Compute kinship and find related individuals

## 12.1 KING
plink --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5_inds_hwe_maf.recode.vcf --make-bed --out ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5_inds_hwe_maf

join -1 1 -2 1 ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5_inds_hwe_maf.fam_old ind_pop.txt | awk '{$1=$7; $7=""; print}' > ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5_inds_hwe_maf.fam

KING/king -b ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5_inds_hwe_maf.bed --kinship --prefix RELATEDNESS/REL_${VCF_filename}_mis10_DP5_GQ20_mis5_inds_hwe_maf

awk '{if ($9>0.0442) print }' REL_${VCF_filename}_mis10_DP5_GQ20_mis5_inds_hwe_maf.kin

## 12.2 PLINK
plink --bfile ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5_inds_hwe_maf --genome --out RELATEDNESS/REL_${VCF_filename}_mis10_DP5_GQ20_mis5_inds_hwe_maf

awk '{if ($10>0.125) print }' REL_${VCF_filename}_mis10_DP5_GQ20_mis5_inds_hwe_maf.genome

##--> Remove the ind from the related pair where missigness is lower and save inds ID in: remove_samples.txt

vcftools --vcf ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5_inds_hwe.recode.vcf --remove ${PATH_VCF_OUT}/remove_samples.txt --recode --out ${PATH_VCF_OUT}/${VCF_filename}_mis10_DP5_GQ20_mis5_inds_hwe_norel


