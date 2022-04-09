# WES Processing and Quality Control
Helper script and documentation to process "whole-exome sequences" (WES) from FASTQ to VCF files, filtering, and QC assessment for population genetics analyses
Pipeline applied in [Font-Porterias, N. et al. *MBE*, 2021](https://doi.org/10.1093/molbev/msab070) on  WES data ([Agilent SureSelect Human All Exon V6 capture kit](https://www.agilent.com/cs/library/datasheets/public/SureSelect%20V6%20DataSheet%205991-5572EN.pdf)), pair-end sequencing at 30X. This repo is not under active development, consider code as it is and check the software versions specified below. 

This pipeline is adapted from: [Van der Auwera, GA., et al. GATK Best Practices pipeline, 2013](https://doi.org/10.1002/0471250953.bi1110s43)

## Citation

If you use these pipeline, please cite: 

> **Font-Porterias N**, Caro-Consuegra R, Lucas-Sánchez M, Lopez M, Giménez A, Carballo-Mesa A, Bosch E, Calafell F, Quintana-Murci L, Comas D. (2021) The Counteracting Effects of Demography on Functional Genomic Variation: The Roma Paradigm. *MBE*, 38(7), 2804-2817, https://doi.org/10.1093/molbev/msab070

> Van der Auwera GA , Carneiro MO , Hartl C , Poplin R , del Angel G , Levy-Moonshine A , Jordan T , Shakir K , Roazen D , Thibault J , et al.  2013. From fastQ data to high-confidence variant calls: the Genome Analysis Toolkit best practices pipeline. *Curr Protoc Bioinformatics* 43:11.10.1–11.10.33, https://doi.org/10.1002/0471250953.bi1110s43



## Pipeline
### 1. Preprocessing: Mapping and variant calling
*Required softwares and version:* ; ; ;
*Required files:* ; ; ;

Brielfy, it performs: (i) a QC assessment of the FASTQ files, (ii) adapters trimming, (iii) mapping to the GRCh37 human reference, (iv) coverage and mapping statistics, and (v) variant calling. 

See documentation [here](Documentation/Pipeline1_WES_preprocessing_pipeline.pdf), for a complete description of the preprocessing pipeline. 

Edit the paths and run:
```bash
sh Script1_mapping_varcalling_exomes.sh
```

### 2. Variant Filtering and Quality Control
*Required softwares and version:* ; ; ;

Brielfy, it performs: (i) Variant and samle filtering, (ii) VCF QC assessment, and (iii) relatedness estimation. 

See documentation [here](Documentation/Pipeline2.VariantFiltering.pdf), for a complete description of the filtering and QC pipeline. 

Edit the paths and run:
```bash
sh Script2_FilteringVCF.sh
```
