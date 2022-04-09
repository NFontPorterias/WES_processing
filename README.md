# WES Processing and Quality Control
Helper script and documentation to process "whole-exome sequences" (WES) from FASTQ to VCF files, filtering, and QC assessment for population genetics analyses
Pipeline applied in [Font-Porterias, N. et al. *MBE*, 2021](https://doi.org/10.1093/molbev/msab070) on  WES data ([Agilent SureSelect Human All Exon V6 capture kit](https://www.agilent.com/cs/library/datasheets/public/SureSelect%20V6%20DataSheet%205991-5572EN.pdf)) sequenced at 30X. This repo is not under active development, consider code as it is and check the software versions specified below. 



## Citation

If you use these pipeline, please cite: 

> **Font-Porterias N**, Caro-Consuegra R, Lucas-Sánchez M, Lopez M, Giménez A, Carballo-Mesa A, Bosch E, Calafell F, Quintana-Murci L, Comas D. (2021) The Counteracting Effects of Demography on Functional Genomic Variation: The Roma Paradigm. *MBE*, 38(7), 2804-2817, https://doi.org/10.1093/molbev/msab070


## Pipeline
### 1. Preprocessing: Mapping and variant calling
*Required softwares and version:* ; ; ;

pair end reads, GRCh37, data paths

See documentation [here](Documentation/Pipeline1_WES_preprocessing_pipeline.pdf), for a complete description of the preprocessing pipeline. 

Edit the paths and run:
```bash
sh Script1_mapping_varcalling_exomes.sh
```

### 2. Variant Filtering and Quality Control
*Required softwares and version:* ; ; ;

See documentation [here](Documentation/Pipeline2.VariantFiltering.pdf), for a complete description of the filtering and QC pipeline. 

Edit the paths and run:
```bash
sh Script2_FilteringVCF.sh
```
