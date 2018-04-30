Store each analyses-specific script and outputs
===============================================
# Directory set-up #
> this analyses directory should be set up to live parallelly to the TCGA_data directory to allow proper data sourcing. Due to the protected nature of some of those files they are not shared publicly. 

# Dependency files #
>global_aes_out.R: used to define plotting functions and color codes for other analysis R codes

>dependency_files.R: load and prep files used for downstream analysis

# Analysis workflow #

# Genotype data analysis #
## 1. process_files/:
> Scripts in process_files/genotype/ are used to preprocess germline variants calls to classified variants.

## 2. PCA_IBD_dist/:
> Scripts used for PCA and relatedness analysis

# Exomic germline variants analysis #
## 1. process_files/:
> Scripts in process_files/germline/ are used to preprocess germline variants calls to classified variants. Steps downstream of GATK/VarScan/Pindel merged calls are described in process_files/germline/local/work.log.sh. 

## 2. sample_listing/:
> Compare PCA germline sample list (based on ISB-CGC manifest) against PCA clinical sample and MC3 sample. Generate population level summary for the 10,467 final samples.

## 3. pathogenic_variants/: 
> Scripts in pathogenic_variants/ are used to further filter variants based on readcount and cancer relevance. 

## 4. association_test/: 
> Conduct association test using ExAC and ExAC-nonTCGA data.

## 5. LOH/:
> Conduct LOH analysis using readcount data. 

## 6. expression_effect/:
> Calculate cohort level expression quantile and retrieve RSEM expression value for each sample within cancer types.

## 7. data_integration/:
> Combine multi-level data for pathogenic variants. 

## 8. hotspot3d/: 
> Co-clustering with somatic mutations using HotSpot3D. 

## 9. data_intergration/:
> Integrate pathogenic variant data with other molecular and omics data
> plot results from data_integration

## All subsequent analysis require integrated file ##
> source dependency_files.R to read in these additional data

## 10. integrated_analysis/:
> generate plots and stats for pathogenic variants that required integrated datap

## 