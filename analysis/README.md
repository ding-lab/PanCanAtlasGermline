Store each analyses-specific script and outputs
===============================================

# dependency files #
>global_aes_out.R: used to define plotting functions and color codes for other analysis R codes

>dependency_files.R: load and prep files used for downstream analysis

# Analysis workflow #

## 1. process_files/:
> Scripts in process_files/germline/ are used to preprocess germline variants calls to classified variants. Steps downstream of GATK/VarScan/Pindel merged calls are described in process_files/germline/local/work.log.sh. 

## 2. pathogenic_variants/: 
> Scripts in pathogenic_variants/ are used to further filter variants based on readcount and cancer relevance. 