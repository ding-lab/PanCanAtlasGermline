# PanCanAtlasGermline #

## Pathogenic Germline Variants in 10,389 Adult Cancers ## 
>The TCGA PanCanAtlas germline analysis working group is investigating germline variants in the largest sequencing cohort of cancer to date: 10,389 cases in 33 cancer types.     
```Pathogenic Germline Variants in 10,389 Adult Cancers.

Huang KL, Mashl RJ, Wu Y, Ritter DI, Wang J, Oh C, Paczkowska M, Reynolds S, Wyczalkowski MA, Oak N, Scott AD, Krassowski M, Cherniack AD, Houlahan KE, Jayasinghe R, Wang LB, Zhou DC, Liu D, Cao S, Kim YW, Koire A, McMichael JF, Hucthagowder V, Kim TB, Hahn A, Wang C, McLellan MD, Al-Mulla F, Johnson KJ; Cancer Genome Atlas Research Network, Lichtarge O, Boutros PC, Raphael B, Lazar AJ, Zhang W, Wendl MC, Govindan R, Jain S, Wheeler D, Kulkarni S, Dipersio JF, Reimand J, Meric-Bernstam F, Chen K, Shmulevich I, Plon SE, Chen F, Ding L.

Cell. 2018 Apr 5;173(2):355-370.e14. doi: 10.1016/j.cell.2018.03.039.

PMID: 29625052
```

## Data Access ##
>De-identified variant-level data for prioritized VUS and pathogenic variants: Table S2 of the publication.

>The protected variants+sample ID and the full callset (Authorized User only):

>GDC link: https://gdc.cancer.gov/about-data/publications/PanCanAtlas-Germline-AWG

>Compressed VCF file of the combined, filtered variant calls using GATK, VarScan2, and Pindel on WES data of the 10,389 final passed-QC samples. - PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz

>Tabix file of the compressed VCF file of the combined, filtered variant calls using GATK, VarScan2, and Pindel on WES data of the 10,389 final passed-QC samples. - PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz.tbi

>Prioritized, cancer related variants discovered in 10,389 cases. Please use "Overall_Classification" column to distinguish between Pathogenic, Likely Pathogenic and Priortizied VUSs. - PCA_pathVar_integrated_filtered_adjusted.tsv

## PanCanAtlas Germline Working Group info (for members) ##
>Wiki: https://wiki.nci.nih.gov/display/TCGAM/PanCanAtlas+Germline+AWG 
>Synapse: syn4602499  

## Directory set-up ##
> The analyses directory should be set up to live parallelly to the TCGA_data directory to allow proper data sourcing. Due to the protected nature of some of those files they are not shared publicly. 
> Analyses scripts are in the "analysis" folder

## Getting started on ISB-CGC [For TCGA PanCanAtlas Germline AWG members Only]: ##  
>1) Install Google Cloud SDK (https://cloud.google.com/sdk/docs/quickstarts) and read through at least the basic gsutil command 
>2) Make sure you can access the project on Google cloud (https://console.cloud.google.com/home/dashboard?project=isb-cgc-06-0004); if not here are some relevant steps (http://isb-cancer-genomics-cloud.readthedocs.io/en/latest/sections/webapp/Gaining-Access-To-Contolled-Access-Data.html)  
>3) Read this one page quick start to google cloud compute engine (https://cloud.google.com/compute/docs/quickstarts) . 
>4) If you want to download the data to your local computer/cluster, you need to make sure with your system administrator that the environment you plan to work on have appropriate access authority and security.  
>5) Read Jay's short getting started guide (https://drive.google.com/file/d/0B0aS3CDIQ_RAd01ld0tKX3JHa00/view). And if you are up to it here is a more detailed guide by Sheila (https://docs.google.com/document/d/1f1YBVG1dAhpF-Un5lp70kMI8Yo2rD3NNQdpl0FDAu-c/edit#heading=h.gv4f8tqq5731)  

>The PCAGermline AWG google doc of shared data (https://docs.google.com/document/d/1ymdfAnRR4o4-20bwHI3vPaRPRuoqtqc0pNUVYO2oiPc/edit) 


> Contact: Kuan-lin Huang [kuan-lin.huang@wustl.edu]
