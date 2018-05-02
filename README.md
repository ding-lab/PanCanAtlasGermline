# PanCanAtlasGermline #

## Pathogenic Germline Variants in 10,389 Adult Cancers ## 
>The TCGA PanCanAtlas germline working group is investigating germline variants in the largest sequencing cohort of cancer to date: 10,389 cases in 33 cancer types.     
>Publication DOI: https://doi.org/10.1016/j.cell.2018.03.039

## Data Access ##
>De-identified variant-level data for prioritized VUS and pathogenic variants: Table S2 of the publication.

>The protected variants+sample ID and the full callset (Authorized User only):
>1) Available on ISB-CGC. The steps to get access and use the cloud resource is documented below (Getting started).
>2) GDC link will likely be available soon. Please keep checking back here.

## working group info ##
>Wiki: https://wiki.nci.nih.gov/display/TCGAM/PanCanAtlas+Germline+AWG 
>Synapse: syn4602499  

## Directory set-up ##
> The analyses directory should be set up to live parallelly to the TCGA_data directory to allow proper data sourcing. Due to the protected nature of some of those files they are not shared publicly. 

## Getting started: ##  
>1) Install Google Cloud SDK (https://cloud.google.com/sdk/docs/quickstarts) and read through at least the basic gsutil command 
>2) Make sure you can access the project on Google cloud (https://console.cloud.google.com/home/dashboard?project=isb-cgc-06-0004); if not here are some relevant steps (http://isb-cancer-genomics-cloud.readthedocs.io/en/latest/sections/webapp/Gaining-Access-To-Contolled-Access-Data.html)  
>3) Read this one page quick start to google cloud compute engine (https://cloud.google.com/compute/docs/quickstarts) . 
>4) If you want to download the data to your local computer/cluster, you need to make sure with your system administrator that the environment you plan to work on have appropriate access authority and security.  
>5) Read Jay's short getting started guide (https://drive.google.com/file/d/0B0aS3CDIQ_RAd01ld0tKX3JHa00/view). And if you are up to it here is a more detailed guide by Sheila (https://docs.google.com/document/d/1f1YBVG1dAhpF-Un5lp70kMI8Yo2rD3NNQdpl0FDAu-c/edit#heading=h.gv4f8tqq5731)  

>The PCAGermline AWG google doc of shared data (https://docs.google.com/document/d/1ymdfAnRR4o4-20bwHI3vPaRPRuoqtqc0pNUVYO2oiPc/edit) 

> Contact: Kuan-lin Huang [kuan-lin.huang@wustl.edu]
