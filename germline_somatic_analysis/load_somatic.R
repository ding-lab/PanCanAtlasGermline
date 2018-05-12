##### load_somatic.R #####
# Kuan-lin Huang 2018
# load somatic mutation/driver/gene list files

### MAIN ###
somaticDriver299_f = "../../TCGA_data/somatic/Driver_BaileyCell2018/299driverGene.txt"
somaticDriver299 = as.vector(t(read.table(header=F, quote = "", sep="\t", file = somaticDriver299_f, stringsAsFactors=FALSE)))

somatic_f = "../../TCGA_data/somatic/mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample.gz"
somatic = read.table(header=T, quote = "", sep="\t", file = gzfile(somatic_f), stringsAsFactors=FALSE)
somatic$bcr_patient_barcode = gsub("(^TCGA-[A-Z0-9][A-Z0-9]-[A-Z0-9][A-Z0-9][A-Z0-9][A-Z0-9])-.*","\\1",somatic$Tumor_Sample_Barcode)

clin_cmap = clin[,c("bcr_patient_barcode","type"),]
colnames(clin_cmap)[2] = "cancer"
somatic = merge(somatic,clin_cmap,by="bcr_patient_barcode")
# somatic_mut_count = data.frame(table(somatic$bcr_patient_barcode))
# colnames(somatic_mut_count) = c("bcr_patient_barcode","MutationCount")

table(somatic$Variant_Classification)
likelyFunctionalTypes = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation",
                          "Nonsense_Mutation","Splice_Site","Translation_Start_Site")
somatic_likelyfunctional = somatic[somatic$Hugo_Symbol %in% c(pathVar$HUGO_Symbol,somaticDriver299) & somatic$Variant_Classification %in% likelyFunctionalTypes,]

# driver mutation
driver_f = "../../TCGA_data/somatic/Driver_BaileyCell2018/Mutation.CTAT.3D.Scores.txt.gz"
driver = read.table(header=T, quote = "", sep="\t", file = gzfile(driver_f), stringsAsFactors=FALSE)
colnames(driver) = gsub("\\.","_",colnames(driver))
driver$numOfEvidence = driver$New_Linear__cancer_focused__flag + driver$New_Linear__functional__flag + driver$New_3D_mutational_hotspot_flag
table(driver$numOfEvidence)
driver_func = driver[driver$numOfEvidence > 1,]
somatic_likelyfunctional_driver = somatic_likelyfunctional[somatic_likelyfunctional$Variant_Classification != "Missense_Mutation" |
                                                             paste(somatic_likelyfunctional$Hugo_Symbol,somatic_likelyfunctional$HGVSp_Short) %in% paste(driver_func$gene,driver_func$protein_change),]
