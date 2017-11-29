#1. get variants
#get somatic mafs with genes of interest
cut -f1,5-7,9,11,13,16,37,38 /gscmnt/gc2741/ding/Drivers/Data/mc3.v0.2.8.PUBLIC.code.filtered.maf | grepList PCA.all.genes.txt 0 > mc3.v0.2.8.PUBLIC.code.filtered.PCAgenes.hotspot.maf
# my @mafcols = ( $mafcols{"Hugo_Symbol"},
# 					$mafcols{"Chromosome"},
# 					$mafcols{"Start_Position"},
# 					$mafcols{"End_Position"},
# 					$mafcols{"Variant_Classification"},
# 					$mafcols{"Reference_Allele"},
# 					$mafcols{"Tumor_Seq_Allele2"},
# 					$mafcols{"Tumor_Sample_Barcode"},
# 					$mafcols{$this->{"transcript_id_header"}},
#Hugo_Symbol     Chromosome      Start_Position  End_Position    Variant_Classification  Reference_Allele        Tumor_Seq_Allele2       Tumor_Sample_Barcode    HGVSp_Short     Transcript_ID

# get germline file
awkt '$13 == "missense_variant" {print $2,$9,$5,$10,"Missense_Mutation",$11,$12,$1,$116,$6}' PCA_pathVar_integrated_filtered.tsv > PCA_pathVar_integrated_filtered_hotspot3d.tsv

#2. combine the somatic file with germline file


cat mc3.v0.2.8.PUBLIC.code.filtered.PCAgenes.hotspot.maf PCA_pathVar_integrated_filtered_hotspot3d.tsv > PCA_somatic_germline_combined.maf
awk 'NR==1 || $9 ~ "p."' PCA_somatic_germline_combined.maf > PCA_somatic_germline_combined_filtered.maf
awk 'NR==1 || $5=="Missense_Mutation" {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$6"\t"$7"\t"$8"\t"$10"\t"$9}' PCA_somatic_germline_combined_filtered.maf > PCA_somatic_germline_combined_missense.maf

#3. hotspot run; following steps on the github
bsubl -oo proximity_search.log 'hotspot3d search --maf-file=PCA_somatic_germline_combined_missense.maf --prep-dir=/gscmnt/gc2706/dinglab/medseq/Structure_Projects/Preprocessing_Output_20141023/'

bsubl -oo post.log 'hotspot3d post --maf-file=PCA_somatic_germline_combined_missense.maf'

bsubl -oo cluster.log 'hotspot3d cluster --pairwise-file=3D_Proximity.pairwise --maf-file=PCA_somatic_germline_combined_missense.maf --vertex-type=recurrence'

bsubl -oo summary.log 'hotspot3d summary --clusters-file=PCA_somatic_germline_combined_missense.maf.3D_Proximity.pairwise.recurrence.l0.r10.clusters'

#4. post processing: plotting
# get presence
perl ~/bin/hotspot3d_KH/scripts/clusterPDBPresence.pl 3D_Proximity.pairwise PCA_somatic_germline_combined_missense.maf.3D_Proximity.pairwise.recurrence.l0.r10.clusters PCA

# plot some of the previous top candidates
musite=canonical.combined.mc3.musites
cluster=PCA_somatic_germline_combined_missense.maf.3D_Proximity.pairwise.recurrence.l0.r10.clusters

# limit cluster file to clusters of interest first

# visualize
grep1 245.0 ${cluster} > tmp.clusters
hotspot3d visual --pairwise-file=3D_Proximity.pairwise --clusters-file=tmp.clusters --pdb=2IVT --output-file=pml_scripts/PCA.2IVT.RET.pml --script-only 
grep1 245.0 ${cluster} > tmp.clusters
hotspot3d visual --pairwise-file=3D_Proximity.pairwise --clusters-file=tmp.clusters --pdb=2X2M --output-file=pml_scripts/PCA.2X2M.RET.pml --script-only 
