#0. run charger
bash run_charger.sh

#1. append results
python liftover_CharGer_result_PCGP.py charged.2015_stJude_germline_nejm_S4_AD_varOnly.vep.tsv /Users/khuang/Box\ Sync/PhD/germline/pan8000_germline_clinical/variant_files/201604_PCGP_variants/2015_stJude_germline_nejm_S4_AD.txt > 2015_stJude_germline_nejm_S4_AD_charger.txt
python liftover_CharGer_result_PCGP.py charged.2015_stJude_germline_nejm_S4_AR_varOnly.vep.tsv /Users/khuang/Box\ Sync/PhD/germline/pan8000_germline_clinical/variant_files/201604_PCGP_variants/2015_stJude_germline_nejm_S4_AR.txt > 2015_stJude_germline_nejm_S4_AR_charger.txt

#2. Run analysis and plotting
Rscript plot_CharGer_summary_PCGP.R
