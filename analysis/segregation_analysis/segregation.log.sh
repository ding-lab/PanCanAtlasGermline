# 201708 Kuan-lin Huang @ WashU

# transfer needed files to VM
gcloud compute scp find_shared_var_relatives.py huangkuanlin@kuan-8cpu-30gb:~/ --zone us-central1-c
gcloud compute scp out/TCGA_z1_z2_relatives.tsv huangkuanlin@kuan-8cpu-30gb:~/ --zone us-central1-c
gcloud compute scp *.sh huangkuanlin@kuan-8cpu-30gb:~/ --zone us-central1-c

# run script to find segregating variants