# data location: /gscmnt/gc3014/dinglab/ExAC/VCF

# normalize exac vcf and subset to only rare variants
bsubl -oo ExAC.r1.multi2biallelic.log '~/bin/bcftools-1.5/bcftools norm -m - ExAC.r1.sites.vep.vcf.gz | bgzip -c > ExAC.r1.sites.vep.biallelic.vcf.gz'
bsubl -oo ExAC_nonTCGA.r1.multi2biallelic.log '~/bin/bcftools-1.5/bcftools norm -m - ExAC_nonTCGA.r1.sites.vep.vcf.gz | bgzip -c > ExAC_nonTCGA.r1.sites.vep.biallelic.vcf.gz'
bsubl -oo tabix.ExAC.r1.sites.vep.biallelic.log 'tabix -p vcf ExAC.r1.sites.vep.biallelic.vcf.gz'
bsubl -oo tabix.ExAC_nonTCGA.r1.sites.vep.biallelic.log 'tabix -p vcf ExAC_nonTCGA.r1.sites.vep.biallelic.vcf.gz'
# annotate nonTCGA frequency using vcfanno
bsub -q bigmem -R"select[mem>80000] rusage[mem=80000]" -M 80000000 -oo ExAC.r1.vcfanno.log '~/bin/vcfanno_linux64.1 -p 8 ExAC_nonTCGA_config.toml ExAC.r1.sites.vep.biallelic.vcf.gz | bgzip -c > ExAC.r1.sites.vep.biallelic.combine.vcf.gz'
bsubl -oo tabix.ExAC.r1.sites.vep.biallelic.combine.log 'tabix -p vcf ExAC.r1.sites.vep.biallelic.combine.vcf.gz'
# get rare variants
bsubl -oo ExAC.r1.getRareAllele.log '~/bin/bcftools-1.5/bcftools view --max-af 0.0005 ExAC.r1.sites.vep.biallelic.combine.vcf.gz | bgzip -c > ExAC.r1.sites.vep.biallelic.combine.maxAF0.01.vcf.gz'
mv ExAC.r1.sites.vep.biallelic.combine.maxAF0.01.vcf.gz ExAC.r1.sites.vep.biallelic.combine.maxAF0.0005.vcf.gz
# Each ethnicity only sums up to AN_Adj; not AN

# his workflow:
# use vcf anno to combine nonTCGA AC/AN to ExAC vcf
vcfanno -p 8 ExAC_nonTCGA_config.toml ExAC.r0.3.1.sites.vep.vcf.gz | bgzip -c > ExAC.r0.3.1.combine_all.vep.vcf.gz

# conduct single variant assoc test
bsubl -oo single_var_association.log 'python2.7 single_var_association.py VCF/ExAC.r1.sites.vep.biallelic.combine.vcf.gz ExAC.r1.sites.vep.biallelic.combine.fisher.anno.tsv'

### from here:
# filter out the non-rare variants; convert to tab delimited
# may need to limit only to regions with sufficient coverage
# Do burden test on variants with specific consequence and AF bins

gzcat ExAC.r1.sites.vep.biallelic.combine.fisher.anno.v2.tsv.gz | grepList /Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/reference_files/20160713_Rahman_KJ_KH_152_gene_table_list.txt 9 > ExAC.r1.sites.vep.biallelic.combine.fisher.anno.v2.152gene.tsv
