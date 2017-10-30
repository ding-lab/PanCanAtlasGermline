#!/bin/bash

#(ads&rjm) 2016-09-22

# use new perl
. /gscmnt/gc2525/dinglab/rmashl/Software/perl/set_envvars
#which ${PERL_BIN}
#exit
#ads_vep="/gscmnt/gc2706/dinglab/medseq/LabCode/AdamDS/ensembl-vep/vep"
#ads_cachevep="/gscmnt/gc2706/dinglab/medseq/LabCode/AdamDS/VEP/.vep/"
vep_cmd="/gscmnt/gc2525/dinglab/rmashl/Software/perl/perl-5.22.0/bin/perl /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v85/ensembl-tools-release-85/scripts/variant_effect_predictor/variant_effect_predictor.pl"
cachedir="/gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v85/cache/"
reffasta="/gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v85/cache/homo_sapiens/85_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
assembly="GRCh37"
#opts="--plugin ExAC,/gscmnt/gc2706/dinglab/medseq/ExAC/VCF/ExAC.r0.3.1.sites.vep.vcf.gz" #"--everything"
opts="--everything"
results="new_run/AnnotatedVCFs/"
variantBuffer=2000
queue="long"
group="/khuang"
forks=4
mem=20000000
export SAMTOOLSDIR="/gscmnt/gc2525/dinglab/rmashl/Software/bin/samtools/1.2/bin"
expoSAMTOOLS="$SAMTOOLSDIR/samtools"

inputDir="/gscmnt/gc3020/dinglab/medseq/Germline/projects/PanCanAtlasGermline/analysis/process_files/germline/local/VCF"

if [ ! -d ${results} ]; then
	mkdir ${results}
fi

for i in {1..22} X Y
do
	vcf=PCA.r1.TCGAbarcode.merge.exon.chr${i}.vcf.gz
	out=PCA.r1.TCGAbarcode.merge.exon.chr${i}.vcf
	runVEP="perl format.pl $inputDir $vcf; ${vep_cmd} ${opts} --offline --cache --dir ${cachedir} --assembly ${assembly} --format vcf --vcf -i new_run/preVEP/anno.$out -o ${results}anno.${out} --force_overwrite --fasta ${reffasta} --fork ${forks} --buffer_size ${variantBuffer};"
	log="${results}anno.${vcf}.log"
	echo "bsub -g ${group} -q ${queue} -n ${forks} -M ${mem} -oo ${log} \"${runVEP}\""
done

