file=$1
gsutil cp $file* .
        vcfName=${file##*/}
        date
        echo "Examining vcf-stats for ${vcfName}"
        ../plinkseq/build/execs/pseq $vcfName v-stats > ${vcfName}.pseq.vstats.tsv
        ../plinkseq/build/execs/pseq $vcfName i-stats > ${vcfName}.pseq.istats.tsv
        #vcf-stats $vcfName > ${vcfName}.vcfstats.json

        rm -f $vcfName
        rm -f ${vcfName}.tbi
