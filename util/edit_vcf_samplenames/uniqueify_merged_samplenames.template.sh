#!/bin/bash
# Jay Mashl, July 2017

SOURCE=variants.vcf.gz
DEST=${SOURCE/variants/variants.fixedHeader}

gunzip -dc $SOURCE  | ./replace_vcf_header_sample_with_source.pl   > newheader.txt
tabix -r newheader.txt   $SOURCE  >  $DEST
tabix -p vcf $DEST
