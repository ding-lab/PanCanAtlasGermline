#!/bin/python
#03 February 2016 - Kuan-Lin Huang @ WashU - 
import os
import sys
import getopt
import gzip
import scipy.stats as stats

def main():
    def usage():
        print """
    single_association.py : why do I exist?

    USAGE: single_association.py [-h] <vcf file with vcfanno-tated ExAC freq> <germline file> <hotspot3d file> # get frequency to variants
     -h    print this message
     <filename>    input file
        """

    args = sys.argv[1:]
    # if len(args) < 2:
    #     usage(); sys.exit("input file missing")

    

    #open ExAC vcf
    try:
        exacF = gzip.open(args[0],"r")
    except IOError:
        print("ExAC_VCF_File: , args[0], does not exist!")

    foutput = args[1]
    fout = open(foutput, 'w')

    header_string = "Var\tExAC_AN\tExAC_AC\tExAC_nonTCGA_AN\tExAC_nonTCGA_AC\tTCGA_AN\tTCGA_AC\tOR\tP\tgene_symbol\timpact\tHGVSp\n"   
    fout.write(header_string)

    for line in exacF:
        if line.startswith("#"):
            continue   

        line=line.strip()
        F = line.split("\t")

        if len(F) < 7:
            print "weird line" + line
            continue   

        var = ":".join(F[0:5])
        
        info = F[7]
        

        ExAC_AC = 0
        ExAC_AN = 0
        ExAC_AF = 0
        ExAC_nonTCGA_AC = 0
        ExAC_nonTCGA_AN = 0
        ExAC_nonTCGA_AF = 0
        gene_symbol = "NA"
        impact = "NA"
        HGVSp = "NA"

        infos = info.split(";")
        for item in infos:
            if item.startswith("AC_Adj="):
                ACstr = item.replace("AC_Adj=","")
                ExAC_AC = int(ACstr)

            elif item.startswith("AN_Adj="):
                ANstr = item.replace("AN_Adj=","")
                ExAC_AN = int(ANstr)

            elif item.startswith("AF_Adj="):
                AFstr = item.replace("AF_Adj=","")
                ExAC_AF = int(AFstr)

            elif item.startswith("ExAC_nonTCGA_AC_Adj="):
                AC_nonTCGAstr = item.replace("ExAC_nonTCGA_AC_Adj=","")
                ExAC_nonTCGA_AC = int(AC_nonTCGAstr)

            elif item.startswith("ExAC_nonTCGA_AN_Adj="):
                AN_nonTCGAstr = item.replace("ExAC_nonTCGA_AN_Adj=","")
                ExAC_nonTCGA_AN = int(AN_nonTCGAstr)

            elif item.startswith("ExAC_AF_nonTCGA_Adj="):
                AF_nonTCGAstr = item.replace("ExAC_AF_nonTCGA_Adj=","")
                ExAC_nonTCGA_AF = int(AF_nonTCGAstr)

            elif item.startswith("CSQ="):
                CSQstr = item.replace("CSQ=","")
                ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. 
                ## Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids
                CSQ = CSQstr.split("|")
                gene_symbol = CSQ[3]
                consequence = CSQ[1]
                HGVSp = CSQ[11]

      # from the nonTCGA subset readme:
      # 1 This directory contains two commonly requested subsets from ExAC release 1 (60,706 samples) 121412
      # 2 ExAC.r1.nonTCGA.sites.vep.vcf.gz - Excluding TCGA cohorts (53,105 samples) 106210 nonTCGA; 15202 TCGA
      # 3 ExAC.r1.nonpsych.sites.vcf.gz - Excluding Psychiatric cohorts (45,376 samples)

      # Mode over first 1 mil lines; if the variant is not present in nonTCGA vcf: for nonTCGA vcf: AN: 106210
      # we gotta capture variants found in ExAC VCF but not nonTCGA VCF
      # assume the largest TCGA AN if the variant is not seen in nonTCGA VCF
        if ExAC_AN > 15204 and ExAC_nonTCGA_AN == 0:
            ExAC_nonTCGA_AN = ExAC_AN - 15202 

        TCGA_AN = int(ExAC_AN) - int(ExAC_nonTCGA_AN)
        TCGA_AC = int(ExAC_AC) - int(ExAC_nonTCGA_AC)
        if TCGA_AN > 16000: # if for some reason it's missing in the nonTCGA ExAC maf; likely because of non-detection
            TCGA_AN = 15202 # maximum count of other alleles (1634)

        if TCGA_AN > 0:
            TCGA_AF = TCGA_AC/TCGA_AN
        
        [odds, testStat] = ["NA","NA"]
        try:
            [odds, testStat] = stats.fisher_exact([[TCGA_AC, TCGA_AN], [ExAC_nonTCGA_AC, ExAC_nonTCGA_AN]])
        except ValueError:
            pass
        
        out_string = var + "\t" + str(ExAC_AN)+ "\t" + str(ExAC_AC) + "\t"+ str(ExAC_nonTCGA_AN) + "\t" + \
            str(ExAC_nonTCGA_AC) + "\t"+ str(TCGA_AN) + "\t" + str(TCGA_AC) + "\t" + str(odds) + "\t" + str(testStat) + "\t" + \
            gene_symbol + "\t" + consequence + "\t" + HGVSp + "\n"
        
        fout.write(out_string)



    # ExAC_AC=2;ExAC_AN=121412;ExAC_AF=1.647e-05;ExAC_AC_Adj=2;ExAC_AN_Adj=121114;ExAC_AF_Adj=1.6513e-05;ExAC_nonTCGA_AC=1;ExAC_nonTCGA_AN=106210;ExAC_nonTCGA_AF=9.415e-06;ExAC_nonTCGA_AC_Adj=1;ExAC_nonTCGA_AN_Adj=106060;ExAC_nonTCGA_AF_Adj=9.4153e-06
    exacF.close()
    fout.close()


if __name__ == "__main__":
    main()
