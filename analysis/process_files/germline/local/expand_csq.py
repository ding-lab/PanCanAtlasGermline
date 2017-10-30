#!/bin/python
#Oct 2017 - Kuan-Lin Huang @ WashU - 
 
import sys
import getopt
import gzip


def main():
    def usage():
        print """
    combine_CharGer2VCF.py : combine CharGer output to its originating VCF file

    USAGE: liftover_CharGer_result.py [-h] <CharGer results file> <VCF>
     -h    print this message
     <filename>    input file
        """

    if len(sys.argv) >= 3:
	vcfheadFH = sys.argv[1]
        CharGerFH= sys.argv[2]
    else:
        usage()
        sys.exit()
    
    try:
        vcfheadF = open(vcfheadFH,"r")
    except IOError:
        print("File , CharGerFH, does not exist!")
    csq_header = ""
    for line in vcfheadF:
	line = line.strip()
	if line.startswith("##INFO=<ID=CSQ"):
		F = line.split("|")
		F[0] = "Allele"
		csq_header = "\t".join(F)
    
    #open CharGer file
    try:
        charGerF = open(CharGerFH,"r")
    except IOError:
        print("File , CharGerFH, does not exist!")
    
    CharGerHeader = charGerF.readline().strip()
    print CharGerHeader + "\t" + csq_header
    #read input file
    for line in charGerF:
        line=line.strip()
        F = line.split("\t")
	csq = F[32].split(",")
	csqF = csq[0].split("|")
	for i in range(len(csqF)):
		if csqF[i] == "":
			csqF[i] = "NA"
	csqFields = "\t".join(csqF)
	
	print line + "\t" + csqFields

    charGerF.close()

if __name__ == "__main__":
    main()
