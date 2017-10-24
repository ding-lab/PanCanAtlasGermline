#!/bin/python
#Oct 2017 - Kuan-Lin Huang @ WashU - 
 
import sys
import getopt
import gzip
from biomine.variant import variant

def vcf2charger( chrom , pos , ref , alt ):
	ensFields = []
	ensFields.append( chrom )
	ensFields.append( str( pos ) )
	#ensFields.append( str( stop ) )
	ref = variant.nullCheck( ref )
	alt = variant.nullCheck( alt )
	ensFields.append( str( ref ) )
	ensFields.append( str( alt ) )
	var = "_".join(ensFields)
	return var

def main():
    def usage():
        print """
    combine_CharGer2VCF.py : combine CharGer output to its originating VCF file

    USAGE: liftover_CharGer_result.py [-h] <CharGer results file> <VCF>
     -h    print this message
     <filename>    input file
        """

    #use getopt to get inputs
    # try:
    #     opts, args = getopt.getopt(sys.argv[1:], 'h') #:after option meaning required arguments
    # except getopt.GetoptError:
    #     print "liftover_CharGer_result.py <CharGerSummaryFile> <inputFile>"

    # for opt, arg in opts: #store the input options
    #     if opt == '-h': # h means user needs help
    #         usage(); sys.exit()

    if len(sys.argv) >= 3:
        CharGerFH= sys.argv[1]
        vcfFH = sys.argv[2]
    else:
        usage()
        sys.exit()

    # store pathogenic/likely pathogenic variants information

    #open CharGer file
    try:
        charGerF = open(CharGerFH,"r")
    except IOError:
        print("File , CharGerFH, does not exist!")
    
    CharGerHeader = charGerF.readline().strip()
    varCharGer = {}
    #read input file
    for line in charGerF:
        line=line.strip()
        F = line.split("\t")
        # chrom = F[1]
        # start = F[2]
        # stop = F[3]
        # ref = F[4]
        # alt = F[5]
        if len(F) > 20 and "Pathogenic" in F[19]:
            	# if F[4] == "-":
            	#     F[4] = "0"
            	# if F[5] == "-":
            	#     F[5] = "0"
            	var = "_".join(F[1:3]+F[4:6])
            	varCharGer[var]=line
		#print "added charger var: " + var
    charGerF.close()

    try:
        vcfF = open(vcfFH,"r")
    except IOError:
        print("VCF file , vcfFH, does not exist!") 

    vcfHeader = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    print vcfHeader + "\t" + "Sample" + "\t" + "Genotype" + "\t" + CharGerHeader
    samples = []
    
    #read input file
    for line in vcfF:
        line=line.strip()
        F = line.split("\t")

        # store sample order
        if line.startswith("#CHR"):
            for i in range(len(F)):
                samples.append(str(F[i]))
		continue
        elif line.startswith("#"):
            continue
	
	if len(F) > 4:
            chrom = F[0]
            pos = F[1]
            #stop = F[2]
            ref = F[3]
            alt = F[4]
            var = vcf2charger(chrom , pos , ref , alt)
	    print "VCF var: " + var
            CharGerAnno = ""
            # print only if it's in CharGer
            if var in varCharGer:
                CharGerAnno = varCharGer[var]
                # check samples have this var
                for i in range(9,len(F)) :
                    sample = samples[i]
		    genotype = F[i]
                    ##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  TCGA-OR-A5K0-10B-01D-A29L-10 
                    if genotype.startswith("./."): # not called in variant files we have
                        continue
                    else: 
                        vcfVar = "\t".join(F[0:9])
                        #sample = samples[i]
                        #print vcfVar + "\t" + sample + "\t" + genotype + "\t" + CharGerAnno  
    
    vcfF.close()

if __name__ == "__main__":
    main()
