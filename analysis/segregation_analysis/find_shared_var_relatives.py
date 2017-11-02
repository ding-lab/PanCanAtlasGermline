#!/bin/python
#Aug 2017 - Kuan-Lin Huang @ WashU - 
# find_shared variants in relatives identified in TCGA samples

import sys
import getopt
import gzip

def main():
    def usage():
        print """
    find_shared_var_relatives.py : why do I exist?

    USAGE: find_shared_var_relatives.py [-h] <file with relative pairs> <vcf.gz file to be scanned>
     -h    print this message
     <filename>    input file
        """

    #use getopt to get inputs
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'h') #:after option meaning required arguments
    except getopt.GetoptError:
        print "find_shared_var_relatives.py <file with relative pairs> <vcf.gz file to be scanned>"

    for opt, arg in opts: #store the input options
        if opt == '-h': # h means user needs help
            usage(); sys.exit()

    args = sys.argv[1:]
    if len(args) < 1:
        usage(); sys.exit("input file missing")

    #open input file
    try:
        fn = args[0]
        relativeF = open(fn ,"r")
    except IOError:
        print("File , args[0], does not exist!")

    # dictionaries
    relative2relative = {}
    existing_relatives = {}
    sample2colID = {}

    # output file
    outFstring = args[0] + "." + args[1]
    #outF = outFstring.replace("vcf","segragatingVar.tsv")
    outF = outFstring.replace("vcf.gz","segragatingVar.tsv")
    outFH = open(outF, "w")

    #read input file
    for line in relativeF:
        line=line.strip()
        F = line.split("\t")
        #print str(len(F)) + "\n"
        if len(F)==2:
            print "Searching for segregating variants in relative pairs: "+ F[0] + ":" + F[1] + "\n"
            relative2relative[F[0]] = F[1]
    relativeF.close()


    try:
        fn = args[1]
        if fn.endswith(".gz"):
            vcfF = gzip.open(fn,"r")
        elif fn.endswith(".vcf"): 
            vcfF = open(fn,"r")
    except IOError:
        print("File , args[1], does not exist or is not a valid vcf!")

    colnames = "sample1\tsample2\tsample1GENO\tsample2GENO\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"
    outFH.write(colnames)
    #read input file
    for line in vcfF:
        line=line.strip()

        # print headers
        if line.startswith("##"):
            #print line
            continue      
        
        F = line.split("\t")
        colNum = len(F)
        
        if line.startswith("#CHR"): # row with column names
            #print line          
            for i in range(0,colNum):
                TCGA_barcode = F[i][0:12]
                sample2colID[TCGA_barcode] = i
            for sample in relative2relative:
                relative = relative2relative[sample]
                # only keep the existing pairs to save iteration time later on
                if sample in sample2colID and relative in sample2colID:
                    existing_relatives[sample] = relative

        else:
            INFO = F[7]
            INFOsplit = INFO.split(";AC=")
            AC = int( INFOsplit[1] )
            if AC > 10:
                continue
            # loop through the relatives; if one has the var, check the other; only print if both has it
            for sample in existing_relatives:
                relative = relative2relative[sample]
                sampleColID = sample2colID[sample]
                relativeColID = sample2colID[relative]
                sampleGeno = F[sampleColID]
                relativeGeno = F[relativeColID]
                # print line if segragating var found! ( BOTH NOT WT )
                if not sampleGeno.startswith("./.") and not relativeGeno.startswith("./."): 
                    outLine = sample + "\t" + relative + "\t" + sampleGeno + "\t" + relativeGeno + "\t".join(F[0:8]) + "\n"
                    outFH.write(outLine)
                    #print outLine
                    
    
    vcfF.close()

if __name__ == "__main__":
    main()
