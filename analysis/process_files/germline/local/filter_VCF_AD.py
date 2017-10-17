#!/bin/python
#03 February 2016 - Kuan-Lin Huang @ WashU - 
 
import sys
import getopt
import gzip

class autovivification(dict):
    '''Implementation of perl's autovivification feature.'''
    def __init__( self , *args , **kwargs ):
        super( autovivification , self ).__init__( *args , **kwargs )
        self.itemlist = super( autovivification , self ).keys()
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

def main():
    def usage():
        print """
    filter_AD_VCF.py : why do I exist?

    USAGE: filter_AD_VCF.py  [-h] <VCF filename> <Allelic Depth (AD) threshold> 
     -h    print this message
     <VCF filename>    input file
        """

    if len(sys.argv) == 3:
        vcfFH= sys.argv[1]
        AD_thres = int(sys.argv[2])
    else:
        usage()
        sys.exit()

    try:
        vcfF = gzip.open(vcfFH,"r")
    except IOError:
        print("VCF file does not exist!")  


    for line in vcfF:
        line=line.strip()
        # print the info lines
        if line.startswith("#"):
            print line
        else:
            F = line.split("\t")

            ref = str(F[3])
            info_f = str(F[7]).split(";")
            format_f = str(F[8]).split(":")
            geno_f = str(F[9]).split(":")
            AD_index = -1

            ### reference filter
            nonpass_ref = False
            if ref == "N":
                nonpass_ref = True

            ### AD filter
            nonpass_AD = False
            for i in range(0,len(format_f)):
                if str(format_f[i]) == "AD":
                    AD_index = i

            genotype = str(geno_f[AD_index])
            # GATK and Pindel calls
            # second int for alt allele
            if "," in genotype: 
                genotypes = genotype.split(",")
                if int(genotypes[1]) < AD_thres:
                    nonpass_AD = True
            # varscan calls
            else:
                if int(genotype) < AD_thres:
                    nonpass_AD = True
            
            if not nonpass_ref and not nonpass_AD:
                #outFH.write(line + "\n")
                print line


if __name__ == "__main__":
    main()
