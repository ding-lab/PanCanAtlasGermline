#!/bin/python
#Aug 2017 - Kuan-Lin Huang @ WashU - 
 
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

    USAGE: filter_AD_VCF.py  [-h] <VCF filename> <ExAC Allelic Frequency (AF) threshold> <Allelic Depth (AD) threshold> 
     -h    print this message
     <VCF filename>    input file
        """

    if len(sys.argv) == 4:
        vcfFH= sys.argv[1]
        MAF_thres = float(sys.argv[2])
        AD_thres = int(sys.argv[3])
    else:
        usage()
        sys.exit()

    try:
        vcfF = open(vcfFH,"r")
    except IOError:
        print("VCF file does not exist!")  

    outFstring = "ExAC_AF." + str(MAF_thres) + ".ExAConly.AD." + str(AD_thres) + ".vcf"
    outF = vcfFH.replace("vcf",outFstring)
    outFH = open(outF, "w")

    all_var = 0
    nonpass_MAF_var = 0
    nonpass_AD_var = 0 
    pass_var = 0 

    for line in vcfF:
        line=line.strip()
        # print the info lines
        if line.startswith("#"):
            outFH.write(line + "\n")
        else:
            F = line.split("\t")
            all_var = all_var + 1

            info_f = str(F[7]).split(";")
            format_f = str(F[8]).split(":")
            geno_f = str(F[9]).split(":")
            AD_index = -1

            # only keep rare, ExAC variants
            ### MAF filter
            nonpass_MAF = True
            for info in info_f:
                # find the cases with annotated ExAC frequency
                if info.startswith("ExAC_AF_Adj"):
                    ExAC_AF = info.replace("ExAC_AF_Adj=","")
                    if "," in ExAC_AF:
                        ExAC_AFs = ExAC_AF.split(",")
                        if float(ExAC_AFs[0]) < MAF_thres: # need to assume it's the first allele
                            nonpass_MAF = False
                    else:
                        if float(ExAC_AF) < MAF_thres:
                            nonpass_MAF = False
            if nonpass_MAF:
                nonpass_MAF_var = nonpass_MAF_var + 1
                continue


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
                    nonpass_AD_var = nonpass_AD_var + 1
            # varscan calls
            else:
                if int(genotype) < AD_thres:
                    nonpass_AD = True
                    nonpass_AD_var = nonpass_AD_var + 1
            
            if not nonpass_MAF and not nonpass_AD:
                pass_var = pass_var + 1
                outFH.write(line + "\n")

    # filter summary
    print "number of total variants:", all_var
    print "number of variants failing MAF filter of", MAF_thres,":", nonpass_MAF_var
    print "number of variants failing AD filter of", AD_thres,":", nonpass_AD_var
    print "number of total passed variants:", pass_var
            # if len(genotypes) == 1 & int(genotypes[0])>=AD_thres:
            #     print line
            # elif len(genotypes) == 2 & int(genotypes[1])>=AD_thres:
            #     print line
    outFH.close()

if __name__ == "__main__":
    main()
