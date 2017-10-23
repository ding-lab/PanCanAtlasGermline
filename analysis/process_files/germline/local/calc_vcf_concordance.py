#!/bin/python
# 2017 - Kuan-Lin Huang @ WashU -
# some unused functions adopted from dream challenge evaluation code (see below)
# mostly used to evaluate against genotype file
 
import sys, os
import getopt
import gzip
import vcf

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

#!/usr/bin/env python


# '''
# Submission evaluation code for TCGA/ICGC/DREAM SMC
# Adam Ewing, ewingad@soe.ucsc.edu
# Requires PyVCF (https://github.com/jamescasbon/PyVCF)
# '''

# def match(subrec, trurec, vtype='SNV'):
#     assert vtype in ('SNV', 'SV', 'INDEL')

#     if vtype == 'SNV' and subrec.is_snp and trurec.is_snp:
#         if subrec.POS == trurec.POS and subrec.REF == trurec.REF and subrec.ALT == trurec.ALT:
#             return True

#     if vtype == 'INDEL' and subrec.is_indel and trurec.is_indel:
#         if subrec.POS == trurec.POS and subrec.REF == trurec.REF and subrec.ALT == trurec.ALT:
#             return True

#     if vtype == 'SV' and subrec.is_sv and trurec.is_sv:
#         trustart, truend = expand_sv_ends(trurec, useCIs=False)
#         substart, subend = expand_sv_ends(subrec, useCIs=False)

#         # check for overlap
#         if min(truend, subend) - max(trustart, substart) > 0:
#             return True

#     return False


# def expand_sv_ends(rec, useCIs=True):
#     '''
#     assign start and end positions to SV calls
#     using conf. intervals if present and useCIs=True
#     '''
#     startpos, endpos = rec.start, rec.end
#     assert rec.is_sv

#     try:
#         if rec.INFO.get('END'): # sometimes this is a list, sometimes it's an int
#             if isinstance(rec.INFO.get('END'), list):
#                 endpos = int(rec.INFO.get('END')[0])
#             if isinstance(rec.INFO.get('END'), int):
#                 endpos = int(rec.INFO.get('END'))

#         if useCIs:
#             if rec.INFO.get('CIPOS'):
#                 ci = map(int, rec.INFO.get('CIPOS'))
#                 if ci[0] < 0:
#                     startpos += ci[0]

#             if rec.INFO.get('CIEND'):
#                 ci = map(int, rec.INFO.get('CIEND')) 
#                 if ci[0] > 0:
#                     endpos += ci[0]

#     except TypeError as e:
#         sys.stderr.write("error expanding sv interval: " + str(e) + " for record: " + str(rec) + "\n")

#     if startpos > endpos:
#         endpos, startpos = startpos, endpos

#     return startpos, endpos
    
    
def evaluate(submission, truth, sampleMatch=False):
    ''' return TP, FP and FN counts '''

    subvcfh = vcf.Reader(filename=submission)
    truvcfh = vcf.Reader(filename=truth)

    outF = submission.replace("vcf.gz","validated.tsv")
    outFH = open(outF, "w")

    truGeno = autovivification()
    sample_validated = {}
    sample_unvalidated = autovivification()
    validated = 0
    unvalidated = 0
    unvalidated_by_type = {}
    mistmatch_types = ["CALLED_0/1;GENOTYPE_0/0" ,"CALLED_0/1;GENOTYPE_1/1", "CALLED_1/1;GENOTYPE_0/0","CALLED_1/1;GENOTYPE_0/1"]

    ''' store list of truth records, otherwise the iterator needs to be reset '''
    for trurec in truvcfh:
        var = str(trurec.CHROM) + ":" + str(trurec.POS) + ":" + str(trurec.REF)
        alt = str(trurec.ALT)
        #print var

        for sample in trurec.samples:
            # get the sample name
            called_name = (str(sample).split(' ')[0]).split('=')[1].strip(',')
            # TCGA barcode name: may extract more down to sample?
            sampleID = called_name[0:12]
            #print sampleID
            genotype = format(sample['GT'])
            #print genotype
            truGeno[var][alt][sampleID] = genotype
            # ++trurecs 
            # #print var + ":" + str(genotype)
            # if genotype == "0/0": #WT
            #     ++truWtCount
    
    for subrec in subvcfh:
        var = str(subrec.CHROM) + ":" + str(subrec.POS) + ":" + str(subrec.REF)
        alt = str(subrec.ALT)
        #print var
        if not var in truGeno:
            continue
        if not alt in truGeno[var]:
            continue

        for sample in subrec.samples:
            called_name = (str(sample).split(' ')[0]).split('=')[1].strip(',')
            sampleID = called_name[0:12]

            if sampleID not in truGeno[var][alt]: # not in the GENOTYPE file
                continue

            genotype = format(sample['GT'])
            if genotype == "./.": # not called in variant files we have
                continue

            truGenotype = truGeno[var][alt][sampleID]

            if genotype == truGenotype:
                validated = validated + 1
                if sampleID in sample_validated:
                    sample_validated[sampleID] = sample_validated[sampleID] + 1
                else: 
                    sample_validated[sampleID] = 1
                #print "Matched: " + sampleID + " CALLED: " + genotype + " vs. GENOTYPE: " + truGenotype
            #elif:
            else:
                
                unvalidated = unvalidated + 1

                mismatch_GT = "CALLED_" + genotype + ";GENOTYPE_" + truGenotype
                if mismatch_GT in unvalidated_by_type:
                    unvalidated_by_type[mismatch_GT] = unvalidated_by_type[mismatch_GT] + 1
                else: 
                    unvalidated_by_type[mismatch_GT] = 1
                
                if sampleID in sample_unvalidated:
                    if mismatch_GT in sample_unvalidated[sampleID]:
                        sample_unvalidated[sampleID][mismatch_GT] = sample_unvalidated[sampleID][mismatch_GT] + 1
                    else:
                        sample_unvalidated[sampleID][mismatch_GT] = 1
                else: 
                    sample_unvalidated[sampleID][mismatch_GT] = 1
                #print "Mistmatched: " + sampleID + " CALLED: " + genotype + " vs. GENOTYPE: " + truGenotype

    outFH.write("Sample\tValidated_count\tCALLED_0/1;GENOTYPE_0/0\tCALLED_0/1;GENOTYPE_1/1\tCALLED_1/1;GENOTYPE_0/0\tCALLED_1/1;GENOTYPE_0/1\n")
    # iterate through samples and print their validated; unvalidated; chr
    for sampleID in sample_validated:
            outString = sampleID + "\t" + str(sample_validated[sampleID])
            outFH.write(outString)  
            for mistmatch_type in mistmatch_types:
                mismatch_string = "\t"  
                unvalidated_count = str(0)
                if sampleID in sample_unvalidated and mistmatch_type in sample_unvalidated[sampleID]:
                    unvalidated_count = str(sample_unvalidated[sampleID][mismatch_type])
                mismatch_string = mismatch_string + unvalidated_count
                outFH.write(mismatch_string)

            outFH.write("\n")


    # for mismatch_GT in unvalidated_by_type:
    #     print mismatch_GT + '\t' +  str(unvalidated_by_type[mismatch_GT]) + '\t' +  outFstring

    outFH.close()

    

if __name__ == '__main__':
    if len(sys.argv) == 3:
        subvcf, truvcf = sys.argv[1:3]

        if not subvcf.endswith('.vcf') and not subvcf.endswith('.vcf.gz'):
            sys.stderr.write("submission VCF filename does not enc in .vcf or .vcf.gz\n")
            sys.exit(1)

        if not os.path.exists(truvcf + '.tbi'):
            sys.stderr.write("truth VCF does not appear to be indexed. bgzip + tabix index required.\n")
            sys.exit(1)

        evaluate(subvcf, truvcf)

    else:
        print "standalone usage for testing:", sys.argv[0], "<submission VCF> <truth VCF (tabix-indexed)>"