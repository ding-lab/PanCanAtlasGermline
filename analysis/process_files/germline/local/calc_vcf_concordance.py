#!/bin/python
#03 February 2016 - Kuan-Lin Huang @ WashU - 
 
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


'''
Submission evaluation code for TCGA/ICGC/DREAM SMC
Adam Ewing, ewingad@soe.ucsc.edu
Requires PyVCF (https://github.com/jamescasbon/PyVCF)
'''

def match(subrec, trurec, vtype='SNV'):
    assert vtype in ('SNV', 'SV', 'INDEL')

    if vtype == 'SNV' and subrec.is_snp and trurec.is_snp:
        if subrec.POS == trurec.POS and subrec.REF == trurec.REF and subrec.ALT == trurec.ALT:
            return True

    if vtype == 'INDEL' and subrec.is_indel and trurec.is_indel:
        if subrec.POS == trurec.POS and subrec.REF == trurec.REF and subrec.ALT == trurec.ALT:
            return True

    if vtype == 'SV' and subrec.is_sv and trurec.is_sv:
        trustart, truend = expand_sv_ends(trurec, useCIs=False)
        substart, subend = expand_sv_ends(subrec, useCIs=False)

        # check for overlap
        if min(truend, subend) - max(trustart, substart) > 0:
            return True

    return False


def expand_sv_ends(rec, useCIs=True):
    '''
    assign start and end positions to SV calls
    using conf. intervals if present and useCIs=True
    '''
    startpos, endpos = rec.start, rec.end
    assert rec.is_sv

    try:
        if rec.INFO.get('END'): # sometimes this is a list, sometimes it's an int
            if isinstance(rec.INFO.get('END'), list):
                endpos = int(rec.INFO.get('END')[0])
            if isinstance(rec.INFO.get('END'), int):
                endpos = int(rec.INFO.get('END'))

        if useCIs:
            if rec.INFO.get('CIPOS'):
                ci = map(int, rec.INFO.get('CIPOS'))
                if ci[0] < 0:
                    startpos += ci[0]

            if rec.INFO.get('CIEND'):
                ci = map(int, rec.INFO.get('CIEND')) 
                if ci[0] > 0:
                    endpos += ci[0]

    except TypeError as e:
        sys.stderr.write("error expanding sv interval: " + str(e) + " for record: " + str(rec) + "\n")

    if startpos > endpos:
        endpos, startpos = startpos, endpos

    return startpos, endpos


def passfilter(rec):
    ''' Return true if a record is unfiltered or has 'PASS' in the filter field (pyvcf sets FILTER to None) '''
    if rec.FILTER is None or rec.FILTER == '.' or not rec.FILTER:
        return True
    return False



def countrecs(result):
    ''' return number of counted mutation calls in submission '''
    assert 'tp' in result and 'fp' in result, "invalid result dictionary!"
    
    ncalls = result['tp'] + result['fp']
    
    return ncalls


def prefix(rec, usechr):
    ''' adjust presence/absence of "chr" prefix according to whether usechr is True or False '''
    if usechr and not rec.CHROM.startswith('chr'):
        rec.CHROM = 'chr' + rec.CHROM
    
    if not usechr and rec.CHROM.startswith('chr'):
        rec.CHROM = rec.CHROM.replace('chr', '')
    
    return rec
    
    
def evaluate(submission, truth, sampleMatch=False):
    ''' return TP, FP and FN counts '''

    subvcfh = vcf.Reader(filename=submission)
    truvcfh = vcf.Reader(filename=truth)
    truWtCount = 0

    tpcount = 0 
    fpcount = 0 
    tncount = 0
    fncount = 0

    trurecs = 0
    subrecs = 0 # counts all predicted variants

    truGeno = autovivification()
    subGeno = autovivification()

    ''' store list of truth records, otherwise the iterator needs to be reset '''
    
    subSample = subvcfh.samples[0][5:16]
    truSample = "none"
    for sample in truvcfh.samples:
        if subSample in sample:
            print "Sample matched! Submission sample: " + subvcfh.samples[0] + " Genotype sample: " + sample
            truSample = sample
    if truSample == "none":
        sys.exit('No matched sample in genotype file.')

    for trurec in truvcfh:
        var = str(trurec.CHROM) + ":" + str(trurec.POS) + ":" + str(trurec.REF)
        alt = str(trurec.ALT)
        #print var

        genotype = format(trurec.genotype(truSample)['GT'])
        truGeno[var][alt] = genotype
        trurecs += 1

        # for sample in trurec.samples:
        #     if sample in # check sample match 
        #     genotype = format(sample['GT'])
        #     #print genotype
        #     truGeno[var][alt] = genotype
        #     trurecs += 1
    
    for subrec in subvcfh:
        var = str(subrec.CHROM) + ":" + str(subrec.POS) + ":" + str(subrec.REF) 
        alt = str(subrec.ALT)
        #print var

        for sample in subrec.samples:
            genotype = format(sample['GT'])
            #print genotype
            subGeno[var][alt] = genotype
            subrecs += 1

    for var in truGeno:
        #print var
        for alt in truGeno[var]:
            genotype = truGeno[var][alt]
            #print genotype
            if genotype == "0/0":
                if var in subGeno:
                    if alt in subGeno[var]:
                        if subGeno[var][alt] == genotype:
                            tncount += 1
                        else:
                            fpcount += 1
                    else: 
                        fpcount += 1
                else:
                    tncount += 1

            elif (genotype == "0/1") or (genotype == "1/0") or (genotype == "1/1"):
                if var in subGeno:
                    if alt in subGeno[var]:
                        tpcount += 1
                        #if subGeno[var][alt] == genotype:
                else:
                    fncount += 1

    print "# of variants in genotype file: " + str(trurecs)
    print "# of variants in variant file: " + str(subrecs)

    print "TP: " + str(tpcount)
    print "TN: " + str(tncount)
    print "FP: " + str(fpcount)
    print "FN: " + str(fncount)
    
    result = { 'tp' : float(tpcount),
               'fp' : float(fpcount),
               'tn' : float(tncount),
               'fn' : float(fncount),
               'sample' : str(subSample) }
    
    return result


def stats(result):
    ''' calculate precision, recall, fscore  from result dictionary '''
    assert 'tp' in result and 'fp' in result and 'fn' in result and 'tn' in result and 'sample' in result, "invalid result dictionary!"
    sample = result['sample']

    sensitivity = float(1)
    if result['tp'] + result['fp'] > 0:
        sensitivity = result['tp'] / (result['tp'] + result['fp'])
    
    specificity = float(1)
    if result['tn'] + result['fp'] > 0:
        specificity = result['tn'] / (result['tn'] + result['fp'])

    recall = float(1)
    if result['tp'] + result['fn'] > 0:
        recall = result['tp'] / (result['tp'] + result['fn'])
    
    precision = float(1)
    if result['tp'] + result['fp'] > 0:
        precision = result['tp'] / (result['tp'] + result['fp'])
    
    fscore = float(0)
    if precision > 0 or recall > 0:
        fscore = 2*((precision*recall) / (precision+recall))
    
    return sample, sensitivity, specificity, recall, precision, fscore
    

if __name__ == '__main__':
    if len(sys.argv) == 4:
        subvcf, truvcf, outF = sys.argv[1:4]

        if not subvcf.endswith('.vcf') and not subvcf.endswith('.vcf.gz'):
            sys.stderr.write("submission VCF filename does not enc in .vcf or .vcf.gz\n")
            sys.exit(1)

        if not os.path.exists(truvcf + '.tbi'):
            sys.stderr.write("truth VCF does not appear to be indexed. bgzip + tabix index required.\n")
            sys.exit(1)

        #print "\nCalculating rates for sample:"
        counts = evaluate(subvcf, truvcf)
        statresults = stats(counts)
        print "sample, sensitivity, specificity, recall, precision, F1-score: " + ','.join(map(str, statresults))
        
        countLine = str(counts['tp']) + '\t' + str(counts['tn']) + '\t' + str(counts['fp'])+ '\t' + str(counts['fn'])
        resultLine = '\t'.join(map(str, statresults))
        outFH = open(outF, "a")
        outFH.write( resultLine + '\t' + countLine)
        outFH.write( '\n' )
        outFH.close()
        #print "number of counted mutations in submission: " + str(ncalls)

    else:
        print "standalone usage for testing:", sys.argv[0], "<submission VCF> <truth VCF (tabix-indexed)> <output file>"