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

    tpcount = 0 # counts TPs that are not masked
    fpcount = 0 # counts FPs that not masked

    tpcountmasked = 0 # counts masked TPs
    fpcountmasked = 0 # counts masked FPs
    trurecs = 0

    subrecs = 0 # counts all predicted variant records
    subuncount = 0 # counts predicted variant records that are not counted as TPs/FPs
    trunotmasked = 0 # counts true variants that are not masked

    truGeno = {}

    ''' store list of truth records, otherwise the iterator needs to be reset '''
    for trurec in truvcfh:
        var = str(trurec.CHROM) + ":" + str(trurec.POS) + ":" + trurec.REF
        alt = trurec.ALT
        #print var

        for sample in trurec.samples:
            genotype = format(sample['GT'])
            #print genotype
            truGeno[var][alt] = genotype
            ++trurecs 
            #print var + ":" + str(genotype)
            if genotype == "0/0": #WT
                ++truWtCount
    
    for subrec in subvcfh:
        var = str(subrec.CHROM) + ":" + str(subrec.POS) + ":" + subrec.REF 
        alt = subrec.ALT
        #print var

        for sample in trurec.samples:
            genotype = format(sample['GT'])
            if alt in truGeno[var]:
                truGenotype = truGeno[var][alt]
                if genotype == truGenotype:
                    ++tpcount
                else:
                    ++fpcount

        else: 
            #print genotype
            truGeno[var] = genotype
            ++trurecs 
    # trulist = [trurec for trurec in truvcfh]
    
#     ''' track whether the truth uses the "chr" prefix (all truth entries are assumed to use the same reference) '''
#     usechr = trulist[0].CHROM.startswith('chr')

#     ''' count records in truth vcf, track contigs/chromosomes '''
#     for trurec in trulist:
#         # check whether these record have concordant genotypes in the evaluated VCF



#         if relevant(trurec, vtype, ignorechroms):
#             truchroms[trurec.CHROM] = True
#             if not mask(trurec, truvcfh, truchroms, active=truthmask):
#                 trunotmasked += 1
    
#     # sanity check
#     if trunotmasked == 0:
#         raise Exception("No unmasked records found in truth file!\n")


#     '''
#     keep track of 'truth' sites used, they should only be usable once
#     if toCount=True, the true variant will be counted; if False, it will not
#     '''
#     used_truth = {}
    
#     '''
#     if submitters use MATEID in their BND calls we can 'tie' them together,
#     indexed by one mate, contains info on other mate in pair
#     '''
#     used_bnd_mates = {}

#     ''' parse submission vcf, compare to truth '''
#     for subrec in subvcfh:
#         subrec = prefix(subrec, usechr)
#         if relevant(subrec, vtype, ignorechroms) and passfilter(subrec):
#             subrecs += 1
            
#             matched = 'UN'

#             startpos, endpos = subrec.start, subrec.end

#             if vtype == 'SV' and subrec.is_sv:
#                 startpos, endpos = expand_sv_ends(subrec, useCIs=False)
            
#             sub_is_masked = mask(subrec, truvcfh, truchroms, active=truthmask)
            
#             if subrec.CHROM in truchroms:
#                 truoverlaplist = [trurec for trurec in truvcfh.fetch(subrec.CHROM, startpos, end=endpos)]
#                 for trurec in truoverlaplist:
#                     if relevant(trurec, vtype, ignorechroms) and match(subrec, trurec, vtype=vtype):
#                         # subrec matches a true variant
#                         tru_is_masked = mask(trurec, truvcfh, truchroms, active=truthmask)
                        
#                         if not sub_is_masked and not tru_is_masked:
#                             if str(trurec) not in used_truth and matched == 'TP':
#                                 # subrec already matches another true variant
#                                 used_truth[str(trurec)] = { 'toCount': False, 'masked': tru_is_masked }
                            
#                             elif matched != 'TP' and (str(trurec) not in used_truth or not used_truth[str(trurec)]['toCount']):
#                                 # subrec matches an unused true variant
#                                 matched = 'TP'
#                                 used_truth[str(trurec)] = { 'toCount': True, 'masked': tru_is_masked }
                        
#                         elif str(trurec) not in used_truth:
#                             # subrec matches an unused true variant but at least one is masked
#                             used_truth[str(trurec)] = { 'toCount': False, 'masked': tru_is_masked }
            
#                         if matched != 'TP':
#                             '''
#                             matched a true variant that was already counted
#                             and/or has a different mask status/both are masked
#                             only note if a true variant wasn't already IDed for subrec
#                             '''
#                             matched = 'T'

#             if matched == 'TP' or matched == 'T':
#                 if subrec.ID in used_bnd_mates and not used_bnd_mates[subrec.ID]['positive']:
#                     # BND mate call was a false positive, remove conflict
#                     if used_bnd_mates[subrec.ID]['masked']:
#                         fpcountmasked -= 1
#                     else:
#                         fpcount -= 1
#                     subuncount += 1
                
#                 elif subrec.INFO.get('MATEID'):
#                     # keep track of the mate info
#                     if isinstance(subrec.INFO.get('MATEID'), list):
#                         mateID = subrec.INFO.get('MATEID')[0]
#                     else:
#                         mateID = subrec.INFO.get('MATEID')
#                     used_bnd_mates[mateID] = { 'positive': True, 'masked': sub_is_masked }
            
#             elif matched == 'UN' and subrec.ID not in used_bnd_mates:
#                 # subrec does not match a true variant and it is not tied to a previous match through a mate
#                 matched = 'FP'
                
#                 if sub_is_masked:
#                     fpcountmasked += 1
#                 else:
#                     fpcount += 1
                    
#                 if subrec.INFO.get('MATEID'):
#                     # keep track of the mate info
#                     if isinstance(subrec.INFO.get('MATEID'), list):
#                         mateID = subrec.INFO.get('MATEID')[0]
#                     else:
#                         mateID = subrec.INFO.get('MATEID')
#                     used_bnd_mates[mateID] = { 'positive': False, 'masked': sub_is_masked }
                
#             if matched == 'TP':
#                 if sub_is_masked:
#                     tpcountmasked += 1
#                 else:
#                     tpcount += 1
#             elif matched == 'UN' or matched == 'T':
#                 subuncount += 1

#     # sanity checks
#     assert (tpcount + fpcount + tpcountmasked + fpcountmasked + subuncount == subrecs)
    
#     if subrecs == 0:
#         raise Exception("No filter-passing variants in submission! Are you sure you selected the correct variant type (SNV/INDEL/SV)?\n")

#     # count the true variants that should not be counted as TPs/FNs
#     truuncountnotmasked = 0
#     for k, v in used_truth.iteritems():
#         if not v['toCount'] and not v['masked']:
#             truuncountnotmasked += 1
    
#     result = { 'tp' : float(tpcount),
#                'fp' : float(fpcount),
#                'fn' : float(trunotmasked - tpcount - truuncountnotmasked) }
    
#     return result


# def stats(result):
#     ''' calculate precision, recall, fscore  from result dictionary '''
#     assert 'tp' in result and 'fp' in result and 'fn' in result, "invalid result dictionary!"
    
#     recall = float(1)
#     if result['tp'] + result['fn'] > 0:
#         recall = result['tp'] / (result['tp'] + result['fn'])
    
#     precision = float(1)
#     if result['tp'] + result['fp'] > 0:
#         precision = result['tp'] / (result['tp'] + result['fp'])
    
#     fscore = float(0)
#     if precision > 0 or recall > 0:
#         fscore = 2*((precision*recall) / (precision+recall))
    
#     return recall, precision, fscore
    

if __name__ == '__main__':
    if len(sys.argv) == 3:
        subvcf, truvcf = sys.argv[1:3]

        if not subvcf.endswith('.vcf') and not subvcf.endswith('.vcf.gz'):
            sys.stderr.write("submission VCF filename does not enc in .vcf or .vcf.gz\n")
            sys.exit(1)

        if not os.path.exists(truvcf + '.tbi'):
            sys.stderr.write("truth VCF does not appear to be indexed. bgzip + tabix index required.\n")
            sys.exit(1)

        print "\nmasked:"
        counts = evaluate(subvcf, truvcf)
        statresults = stats(counts)
        ncalls  = countrecs(counts)
        print "recall, precision, F1-score: " + ','.join(map(str, statresults))
        print "number of counted mutations in submission: " + str(ncalls)

    else:
        print "standalone usage for testing:", sys.argv[0], "<submission VCF> <truth VCF (tabix-indexed)>"