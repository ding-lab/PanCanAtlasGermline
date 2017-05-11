# Sums up mutation data for each individual gene
def getSums(filename):
    ind = [1, 7, 38, 39, 41, 42, 68, 69] # chrom., gene, t-ref, t-var, n-ref, n-var, FClass, bin-type, num-mutations
    
    with open(filename, 'r') as file:
        geneSums = []
        next(file)
        current = ['chr', 'gene', 't-ref', 't-var', 'n-ref', 'n-var', 'FClass', 'bin-type', 'num']
        num = 0
        # Iterating through each line
        for line in file:
            parsed = line.split('\t')
            
            if(parsed[ind[1]] != current[1]):
                  
                if(current[0] != 'chr'):
                    geneSums.append([])
                    for i in current:
                        geneSums[-1].append(i)

                current[0] = int(parsed[ind[0]])
                current[1] = parsed[ind[1]]
                for i in range(2, 6):
                    current[i] = int(parsed[ind[i]])
                if(parsed[ind[6]] == "Uncertain Significance"):
                    current[6] = 0
                else:
                    current[6] = 1
                current[7] = parsed[ind[7]]
                current[8] = 1 
            
            else:
                for i in range(2, 6):
                    current[i] += int(parsed[ind[i]])
                if(parsed[ind[6]] != "Uncertain Significance"):
                    current[6] += 1
                current[8] += 1
        else:
            geneSums.append([])
            for i in current:
                geneSums[-1].append(i)

    return geneSums

# Does chi^2 contingency test based on independence of reference/variant with tumor/normal
# Returns significant genes with two extra arguments: the p-value, and whether tumors see more variant than normal.
import scipy.stats
from scipy.stats import chi2_contingency

# data should be in the form of geneSums, sig = significance level.
def sigGenes(data, sig = 0.00001):
    sigList = []
    for i in data:
        table = [[i[2], i[3]], [i[4], i[5]]]
        pvalue = chi2_contingency(table)[1]
        # print(pvalue)
        if(pvalue < sig):
            sigList.append([])
            for j in i:
                sigList[-1].append(j)
            sigList[-1].append(pvalue)
            if(float(i[3])/i[2] > float(i[5])/i[4]):
                sigList[-1].append("More Variant")
            else:
                sigList[-1].append("Less Variant")
    return sigList
	
# Check if this method has any correlation with the classifications
def check(data, sigList):
    total = 0
    for i in data:
        if(i[6] > 0):
            total += 1
    sig = 0
    for i in sigList:
        if(i[6] > 0):
            sig += 1
    table = [[total, len(data)-total], [sig, len(sigList)-sig]]
    pvalue = chi2_contingency(table)[1]
    print("Genes classified as at least Likely Pathogenic for some mutation: " + str(total) + "/" + str(len(data)) + "=" + str(float(total)/len(data)))
    print("Those classified genes out of those that were found to be significant: " + str(sig) + "/" + str(len(sigList)) + "=" + str(float(sig)/len(sigList)))
    print("p-value under a chi^2 test for independence: " + str(pvalue))