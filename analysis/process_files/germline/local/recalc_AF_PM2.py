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

    if len(sys.argv) >= 2:
        CharGerFH= sys.argv[1]
    else:
        usage()
        sys.exit()
    
    #open CharGer file
    try:
        charGerF = open(CharGerFH,"r")
    except IOError:
        print("File , CharGerFH, does not exist!")
    
    CharGerHeader = charGerF.readline().strip()
    print CharGerHeader + "\tExAC_adj_AF"
    #read input file
    for line in charGerF:
        line=line.strip()
        F = line.split("\t")
	AF = 0
	AFstr = F[80]
        AFstr1 = AFstr.split("&")
        AFrelevantStr = ""
        for AFstrings in AFstr1:
        	if AFstrings.startswith(F[7]):
                	AFrelevantStr = AFstrings
        if len(AFrelevantStr.split(":")) >1:
                AF = float( AFrelevantStr.split(":")[1] )
	if len(F) > 14 and "PM2" in F[14]:
		score = F[18]
		if AF > 0.0005: 
			F[18] = int(F[18])-2 # score
			F[14] = F[14].replace("PM2,","")
			if int(F[18]) < 5:
				continue
	F[18] = str(F[18])
	AF = str(AF)
	F.append(AF)
	print "\t".join(F)

    charGerF.close()

if __name__ == "__main__":
    main()
