#!/bin/python
#03 February 2016 - Kuan-Lin Huang @ WashU - 
 
import sys
import getopt

def main():
    def usage():
        print """
    liftover_CharGer_result.py : why do I exist?

    USAGE: liftover_CharGer_result.py [-h] <CharGer results file> <PCGP input file>
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

    args = sys.argv[1:]
    if len(args) < 1:
        usage(); sys.exit("input file missing")

    #open input file
    try:
        charGerF = open(args[0],"r")
    except IOError:
        print("File , args[0], does not exist!")
    
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
        if len(F) > 4:
            # if F[4] == "-":
            #     F[4] = "0"
            # if F[5] == "-":
            #     F[5] = "0"
            var = "_".join(F[1:3]+F[4:6])
            varCharGer[var]=line
    charGerF.close()


    try:
        inputF = open(args[1],"r")
    except IOError:
        print("File , args[1], does not exist!")

    header = inputF.readline().strip()
    headerF = header.split("\t") #arrays mutable, I'm farily sure
    i = 1
    #for headerItem in headerF:
    for k in range(0, len(headerF)):
        if headerF[k] == "":
             headerF[k] = "Missing" + str(i)
             i+=1
        headerF[k] = headerF[k].replace("#","").strip() #this is the same problem as self.userVariant or something
        # print ":" + headerF[k] + ":"
        # print i
        # i +=1
         

    print "\t".join(headerF[0:35]) + "\t" + CharGerHeader

    #read input file
    for line in inputF:
        line=line.strip()
        F = line.split("\t")
        if len(F) > 20:
            F[16] = F[16].upper().replace("CHR", "")
            start = F[17]
            #stop = F[2]
            ref = F[18]
            alt = F[19]
            #sample = F[21]
            var = "_".join(F[16:20])
        
            CharGerAnno = ""
            if var in varCharGer:
                CharGerAnno = varCharGer[var]
        
        print "\t".join(F[0:35]) + "\t" + CharGerAnno  
        #print line + "\t" + CharGerAnno
    
    inputF.close()

if __name__ == "__main__":
    main()
