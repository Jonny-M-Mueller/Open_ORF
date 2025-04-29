#! /bin/ env python3
import re
import itertools
INPUT_FILE=open('INPUT/HTR2A_short.txt')
OUTPUT_FILE=open('OUTPUT.txt', "w")
PRE_INPUT_STRING=INPUT_FILE.read()
INPUT_FILE.close()
MD_FILE=open('OUTPUT.md', "w")
INPUT_STRING=PRE_INPUT_STRING.replace("\n", "")

INDEX=0 #is not a true index but is used to slice the iterable
pattern=re.compile('(?=(ATG([CTGA]{3})*?(TGA|TAA|TAG)))')
ORF_ITER=pattern.finditer(INPUT_STRING)  #pattern.finditer = LIST OF FOUND PATTERNS; ORF = found pattern at loop position in LIST 'pattern.finditer'
for ORF in pattern.finditer(INPUT_STRING):  #BUG RELATED earlier problem is because itteratortools.islice CONSUMES the sliced part (i.e. the orignial iterator advances past the sliced part)
    INDEX=INDEX+1
    ABOVE=INDEX+1
    BELOW=INDEX-1
    length=ORF.end(1) - ORF.start(1)
    length_trip=length//3
    PEP=""
    START_IN_EARLIER="false"
    END_IN_EARLIER="false"
    END_IN_LATER="false"

    EARLIER_ORFS=itertools.islice(ORF_ITER, BELOW)  #BUG PROBLEM IS LIKELY HERE!! I think i misunderstood islice ITERSLICE CREATES AN ITER FROM ELEMENTS not A SLICE FROM AN ITER
    LATER_ORFS=itertools.islice(ORF_ITER, ABOVE, None)  #likely: cosumption of ORF_ITER in both directions prevents the function creating usefull output after the first time
    print("Earlier_start=", EARLIER_ORFS)
    print("later=", LATER_ORFS)
    for e in EARLIER_ORFS:                       # BUG all if here return false
        print("e.start=", e.start(1))
        if e.end(1) > ORF.start(1):
            START_IN_EARLIER="true"
            if e.end(1) > ORF.end(1):
                END_IN_EARLIER="true"
                break

    for l in LATER_ORFS:
        print("l=", l)
        if l.start(1) < ORF.end(1):
            END_IN_LATER="true"
            break
    #search an replace block
    if  START_IN_EARLIER=="false" and END_IN_LATER=="false":    # ORF all alone
        REPLACE_STRING="\*\*" + ORF.group(1) + "\*\*"

    if  START_IN_EARLIER=="true" and END_IN_LATER=="false":      # START of ORF is located within 1 or more earlier ORFS but end is alone
        REPLACE_STRING=ORF.group(1) + "\*\*"
    if  START_IN_EARLIER=="false" and END_IN_LATER=="true":      # END of ORF is located within 1 or more other but START is alone
        REPLACE_STRING="\*\*" + ORF.group(1)
    if  END_IN_EARLIER=="true":                                  # END and START of ORF is located within 1 or more other ORFS
        REPLACE_STRING=ORF.group(1)
    INPUT_STRING.replace(ORF.group(1), REPLACE_STRING)
    print("Se ", START_IN_EARLIER, " Ee ", END_IN_EARLIER, " El ", END_IN_LATER)

    for i in range(length_trip):
        j=i*3
        k=j+3
        print(i,j,k)
        CURRENT_TRIP=ORF.group(1)[j:k]
        print(CURRENT_TRIP )
        match CURRENT_TRIP:
             case 'ATG':
                 PEP=PEP + "M"
             case 'AAA' | 'AAG':
                 PEP=PEP + "K"
             case 'AAC' | 'AAT':
                 PEP=PEP + "N"
             case 'ACA' | 'ACC' | 'ACG' | 'ACT' :
                 PEP=PEP + "T"
             case 'ATA' | 'ATC' | 'ATT':
                 PEP=PEP + "I"
             case 'AGG' | 'AGA' | 'CGA' | 'CGC' | 'CGG' | 'CGT':
                 PEP=PEP + "R"
             case 'AGC' | 'AGT' | 'TCA' | 'TCC' | 'TCG' | 'TCT':
                 PEP=PEP + "S"
             case 'CAG' | 'CAA':
                 PEP=PEP + "Q"
             case 'CAC' | 'CAT':
                 PEP=PEP + "H"
             case 'CCA' | 'CCC' | 'CCG' | 'CGT':
                 PEP=PEP + "P"
             case 'CTA' | 'CTC' | 'CTG' | 'CTT' | "TTA" | "TTG":
                 PEP=PEP + "L"
             case 'TTT' | 'TTC':
                 PEP=PEP + "F"
             case 'TAT' | 'TAC':
                 PEP=PEP + "Y"
             case 'TGT' | 'TGC':
                 PEP=PEP + "C"
             case 'TGG':
                 PEP=PEP + "W"
             case 'GTA' | 'GTC' | 'GTG' | 'GTT':
                 PEP=PEP + "V"
             case 'GCA' | 'GCC' | 'GCG' | 'GCT':
                 PEP=PEP + "A"
             case 'GAT' | 'GAC':
                 PEP=PEP + "D"
             case 'GAA'|'GAG':
                 PEP=PEP + "E"
             case 'GGA'| 'GGC' | 'GGG' | 'GGT':
                 PEP=PEP + "G"
             case 'TAA' | 'TGA' | 'TAG':
                 PEP=PEP

    print('length',length,'   #triplets', length_trip)
    print('======================================',file=OUTPUT_FILE )
    print('START', ORF.start(1), '    END', ORF.end(1),'LENGTH', length, file=OUTPUT_FILE )
    print( ORF.group(1),file=OUTPUT_FILE)
    print(PEP,file=OUTPUT_FILE)
print(INPUT_STRING, file=MD_FILE)
OUTPUT_FILE.close()
MD_FILE.close()
