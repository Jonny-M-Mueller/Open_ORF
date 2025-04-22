#! /bin/ env python3
import re
INPUT_FILE=open('INPUT/HTR2A.txt')
OUTPUT_FILE=open('OUTPUT.txt', "w")
PRE_INPUT_STRING=INPUT_FILE.read()
INPUT_FILE.close()
INPUT_STRING=PRE_INPUT_STRING.replace("\n", "")

pattern=re.compile('(?=(ATG([CTGA]{3})*?(TGA|TAA|TAG)))')
for ORF in pattern.finditer(INPUT_STRING):
    length=ORF.end(1) - ORF.start(1)
    length_trip=length//3
    PEP=""
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
   # print(ORF.group(1), ORF.span(1))
