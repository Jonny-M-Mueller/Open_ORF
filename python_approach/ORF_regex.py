#! /bin/ env python3
import re
import itertools
INPUT_FILE=open('test.txt')
OUTPUT_FILE=open('OUTPUT.txt', "w")
PRE_INPUT_STRING=INPUT_FILE.read()
INPUT_FILE.close()
MD_FILE=open('OUTPUT.md', "w")
INPUT_STRING=PRE_INPUT_STRING.replace("\n", "")

pattern=re.compile('(?=(ATG([CTGA]{3})*?(TGA|TAA|TAG)))')
ORF_ITER=pattern.finditer(INPUT_STRING)
STARCOUNTER=0
ORF_LIST=[]
L_STAR_TEST=[]
for TOUPLE in ORF_ITER:
    ORF_LIST.append(TOUPLE)

for ORF in ORF_LIST:        #ORF is type match --> match.end([group]) match.start() is possible
    length=ORF.end(1) - ORF.start(1)
    length_trip=length//3
    PEP=""
    START_IN_EARLIER=False
    END_IN_EARLIER=False
    END_IN_LATER=False
    PREVENT_END_SHIFT=False
    END_EQUAL=False
    L_NO_STAR=False
    EARLIER_ORFS=ORF_LIST[0 : ORF_LIST.index(ORF) ]
    LATER_ORFS=ORF_LIST[ORF_LIST.index(ORF) +1 : len(ORF_LIST) ]
    for e in EARLIER_ORFS:
        if e.end(1) > ORF.start(1):
            START_IN_EARLIER=True
            if e.end(1) > ORF.end(1):
                END_IN_EARLIER=True
                break

    for l in LATER_ORFS:
        if not L_STAR_TEST==[]:
            for x in L_STAR_TEST:
                if (l.start(1)<x.end(1) and l.end(1) > x.end(1)) or l.end(1)==x.end(1) :
                    L_NO_STAR=True
        if l.start(1)<ORF.end(1) and l.end(1) > ORF.end(1): #BUGfix  later Start< current End could also mean OS- LS-LE -OE --> ORF ENd is outside or equal to later END
            END_IN_LATER=True
            L_STAR_TEST.append(l)
        if l.end(1)==ORF.end(1): # only true for the earlier starting orf --> when L is the current ORF it does not get flagged as equal --> is marked with stars correctly
            END_EQUAL=True
            L_STAR_TEST.append(l)

    #search an replace block
    if  START_IN_EARLIER==False and END_IN_LATER==False and END_EQUAL==False:                            # ORF all alone
        REPLACE_STRING="**" + ORF.group(1) + "**"
        print("**ORF**")
        STARADD=4
    if  START_IN_EARLIER==True and END_IN_LATER==False and END_IN_EARLIER==False and END_EQUAL==False:    # START of ORF is located within 1 or more earlier ORFS but end is alone
        REPLACE_STRING=ORF.group(1) + "**"
        print("ORF**")
        STARADD=2
        STAR_POS="after"
    if  START_IN_EARLIER==False and (END_IN_LATER==True or END_EQUAL==True):                             # END of ORF is located within 1 or more other but START is alone
        REPLACE_STRING="**" + ORF.group(1)
        print("**ORF")
        STARADD=2
        STAR_POS="before"
    if  END_IN_EARLIER==True or (START_IN_EARLIER==True and (END_IN_LATER==True or END_EQUAL==True)):    # END and START of ORF is located within 1 or more other ORFS
        REPLACE_STRING=ORF.group(1)
        print ("ORF")
        STARADD=0
        if (END_IN_EARLIER==True or END_EQUAL==True) and L_NO_STAR==-False:
           PREVENT_END_SHIFT=True
           print("Shift prevented")

    print("L_NO_STAR", L_NO_STAR)
    if PREVENT_END_SHIFT==True:
        print("ORF.start", ORF.start(1), "ORF.end(1)", ORF.end(1), "STARCOUNT= ", STARCOUNTER, "||ORF.start+STARCOUNT", ORF.start(1)+ STARCOUNTER, "ORF.End+STARCOUNT", ORF.end(1)+ STARCOUNTER)
        INPUT_STRING=INPUT_STRING[:ORF.start(1) + STARCOUNTER-2] + REPLACE_STRING + INPUT_STRING[ORF.end(1) + STARCOUNTER-2:] # WORKS only one orf DEEP

    else:
        print("ORF.start", ORF.start(1), "ORF.end(1)", ORF.end(1), "STARCOUNT= ", STARCOUNTER, "||ORF.start+STARCOUNT", ORF.start(1)+ STARCOUNTER, "ORF.End+STARCOUNT", ORF.end(1)+ STARCOUNTER )
        INPUT_STRING=INPUT_STRING[:ORF.start(1) + STARCOUNTER] + REPLACE_STRING + INPUT_STRING[ORF.end(1) + STARCOUNTER:]   #IM ORF startende und endende Sequenzen werden auch durch abschließende

    STARCOUNTER=STARCOUNTER+STARADD
    print(REPLACE_STRING)
    print("Starts in earlier orf ", START_IN_EARLIER, " Ends in earlier ORF ", END_IN_EARLIER, " Ends in later ORF ", END_IN_LATER, "EQUAL ENDS", END_EQUAL)
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
    print("==============")

print(INPUT_STRING, file=MD_FILE)
#print(INPUT_STRING)
OUTPUT_FILE.close()
MD_FILE.close()
