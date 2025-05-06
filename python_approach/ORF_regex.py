#! /bin/ env python3
import re
import itertools
INPUT_FILE=open('INPUT/HTR2A_short.txt')
OUTPUT_FILE=open('OUTPUT.txt', "w")
PRE_INPUT_STRING=INPUT_FILE.read()
INPUT_FILE.close()

INPUT_STRING=PRE_INPUT_STRING.replace("\n", "")

REVERSE_STRING=INPUT_STRING[::-1]
COMPLEMENTARY=REVERSE_STRING.replace("A","x")
COMPLEMENTARY=COMPLEMENTARY.replace("T","A")
COMPLEMENTARY=COMPLEMENTARY.replace("x","T")
COMPLEMENTARY=COMPLEMENTARY.replace("G","y")
COMPLEMENTARY=COMPLEMENTARY.replace("C","G")
COMPLEMENTARY=COMPLEMENTARY.replace("y","C")

pattern=re.compile('(?=(ATG([CTGA]{3})*?(TGA|TAA|TAG)))')
CONTINUED_ORF_COUNT=0

for SEQUENCE in (INPUT_STRING, COMPLEMENTARY):
    if SEQUENCE==INPUT_STRING:
        STRAND="+"
    if SEQUENCE==COMPLEMENTARY:
        STRAND="-"
    if SEQUENCE==INPUT_STRING:
        MD_FILE=open('forward.md', "w")
    if SEQUENCE==COMPLEMENTARY:
        MD_FILE=open('COMPLEMENTARY.md', "w")

    ORF_ITER=pattern.finditer(SEQUENCE)
    STARCOUNTER=0
    ORF_LIST=[]
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
        E_FREE_END_TEST=[]
        E_FREE_END=True
        EARLIER_ORFS=ORF_LIST[0 : ORF_LIST.index(ORF) ]
        LATER_ORFS=ORF_LIST[ORF_LIST.index(ORF) +1 : len(ORF_LIST) ]


        for e in EARLIER_ORFS:
            if e.end(1) > ORF.start(1):
                START_IN_EARLIER=True
                if e.end(1) > ORF.end(1):
                    END_IN_EARLIER=True
                    E_FREE_END_TEST.append(e)
                    print(e.start(1), e.end(1))

        for l in LATER_ORFS:
            if not E_FREE_END_TEST==[]:
                for x in E_FREE_END_TEST:
                    if x.end(1)> l.start(1) and (x.end(1)<l.end(1) or x.end(1)==l.end(1)):
                        E_FREE_END=False
                    print( "X_ENDE", x.end(1), "l_start", l.start(1), "l_ende",l.end(1))

            if l.start(1)<ORF.end(1) and l.end(1) > ORF.end(1): #BUGfix  later Start< current End could also mean OS- LS-LE -OE --> ORF ENd is outside or equal to later END
                END_IN_LATER=True
            if l.end(1)==ORF.end(1): # only true for the earlier starting orf --> when L is the current ORF it does not get flagged as equal --> is marked with stars correctly
                END_EQUAL=True


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
            if (START_IN_EARLIER==True and END_IN_EARLIER==True) and E_FREE_END==True:
                PREVENT_END_SHIFT=True
                print("Shift prevented")

        print("E_FREE_END", E_FREE_END)
        if PREVENT_END_SHIFT==True:
            print("ORF.start", ORF.start(1), "ORF.end(1)", ORF.end(1), "STARCOUNT= ", STARCOUNTER, "||ORF.start+STARCOUNT", ORF.start(1)+ STARCOUNTER, "ORF.End+STARCOUNT", ORF.end(1)+ STARCOUNTER)
            SEQUENCE=SEQUENCE[:ORF.start(1) + STARCOUNTER-2] + REPLACE_STRING + SEQUENCE[ORF.end(1) + STARCOUNTER-2:]

        else:
            print("ORF.start", ORF.start(1), "ORF.end(1)", ORF.end(1), "STARCOUNT= ", STARCOUNTER, "||ORF.start+STARCOUNT", ORF.start(1)+ STARCOUNTER, "ORF.End+STARCOUNT", ORF.end(1)+ STARCOUNTER )
            SEQUENCE=SEQUENCE[:ORF.start(1) + STARCOUNTER] + REPLACE_STRING + SEQUENCE[ORF.end(1) + STARCOUNTER:]

        STARCOUNTER=STARCOUNTER+STARADD
        print(REPLACE_STRING)
        print("Starts in earlier orf ", START_IN_EARLIER, " Ends in earlier ORF ", END_IN_EARLIER, " Ends in later ORF ", END_IN_LATER, "EQUAL ENDS", END_EQUAL)
        for i in range(length_trip):
            j=i*3
            k=j+3
        # print(i,j,k)
            CURRENT_TRIP=ORF.group(1)[j:k]
            #print(CURRENT_TRIP )
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

        print('length',length,'   #triplets', length_trip,)
        print('======================================',file=OUTPUT_FILE )
        if STRAND=="+":
            print("ORF:",ORF_LIST.index(ORF), file=OUTPUT_FILE)
            print('START:', ORF.start(1), '    END:', ORF.end(1),'    LENGTH:', length, 'STRAND:  +',  file=OUTPUT_FILE )
            CONTINUED_ORF_COUNT=CONTINUED_ORF_COUNT +1 # ORF_LIST is reset for complimentary strand; CONTINUED_ORF_COUNT keeps track of ORF number
        if STRAND=="-":
            print("ORF:",ORF_LIST.index(ORF) + CONTINUED_ORF_COUNT, file=OUTPUT_FILE)
            print('START:', ORF.start(1), '    END:', ORF.end(1),'    LENGTH:', length, 'STRAND:  -',  file=OUTPUT_FILE )
            print("Matching positions on forward strand:",file=OUTPUT_FILE)
            print(' START:', len(INPUT_STRING)- ORF.start(1),'   END:',len(INPUT_STRING)- ORF.end(1), file=OUTPUT_FILE)
        print( ORF.group(1),file=OUTPUT_FILE)
        print(PEP,file=OUTPUT_FILE)
        print("==============")
    print(SEQUENCE, file=MD_FILE)
    #print(SEQUENCE)
    MD_FILE.close()
OUTPUT_FILE.close()
