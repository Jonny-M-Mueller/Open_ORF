#! /bin/ env python3
import re
import itertools
import argparse

def ORF_OUTPUT():
        global  CONTINUED_ORF_COUNT
        print('======================================',file=OUTPUT_FILE )
        if STRAND=="+":
            STRAND_FOR_POS=INPUT_STRING
            print("ORF:",ORF_LIST.index(ORF), file=OUTPUT_FILE)
            print('START:', ORF.start(1), '    END:', ORF.end(1),'    LENGTH:', length, 'STRAND:  +',  file=OUTPUT_FILE )
            CONTINUED_ORF_COUNT=+1 # ORF_LIST is reset for complementary strand; CONTINUED_ORF_COUNT keeps track of ORF number
        if STRAND=="-":
            STRAND_FOR_POS=COMPLEMENTARY
            print("ORF:",ORF_LIST.index(ORF) + CONTINUED_ORF_COUNT, file=OUTPUT_FILE)
            print('START:', ORF.start(1), '    END:', ORF.end(1),'    LENGTH:', length, 'STRAND:  -',  file=OUTPUT_FILE )
            print("Matching positions on forward strand:",file=OUTPUT_FILE)
            print(' START:', len(INPUT_STRING)- ORF.start(1),'   END:',len(INPUT_STRING)- ORF.end(1), file=OUTPUT_FILE)
        print( ORF.group(1),file=OUTPUT_FILE)
        print(PEP,file=OUTPUT_FILE)

        if args.tata==True:
            for TATA in Tata_RANGE_LIST:
                print("TATA-Box found within", TATA_DISTANCE, "bp", file=OUTPUT_FILE)
                print("TATA-Box start:", TATA.start(1), "relative to ORF:", ORF.start(1)-TATA.start(1), file=OUTPUT_FILE)
                print("Sequence:", STRAND_FOR_POS[TATA.start(1): ORF.start(1)+3],"...", file=OUTPUT_FILE)
        if args.kozak==True and ORF.start(1) >= 6:
            if KOZAK_COMPARE[0]=="1":
                print("Strong Kozak-Sequence found (-3 bp = G/A AND +4 bp = G )", file=OUTPUT_FILE)
                print("STROOOOOOOOOOOOOOOOOONG")
            if KOZAK_COMPARE[1]=="1":
                print("Medium strength Kozak-Sequence found (-3 bp = G/A OR +4 bp = G, but not both )", file=OUTPUT_FILE)
            if KOZAK_COMPARE[2]=="1":
                print("Weak Kozak-Sequence found (-3 bp = T/C AND +4 bp = T/C/A )", file=OUTPUT_FILE)
            print("Sequence:", STRAND_FOR_POS[KOZAK_SEQUENCE.start(1): KOZAK_SEQUENCE.end(1)], file=OUTPUT_FILE)



parser = argparse.ArgumentParser(description='A script to locate and translate all Open Reading Frames (ORFs) in a Sequence. Optional analysis of complementary strand and Markdown Output')
parser.add_argument('-m', '--markdown',action='store_true', help="enables output of markdown files 'forward.md' (and 'reverse.md' if '-c' flag is set)")
parser.add_argument('-c', '--complementary',action='store_true', help="complementary strand is generated and analysed")
parser.add_argument('-t', '--tata', action= 'store_true', help="checks for presence of TATA-Box ('TATA[TA]A[TA]')")
parser.add_argument('-d', '--distance', type=int, help="sets maximum distance of TATA-Box (counted from first T) to start of ORF (to A in ATG), defaults to 35 bp if not set, ignored when -t flag is not set")
parser.add_argument('-k', '--kozak', action= 'store_true', help="evaluates strength of Kozak_sequence by checking the -3 bp and +4 bp site relative to the ORf's start, assumes eukariotic Sequence, DOES NOT evaluate less conserved bases in the putative Kozak-Sequence")
parser.add_argument('-f', '--Filter', type=str, help="Filters ORFs depending on char: 'S': ORFs with a strong Kozak sequence, 'M': ORFs with strong or medium Kozak-sequence, 'T': ORFs [distance] bp downstream of TATA-sequence. Logic: Kozak filters are linked with 'OR' (-f SM --> ORF is filtered if it has atleast a medium strength Kozak sequence), while TATA_obx filter is linked with 'AND' (-f TS --> ORF is only filtered if it has a Strong Kozak Sequence and a TATA-Box) by default. 'OR' behaviour can be set by including the lettter 'O'. 'XOR' behaviour is forced by includign the letter 'X'  ")
parser.add_argument('file_name', help="location of file containing the Input Sequence")
args = parser.parse_args()

if not args.Filter==None:
    Filter_ONLY="Both" # used later in case only TATA or Kozak is filtered for
    Filter_TATA=False
    Filter_KOZAK="000"
    if 'T' in args.Filter:
        Filter_TATA=True
    if 'S' in args.Filter:
        Filter_KOZAK="1"+Filter_KOZAK[1:]
    if 'M' in args.Filter:
        Filter_KOZAK=Filter_KOZAK[:1]+"1"+Filter_KOZAK[2:]
    if 'W' in args.Filter:
        Filter_KOZAK=Filter_KOZAK[:2]+"1"
    if Filter_TATA==True and not Filter_KOZAK=="000":
        if 'O' in args.Filter:
            Filter_Operator='Or'
        elif 'X' in args.Filter:
            Filter_Operator='Xor'
        elif 'N' in args.Filter:
            Filter_Operator='Nand'
        else:
            Filter_Operator='And'
    elif Filter_TATA==False and not Filter_KOZAK=="000":
        Filter_ONLY="KOZAK"
    elif Filter_TATA==True and Filter_KOZAK=="000":
        Filter_ONLY="TATA"
    elif Filter_TATA==False and Filter_KOZAK=="000":
        print("Error, selected no filter conditions!")
        exit()


if args.complementary:
    print("Complementrary=true")
print(args.file_name)
INPUT_FILE=open(args.file_name)
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

if args.tata==True:
    print("===============================================================================================================")
    print("NOTE:")
    print("Presence of TATA-Box is checked via pattern 'TATA[TA]A[TA]', bases in [brackets] mean one of the listed bases")
    print("===============================================================================================================")
    TATA_pattern=re.compile('(?=(TATA[TA]A[TA]))')
if args.kozak==True:
    print("===============================================================================================================")
    print("NOTE:")
    print("Does not check if there is a Kozak Sequence but checks it's putative strength")
    print("Highly conserved sites -3 bp (A or G) and + 4 bp (G) relative to the Startcodons 'A'")
    print("Strong Kozak: both conserved bases present (-3 bp is A/G and  +4 is G)")
    print("Medium Kozak: one of the conserved bases present (either -3 bp is A/G or +4 bp is G but not both)")
    print("Weak Kozak: none of the conserved bases present (neither -3 bp is A/G nor +4 bp is G)")
    print("===============================================================================================================")
    KOZAK_strong_pattern=re.compile('(?=([ATGC]{3}[AG][ATGC]{2}ATGG))') # Kozak_sequence in which the -3 Base being A/G AND the +4 Base being G
    KOZAK_medium_pattern=re.compile('(?=(([ATGC]{3}[TC][ATGC]{2}ATGG)|([ATGC]{3}[AG][ATGC]{2}ATG[ATC])))') #Kozak_sequence in which the -3 Base being A/G XOR the +4 Base being G (but not both)
    KOZAK_weak_pattern=re.compile('(?=([ATGC]{3}[TC][ATGC]{2}ATG[ATC]))') #KOZAK_sequence in which the -3 base is T/C AND the +3 Base is NOT G

ORF_pattern=re.compile('(?=(ATG([CTGA]{3})*?(TGA|TAA|TAG)))')

CONTINUED_ORF_COUNT=0

if args.complementary==True:
    print("Analysing forward (+) and complimentary (-) strand")
    SEQUENCES_LIST=[INPUT_STRING, COMPLEMENTARY]

else:
    SEQUENCES_LIST=[INPUT_STRING]
    print("Analysing only forward (+) strand")

for SEQUENCE in (SEQUENCES_LIST):
    if SEQUENCE==INPUT_STRING:
        STRAND="+"
        if args.markdown==True:
            MD_FILE=open('forward.md', "w")
    if SEQUENCE==COMPLEMENTARY:
        STRAND="-"
        if args.markdown==True:
            MD_FILE=open('complementary.md', "w")

    #finds all putative TATA-Boxes and all putative Kozak-Sequences. For TATA-Boxes: handles max distance to ORF
    if args.tata==True:
        TATA_ITER=TATA_pattern.finditer(SEQUENCE)
        TATA_LIST=[]
        if args.distance== None:
            print("no distance set, defaulting to 35 bp")
            TATA_DISTANCE=35
        else:
            TATA_DISTANCE=args.distance
            print("distance set to", TATA_DISTANCE)
        for TOUPLE in TATA_ITER:
            TATA_LIST.append(TOUPLE)
    if args.kozak==True:
        K_STRONG_ITER=KOZAK_strong_pattern.finditer(SEQUENCE)
        K_STRONG_LIST=[]
        for TOUPLE in K_STRONG_ITER:
            K_STRONG_LIST.append(TOUPLE)
        K_MEDIUM_ITER=KOZAK_medium_pattern.finditer(SEQUENCE)
        K_MEDIUM_LIST=[]
        for TOUPLE in K_MEDIUM_ITER:
            K_MEDIUM_LIST.append(TOUPLE)
        K_WEAK_ITER=KOZAK_weak_pattern.finditer(SEQUENCE)
        K_WEAK_LIST=[]
        for TOUPLE in K_WEAK_ITER:
            K_WEAK_LIST.append(TOUPLE)

    ORF_ITER=ORF_pattern.finditer(SEQUENCE)
    STARCOUNTER=0
    ORF_LIST=[]
    for TOUPLE in ORF_ITER:
        ORF_LIST.append(TOUPLE)

    for ORF in ORF_LIST:        #ORF is type match --> match.end([group]) match.start() is possible
        length=ORF.end(1) - ORF.start(1)
        length_trip=length//3
        PEP=""

        if args.markdown==True:
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

            for l in LATER_ORFS:
                if not E_FREE_END_TEST==[]:
                    for x in E_FREE_END_TEST:
                        if x.end(1)> l.start(1) and (x.end(1)<l.end(1) or x.end(1)==l.end(1)):
                            E_FREE_END=False

                if l.start(1)<ORF.end(1) and l.end(1) > ORF.end(1): #BUGfix  later Start< current End could also mean OS- LS-LE -OE --> ORF ENd is outside or equal to later END
                    END_IN_LATER=True
                if l.end(1)==ORF.end(1): # only true for the earlier starting orf --> when L is the current ORF it does not get flagged as equal --> is marked with stars correctly
                    END_EQUAL=True


            #search an replace block
            if  START_IN_EARLIER==False and END_IN_LATER==False and END_EQUAL==False:                            # ORF all alone
                REPLACE_STRING="**" + ORF.group(1) + "**"
                STARADD=4
            if  START_IN_EARLIER==True and END_IN_LATER==False and END_IN_EARLIER==False and END_EQUAL==False:    # START of ORF is located within 1 or more earlier ORFS but end is alone
                REPLACE_STRING=ORF.group(1) + "**"
                STARADD=2
                STAR_POS="after"
            if  START_IN_EARLIER==False and (END_IN_LATER==True or END_EQUAL==True):                             # END of ORF is located within 1 or more other but START is alone
                REPLACE_STRING="**" + ORF.group(1)
                STARADD=2
                STAR_POS="before"
            if  END_IN_EARLIER==True or (START_IN_EARLIER==True and (END_IN_LATER==True or END_EQUAL==True)):    # END and START of ORF is located within 1 or more other ORFS
                REPLACE_STRING=ORF.group(1)
                STARADD=0
                if (START_IN_EARLIER==True and END_IN_EARLIER==True) and E_FREE_END==True:
                    PREVENT_END_SHIFT=True

            if PREVENT_END_SHIFT==True:
                SEQUENCE=SEQUENCE[:ORF.start(1) + STARCOUNTER-2] + REPLACE_STRING + SEQUENCE[ORF.end(1) + STARCOUNTER-2:]

            else:
                SEQUENCE=SEQUENCE[:ORF.start(1) + STARCOUNTER] + REPLACE_STRING + SEQUENCE[ORF.end(1) + STARCOUNTER:]
            STARCOUNTER=STARCOUNTER+STARADD

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

        #resets List of potential TATA-Boxes in range and resets the Kozak_sequence of the last ORF before processing next ORF
        #also resets the flags for Found Tatabox/ flag for Kozak-Sequence strength

        #Checks wether (a) TATA-Box(es) is/are in range and the strength of the Kozak-Sequence of the current ORF (measured at -3/+4 bp conservation relative to ORF start)
        #Sequence and location of TATA-BOX(es) and the KOZAK_sequence is saved
        if args.tata==True:
            TATA_IN_RANGE=False
            Tata_RANGE_LIST=[]
            for TATA in TATA_LIST:
                if (ORF.start(1)-TATA.start(1) <= TATA_DISTANCE) and (ORF.start(1) - TATA.start(1) >0):
                    TATA_IN_RANGE=True
                    Tata_RANGE_LIST.append(TATA)

        if args.kozak==True:
            KOZAK_COMPARE="000"
            KOZAK_SEQUENCE=""
            for K_STRONG in K_STRONG_LIST:
                if ORF.start(1)-K_STRONG.start(1)==6:
                    KOZAK_COMPARE="1"+KOZAK_COMPARE[1:]
                    KOZAK_SEQUENCE=K_STRONG
                    print("K_COMPARE:", KOZAK_COMPARE)

            for K_MEDIUM in K_MEDIUM_LIST:
                if ORF.start(1)-K_MEDIUM.start(1)==6:
                    KOZAK_COMPARE=KOZAK_COMPARE[:1] +"1"+KOZAK_COMPARE[2:]
                    KOZAK_SEQUENCE=K_MEDIUM
                    print("K_COMPARE:", KOZAK_COMPARE)

            for K_WEAK in K_WEAK_LIST:
                if ORF.start(1)-K_WEAK.start(1)==6:
                    KOZAK_COMPARE=KOZAK_COMPARE[:2] +"1"
                    KOZAK_SEQUENCE=K_WEAK
                    print("K_COMPARE:", KOZAK_COMPARE)

        #matches Filter Logic Operator from input, then checks the flags set for the current ORF with that Operator. If a true statement results the ORF and related Output is appended to the OUTPUT.txt file
        if not args.Filter==None: # TODO: have to filter for NOT!!!!
            if Filter_ONLY=="TATA":
                if TATA_IN_RANGE==Filter_TATA:
                    ORF_OUTPUT()
            if Filter_ONLY=="KOZAK":
                if Filter_KOZAK==KOZAK_COMPARE:
                    ORF_OUTPUT()
            if Filter_ONLY=="Both":
                match Filter_Operator:
                    case "Or":
                        if TATA_IN_RANGE==Filter_TATA or (Filter_KOZAK[0]==KOZAK_COMPARE[0]=="1" or Filter_KOZAK[1]== KOZAK_COMPARE[1]=="1" or Filter_KOZAK[2]==KOZAK_COMPARE[2]=="1"):
                            ORF_OUTPUT()
                    case "Xor":
                        if TATA_IN_RANGE==Filter_TATA ^ (Filter_KOZAK[0]==KOZAK_COMPARE[0]=="1" or Filter_KOZAK[1]==KOZAK_COMPARE[1]=="1" or Filter_KOZAK[2]==KOZAK_COMPARE[2]=="1"):
                            ORF_OUTPUT()
                    case "Nand":
                        if not (TATA_IN_RANGE==Filter_TATA and (Filter_KOZAK[0]==KOZAK_COMPARE[0]=="1" or Filter_KOZAK[1]==KOZAK_COMPARE[1]=="1" or Filter_KOZAK[2]==KOZAK_COMPARE[2]=="1")):
                            ORF_OUTPUT()
                    case "And":
                        if TATA_IN_RANGE==Filter_TATA and (Filter_KOZAK[0]==KOZAK_COMPARE[0]=="1" or Filter_KOZAK[1]==KOZAK_COMPARE[1]=="1" or Filter_KOZAK[2]==KOZAK_COMPARE[2]=="1"):
                            ORF_OUTPUT()
        #if no filter was set, simply print the output to the OUTPUT.txt file
        if args.Filter==None:
            ORF_OUTPUT()


if args.markdown==True:
    print(SEQUENCE, file=MD_FILE)
    MD_FILE.close()
OUTPUT_FILE.close()
