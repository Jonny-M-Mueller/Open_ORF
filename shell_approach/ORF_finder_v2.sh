#! /bin/bash
#PROBLEM MAXIMALE ZEILENLÄNGE 10 000 CHarakters
INPUT_FILE=$1
INPUT_STRING=$(cat $INPUT_FILE)

#Appends to Resultfile
  #clear resultsfile
  echo "" |cat > OUTPUT/ORF_results.txt
echo "==========================="|cat >> OUTPUT/ORF_results.txt
echo "$INPUT_FILE"|cat >> OUTPUT/ORF_results.txt
echo "==========================="|cat >> OUTPUT/ORF_results.txt
INPUT_STRING=$(echo "$INPUT_STRING"|sed -z 's/\n//g')
INPUT_STRING=$(echo "$INPUT_STRING"|sed 's/ //g')
echo "$INPUT_STRING"|cat > temp/cleanup.txt #> overwrites >> apends to file

#empty fragment files from last run
#clear temp_fragments
> temp/temp_fragments1.txt
> temp/temp_fragments2.txt
> temp/temp_fragments3.txt
> temp/ATGPOS.txt

#Fragmentizes first, creates 3 reading frames starting at the 1rst, 2nd and 3rd character
bash Fragmentize.sh temp/cleanup.txt 1
bash Fragmentize.sh temp/cleanup.txt 2
bash Fragmentize.sh temp/cleanup.txt 3

READING_FRAME[1]=temp/temp_fragments1.txt
READING_FRAME[2]=temp/temp_fragments2.txt
READING_FRAME[3]=temp/temp_fragments3.txt

for i in {1..3}
    do
        echo "====="
        echo "READINGFRAME[$i]"

        ATG_POS=$(cat "${READING_FRAME[$i]}"|grep -no "ATG") #gets only line number and the expression "ATG" = position of ATG-triplett (unit= triplets) in the current global frame
        ATG_POS=$(echo $ATG_POS|sed 's/:ATG//g') #removes everything but the line number
        echo "$ATG_POS"|cat >>temp/ATGPOS.txt
        echo "=================================="
        read -r -a START_POS_ARRAY <<< "$ATG_POS"
        Number_of_Starts=${#START_POS_ARRAY[@]}

        ORF_STRING=""
        for k in $(seq 0 $((Number_of_Starts-1)))
        do
            TRIP="-"
            PEP=""
            CURRENT_TRIP_POS=${START_POS_ARRAY[k]} # P
            while [ "$TRIP" != "TAA" ] && [ "$TRIP" != "TAG" ] && [ "$TRIP" != "TGA" ] && [ "$TRIP" != "" ]
            do
                echo "Current pos "$CURRENT_TRIP_POS
                TRIP=$(sed -n "${CURRENT_TRIP_POS}p" < ${READING_FRAME[$i]})
                ORF+=" $TRIP"
                case $TRIP in
                  ATG)
                  PEP="${PEP}M";;
                  AAA|AAG)
                  PEP="${PEP}K";;
                  AAC|AAT)
                  PEP="${PEP}N";;
                  AC*)
                  PEP="${PEP}T";;
                  ATA|ATC|ATT)
                  PEP="${PEP}I";;
                  AGG|AGA|CG*)
                  PEP="${PEP}R";;
                  AGC|AGT|TC*)
                  PEP="${PEP}S";;
                  CAG|CAA)
                  PEP="${PEP}Q";;
                  CAC|CAT)
                  PEP="${PEP}H";;
                  CC*)
                  PEP="${PEP}P";;
                  CT*)
                  PEP="${PEP}L";;
                  TTT|TTC)
                  PEP="${PEP}F";;
                  TAT|TAC)
                  PEP="${PEP}Y";;
                  TGT|TGC)
                  PEP="${PEP}C";;
                  TGG)
                  PEP="${PEP}W";;
                  GT*)
                  PEP="${PEP}V";;
                  GC*)
                  PEP="${PEP}A";;
                  GAT|GAC)
                  PEP="${PEP}D";;
                  GAA|GAG)
                  PEP="${PEP}E";;
                  GG*)
                  PEP="${PEP}G";;
                  TAA|TGA|TAG)
                  break;;
                esac
                echo "POS $CURRENT_TRIP_POS"
                CURRENT_TRIP_POS=$((CURRENT_TRIP_POS + 1)) # TRIP_POS increases by +1 anyway --> end position is -1 actually
            done
           #prints information about current orf to resultsfile $j and $l accounts for the start and end of an ORF being shifted 1 or 2 positions to the right due to the global reading frame
          if [ "$i" == "1" ] # PROBLEM: POS currently means triplet, not base pairs
           then
            GLOBAL_START=$(((START_POS_ARRAY[k]-1)*3+1)) #
            GLOBAL_END=$(((CURRENT_TRIP_POS-1)*3))
          fi
          if [ "$i" == "2" ]
            then
            GLOBAL_START=$(((START_POS_ARRAY[k]-1)*3+2))
            GLOBAL_END=$(((CURRENT_TRIP_POS-1)*3+1))
          fi
          if [ "$i" == "3" ]
           then
            GLOBAL_START=$(((START_POS_ARRAY[k])*3))
            GLOBAL_END=$(((CURRENT_TRIP_POS-1)*3+2))
          fi
            echo "GLOBAL FRAME: RF$i Frame Number: $((k+1)) Start: $GLOBAL_START End:$GLOBAL_END LENGTH: $((GLOBAL_END - $GLOBAL_START+1 )) bp"|cat >> OUTPUT/ORF_results.txt # current results are bullshit (triplets instead bp)
            echo "$ORF" |cat  >> OUTPUT/ORF_results.txt
            echo "$PEP" |cat  >> OUTPUT/ORF_results.txt
            echo "" |cat >> OUTPUT/ORF_results.txt
            echo "==========================" |cat >> OUTPUT/ORF_results.txt
            ORF=""
        done
    done
