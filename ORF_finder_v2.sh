#! /bin/bash
#PROBLEM MAXIMALE ZEILENLÄNGE 10 000 CHarakters
INPUT_FILE=$1
INPUT_STRING=$(cat $INPUT_FILE)

#Appends to Resultfile
echo "==========================="|cat >> OUTPUT/ORF_results.txt
echo "$INPUT_FILE"|cat >> OUTPUT/ORF_results.txt
echo "==========================="|cat >> OUTPUT/ORF_results.txt
INPUT_STRING=$(echo "$INPUT_STRING"|sed -z "s/\n//g")
INPUT_STRING=$(echo "$INPUT_STRING"|sed 's/ //g')
echo "$INPUT_STRING"|cat > temp/cleanup.txt #> overwrites >> apends to file

#empty fragment files from last run
#clear temp_fragments1
echo ""|cat > temp/temp_fragments1.txt
echo ""|cat > temp/temp_fragments2.txt
echo ""|cat > temp/temp_fragments3.txt

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

        ATG_POS=$(cat "${READING_FRAME[$i]}"|grep -no "ATG")
        ATG_POS=$(echo $ATG_POS|sed 's/:ATG//g')
        echo "=================================="
       # echo "ATG POS $ATG_POS"

      #  TGA_POS=$(cat "${READING_FRAME[$i]}"|grep -no "TGA")
       # TGA_POS=$(echo $TGA_POS|sed 's/:TGA//g')
        #echo "=================================="
        #echo "TGA POS $TGA_POS"

        #TAG_POS=$(cat "${READING_FRAME[$i]}"|grep -no "TAG")
        #TAG_POS=$(echo $TAG_POS|sed 's/:TAG//g')
        #echo "=================================="
        #echo "TAG POS $TAG_POS"

        #TAA_POS=$(cat "${READING_FRAME[$i]}"|grep -no "TAA")
        #TAA_POS=$(echo $TAA_POS|sed 's/:TAA//g')
        #echo "=================================="
        #echo "TAA POS $TAA_POS"

        read -r -a START_POS_ARRAY <<< "$ATG_POS"
        Number_of_Starts=${#START_POS_ARRAY[@]}

        ORF_STRING=""
        for k in $(seq 0 $((Number_of_Starts-1)))
        do
            TRIP="-"
            CURRENT_TRIP_POS=${START_POS_ARRAY[k]} # P
            while [ "$TRIP" != "TAA" ] && [ "$TRIP" != "TAG" ] && [ "$TRIP" != "TGA" ] && [ "$TRIP" != "" ]
            do
                echo "Current pos "$CURRENT_TRIP_POS
                TRIP=$(sed -n "${CURRENT_TRIP_POS}p" < ${READING_FRAME[$i]})
                ORF+=" $TRIP"
                echo "POS $CURRENT_TRIP_POS"
              #  echo "ORF $ORF"
                CURRENT_TRIP_POS=$((CURRENT_TRIP_POS + 1))
            done
            echo "$ORF" |cat  >> OUTPUT/ORF_results.txt
            echo "" |cat >> OUTPUT/ORF_results.txt
            ORF=""
        done
      #  echo "====== NEXT FRAME ====" |cat >> OUTPUT/ORF_results.txt
    done
