#! /bin/bash

INPUT_FILE=$1
INPUT_STRING=$(cat $INPUT_FILE)
echo ""|cat > temp/reverse.txt
echo ""|cat > OUTPUT/min_strand.txt

#currently also in ORF_finder --> after integration remove one instance
INPUT_STRING=$(echo "$INPUT_STRING"|sed -z "s/\n//g")
INPUT_STRING=$(echo "$INPUT_STRING"|sed 's/ //g')
echo "$INPUT_STRING"|cat > temp/cleanup.txt #> overwrites >> apends to file
#=== end of duplicate code===
echo $INPUT_STRING
INPUT_STRING_LENGTH=${#INPUT_STRING}
echo $INPUT_STRING_LENGTH

for i in $(seq $INPUT_STRING_LENGTH -1 1)
    do
    echo $i
 echo "$INPUT_STRING"|cut -c $i >>temp/reverse.txt
    done
OUTPUT=$(tr -d '\n' <temp/reverse.txt)
OUTPUT=$(echo "$OUTPUT"|sed 's/ //g')
# A/T exchange
OUTPUT=$(echo "$OUTPUT"|sed 's/T/y/g' )
OUTPUT=$(echo "$OUTPUT"|sed 's/A/T/g')
OUTPUT=$(echo "$OUTPUT"|sed 's/y/A/g' )

#G/C exchange
OUTPUT=$(echo "$OUTPUT"|sed 's/G/x/g')
OUTPUT=$(echo "$OUTPUT"|sed 's/C/G/g' )
OUTPUT=$(echo "$OUTPUT"|sed 's/x/C/g')
echo "$OUTPUT"
echo "$OUTPUT" > OUTPUT/min_strand.txt
