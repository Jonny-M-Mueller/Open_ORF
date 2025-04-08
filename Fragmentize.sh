#! /bin/bash
echo "========START FRAGMENTIZE $2"
INPUT=$1  #either .txt-file or string
echo "FRAGMENTIZER INPUT" $INPUT
#FRAGMENTFILE=$2
FRAGMENT_LENGHT=3
if [ -z $FRAGMENT_LENGHT ]      #FRAGMENT_LENGHT=3 is the default case if no FRAGMENT_LENGHT has been declared (DNA tripplets)
then
echo "FRAGMENT LENGTH=3"
FRAGMENT_LENGHT=3
fi
STARTPOSITION=$2
echo "STARTPOSITION $STARTPOSITION"
if [ -z $STARTPOSITION ]        # STARTPOSITION=1 is the default case if no different position is declared
    then
        STARTPOSITION=1
    fi
FRAGSTRING=''
declare -a FRAGMENT_ARRAY=()
#Defining INPUT_STRING and its Length; Checks if INPUT_STRING is a file or not (If not -> assumes String)
if [ -f $INPUT ]
    then
        echo "INPUT==FILE: TRUE"
        INPUT_STRING=$(cat "$INPUT")
    else
        echo "INPUT==FILE: FALSE"
        INPUT_STRING=$INPUT
fi
echo "Inputstring" $INPUT_STRING
INPUT_STRING=${INPUT_STRING:($STARTPOSITION - 1)}
echo echo "Inputstring reduced" $INPUT_STRING
INPUT_STRING_LENGTH=${#INPUT_STRING}
echo "LÄNGE="$INPUT_STRING_LENGTH

#Making sure input string is divisible by fragment lenghth, calculating total number of fragments
MODULO_LENGTH=$((INPUT_STRING_LENGTH % FRAGMENT_LENGHT))
echo "MODULO"$MODULO_LENGTH
if [ $MODULO_LENGTH != "0" ]
then
    echo "NICHT TEILBAR"
    MISSING_CHARS=$(( FRAGMENT_LENGHT - MODULO_LENGTH))     #Missing_Chars = characters that would be missing in the final fragment
    echo "MISSING $MISSING_CHARS"
       for i in $(seq 1 $MISSING_CHARS)
            do
             INPUT_STRING="${INPUT_STRING}-"               # Final fragment is filled up with '-' will be removed at end of program
            done
fi
NR_OF_FRAGMENTS=$(( INPUT_STRING_LENGTH / FRAGMENT_LENGHT))

#TESTBLOCK
echo 'STARTPOSITION=' $STARTPOSITION
echo 'FRAGMENT_LENGHT=' $FRAGMENT_LENGHT
echo 'Number of Fragments=' $NR_OF_FRAGMENTS

 echo ">"|cat >> "temp/temp_fragments${2}.txt" #ads ">" as a seperator before ORFS
for j in $(seq 1 $NR_OF_FRAGMENTS)
    do
        for k in $(seq 1 $FRAGMENT_LENGHT)
            do
            CHAR=$(echo $INPUT_STRING|cut -c 1)
            FRAG_STRING+=$CHAR
          #  echo "INPUT=" $INPUT_STRING #horrible for longer input
            echo "CHAR=" $CHAR
            echo "FragString=" $FRAG_STRING
            INPUT_STRING=${INPUT_STRING:1}
            done
        FRAGMENT_ARRAY+=($FRAG_STRING)
        FRAG_STRING=""
        k=$((j-1))
        echo "k= $k"
        #echo "ARRAY[k]= "${FRAGMENT_ARRAY[$k]}
       # echo "ARRAY= "${FRAGMENT_ARRAY[*]}
        echo "${FRAGMENT_ARRAY[$k]}"|cat >> "temp/temp_fragments${2}.txt"
    done
