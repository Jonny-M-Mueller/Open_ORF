# Open_ORF
bash scripts to extract all open reading frames (ORFs) from a DNA sequence, aditionall features comming later

Goal of the project is to provide a fast terminal based collection of scripts to quickly analyse a DNA Sequence, find possible open reading frames, translating them to peptides and comparing different sequences for matching strings without having to use multiple tools or webservices. The output will be kept simple to be readable by humans but also so it can be used by another program.


Planned UPDATE
=

ORF_regex.py: compatibility with RNA
ORF_regex.py: variable In- and Outputlocations
Organization: seperate project into different branches instead of directories
README.md:    make file more readable

Current version: April 22th 2025
=

ORF_regex.py: redid project in python using regular expressions
ORF_regex.py: currently lacks variable input and output location
Organisation: seperated project into python_approach and shell_approach and reorganized directories

April 18th 2025
=

added script 'min_strander.sh'
BUGFIX: ORF_finder now creates empty temporary files instead of writing an empty first line to them
ORF_finder now prints information about the ORFs (Start and End position, Global Reading frame, ORF number in current global frame, lenght in bp)
ORF_finder now also prints the peptides resulting from the ORFs in single letter code to the results file
EARLIER UPDATES

INITIAL VERSION
=
*extracts ORFs from + strand only
*works for short and longer sequences
*removes spaces and line breaks automatically

NOTE:
I am a Biochemist and a self taught coder. I have tested the scripts and compared the results whith NCBIs ORF_finder (https://www.ncbi.nlm.nih.gov/orffinder/) but there might still be some bugs I have not found yet, espeacially when working with longer sequences.
The provided example sequence was also sourced from the NCBI https://www.ncbi.nlm.nih.gov/gene/3356.

