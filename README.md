# Open_ORF
bash based bash scripts to extract all open reading frames (ORFs) from a DNA sequence, aditionall features comming later

Goal of the project is to provide a fast terminal based collection of scripts to quickly analyse a DNA Sequence, find possible open reading frames, translating them to peptides and comparing different sequences for matching strings without having to use multiple tools or webservices. The output will be kept simple to be readable by humans but also so it can be used by another program.

Use Instruction:
The scripts are intended for use with a unix shell.
Keep the files 'Fragmentize.sh' and the file 'ORF_finder_v2' in the same directory.
Use a FASTA file as input but remove the initial line manually, line breaks and empty spaces are removed automatically.
To generate all ORFs of the + strand type in your terminal:

./ORF_finder_v2.sh 'input_file.txt'

or 

./ORF_finder_v2.sh 'path_to/input_file.txt'

for the provided example file this would be:
./ORF_finder_v2.sh INPUT/HTR2A.txt

before using the scripts for the first time you need to make them executable
chmod +x ORF_finder_v2.sh; chmod +x Fragmentize.sh


Current version:
=
*extracts ORFs from + strand only
*works for short and longer sequences
*removes spaces and line breaks automatically

Planned features:
*improved input sanitation
*finding ORFs on the - strand
*annotation of found ORFS (position in original sequence, + or - strand)
*translating ORFs to peptides
*rewriting the scripts in python for platform independent work

NOTE:
I am a Biochemist and a self taught coder. I have tested the scripts and compared the results whith NCBIs ORF_finder (https://www.ncbi.nlm.nih.gov/orffinder/) but there might still be some bugs I have not found yet, espeacially when working with longer sequences.
The provided example sequence was also sourced from the NCBI https://www.ncbi.nlm.nih.gov/gene/3356.

