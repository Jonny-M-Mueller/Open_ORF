# Open_ORF
A Python script to extract all open reading frames (ORFs) from a DNA sequence, aditionall features comming later.  
earlier (much slower) shell scripts with less features can be found in the earlier commits but are not maintained anymore  



Planned UPDATEs
=
ORF_regex.py: compatibility with RNA  
ORF_regex.py: detection of further important sequences (e.g. TATA-box, Kozak sequence, etc)  


**Current version: May 06 2025**
=

removed all shell scripts (still available via earlier commits)  
merged up-to-date python branch into main branch  

ORF_regex: input-file is now determined directly from the terminal/command line  
ORF_regex: a help page with syntax and a short description is now available (OPTION_ -h/--help)  
ORF_regex: now supports creation of markdown files highlighting the coding region of the input sequence (OPTION_ -m/--markdown)  
ORF_regex: can now analyse the complementary strand (OPTION_ -c/--complementary)


April 22th 2025
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

