# Open_ORF
A Python script to extract all open reading frames (ORFs) from a DNA sequence, aditionall features comming later.  
earlier (much slower) shell scripts with less features can be found in the earlier commits but are not maintained anymore  



Planned UPDATEs
=

- adding the abillity to filter for the Absence of TATA-Boxes (I forgot about that)

**Current commit May 16th 2025**
=

ORF_regex.py Output can now include information about TATA-Boxes in Range and the strength of it's Kozak-Sequence  
- TATA_Box  
    - set Flag '-T' to detect TATA-Boxes in range of your ORF  
    - Default max Distance between TATA-Box and ORF is 35 bp but can be changed by setting the '-d' flag followed by an integer  
- Kozak_Sequence  
    - set flag '-k' to detect Kozak-Sequences  
    - Kozak Sequence "strength" (capacity of stalling the pre-initiation complex due to eIF2 interaction with the Base at -3/+4 bp from the start codon) is defined by conservation of the -3/+4 Base  
        - G or A at the -3 site  
        - G at the +4 site  
        - 2 conserved bases is classfied as strong, 1 conserved base is classified as medium, 0 conserved basses is classified as weak  
        - other less conserved bases are not factored in --> if the ORF starts atleast +6 bases after the beginning it, a Kosak sequence is always determined  
        
- Filtering  
    - can be activated by setting the -f flag followed by a combination of the Following characters (characters MUST be capitalized)
        - 'T' only includes ORFs with TATA_Boxes in range
        - 'S' only includes ORFs with a strong Kozak-Sequence in their resulting mRNA
        - 'M' only includes ORFs with a medium Kozak-Sequence in their resulting mRNA
        - 'W' only includes ORFs with a weak Kozak-Sequence in their resulting mRNA
        - 'O' sets the logic operator to 'OR'
        - 'X' sets the logic operator to 'XOR'
        - 'N' sets the logic operator to 'NAND'
    - Logic
        - Flags set regarding the Kozak strength function as or (i.e. '-f SM' would Output all ORFs which start with a Strong or Medium Kozak Sequence but include weak Kozak Sequences)
        - If you filter for BOTH TATA-Box and any Kozak-Sequence strength the AND operator is set as a default (i.e. a TATA-Box must be present and a Kozak-Sequence of the selected strength)
            -The filtered relationship between the TATA-Box and the Kozak-sequence can be changed by the letters 'O'/'X'/'N' (resulting in OR/XOR/NAND respectively. meaning at least TATABOX or Kozak criterion fullfilled/ exactly one of both fullfilled/ neither both TATA-Box and  selected strength Kozak-Sequence present or absent at the same time)
    

**May 06th 2025: 1.0**
=

removed all shell scripts (still available via earlier commits)  
merged up-to-date python branch into main branch  

ORF_regex: input-file is now determined directly from the terminal/command line  
ORF_regex: a help page with syntax and a short description is now available (OPTION: -h/--help)  
ORF_regex: now supports creation of markdown files highlighting the coding region of the input sequence (OPTION_ -m/--markdown)  
ORF_regex: can now analyse the complementary strand (OPTION: -c/--complementary)


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

