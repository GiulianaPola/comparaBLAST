# comparaBLAST
Program for comparing sequence BLAST results against two databases

## Input file (-a or -b)
The results are in tabular format files and were generated with BLAST used e-value filter and the parameter -max_target_seqs 1 (only the best hit for each sequence). 
Columns:
qseqid,sseqid,length,pident,bitscore,score,evalue,stitle

## Versions

### 10/10/2021 (1.0.0)
Redundancies were reduced by eliminating results with a lower e-value or bitscore since a sequence may encounter numerous alignment blocks. It was established if the affirmative sequences are present exclusively in databases A and B, or in both. It was determined which bank had the most significant result for the set of positive sequences in both databases based on the e-value (the more negative the e-value, the better). Both banks A and B's positive end sequences were listed. The positive end sequences for banks A, B, or both were saved to output file:

ID seq columns for A, evalue for A, ID seq for B, e-value, ID seq for A/B, e-value in A, evalue in B, decision (A or B)

### 18/10/2021 (1.0.1)
- Inserting the Python shebang line: 

!/usr/bin/python

- Incorporating help with -h or -help or no parameter. 

### 25/04/2022 (1.1.0)
- Incorporating the stitle column of databases A and B:

ID seq columns for A, evalue for A, stitle for A, ID seq for B, e-value, stitle for B, ID seq for A/B, e-value in A, stitle for A, evalue in B, stitle for B, decision (A or B)

### 28/04/2022
#### (1.1.1)
- Corrected column formatting in the output table;
- Changed "Escolha" column to "Selected_choice.

#### (1.1.2)
- Corrected error when extracting columns

### 29/04/2022 (1.1.3)
- Sorting the output table to align rows of the same qseqid from database A, B, and those in both