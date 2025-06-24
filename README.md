# oligo-aligner
align short DNA sequences to a reference DNA sequence

This scripts aligns several short DNA sequences to a single reference DNA sequence. This script is primarily for use in checking the specificity and coverage of primers, in particular for whole genome amplification and overlap-extension PCR. 

**#Input**

The script requires 2 inputs:

-reference sequence in fasta format

-oligo sequence in fasta format

**#Running the script:**

python3 oligo-aligner.py -ref "reference sequence in fasta format" -oligo "oligo sequence in fasta format" -output "output file name in fasta format"

**#output:**

-alignment file in fasta format

-mismatch.err: this file contains information on which oligos were not aligned

The output alignment file can then be loaded into any sequence alignment software (e.g. MEGA7, MEGAX, etc)

**#limitations:**

-The script only aligns oligos that have 100% match to the reference sequence. Degenerate sequences are not recognized. 
