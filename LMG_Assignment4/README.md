# CODE FOR ASSIGNMENT #4: Searching for Orthologues

```
This code was created in collaboration with my colleague Miguel La Iglesia Mirones
Lucía Muñoz Gil, MSc in Computational Biology, Bioinformatics Programming Challenges course
December 2023
```
**DESCRIPTION**

This code would perform a "reciprocal best hits BLAST" between the proteomes of Arabidopsis and S.pombe with the aim to find pairs of putative orthologues candidates. Arabidopsis proteome, as provide as a CDS file (containing nucleotide sequences) will be translated into protein sequences. Databases for both species will be created. Code can recognize the type of sequence to create the database, in this case both are protein type. Then, as seen in other similar projects, a double blastp would ne performed. Each sequence in the proteome of S.pombe will be blasted against the Arabidopsis proteome database. The best hit, if any, will then be blasted agaist the S.pombe proteome database. If second blast's best hit is the original S.pombe sequence, then we assume this is a pair of candidates 
orthologue protein candidates.

**CODE USAGE**
```
ruby main.rb ArabidopsisSubNetwork_GeneList.txt
```

**OUTPUT**
- AT_repeats_chromosomal.gff: GFF3 file for CTTCTT annotated motifs with coordinates respect to the whole chromosome
- AT_repeats_sequence.gff: GFF3 file for the CTTCTT annotated motifs with coordinates respect to the specific gene's sequence
- gff_report.txt: Report on those genes that do not show the CTTCTT motif on their exonic regions

**"SENSIBLE" BLAST PARAMETERS AND REFERENCES**
- Maximum E-value threshold of  ($evalue^1e-6$)

**FUTURE PERSPECTIVES**
