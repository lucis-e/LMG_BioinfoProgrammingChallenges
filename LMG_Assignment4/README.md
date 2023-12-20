# CODE FOR ASSIGNMENT #4: Searching for Orthologues

```
This code was created in collaboration with my colleague Miguel La Iglesia Mirones
Lucía Muñoz Gil, MSc in Computational Biology, Bioinformatics Programming Challenges course
December 2023
```
**DESCRIPTION**

This code would perform a "reciprocal best hits BLAST" between the proteomes of Arabidopsis and S.pombe with the aim to find pairs of putative orthologues candidates. Arabidopsis proteome, as provide as a CDS file (containing nucleotide sequences) will be translated into protein sequences. Databases for both species will be created (and storaged at Databases folder). Code can recognize the type of sequence to create the database, in this case both are protein type. Then, as seen in other similar projects, a double blastp would ne performed. Each sequence in the proteome of S.pombe will be blasted against the Arabidopsis proteome database. The best hit, if any, will then be blasted agaist the S.pombe proteome database. If second blast's best hit is the original S.pombe sequence, then we assume this is a pair of candidates 
orthologue protein candidates.

**CODE USAGE**
```
ruby main.rb TAIR10_cds.fa proteome_pombe.fa
```
or
```
ruby main.rb proteome_pombe.fa TAIR10.cds
```

**OUTPUT**
- TAIR10_translated.fa: translated proteome of Arabidopsis (from CDS nucleotide sequences to protein)
- Orthologue_candidates_report: report with the total number of pairs of orthologue candidates, the total number of candidates and all names and belonging species for the pairs of candidates.

**"SENSIBLE" BLAST PARAMETERS AND REFERENCES**
- Maximum E-value threshold of $evalue = 1e^{-6}$

**FUTURE PERSPECTIVES**
