# CODE FOR ASSIGNMENT #4: Searching for Orthologues

```
This code was created in collaboration with my colleague Miguel La Iglesia Mirones
Lucía Muñoz Gil, MSc in Computational Biology, Bioinformatics Programming Challenges course
December 2023
```
**DESCRIPTION**

This code would perform a "reciprocal best hits BLAST" between the proteomes of Arabidopsis and S.pombe with the aim to find pairs of putative orthologues candidates. Arabidopsis proteome, as provide as a CDS file (containing nucleotide sequences) will be translated into protein sequences. Databases for both species will be created (and storaged at Databases folder). Our code can recognize the type of sequence to create the database, in this case both are protein type. Then, as seen in other similar projects, a double blastp would be performed. We decided to perform Best Reciprocal Hits by **two blastp instead of blastx and tblastn**, considering blastp between two proteomes is the most common way to search for orthologues (Camacho et al., 2009). Each sequence in the proteome of S.pombe will be blasted against the Arabidopsis proteome database. The best hit, if any, will then be blasted agaist the S.pombe proteome database. If second blast's best hit is the original S.pombe sequence, then we assume this is a pair of candidates 
orthologue protein candidates.

References:

[1] Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., &amp; Madden, T. L. (2009). BLAST+: Architecture and applications. BMC Bioinformatics, 10(1). https://doi.org/10.1186/1471-2105-10-421 

**CODE USAGE**
```
ruby main.rb TAIR10_cds.fa proteome_pombe.fa
```
or
```
ruby main.rb proteome_pombe.fa TAIR10.cds
```

**OUTPUT**
- **TAIR10_translated.fa**: translated proteome of Arabidopsis (from CDS nucleotide sequences to protein)
- **Orthologue_candidates_report.txt**: report with the total number of pairs of orthologue candidates, the total number of candidates and all names and belonging species for the pairs of candidates.

**"SENSIBLE" BLAST PARAMETERS AND REFERENCES**
- Maximum E-value threshold of $evalue = 1e^{-6}$
- Sorting from minimum to maximum E-value, and best hit is the one with lower E-value.
- Required coverage of at least 50% in the alignments (overlap)

References:

[2] Ward, N., &amp; Moreno-Hagelsieb, G. (2014). Quickly finding orthologs as reciprocal best hits with Blat, last, and UBLAST: How much do we miss? PLoS ONE, 9(7). https://doi.org/10.1371/journal.pone.0101850  
[3] Moreno-Hagelsieb, G., &amp; Latimer, K. (2007). Choosing blast options for better detection of orthologs as Reciprocal Best Hits. Bioinformatics, 24(3), 319–324. https://doi.org/10.1093/bioinformatics/btm585  

**POSSIBLE FUTURE STEPS**  
In order to continue our analysis, based on similar projects we have seen, we propose the following
- Employing other well-stablished methods including tree-methods or other graph-methods (Pairwise sequence similarity comparisons) appart from Best Reciprocal Hits) like Reciprocal Smallest Distance (RSD).
- Employing stardard databases of putative orthologues as EggNOG, Homologene or InParanoid
- Confirming Best Reciprocal Hits using fast BRH analysis software: lastal, diamond or MM2seqs2. (Hernández-Salmerón et al., 2020)
- Phylogenetic Analysis: construct phylogenetic trees to analyse the evolutionary relationship between the sequences.

References:

[4] Altenhoff, A. M., Boeckmann, B., Capella-Gutierrez, S., Dalquen, D. A., DeLuca, T., Forslund, K., Huerta-Cepas, J., Linard, B., Pereira, C., Pryszcz, L. P., Schreiber, F., da Silva, A. S., Szklarczyk, D., Train, C.-M., Bork, P., Lecompte, O., von Mering, C., Xenarios, I., Sjölander, K., … Dessimoz, C. (2016). Standardized benchmarking in the quest for orthologs. Nature Methods, 13(5), 425–430. https://doi.org/10.1038/nmeth.3830 
[5] https://www.flyrnai.org/RNAi_orthology.html
[6] Hernández-Salmerón, J.E., Moreno-Hagelsieb, G. (2020). Progress in quickly finding orthologs as reciprocal best hits: comparing blast, last, diamond and MMseqs2. BMC Genomics 21, 741. https://doi.org/10.1186/s12864-020-07132-6
[7] Wall, D. P., &amp; DeLuca, T. (2007). Ortholog detection using the reciprocal smallest distance algorithm. Comparative Genomics, 95–110. https://doi.org/10.1007/978-1-59745-515-2_7 

