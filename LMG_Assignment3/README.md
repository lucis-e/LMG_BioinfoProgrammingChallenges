# CODE FOR ASSIGNMENT #3: GFF feature files and visualization

```
This code was created in collaboration with my colleague Miguel La Iglesia Mirones
Lucía Muñoz Gil, MSc in Computational Biology, Bioinformatics Programming Challenges course
November 2023
```
**DESCRIPTION**

This code would take a list of gene locus names and make a get request to retrieve information about each gene's sequences, constructing a final EMBL file. Then, employing the BioRuby library it would seach for the repetitive motif CTTCTT in all exonic regions of the given genes in both positive and complemnt strand. Finally, found motifs and its specific chromosomic / gene sequence coordinates would be annotated in two GFF3 files. Also, a report with gene information and count of those genes without CTTCTT motifs in its exons is created. This repository also includes an screeshot of the results of uploading the GFF file with chromosomal coordiantes to Ensembl web beside de **AT2G46340** gene. Documentation for this code is also available at **doc directory**.

**CODE USAGE**
```
ruby main.rb ArabidopsisSubNetwork_GeneList.txt
```

**OUTPUT**
- AT_repeats_chromosomal.gff: GFF3 file for CTTCTT annotated motifs with coordinates respect to the whole chromosome
- AT_repeats_sequence.gff: GFF3 file for the CTTCTT annotated motifs with coordinates respect to the specific gene's sequence
- gff_report.txt: Report on those genes that do not show the CTTCTT motif on their exonic regions

**ADDITIONAL**
- Arabidopsis_thaliana_219022016_19028298.pdf: Screeshot of the results of uploading the GFF file with chromosomal coordiantes to Ensembl web beside de **AT2G46340** gene
