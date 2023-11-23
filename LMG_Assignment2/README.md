# CODE FOR ASSIGNMENT #2: INTENSIVE INTEGRATION USING WEB APIs

```
This code was created in collaboration with my colleague Miguel La Iglesia Mirones
Lucía Muñoz Gil, MSc in Computational Biology, Bioinformatics Programming Challenges course
November 2023
```

**PROJECT DESCRIPTION:**

The goal of this project was to programmatically access several databases and extract available protein-protein interaction inormation to build regulatory networks to determine if a set of predicted co-expressed genes are known to bind one another.

For this aim, we have develop code that would iterativelly request protein-protein interaction data for a list of genes defined by its locus name and Uniprot id. Recursive search with a depth of 2 was implemented by enabling our code to search for direct interactors of the direct interactors of every initial gene specified on the input file. However, further depth it is allowed by changing the value of DEPTH constant. We considerated interactions trustworthy if the quality score is higher than 0.5, were not determined by two-hybrid interaction detection methods (high rate of false positives) and are physical interactions.

Found interactors and initial gene would be integrated as part of a network. We assume that for a network to exist should have at least 2 members. Nodes without interactions are not included in any network. Networks with commmon members are merged into one. All networks are annotated with GO Terms and KEGG Pathways its members are part of / associated with. 

**CODE USAGE:**

```
ruby main.rb ArabidopsisSubNetwork_GeneList.txt Final_report.txt
```

**OUTPUT:**

A final report including:
-  A GLOBAL report for all genes: total number of networks built, total number of nodes (initial genes and interactors) and number of genes not included in any network
-  A report of the FEATURES of every netowork: network ID, number of nodes, gene locus names from file included in the network and KEGG and GO annotations.
  
**CONCLUSION:**

Only 7 interaction netowrks were found, out of which 6 included only 1 gene of the set of "co-expressed" genes. The remaining network included just 12 out of the 168 genes (162 if we exclude those in individual networks) genes inputed as part of the co-expression genes. This results demonstrate than 150 are not part of any network, or at least, there is no available and trustworthy information about its direct binding to any of the genes (or its interactors) in he IntAct database, which we employed for this assignment. By looking at these results we suspect this genes might not be coexpressed or part of the same coexpression set.
