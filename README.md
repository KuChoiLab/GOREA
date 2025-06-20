# GOREA
### A Framework for Systematic and Unbiased Interpretation of Gene Ontology Enrichment
## Key algorithm for GOREA
<img width="830" alt="image" src="https://github.com/user-attachments/assets/17632e57-f23a-41e7-b8ee-aa371045a855" />

To conduct GOREA, a dataframe containing significant GOBP terms with proportion of overlapping genes relative to the total number of genes in each GOBP term from Over-Represenation Analysis (ORA) or Normalized enrichment score (NES) from Gene Set Enrichment Analysis (GSEA) must be assigned as input data (these are indicated as *Score in the figure). First, using the input data, clustering step is conducted, wherein combined method that we devised to apply the positive aspects of binary cut and hierarchical clustering is utilized. To define representative terms, information about ancestor terms and levels of GOBP terms was used. Specifically, the following steps were applied iteratively within each cluster: (1) Identify the common ancestor term at the highest level that covers a subset of the input GOBP terms. (2) For any remaining GOBP terms not explained by the representative term from step (1), the process is repeated on the remaining terms.

Note:
1. This algorithm can be applied to human and mouse.
2. When performing GOREA based on the results from GSEA, it is important to account for the directionality of enrichment by creating input sets that include only genes with either postiive or negative NES values.

## The result of GOREA
<img width="913" alt="image" src="https://github.com/user-attachments/assets/8a36e81f-62d6-4c69-884d-c987c90faaf9" />

The clustering results are displayed as a heatmap, with the representative terms shown on the right. On the right side of the heatmap, each cluster with its representative terms is ordered by singificance, based on either the average proportion of overlapping genes or the average absoulte value of the NES. To make more general observations about the resulting GOBP terms, we created a panel above the heatmap. First, we defined broad GOBP terms. Broad GOBP terms including input GOBP terms as their child term are clustered. For each cluster, the broad GOBP terms that contain the highest number of significant GOBP terms as their child term are selected and displayed. On the right side of the broad GOBP terms’ panel, the percentage of GOBP terms that each broad GOBP terms encompasses as child term are indicated.

## Run GOREA (human)
The following code is used to perform GOREA analysis.

```R
### human tutorial ####
library(dplyr)
library(plyr)
library(fgsea)
library(tibble)
library(ggplot2)
library(GOSemSim)
library(WriteXLS)
library(colorRamp2)
library(simplifyEnrichment)
library(ComplexHeatmap)
library(org.Hs.eg.db) # human

setwd("output path")
source("path to GOREA/GOREA/human/20250603_gorea_function_human_hj.R")

# 1. Setting environment for analysis ----

localdir <- "/Users/hojin/Dropbox/project/GOREA/20250603/" # this directory has to be parents directory of GeneOntology directory
gorea_enviromnet(localdir)

# 2. example ---

## 2.1 make test data ----
test <- sample(GOID_TERM$GOID, 500, replace = F) 
input_df <- data.frame(GOID = test)
input_df$NES <- sample(seq(0.1, 4, 0.01), replace = T, size = nrow(input_df))

head(input_df)

## example ##
# GOID  NES
# 1 GO:0071674 1.29
# 2 GO:0031110 3.55
# 3 GO:0007266 3.92
# 4 GO:0046112 2.65
# 5 GO:1903911 2.08
# 6 GO:0048670 1.93

## 2.2 outlier plot (additional step) ----
# before starting clustering steps, you can check a plot for the number of small clusters depending on cutoff (the cutoff is the value that you can assign in the gorea function as a parameter)
w <- gorea_sim_mat(input = input_df, godata_GO = godata_GO)
gorea_outlier_plot(w = w)

## 2.3 GOREA main function ----
# input_df; input dataframe from tools such as fgsea and enrichGO.
# k_val; when increasing this value, the number of cluster is increased.
# cutoff; if you want to remove broad amount of small clusters, decrease this cutoff. but, according to simplifyenrichment, 0.85 is a default value. 
# top_ancestor_annotation; for general description, top panel for broad GOBP terms is implemented.
# top_ancestor_annotation_number; clustered broad GO terms will be split according to this number. and the high ranked broad GO terms based on the number of input GOBP terms as child terms for each cluster can be illustrated in the top panel.
# filename1; this file includes total input GOBP terms and assigned ACNESTOR terms for each GOBP term.
# filename2; ㅅhis file includes representative terms (referred to as ANCESTOR terms), and among them, the ones associated with the largest GOBP terms are displayed in the result plot.
res <- gorea(input = input_df,
             k_val = 20, # considering your total number of GOBP terms.
             godata_GO = godata_GO,
             cutoff = 0.85, # default (you can change this value, according to the results from gorea_outlier_plot function)
             outlier_detect = T,
             min_cluster = 3,
             representative_term_level_cutoff = 1, GO_explain = 3,
             score = "NES", # "NES" or "Overlap_freq"
             filename1 = "testfile1.xlsx",
             filename2 = "testfile2.xlsx",
             heatmap_filename = "testplot.png",
             plot = T,
             heatmap_width = 40, heatmap_height = 30,
             ancestor_annotation = T,
             right_annotation_font_size = 10,
             cluster_font_size = 4,
             top_ancestor_annotation = T,
             top_ancestor_annotation_number = 3,
             color = c("gold"))
```

## Real example for GOREA analysis
<img width="1249" alt="image" src="https://github.com/user-attachments/assets/3dd9eef9-8459-44e1-85a5-131c88a5a369" />


## Reference
1. Gu, Zuguang, and Daniel Hübschmann. "simplifyEnrichment: a Bioconductor package for clustering and visualizing functional enrichment results." Genomics, Proteomics & Bioinformatics 21.1 (2023): 190-202.
2. Liberzon, Arthur, et al. "Molecular signatures database (MSigDB) 3.0." Bioinformatics 27.12 (2011): 1739-1740.
