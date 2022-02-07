# Synovium Single Cell Analysis

Code and primary data for expression and BCR for single cells from peripheral blood and synovium from RA134, RA172, RA195, RA221, e.g. figure 1A-1D, 2A-2D.
Uses R 3.5.1

[02vdj_explore.Rmd](02vdj_explore.md) Filter and QC of VDJ data
[03cluster_cdr3.Rmd](03cluster_cdr3.Rmd)  Clonotype definitions, etc of VDJ data.
[04expression_cluster_de.Rmd](04expression_cluster_de.md) Unsupervised clustering and differential expression using Seurat. Also some analysis connecting VDJ props to clustering.
[05slingshot.Rmd](05slingshot.md) Trajectory analysis using slingshot.
[06rep_and_expression.Rmd](06rep_and_expression.md) Further analysis of clonotype properties (e.g. permutation tests) compared to gene expression defined subpopulations defined in 04expression_cluster_de.md.