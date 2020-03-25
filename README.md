
<!-- README.md is generated from README.Rmd. Please edit that file -->

# COmprehend fuNctional conseQUencEs R <img src="inst/logo/CONQUER.png" align="right" width="120" />

*Gerard Bouland, Joline Beulens, Joey Nap, Arno van der Slik, Arnaud
Zaldumbide, Leen ’t Hart and Roderick Slieker*

## Installation

Install the the development version from GitLab:  
`install.packages("devtools")`

``` r
devtools::install_github("roderickslieker/CONQUER.DB")
devtools::install_github("roderickslieker/CONQUER.d3")
devtools::install_github("roderickslieker/CONQUER")
```

## Overview

With the use of two functions, SNPs are summarised and visualised,
namely: `summarise()` and `visualise()`.

  - The `summarise()` function is used to collect all data related to
    SNPs.
  - The `visualise()` function initiates a RStudio Shiny-based dashboard
    that visualises all relevant plots.

Note: We use the LD data from the API of NIH. You will need to register
on the site to obtain a token. Please see:

<https://ldlink.nci.nih.gov/?tab=apiaccess>

The token is send by email and can be provided as character string.

``` r
library(CONQUER)

summarise(
  variants = "rs184660829",
  directory = "somedirectory",
  token = "sometoken"
)
```

``` r
visualise(variant = "rs184660829", 
          directory = "somedirectory")
```

       

## Data resources

| Source              | Variable                           | Details                                                                                                                                                                                                                                                                                       | Genome build      |
| :------------------ | :--------------------------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :---------------- |
| ENSEMBL REST API    | Linkage disequilibrium             | 1000 genomes phase 3, default setting CEU                                                                                                                                                                                                                                                     | \-                |
|                     | Variant location                   |                                                                                                                                                                                                                                                                                               | hg38              |
|                     | Recombination                      | 1000 genomes phase 3, default setting CEU                                                                                                                                                                                                                                                     | \-                |
|                     | Genes                              | Ensembl HAVANA                                                                                                                                                                                                                                                                                | hg38              |
| 4D genomes          | 3d interactions                    | IM-PET of A549, Brain hippocampus middle, CD34+, CD4+ Memory, CD4+ Naive, CD4+ T, CD8+ Naive, Foreskin fibroblast, Foreskin keratinocyte, Foreskin melanocyte, GM12878, H1 derived mesenchymal stem cell, H1 derived mesendoderm cell, H1 derived neural progenitors, H1 derived trophoblast, | Lift-over to hg38 |
|                     |                                    | H1ESC, HCC1954, HCT116, HELA, HEP,                                                                                                                                                                                                                                                            |                   |
|                     |                                    | HMEC, HSMM, HUVEC, IMR90, IPS19.11,                                                                                                                                                                                                                                                           |                   |
|                     |                                    | IPS6.9, K562, MCF7, NHEK, NHLF, PANC                                                                                                                                                                                                                                                          |                   |
| Epigenomics roadmap | Chromatin state segmentations      | 127 cell types                                                                                                                                                                                                                                                                                | Lift-over hg38    |
| ENCODE              | Transcription factor binding sites |                                                                                                                                                                                                                                                                                               | Lift-over hg38    |
| Multiple sources    | Experimental miQTLs                |                                                                                                                                                                                                                                                                                               |                   |
|                     | Predicted miQTLs                   |                                                                                                                                                                                                                                                                                               |                   |
|                     | meQTLs                             | Blood                                                                                                                                                                                                                                                                                         |                   |
|                     | pQTLs                              |                                                                                                                                                                                                                                                                                               |                   |
| GWAS catalogue      | GWAS catalogue                     | Latest version available                                                                                                                                                                                                                                                                      |                   |
