
<!-- README.md is generated from README.Rmd. Please edit that file -->

# COmprehend fuNctional conseQUencEs R <img src="https://github.com/roderickslieker/CONQUER.test/blob/master/CONQUER.png" align="right" width="120" />

*Gerard Bouland, Joline Beulens, Joey Nap, Arno van der Slik, Arnaud
Zaldumbide, Leen ’t Hart and Roderick Slieker*

## Installation

Install the the development version from GitLab:  
`install.packages("devtools")`

``` r
devtools::install_github("roderickslieker/CONQUER.db")
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

## Example without multianalyze

``` r
DIR <- "somedirectory"

library(CONQUER)

summarize(variants = c("rs878521","rs10830963"),
          directory=DIR,
          multiAnalyze=FALSE,
          token="sometoken",
          tissues=NULL)
```

## Example with multianalyze

The available tissues can be viewed with the following command:

``` r
tissues <- conquer.db::gtexTissuesV8
```

The summary files from the example below can also be obtained from
<https://github.com/roderickslieker/CONQUER.test/tree/master/Test>

``` r
library(CONQUER)
snps <- c("rs11642430","rs11820019","rs11842871","rs13426680","rs1377807","rs1783541",
"rs1801212","rs1801645","rs2268078","rs2581787","rs34855406","rs3802177",
"rs3810291","rs4148856","rs5213","rs6011155","rs601945","rs75423501",
"rs8010382","rs8046545")

CONQUER::summarize(variants = snps,
          directory=DIR,
          multiAnalyze=TRUE,
          token=NULL,
          tissues=c("Pancreas","Muscle_Skeletal","Liver"))


visualize(directory = "somedirectory", SNPs = snps)
```

## Modules

<img src="https://github.com/roderickslieker/CONQUER.test/blob/master/SunBurst_pathways.png" align="center" width="500" />

 

## Enrichment

<img src="https://github.com/roderickslieker/CONQUER.test/blob/master/Modules_pathways.png" align="center" width="500" />

 

## Pathways shared by tissues

<img src="https://github.com/roderickslieker/CONQUER.test/blob/master/Pathways_acrossTissue.png" align="center" width="500" />

 

## LD

<img src="https://github.com/roderickslieker/CONQUER.test/blob/master/LD_rs11642430.png" align="center" width="500" />

 

## Chromosomal interactions

<img src="https://github.com/roderickslieker/CONQUER.test/blob/master/Chromosomal_interactions.png" align="center" width="500" />

 

## Chromatin states

<img src="https://github.com/roderickslieker/CONQUER.test/blob/master/Chromatin_states.png" align="center" width="500" />

 

## eQTLs

<img src="https://github.com/roderickslieker/CONQUER.test/blob/master/eQTLs.png" align="center" width="500" />
 

## Gene expression

<img src="https://github.com/roderickslieker/CONQUER.test/blob/master/Expression.png" align="center" width="500" />
