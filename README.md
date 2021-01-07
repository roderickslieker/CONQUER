CONQUER
================
Gerard Bouland, Joline Beulens, Joey Nap, Arno van der Slik, Arnaud
Zaldumbide, Leen ’t Hart and Roderick Slieker
06 januari, 2021

  - [0.1 Change log](#change-log)
  - [0.2 Installation](#installation)
      - [0.2.1 Install the depencies](#install-the-depencies)
      - [0.2.2 Install CONQUER](#install-conquer)
  - [0.3 Overview](#overview)
  - [0.4 Note on `precalculated`
    argument](#note-on-precalculated-argument)
  - [0.5 Note on `pcutoff` argument](#note-on-pcutoff-argument)
  - [0.6 Note on multianalyze](#note-on-multianalyze)
  - [0.7 Citation](#citation)
  - [0.8 Example without multianalyze](#example-without-multianalyze)
  - [0.9 Example with multianalyze](#example-with-multianalyze)
  - [0.10 Figure examples](#figure-examples)
      - [0.10.1 Modules](#modules)
      - [0.10.2 Enrichment](#enrichment)
      - [0.10.3 Pathways shared by tissues](#pathways-shared-by-tissues)
      - [0.10.4 LD](#ld)
      - [0.10.5 Chromosomal interactions](#chromosomal-interactions)
      - [0.10.6 Chromatin states](#chromatin-states)
      - [0.10.7 eQTLs](#eqtls)
      - [0.10.8 Gene expression](#gene-expression)

<img src="https://github.com/roderickslieker/CONQUER.test/blob/master/CONQUER.png" align="right" width="120" />

## 0.1 Change log

*v.1.1.3*

  - Added support for different P-value cut-offs in multi-analyze (see
    below)

*v.1.1.2*

  - Allow users to use precalculated SNPs in GTEx instead of calculating
    them all.

  - Added overlap tab of QTL to show overlap between types of QTLs

## 0.2 Installation

### 0.2.1 Install the depencies

``` r
depp <- c("BioCircos","cluster","ggplot2","enrichR","htmlwidgets",
          "rio","shiny","shinycssloaders","stringr","viridis",
          "DT","coloc","curl","dplyr","grDevices","jsonlite","plotly",
          "shinyjs","reshape2","shinythemes","stats","purrr","readr","UpSetR")

BioDepp <- c("IRanges","BiocGenerics","clusterProfiler","GenomicRanges")


# Check present packages
depp.new<-depp[!(depp%in%installed.packages())]
if (length(depp.new)) {
  install.packages(depp.new)
}
# Bioconductor
BioDepp.new<-BioDepp[!(BioDepp%in%installed.packages())]
if (length(BioDepp.new)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(BioDepp, type="source")
}
# load packages
sapply(depp, library, character.only = TRUE)
sapply(BioDepp, library, character.only = TRUE)
```

### 0.2.2 Install CONQUER

Install the the development version from GitLab:

`install.packages("devtools")`

``` r
# Install CoNQUER.db
install.packages("https://github.com/roderickslieker/CONQUER.db/releases/download/v0.1.2/conquer.db_0.1.2.tar.gz", type="source", repos=NULL)

#Install other two dependencies
devtools::install_github("roderickslieker/CONQUER.d3")
devtools::install_github("roderickslieker/CONQUER")
```

## 0.3 Overview

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

## 0.4 Note on `precalculated` argument

To also allow faster pre-process of SNPs, we also allow users to only
take the precalculated eQTLs from GTEx. Pros: much faster to summarize.
Cons: you will miss interesting QTLs from your results because GTEx only
includes genes with TSSs 1Mb from the SNP. So, when precalculated is
`FALSE` the GTEx API will be used to test the lead SNP against genes in
cis and trans.

## 0.5 Note on `pcutoff` argument

Given that the numbers of individuals per tissue vary in GTEx one may be
interested to analyse the data with a more liberal P-value to adjust for
the strong correlation between the number of eQTLs and the number of
individuals in that dataset. As such the package allows to analyze with
three different settings:

  - `stringent` : only allow eQTLs that match the pValue threshold as
    defined by GTEx for SNP

  - `liberal` : P\<0.001 and P-value at max 2x higher than GTEx’s
    threshold.

  - `veryliberal` : All eQTLs that are nominally significant P\<0.05.

## 0.6 Note on multianalyze

Multianalyze works best if the number of SNPs is large (for example \>50
SNPs). Lower numbers may result in spurious results. Instead one should
focus on the single SNPs.

## 0.7 Citation

**CONQUER: an interactive toolbox to understand functional consequences
of GWAS hits.**

Gerard A Bouland, Joline W J Beulens, Joey Nap, Arno R van der Slik,
Arnaud Zaldumbide, Leen M ’t Hart, Roderick C Slieker

*NAR Genomics and Bioinformatics, Volume 2, Issue 4, December 2020,
lqaa085, <https://doi.org/10.1093/nargab/lqaa085>*

## 0.8 Example without multianalyze

``` r
DIR <- "somedirectory"

library(CONQUER)

summarize(variants = c("rs878521","rs10830963"),
          directory=DIR,
          precalculated = TRUE, 
          multiAnalyze=FALSE,
          token="sometoken",
          tissues=NULL)
```

## 0.9 Example with multianalyze

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
          precalculated = TRUE, 
          token=NULL,
          tissues=c("Pancreas","Muscle_Skeletal","Liver"))


visualize(directory = "somedirectory", SNPs = snps)
```

## 0.10 Figure examples

### 0.10.1 Modules

<img src="https://github.com/roderickslieker/CONQUER.test/blob/master/SunBurst_pathways.png" align="center" width="500" />

 

### 0.10.2 Enrichment

<img src="https://github.com/roderickslieker/CONQUER.test/blob/master/Modules_pathways.png" align="center" width="500" />

 

### 0.10.3 Pathways shared by tissues

<img src="https://github.com/roderickslieker/CONQUER.test/blob/master/Pathways_acrossTissue.png" align="center" width="500" />

 

### 0.10.4 LD

<img src="https://github.com/roderickslieker/CONQUER.test/blob/master/LD_rs11642430.png" align="center" width="500" />

 

### 0.10.5 Chromosomal interactions

<img src="https://github.com/roderickslieker/CONQUER.test/blob/master/Chromosomal_interactions.png" align="center" width="500" />

 

### 0.10.6 Chromatin states

<img src="https://github.com/roderickslieker/CONQUER.test/blob/master/Chromatin_states.png" align="center" width="500" />

 

### 0.10.7 eQTLs

<img src="https://github.com/roderickslieker/CONQUER.test/blob/master/eQTLs.png" align="center" width="500" />
 

### 0.10.8 Gene expression

<img src="https://github.com/roderickslieker/CONQUER.test/blob/master/Expression.png" align="center" width="500" />
