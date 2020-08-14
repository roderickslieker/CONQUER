*Gerard Bouland, Joline Beulens, Joey Nap, Arno van der Slik, Arnaud
Zaldumbide, Leen ’t Hart and Roderick Slieker*    

CONQUER collects data from multiple resources. This page gives the
overview of the used resources.    

Linkage disequilibrium (LD)
---------------------------

 

All SNPs in LD are retrieved from ENSEMBL REST API. The lead and LD
(r2&gt;0.8) SNPs are used in all subsequent analyses and lookups.    

QTL data
--------

  \#\#\# eQTLs   Expression QTLs were obtained from GTEx. For all tested
SNPs, CONQUER tests the lead SNP against all genes in cis and trans
(those with an interaction) for all 44 tissues in GTEx. The
precalculated eQTLs are used in the colocalization part.
