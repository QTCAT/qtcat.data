
# QTCAT workflow

***Attention: the package contains a large genetic data set, in order to analyse it, make sure you have at least 8 GB of free RAM.***

The purpose of this package is to illustrate [QTCAT](https://github.com/QTCAT/qtcat)s functionality (Klasen et al. 2016), on the basis of real genetic data and simulated phenotypes.  This approach allows a comparison of QTCATs results to the true (simulated) gene effects.


The package contains a data set of 1,307 Arabidopsis accessions which are genotyped for 214,051 SNPs (Horton et al. 2012) and functions to simulate a phenotype on basis of this data.

________________________________________________________________________________

## QTCAT in action
### Installation:
The package can be installed from an R console via ``devtools``.

```{R}
################################################################################
# installation
# install.packages("devtools")
devtools::install_github("QTCAT/qtcat") 
devtools::install_github("QTCAT/qtcat.data") 

```

### Simulating a phenotype:

```{R}
require(qtcat)
require(qtcat.data)

#------------------------------------------------------------------------------#
# load data
data("snp")
data("genepos")

#------------------------------------------------------------------------------#
# Simulate a phenotype

# gene density to SNP probabilities
snpp <- snpProb(gene.pos = genepos, snp.pos = snpInfo(snp)[, 1:2])

# phenotype with 50 loci
pdat <- normalPheno(snp = snp,             # SNP data object
                    n.rep = 1,             # No. of replications per individual
                    n.loci = 50,           # No. of loci affecting the phenotype
                    eff.dist = "gaussian", # Dist. from which the effects are drawn
                    snp.probs = snpp,      # Probability for SNPs to be selected as causal
                    h2 = 0.7)              # Heritability of the simulated phenotype

```

In order to find at least some of the simulated loci, the data are analysed in the following steps using QTCAT.

### QTCAT analysis:
The first step of the analysis is a hierarchical clustering of all SNPs.  This can take some time (it has to be done only once per SNP data set).  In order to allow you to move directly forward to the HIT analysis below, the package includes a clustering result for the example data set.

```{R}
#------------------------------------------------------------------------------#
# clustering of the SNPs
# data(snpclust) # instead of running the clustering, you can load a pre-calculated result 
snpclust <- qtcatClust(snp = snp)

```

The main tests is done using the HIT algorithm. The ``qtcatHit``-function uses a phenotype object, containing the actual phenotype and eventual covariates (e.g. environment, replications, ...).  Furthermore a genotype object is needed which contains the SNP data and SNP clustering information. 

Finally significant Quantitative Trait Clusters (QTCs) at a specific alpha level can be summarized.

```{R}
#------------------------------------------------------------------------------#
# create a phenotype object for the HIT analysis
pheno <- qtcatPheno(names = pdat$pheno$id,
                    pheno = pdat$pheno$pheno,
                    family = "gaussian")

#------------------------------------------------------------------------------#
# create a genotype object for the HIT analysis
geno <- qtcatGeno(snp = snp,
                  snpClust = snpclust)


#------------------------------------------------------------------------------#
# run the core analysis, the HIT
hitfit <- qtcatHit(pheno, geno)
(result <- qtcatQtc(hitfit, alpha = 0.05))

```

QCAT is a model selection method, in contrast to multiple testing methods, it is not giving p-values for non significant SNPs.  p-values are in addition to the common influences, like sample size, effect size, and variance, also dependent at the correlation between SNPs.  Interpretation about the importance of a loci only based at the size of the p-value is even less meaningful than in scenarios without correlated variables.

Visualization of the results can by done with ``plot Qtc``, which shows the detected QTCs at there genomic position.  ``plotSelFreq`` gives the selection frequency per marker in the repeated LASSO selection.

```{R}
#------------------------------------------------------------------------------#
# Plot results
plotQtc(hitfit, alpha = 0.05, pch = 20)
plotSelFreq(hitfit, pch = 20)

################################################################################
```

The full procedure is described in more detail in Klasen et al. (2016).

________________________________________________________________________________

### Literature

**Klasen, J. R. et al. (2016)**. *A multi-marker association method for genome-wide 
association studies without the need for population structure correction*. Nature 
Communications. [Paper](http://www.nature.com/articles/ncomms13299)

**Horton et al. (2012)**. *Genome-wide patterns of genetic variation in worldwide Arabidopsis thaliana 
accessions from the RegMap panel*. Nat. Genet.

________________________________________________________________________________

[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg)](https://www.gnu.org/licenses/gpl-2.0.html)
&copy; 2016 JR Klasen
