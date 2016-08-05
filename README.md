
# QTCAT workflow

***Attention: the package contains a huge genetic data set, if you would like to analyse it, make sure you have at least 8 GB free RAM.***

The purpose of this package is to illustration [QTCAT](https://github.com/QTCAT/qtcat) with real genetic data in combination with simulated phenotypes (similar to Klasen et al.).  It allows a comparison of the results to the true (simulated) gene effects.


This package includes a data set of 1,307 Arabidopsis accessions genotyped for 214,051 SNPs (Horton et al. 2012) and functions to simulate phenotypes on top of this data.

________________________________________________________________________________

## QTCAT in action
### Installation:
The package can be installed from an R console via devtools.

```{R}
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
snpp <- snpProb(gene.pos = t(genepos), snp.pos = getPos(snp))

# phenotype with 50 loci
pdat <- normalPheno(snp = snp,             # SNP data object
                    n.rep = 1,             # No. of replication per individuals
                    n.loci = 50,           # No. of loci affecting phenotype
                    eff.dist = "gaussian", # Dist. from which effects are drawn
                    snp.probs = snpp,      # Probability of SNPs to be selected
                    h2 = 0.7)              # Heritability of the trait

```

The data are in the following steps analysed in order to find at least some of the phenotype effecting loci we have simulated.

### QTCAT analysis:
The first step of the analysis is a hierarchical clustering of all SNPs.  It can take a substantial amount of time, however, it has to be done only once per SNP data set.  The package includes clustering result for the Arabidopsis SNP data, which allows you moving directly forward to the HIT analysis below.

```{R}
#------------------------------------------------------------------------------#
# clustering of the SNPs
# data(snpclust) # instead of running the SNP clustering you can load pre calculated result 
snpclust <- qtcatClust(snp = snp)

```

The actual testing procedure is the HIT algorithm, the function uses a phenotype object, containing the actual phenotype and eventual covariates (e.g. environment).  Furthermore a genotype object is needed which contains the SNPs and indices generated from the clustering. 

Significant QTCs at a specific alpha level can be obtained from the HIT results.  The procedure is described in more detail in Klasen et al.

```{R}
#------------------------------------------------------------------------------#
# create a phenotype object for the HIT analysis
pheno <- qtcatPheno(names = pdat$pheno$id,
                    pheno = pdat$pheno$pheno,
                    family = "gaussian")

#------------------------------------------------------------------------------#
# create a genotype object for the HIT analysis
geno <- qtcatGeno(snp = snp,
                  snpClust = snpclust,
                  min.absCor = 0.3)   # ignore everything less correlated than this value


#------------------------------------------------------------------------------#
## run the core analysis, the HIT
hitfit <- qtcatHit(pheno, geno)
(result <- qtcatQtc(hitfit, alpha = 0.05))

```

QCAT is a model selection method, in contrast to multiple testing methods it is not giving p-values for non significant SNPs.  p-values are in addition to the common influences like sample size, effect size, and variance also dependent at the correlation in the data. Interpretation about the importance of a loci only based at the size of the p-values are in this case even more misleading than usual.

________________________________________________________________________________

__Litriture__

Klasen et al. A multi-marker association method for genome-wide association studies without the need for population structure correction (under review).

Horton et al. Genome-wide patterns of genetic variation in worldwide Arabidopsis thaliana accessions from the RegMap panel. Nat. Genet. 44, 212–216 (2012).

________________________________________________________________________________

[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg)](https://www.gnu.org/licenses/gpl-2.0.html)
&copy; 2015 JR Klasen
