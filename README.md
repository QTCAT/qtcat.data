# QTCAT in action
***This package contains a huge genetic data set of Arabidopsis, only analyse this data if you have at least 6 GB free RAM. Even this can be not enough, if you copy the data around during analysis.***

The purpose of this package is the illustration QTCAT at a real genetic data set with simulated phenotypes.  It allows the user to create phenotype simulation similar to Klasen et al. (under review).

This package includes a data set of 1,307 Arabidopsis accessions genotyped for 214,051 SNP (Horton 2012) and functions to simulate phenotypes on top of this data set.

________________________________________________________________________________

## Example of QTCAT usage

### Simulation of a phenotype
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
                    n.rep = 1,             # No. of replication per individual
                    n.loci = 50,           # No. of loci affecting phenotype
                    eff.dist = "gaussian", # Dist. from which effects are drawn
                    snp.probs = snpp,      # Probability of SNPs to be selected
                    h2 = 0.7)              # Heritability of the trait

```

The SNP data together with phenotype give us a full GWA data set.

### QTCAT analysis
The first step of the analysis is the clustering of all SNP data.  It can take a substantial amount of time, however it has to be done only once per SNP data set.  The package includes clustering result for the Arabidopsis SNP data set which would allows directly moving to the HIT analysis below.

```{R}
#------------------------------------------------------------------------------#
# clustering of the SNPs
# !!! this step can easily take a day or so if done only at one core
# or load the SNP clustering insted
# data(snpclust)
snpclust <- qtcatClust(snp = snp)

```

The actual testing procedure is the hit algorithm. The function uses a phenotype object containing the actual phenotype and eventual covariates (e.g. environment). Second a genotype object is needed containing SNP information and indices generated from the clustering. From this significant QTCs are estimated. The procedure is described in Klasen et al. (under review).

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
                  max.height = .3)     # ignor evrthing less  correlated


#------------------------------------------------------------------------------#
## run the core
hitfit <- qtcatHit(pheno, geno)
(result <- qtcatQtc(hitfit, alpha = 0.05))

```
The results represents all significant QTCs at a given alpha level. As QTCAT is a model selection method it is not giving p-values for the non selected SNPs.

________________________________________________________________________________

__Litriture__

Horton, M. W. et al. Genome-wide patterns of genetic variation in worldwide Arabidopsis thaliana accessions from the RegMap panel. Nat. Genet. 44, 212–216 (2012).
