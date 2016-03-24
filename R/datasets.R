#' An Arabidopsis thaliana SNP data object.
#'
#' A dataset of 1307 Arabidopsis thaliana accesions genotyped and 214051 SNPs. 
#' The data required 285 Mb of memory only to be loaded, so make sure your 
#' memory is big enough. 
#'
#' @format SnpData object:
#' \describe{
#'   \item{snp}{An object of S4 class \linkS4class{snpData}}
#' }
#' @source Horton, Matthew W., Angela M. Hancock, Yu S. Huang, Christopher 
#' Toomajian, Susanna Atwell, Adam Auton, N. Wayan Muliyati, et al. 2012. 
#' Genome-Wide Patterns of Genetic Variation in Worldwide Arabidopsis thaliana 
#' Accessions from the RegMap Panel. Nature Genetics 44 (2): 212-16. 
#' \url{http://doi.org/10.1038/ng.1042.} 
#' 
#' Data are avalebel from: 
#' \url{http://bergelson.uchicago.edu/wp-content/uploads/2015/04/call_method_75.tar.gz}
"snp"


#' Gene positions of Arabidopsis thaliana.
#'
#' A dataset of of gene positions of Arabidopsis thaliana genes. 
#'
#' @format A data frame:
#' \describe{
#'   \item{chr}{Chromosomes were genes are located}
#'   \item{pos}{Position at the chromosomes were genes are located}
#' }
#' 
#' @source Data are avalebel from: 
#' \url{ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff}
"genepos"


#' A SNP cluster object of the Arabidopsis thaliana data.
#'
#' A SNP cluster object of the Arabidopsis thaliana data.
#'
#' @format Clustering object:
#' \describe{
#'   \item{snpclust}{A clustering of all SNPs}
#' }
#' 
"snpclust"
