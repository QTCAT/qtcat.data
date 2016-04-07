#' Estimate sample probability for each SNP from gene density
#' 
#' @param gene.pos a matrix of two rows. The first row contains the 
#' chromosome, the second row contains the position of a gene (e.g. start in bp).
#' @param snp.pos a matrix of two rows. the first row contains the 
#' chromosome, the second row contains the position of a SNP  bp.
#' @param bw see \code{\link[stats]{density}}.
#' 
#' @importFrom stats density approxfun
#' @export
snpProb <- function(gene.pos, snp.pos, bw = 10000) {
  stopifnot(nrow(gene.pos) == 2)
  stopifnot(nrow(snp.pos) == 2)
  stopifnot(all(!is.null(rownames(snp.pos))))
  stopifnot(all(!is.null(rownames(gene.pos))))
  rownames(snp.pos) <- substring(tolower(rownames(snp.pos)), 1, 3)
  rownames(gene.pos) <- substring(tolower(rownames(gene.pos)), 1, 3)
  stopifnot(all(rownames(snp.pos) == c("chr", "pos")))
  stopifnot(all(rownames(gene.pos) == c("chr", "pos")))
  allchr <-  sort(unique(snp.pos["chr", ]))
  stopifnot(all(allchr %in% unique(gene.pos["chr", ])))
  snp.probs <- c()
  for (chr in sort(unique(snp.pos["chr", ]))) {
    gene.chr <- gene.pos[2, gene.pos[1, ] == chr]
    snp.chr <- snp.pos[2, snp.pos[1, ] == chr]
    genedens <- density(gene.chr, bw = bw)
    # a function for linear approxiamtion
    density.approx <- approxfun(genedens$x, genedens$y * length(gene.chr))
    snp.probs <- c(snp.probs, density.approx(snp.chr))
  }
  snp.probs <- snp.probs / sum(snp.probs)
  snp.probs
}


#' Simulates a normal distributed phenotype with effects either normal 
#' or gamma distributed.
#' 
#' @param snp An object of S4 class \linkS4class{snpData}.
#' @param n.rep Number of replication per individual. For the replication no 
#' effects at the phenotype are simulated.
#' @param n.loci Number of loci with effect at the phenotype. 
#' @param h2 Trait heritability.
#' @param snp.probs A vector of probabilities for each SNP to be selected 
#' as a causal loci.
#' @param eff.dist Distribution from which the effects are drown.
#' @param shape See \code{\link[stats]{rgamma}} (only considered if 'eff.dist = 
#' gamma').
#' 
#' @importFrom stats rnorm rgamma
#' @importFrom qtcat getPos as.matrix
#' @export
normalPheno <- function(snp, n.rep = 1, n.loci = 50, h2 = .5, snp.probs = NULL,
                        eff.dist = c("gaussian", "gamma"), shape = .5) {
  eff.dist <- match.arg(eff.dist, c("gaussian", "gamma"))
  if (eff.dist == "gamma") {
    eff.sign <- sample(c(-1, 1), size = n.loci, replace = TRUE) 
    loci.eff <- rgamma(n.loci, shape) * eff.sign
  } else {
    loci.eff <- rnorm(n.loci)
  }
  n.snp <- ncol(snp)
  snp.inx <- sort(sample.int(n.snp, n.loci, prob = snp.probs))
  snp.loci <- as.matrix(snp[, snp.inx])
  geno <- rep(drop(snp.loci %*% loci.eff), each = n.rep)
  res <- rres(length(geno), var(geno), h2)
  id <- rep(rownames(snp.loci), each = n.rep)
  pheno <- data.frame(id = id, pheno = geno + res, geno = geno, 
                      stringsAsFactors = FALSE)
  pos <- t(getPos(snp)[, snp.inx])
  effects <- data.frame(colnames(snp)[snp.inx], pos, loci.eff, 
                        stringsAsFactors = FALSE)
  colnames(effects) <- c("snp", "chr", "pos", "effects")
  rownames(effects) <- NULL
  list(pheno = pheno, effects = effects, h2 = h2)
}


#' Random residuals which match the expected heritability
#' 
#' @param n Number of residuals.
#' @param Va Genetic variance.
#' @param h2 Expected heritability.
#' 
#' @keywords internal
rres <- function(n, Va = 1, h2 = .5) {
  Ve <- Va / h2 - Va
  res <- rnorm(n)
  res <- 1 / sd(res) * sqrt(Ve) * res
  res
}


# A first attempt to implement case control data, but it is still a bit too naive, as it 
# ignores prevalence and thereby the liability scale.
# Tenesa & Haley (2013) The heritability of human disease: estimation, 
# uses and abuses. Nat Rev Genet
#
# #' Simulates a Bernoulli distributed phenotype with effects either normal
# #' or gamma distributed.
# #'
# #' @param snp An object of S4 class \linkS4class{snpData}.
# #' @param n.rep Number of replication per individual. For the replication no
# #' effects at the phenotype are simulated.
# #' @param n.loci Number of loci with effect at the phenotype.
# #' @param h2 Trait heritability.
# #' @param snp.probs A verctor of probertibilities for each SNP to be selected
# #' as a causal loci.
# #' @param eff.dist Distribution from which the effects are drown.
# #' @param shape See \code{\link[stats]{rgamma}} (only considert if 'eff.dist =
# #' gamma').
# #'
# #' @importFrom stats rnorm rgamma rbinom
# #' @importFrom qtcat getPos as.matrix
# #' @export
# bernoulliPheno <- function(snp, n.rep = 1, n.loci = 50,  h2, snp.probs = NULL,
#                         eff.dist = c("gaussian", "gamma"), shape = .5) {
#   eff.dist <- match.arg(eff.dist, c("gaussian", "gamma"))
#  if (eff.dist == "gamma") {
#     eff.sign <- sample(c(-1, 1), size = n.loci, replace = TRUE)
#     loci.eff <- rgamma(n.loci, shape) * eff.sign
#   } else {
#     loci.eff <- rnorm(n.loci)
#   }
#   n.snp <- ncol(snp)
#   snp.inx <- sort(sample.int(n.snp, n.loci, prob = snp.probs))
#   snp.loci <- as.matrix(snp[, snp.inx])
#   geno.raw <- rep(drop(snp.loci %*% loci.eff), each = n.rep)
#   geno.raw <- geno.raw - mean(geno.raw)  # scale intercept gives 50% cases and contols
#   geno.scale <- bernoulliscaling(geno.raw, h2)
#   geno <- 1 / (1 + exp(-(geno.raw * geno.scale)))
#   id <- rep(rownames(snp.loci), each = n.rep)
#   pheno <- data.frame(id = id, pheno = rbinom(length(geno), 1, geno), geno = geno,
#                       stringsAsFactors = FALSE)
#   pos <- t(getPos(snp)[, snp.inx])
#   effects <- data.frame(colnames(snp)[snp.inx], pos, loci.eff * geno.scale,
#                         stringsAsFactors = FALSE)
#   colnames(effects) <- c("snp", "chr", "pos", "effects")
#   rownames(effects) <- NULL
#   list(pheno = pheno, effects = effects, h2 = h2)
# }
# 
# 
# #' Sacling factore for geno inorder to get the expected heritability
# #' from a Bernoulli distribution
# #'
# #' @param geno Vector of genetic values.
# #' @param h2 Expected heritability.
# #'
# #' @importFrom stats optimise rbinom
# #' @keywords internal
# bernoulliscaling <- function(geno, h2 = .5) {
#   out <- optimise(function(scale, geno, h2){
#     geno <- geno * scale
#     geno2 <- 1 / (1 + exp(-geno))
#     y <- rbinom(length(geno2), 1, geno2)
#     abs(var(geno2) / var(y) - h2)
#   }, interval = c(0, 150), geno = geno, h2 = h2)
#   out$minimum
# }


# #' Variance if true mean is known as in simulations 
# #' as we should always have many observationsi this is a bit pedantic 
# Var <- function(x) sum((x - mean(x))^2) / length(x)
