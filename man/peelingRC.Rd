\name{peelingRC}
\alias{peelingRC}
\title{
Peeling interface in R.
}
\description{
peelingRC is based on Elston-Stewart algorithm and it is a low level function has been integrated in \link{lfsproC}, \link{lfsproC.cs} and \link{lfsproC.mpc}. We implemented the function in C++ to improve computation so peelingRC linked the peeling algorithm in C++ version to R.
}
\usage{
peelingRC(allef, LIK, ped, counselee.id, nloci = 1, mRate = 0)
}
\arguments{
  \item{allef}{
List. Allele frequency for each locus/gene. If there is only one gene and two alleles in the gene (allele frequency is 0.1 and 0.9), allef = list(c(0.1,0.9)), If there are two genes,two alleles (allele frequency is 0.1 and 0.9) for gene 1  and three alleles (allele frequncy is 0.2, 0.2 and 0.6) for gene 2. allef = list(c(0.1,0.9),c(0.2,0.2,0.6))
}
  \item{LIK}{
Matrix, likelihood, Pr(D|G), for three genotypes for all individuals in the family. D: healthy status. G: genotype.
}
  \item{ped}{
Pedigree structure. A data frame with four varaibles: ID (individual id), Gender (sex, 0: female, 1: male), FatherID (father id) and MotherID (mother id).
}
  \item{counselee.id}{
Individual id for the counselee. If you want to estimate multiple samples at the same time, just set counselee.id as a vector of IDs for all the samples.
}
  \item{nloci}{
Number of loci/genes in the model.
}
  \item{mRate}{
Mutation rate.
}
}
\details{
One family a time.
}
\value{
The posterior probability (Pr(G|D)) for each counselee
}
\references{
Elston, R. C., Stewart, J. (1971) A general model for the genetic analysis of pedigree data. \emph{Hum Hered.}, \bold{21}, 523-542.
}
\author{
Gang Peng
}

\seealso{
\code{\link{lfsproC}}
}
