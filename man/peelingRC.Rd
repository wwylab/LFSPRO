\name{peelingRC}
\alias{peelingRC}
\title{
Peeling interface in R.
}
\description{
Peeling (Elston-Stewart algorithm) is the key function in LFSpro. We implmented it in C++ to make it fast. peelingRC is used to link the peeling algorithm in C++ version to R.
}
\usage{
peelingRC(allef, LIK, ped, counselee.id, nloci = 1, mRate = 0)
}
\arguments{
  \item{allef}{
List. Allele frequency for each locus/gene. If there is only one gene and two alleles in the gene (allele frequency is 0.1 and 0.9), allef = list(c(0.1,0.9)), If there are two genes,two alleles (allele frequency is 0.1 and 0.9) for gene 1  and three alleles (allele frequncy is 0.2, 0.2 and 0.6) for gene 2. allef = list(c(0.1,0.9),c(0.2,0.2,0.6))
}
  \item{LIK}{
Matrix, likelihood, Pr(D|G), for three kinds of genotype for all the individuals in the family. D: healthy status. G: genotype.
}
  \item{ped}{
Pedigree structure. A data frame with four varaibles: ID(individual id), Gender (gender, 0: female, 1: male), FatherID(father id) and MotherID(mother id).
}
  \item{counselee.id}{
Individual id for the counselee. If you want to estimate multiple samples at the same time, just set counselee.id as a vector of IDs for all the samples you want to estimate.
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
\examples{
  fam.cancer.data <- CombineData(fam.data, cancer.data)
  allef <- allef.g <- list(c(0.9997,0.0003))
  mRate.g <- 6e-05
  FamData <- fam.cancer.data[[1]]
  lik <- calLK(FamData, LFSpenet.2010)
  
  ##################################
  # convert data
  ##################################
  counselee.id <- FamData$id[1]
  id <- as.integer(FamData$id)
  fid <- as.integer(FamData$fid)
  mid <- as.integer(FamData$mid)
  counselee.id <- as.integer(counselee.id)
  
  if(min(id)==0){
    id <- id+1
    fid <- fid+1
    mid <- mid+1
    counselee.id <- counselee.id+1
  }
  fid[is.na(fid)] <- 0
  mid[is.na(mid)] <- 0
  
  ped <- data.frame(ID=id, Gender = FamData$gender,
                    FatherID = fid, 
                    MotherID = mid, 
                    stringsAsFactors=FALSE)
  
  peelingRC(allef, lik, ped, counselee.id, 1, mRate.g)
}

\keyword{model}
