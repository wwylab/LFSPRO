\name{risk.cs}
\alias{risk.cs}
\title{
  Predict the cancer-specific risk
}
\description{
  Predict the future risks of having breast cancer, sarcoma and other cancers, who never had cancers. 
}
\usage{
  risk.cs(fam.cancer.data, penetrance.all, counselee.id, posterior)
}

\arguments{
  \item{fam.cancer.data}{
    Data frame storing family and cancer information. See \link{fam.cancer.data} for details.
  }
  \item{penetrance.all}{
 List. Two lists indicate the cancer-specific penetrances of male and female (fMX: male, fFX: female), respectively.
For each list, it contains 4 of 110x3 penetrance matrices, row: from age 1 to age 110, column: three
kinds of genotype: P530 (wild type, neither of allele mutated), P531 (one allele mutated in P53) and P532 (two allele mutated in P53). Each matrix represents one penetrance matrix of breast cancer, sarcoma, other cancers, and death, respectively. See \link{lfspenet.cs} for details.
  }
  \item{counselee.id}{Individual id for the counselee. If you want to estimate multiple samples at the same time, just set counselee.id as a vector of IDs for all the samples.}
  \item{posterior}{
  Posterior probability of carrying TP53 mutation}
}
\value{
  List. Each matrix of the list indicates the risks of counselees in 5, 10, 15, 20 years of having breast cancer, sarcoma and other cancers, respectively.  
}

\author{
  Jingxiao Chen
}
\examples{
nloci <- 1 
allef <- list(c(1 - 1.0e-4, 1.0e-4))
counselee.id <- data.frame(fam.data[1:5,1:2]) 
mRate <- 1.0e-6
cancer.data <- data.frame(fam.id = c("fam1", "fam1", "fam1"),
                          id = c(0,2,4),
                          cancer.type = c(2,4,1),
                          diag.age = c(12,50,14))
fam.cancer.data <- combinedata(fam.data, cancer.data) 
cid <- c(0,5,10)
pp <- lfsproC.cs(fam.cancer.data[[1]], lfspenet.cs, 
                 cid, allef, nloci, mRate)
risk.cs(fam.cancer.data[[1]], lfspenet.cs, cid, pp)
}

