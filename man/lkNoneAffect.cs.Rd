\name{lkNoneAffect.cs}
\alias{lkNoneAffect.cs}
\title{
Likelihood for un-affected individuals for cancer-specific model
}
\description{
The probability of observing cancer outcome given genotype for each un-affected individual in the family in the cancer-specific model. 
}
\usage{
lkNoneAffect.cs(penetrance, age)
}
\arguments{
  \item{penetrance}{
List. Two lists indicate the cancer-specific penetrances of male and female (fMX: male, fFX: female), respectively.
For each list, it contains 4 of 110x3 penetrance matrices, row: from age 1 to age 110, column: three
kinds of genotype: P530 (wild type, neither of allele mutated), P531 (one allele mutated in P53) and P532 (two alleles mutated in P53). Each matrix represents one penetrance matrix of breast cancer, sarcoma, other cancers, and death, respectively. See \link{lfspenet.cs} for details.
  }
  \item{age}{
    The age when the individual is still healthy.
  }
}

\value{
  The likelihood (Pr(D|G)) for the un-affected individual at 'age' years old. G, genotype, can be TP530 (wild type), TP531 (one allele of TP53 mutated) and TP532 (both alleles of TP53 mutated).
}

\author{
  Seung Jun Shin
}

\seealso{
  \code{\link{calLK.cs}}
}
\examples{
  lkNoneAffect.cs(lfspenet.cs$fFX, 50)
}

