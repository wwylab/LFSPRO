\name{lkNoneAffect}
\alias{lkNoneAffect}
\title{
Likelihood for un-affected individuals
}
\description{
The probability of observing cancer outcome given genotype for each un-affected individual in the family.
}
\usage{
lkNoneAffect(penetrance, age)
}
\arguments{
  \item{penetrance}{
List. Two component of list indicates the penetrance of male and female (fMX: male, fFX: female). For each component, it is a 110*3 matrix containing the penetrance from age 1 to age 110 for three kinds of genotype: P530 (wild type), P531 (one allele mutation) and P532 (two alleles mutation). See \link{lfspenet.2010} for details.
}
  \item{age}{
The age when the individual is still healthy.
}
}

\value{
 The likelihood (Pr(D|G)) for the un-affected individual at 'age' years old. G, genotype, can be TP530 (wild type), TP531 (one allele of TP53 mutated) and PT532 (both alleles of TP53 mutated).
}

\author{
Gang Peng
}

\seealso{
\code{\link{calLK}}
}
\examples{
lkNoneAffect(lfspenet.2010$fFX, 50)
}

