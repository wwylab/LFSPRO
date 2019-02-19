\name{lfspenet.cs.nodeath}
\alias{lfspenet.cs.nodeath}
\docType{data}
\title{
  Penetrance table for cancer-specific mutation with three competing risks (breast cancer, sarcoma and other cancers).   Default penetrance table in predicting cancer-specific risks. 
}
\description{
  Penetrance table
}
\usage{lfspenet.cs.nodeath}
\format{
List. Two lists indicate the cancer-specific penetrances of male and female (fMX: male, fFX: female), respectively.
For each list, it contains 3 of 110x3 penetrance matrices, row: from age 1 to age 110, column: three
kinds of genotype: P530 (wild type, neither of allele mutated), P531 (one allele mutated in P53) and P532 (two allele mutated in P53). Each matrix represents one penetrance matrix of breast cancer, sarcoma and other cancers, respectively.
}
\examples{
  lfspenet.cs.nodeath
}
