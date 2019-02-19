\name{lfspenet.cs.death}
\alias{lfspenet.cs.death}
\docType{data}
\title{
Penetrance table for cancer-specific mutation with four competing risks (breast cancer, sarcoma, other cancers and death).
}
\description{
Penetrance table
}
\usage{lfspenet.cs.death}
\format{
List. Two lists indicate the cancer-specific penetrances of male and female (fMX: male, fFX: female), respectively.
For each list, it contains 4 of 110x3 penetrance matrices, row: from age 1 to age 110, column: three
kinds of genotype: P530 (wild type, neither of allele mutated), P531 (one allele mutated in P53) and P532 (two allele mutated in P53). Each matrix represents one penetrance matrix of breast cancer, sarcoma, other cancers, and death, respectively.
}
\examples{
lfspenet.cs.death
}
