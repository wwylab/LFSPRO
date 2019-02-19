\name{parameter.mpc}
\alias{parameter.mpc}
\docType{data}
\title{
Estimated parameter for semiparametric recurrent event model of multiple primary cancers 
}
\description{
Parameters used in the multiple primary cancer model.  
}
\usage{parameter.mpc}
\format{
  List
  \describe{
    \item{\code{beta}}{Genotype, sex, cancer status, interaction of genotype and sex, interaction of genotype and cancer status.}
    \item{\code{gamma}}{Non-negative flat priors for the baseline intensity.}
}
}
\references{
Shin, S. J., Ning, J., Bojadzieva, J., Strong, L. C., and Wang, W. (2018). Bayesian estimation of a semiparametric recurrent event model with applications to the penetrance estimation of multiple primary cancers in Li-Fraumeni syndrome. \emph{Biostatistics}, 00, 1–16. https://doi.org/10.1093/biostatistics/kxy066
}
\examples{
parameter.mpc
}
