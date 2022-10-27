\name{calLK.mpc}
\alias{calLK.mpc}
\title{
  Calculate the likelihood Pr(D|G) of multiple primary cancers.
}
\description{
Calculate the likelihood (Pr(D|G), the probability of disease status D given the genotype G for each individual in the family.
}
\usage{
  calLK.mpc(data1, data2, parameter)
}
\arguments{
  \item{data1}{sex, father ID and mother ID}
  \item{data2}{age of diagnosis, sex and multiple primary cancers index}
  \item{parameter}{
    Parameter estimates for semiparametric recurrent event model of multiple primary cancers. See \link{parameter.mpc} for details.
  }
}

\value{
Return a n*3 matrix. The likelihood of observing the cancer outcome given genotype for each individual in the family. n denotes the total number of individuals in the family.
}
\references{
Shin, S. J., Ning, J., Bojadzieva, J., Strong, L. C., and Wang, W. (2018). Bayesian estimation of a semiparametric recurrent event model with applications to the penetrance estimation of multiple primary cancers in Li-Fraumeni syndrome. \emph{Biostatistics}, 00, 1–16. https://doi.org/10.1093/biostatistics/kxy066
}
\author{
  Seung Jun Shin, Nam Nguyen
}

\seealso{
  \code{\link{lfsproC.mpc}}
}
\examples{
  fam.cancer.data <- combinedata(fam.data, cancer.data)
  data.obj <- convert.data(fam.cancer.data)
  data.obj1 <- data.obj[[1]]
  data.obj2 <- data.obj[[2]]
  calLK.mpc(data.obj1[[2]], data.obj2[[2]], parameter.mpc)
}

