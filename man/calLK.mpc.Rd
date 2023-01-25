\name{calLK.mpc}
\alias{calLK.mpc}
\title{
  Calculate the likelihood Pr(D|G) of multiple primary cancers.
}
\description{
Calculate the likelihood (Pr(D|G), the probability of disease status D given the genotype G for each individual in the family.
}
\usage{
  calLK.mpc(data1, data2, parameter, mut.info = TRUE)
}
\arguments{
  \item{data1}{
    Family data that describe how the individuals are related, for \emph{ONE FAMILY ONLY}.
  }
  \item{data2}{
    A comprehensive dataset that describes each individual in details, including age, gender (1 for male, 2 for female), mutation status (0 for wildtype, 1 for mutation, NA if unknown), vital status (A if alive, D if dead), and cancer history (cancer type and age at diagnosis). Both data1 and data2 can be generated from family and cancer datasets using \link{combinedata} and \link{convert.data}. See the example below.
  }
  \item{parameter}{
    Parameter estimates for semiparametric recurrent event model of multiple primary cancers. See \link{parameter.mpc} for details.
  }
  \item{mut.info}{
    Set to TRUE (default) if you want to use known genetic testing results of family members in the calculation of mutation probabilities. Set to FALSE otherwise. 
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
calLK.mpc(data.obj1[[1]], data.obj2[[1]], parameter.mpc)
}
