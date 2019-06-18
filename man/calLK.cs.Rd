\name{calLK.cs}
\alias{calLK.cs}
\title{
  Calculate the likelihood Pr(D|G) for cancer specific.
}
\description{
Calculate the likelihood (Pr(D|G), probability of disease status D given the genotype G for each individual in the family.
}
\usage{
  calLK.cs(fam.cancer.data, penetrance.all)
}
\arguments{
  \item{fam.cancer.data}{
    Data including family information and cancer information. See \link{fam.cancer.data} for details.
  }
  \item{penetrance.all}{
    The penetrance of three genotype (wild type, one copy and two copies of TP53 mutated alleles) for male and female in the population. See \link{lfspenet.cs} for details.
  }
}

\value{
Return a n*3 matrix. The likelihood of observing the cancer outcome given genotype for each individual in the family. n denotes the total number of individuals in the family.
}
\references{
Shin, S. J., Yuan, Y., Strong, L. C., Bojadzieva, J., and Wang, W. (2018). Bayesian Semiparametric Estimation of Cancer-Specific Age-at-Onset Penetrance With Application to Li-Fraumeni Syndrome. \emph{Journal of the American Statistical Association}, 1–12. https://doi.org/10.1080/01621459.2018.1482749
}
\author{
  Seung Jun Shin
}

\seealso{
  \code{\link{lkNoneAffect.cs}}, \code{\link{lfsproC.cs}}, \code{\link{peelingRC}}
}
\examples{
cancer.data$cancer.type <- c(2,4,3,7,2,2,6,50,50,2,4,50,50,4,50) #cancer types
fam.cancer.data <- combinedata(fam.data, cancer.data)
calLK.cs(fam.cancer.data[[1]], lfspenet.cs)
}

