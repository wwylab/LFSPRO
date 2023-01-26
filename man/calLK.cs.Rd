\name{calLK.cs}
\alias{calLK.cs}
\title{
  Calculate the likelihood Pr(D|G) for cancer specific.
}
\description{
Calculate the likelihood (Pr(D|G), probability of disease status D given the genotype G for each individual in the family.
}
\usage{
  calLK.cs(ind.data, penetrance.all, mut.info)
}
\arguments{
  \item{ind.data}{
    A comprehensive dataset that describes each individual in details, including age, gender (1 for male, 2 for female), mutation status (0 for wildtype, 1 for mutation, NA if unknown), vital status (A if alive, D if dead), and cancer history (cancer type and age at diagnosis). Both data1 and data2 can be generated from family and cancer datasets using \link{combinedata} and \link{convert.data}. See the example below.
  }
  \item{penetrance.all}{
    The penetrance of three genotype (wild type, one copy and two copies of TP53 mutated alleles) for male and female in the population. See \link{lfspenet.cs} for details.
  }
  \item{mut.info}{
    Set to TRUE (default) if you want to use known genetic testing results of family members in the calculation of mutation probabilities. Set to FALSE otherwise. 
  }
}

\value{
Return a n*3 matrix. The likelihood of observing the cancer outcome given genotype for each individual in the family. n denotes the total number of individuals in the family.
}
\references{
Shin, S. J., Yuan, Y., Strong, L. C., Bojadzieva, J., and Wang, W. (2018). Bayesian Semiparametric Estimation of Cancer-Specific Age-at-Onset Penetrance With Application to Li-Fraumeni Syndrome. \emph{Journal of the American Statistical Association}, 1–12. https://doi.org/10.1080/01621459.2018.1482749
}
\author{
  Seung Jun Shin, Nam Nguyen
}

\seealso{
  \code{\link{lkNoneAffect.cs}}, \code{\link{lfsproC.cs}}, \code{\link{peelingRC}}
}
\examples{
fam.cancer.data <- combinedata(fam.data, cancer.data)
data.obj <- convert.data(fam.cancer.data)
data.obj1 <- data.obj[[1]]
data.obj2 <- data.obj[[2]]
calLK.cs(data.obj2[[1]], lfspenet.cs)
}
