\name{calLK}
\alias{calLK}
\title{
Calculate the likelihood Pr(D|G).
}
\description{
Calculate the likelihood (Pr(D|G), probability of healthy status D given the genotype G) for each individual in the family.
}
\usage{
calLK(fam.cancer.data, penetrance.all)
}
\arguments{
  \item{fam.cancer.data}{
Data including family information and cancer information. See \link{fam.cancer.data} for details.
}
  \item{penetrance.all}{
The penetrance of three kinds of genotype (wild type, one copy TP53 mutation and two copy TP53 mutation) for male and female in the population. See \link{LFSpenet.2010} for details.
}
}

\value{
Return a n*3 matrix. The likelihood of three kinds of genotype for each individual in the family. n denotes the number of individuals in the family.
}

\author{
Gang Peng
}

\seealso{
\code{\link{lkNoneAffect}}, \code{\link{lfsproC}}, \code{\link{peelingRC}}
}
\examples{
fam.cancer.data <- CombineData(fam.data, cancer.data)
calLK(fam.cancer.data[[1]], LFSpenet.2010)
}

\keyword{models}
