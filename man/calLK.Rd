\name{calLK}
\alias{calLK}
\title{
Calculate the likelihood Pr(D|G).
}
\description{
Calculate the likelihood (Pr(D|G), probability of disease status D given the genotype G for each individual in the family.
}
\usage{
calLK(fam.cancer.data, penetrance.all)
}
\arguments{
  \item{fam.cancer.data}{
Data including family information and cancer information. See \link{fam.cancer.data} for details.
}
  \item{penetrance.all}{
The penetrance of three genotype (wild type, one copy and two copies TP53 mutated alleles) for male and female in the population. See \link{lfspenet.2010} for details.
}
}

\value{
Return a n*3 matrix. The likelihood of observing the cancer outcome given genotype for each individual in the family. n denotes the total number of individuals in the family.
}
\references{
Peng, G., Bojadzieva, J., Ballinger, M. L., Li, J., Blackford, A. L., Mai, P. L., Savage, S. A., Thomas, D. M., Strong, L. C., … Wang, W. (2017). Estimating TP53 Mutation Carrier Probability in Families with Li-Fraumeni Syndrome Using LFSPRO. \emph{Cancer Epidemiology, Biomarker and Prevention}, \bold{26(6)}, 837-844. DOI: 10.1158/1055-9965.EPI-16-0695
}

\author{
Gang Peng
}

\seealso{
\code{\link{lkNoneAffect}}, \code{\link{lfsproC}}, \code{\link{peelingRC}}
}
\examples{
cancer.data$cancer.type <- c(2,4,3,7,2,2,6,50,50,2,4,50,50,4,50) #cancer types
fam.cancer.data <- combinedata(fam.data, cancer.data)
calLK(fam.cancer.data[[1]], lfspenet.2010) #likelihood in family 1
}
