\name{combined.risk.mpc}
\alias{combined.risk.mpc}
\title{
Calculate the future cancer risk using the MPC model
}
\description{
Calculate the future cancer risk of each individual who is either healthy or has only been diagnosed with one primary cancer using our MPC model (this model cannot predict the cancer risks for patients who have had multiple primary cancer).
}
\usage{
  combined.risk.mpc(posterior, risk, counselee.id)
}

\arguments{
  \item{posterior}{
    Genotype likelihood, n*3 matrix.
  }
  \item{risk}{
  Matrix, future cancer risks for each genotype, wild-type and mutated, in 5 years, 10 years and 15 years.
  }
  \item{counselee.id}{
Data frame including two variables: fam.id (family ids of counselees) and id (individual ids of counselees).
  }
}
\value{
  Matrix, family ID, personal ID and the future cancer risks in 5, 10 and 15 years.
}

\author{
  Jingxiao Chen
}
\examples{
counselee.id <- data.frame(fam.id = "fam1", cid=1)
pp <- matrix(c(0.365762982, 0.634046770, 1.902481e-04),nrow=1, ncol=3)
risk <- matrix(c("fam1",1,35,0.008601,0.103538,0.020350,0.220971,0.035787,0.3410216),
               nrow=1, ncol=9)
colnames(risk) <- c("fam.id", "ID", "age","5 years(wildtype)", "5 years(mutation)", 
                    "10 years(wildtype)", "10 years(mutation)", "15 years(wildtype)", 
                    "15 years(mutation)")
combined.risk.mpc(pp, risk, counselee.id)
}

