\name{risk.mpc}
\alias{risk.mpc}
\title{
  Predict the risks of developing multiple primary cancers
}
\description{
Calculate the future cancer risk of each counselee who has had one primary cancer before using the MPC model (this model is not built to predict the cancer risks for patients who have had multiple primary cancer) 
}
\usage{
  risk.mpc(fam.cancer.data, cid, parameter)
}

\arguments{
  \item{fam.cancer.data}{
    Data frame storing family and cancer information. See \link{fam.cancer.data} for details.
  }
  \item{cid}{Individual id for the counselee. If you want to estimate multiple samples at the same time, just set counselee.id as a vector of IDs for all the samples.}
  \item{parameter}{See \link{parameter.mpc} for details.}
}
\value{
  Matrix, family ID, personal ID and the cancer risks in 5, 10 and 15 years.
}
\references{
Shin SJ, Li J, Ning J, Bodjadzieva J, Strong LC and Wang W. A Bayesian estimation of semiparametric recurrent event model with applications to the penetrance estimation of multiple primary cancers in Li-Fraumeni syndrome. \emph{Biostatistics}
  
}

\author{
  Nam Nguyen
}
\examples{
fam.id <- unique(fam.data$fam.id) 
famdata <- fam.data[fam.data$fam.id == fam.id[1], c(1:6)] 
cancerdata <- cancer.data[cancer.data$fam.id == fam.id[1],] 
cid <- c(0,1,5,10)
fam.cancer.data <- combinedata(fam.data, cancer.data)
risk.mpc(fam.cancer.data[[1]], cid, parameter.mpc)
}

