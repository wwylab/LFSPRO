\name{risk.mpc}
\alias{risk.mpc}
\title{
  Predict the risks of developing multiple primary cancers
}
\description{
Calculate the future cancer risk of each counselee who is either healthy or has had one primary cancer before using the MPC model (this model is not built to predict the cancer risks for patients who have had multiple primary cancer). This function has been integrated in \link{lfspro.mode} to predict cancer risks when mode is "mpc". 
}
\usage{
  risk.mpc(fam.cancer.data, cancer.data, cid, data2, parameter)
}

\arguments{
  \item{fam.cancer.data}{
    Data frame storing family and cancer information. See \link{fam.cancer.data} for details.
  }
  \item{cancer.data}{
  Data frame storing cancer information. See \link{cancer.data} for details.
  }
  \item{cid}{Individual id for the counselee. If you want to estimate multiple samples at the same time, just set counselee.id as a vector of IDs for all the samples.}
  \item{data2}{Converted family and cancer data in a good format}
  \item{parameter}{See \link{parameter.mpc} for details.}
}
\value{
  Matrix, family ID, personal ID and the cancer risks in 5, 10 and 15 years.
}
\references{
Shin SJ, Li J, Ning J, Bodjadzieva J, Strong LC and Wang W. A Bayesian estimation of semiparametric recurrent event model with applications to the penetrance estimation of multiple primary cancers in Li-Fraumeni syndrome. \emph{Biostatistics}
  
}

\author{
  Jasper Chen
}
\examples{
fam.id <- unique(fam.data$fam.id) 
famdata <- fam.data[fam.data$fam.id == fam.id[1], c(1:6)] 
cancerdata <- cancer.data[cancer.data$fam.id == fam.id[1],] 
cid <- c(0,1,5,10)
fam.cancer.data <- combinedata(fam.data, cancer.data)
data.obj <- convert.data(fam.cancer.data)
data.obj1 <- data.obj[[1]]
data.obj2 <- data.obj[[2]]
risk.mpc(fam.cancer.data[[1]], cancerdata, cid, data.obj2[[1]], parameter.mpc)
}

