\name{LFSChompret2015}
\alias{LFSChompret2015}
\title{
The Chompret criteria for Li-Fraumeni syndrome
}
\description{
Use LFS Chompret criteria 2015 version to find out the TP53 mutation carriers in the family.
}
\usage{
LFSChompret2015(fam.data, cancer.data, counselee.id)
}

\arguments{
  \item{fam.data}{
Family information data. See \link{fam.data} for details.
}
  \item{cancer.data}{
Cancer information data. See \link{cancer.data} for details.
}
  \item{counselee.id}{
Data frame including two variables: fam.id (family id of counselees) and id (individual id of counselees).
}
}
\value{
A data frame of three columns:fam.id, id and result. The result column is a vector of boolean variables for all the counselees indicating whether they are TP53 mutation carriers. TRUE: carrier. FALSE: non-carrier. 
}
\references{
Chompret, A., et al. (2001). Sensitivity and predictive value of criteria for p53 germline mutation screening. \emph{J Med Genet}, \bold{38(1)}: 43-47.

Tinat, J., et al. (2009). 2009 version of the Chompret criteria for Li Fraumeni syndrome. \emph{J Clin Oncol}, \bold{27(26)}: e108-109.

Bougeard, G., et al.(2015). Revisiting Li-Fraumeni Syndrome From TP53 Mutation Carriers. \emph{J Clin Oncol}, \bold{33(21)}: p. 2345-52.
}
\author{
Gang Peng
}

\seealso{
\code{\link{LFSClassic}}, \code{\link{lfspro}}, \code{\link{LFSChompret}}
}
\examples{
options(stringsAsFactors = FALSE)
fam.id <- c("fam1","fam2","fam2","fam2","fam2")
id <- c(0,0,2,100,200)
counselee.id <- data.frame(fam.id, id)
LFSChompret2015(fam.data, cancer.data, counselee.id)
}

\keyword{model}
