\name{lfsClassic}
\alias{lfsClassic}
\title{
The Classic criteria for Li-Fraumeni syndrome
}
\description{
Use the Classic criteria to identify TP53 germline mutation carriers in the family.
}
\usage{
lfsClassic(fam.data, cancer.data, counselee.id)
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
A data frame of three columns:fam.id, id and result. The result column is a vector of TRUE/FALSE. for all counselees indicating whether they are predicted to be TP53 mutation carriers. TRUE: carrier. FALSE: non-carrier.
}
\references{
Li, F. P., et al. (1988). A cancer family syndrome in twenty-four kindreds. \emph{Cancer Res}, \bold{48(18)}: 5358-5362.
}
\author{
Gang Peng
}
\seealso{
\code{\link{lfsChompret2015}}, \code{\link{lfspro}}
}
\examples{
options(stringsAsFactors = FALSE)
fam.id <- c("fam1","fam2","fam2","fam2","fam2")
id <- c(0,0,2,100,200)
counselee.id <- data.frame(fam.id, id)
lfsClassic(fam.data, cancer.data, counselee.id)
}

