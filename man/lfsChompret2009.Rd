\name{lfsChompret2009}
\alias{lfsChompret2009}
\title{
The Chompret criteria for Li-Fraumeni syndrome
}
\description{
Use the Chompret criteria to identify TP53 germline mutation carriers in the family.
}
\usage{
lfsChompret2009(fam.data, cancer.data, counselee.id)
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
A data frame of three columns: fam.id, id and result. The result column is a vector of TRUE/FALSE for all counselees indicating whether they are predicted to be TP53 mutation carriers. TRUE: carrier. FALSE: non-carrier.
}
\references{
Chompret, A., et al. (2001). Sensitivity and predictive value of criteria for p53 germline mutation screening. \emph{J Med Genet}, \bold{38(1)}: 43-47.

Tinat, J., et al. (2009). 2009 version of the Chompret criteria for Li Fraumeni syndrome. \emph{J Clin Oncol}, \bold{27(26)}: e108-109; author reply e110.
}
\author{
Gang Peng
}

\seealso{
\code{\link{lfsClassic}}, \code{\link{lfspro}}
}
\examples{
fam.id <- c("fam1","fam2","fam2","fam2","fam2")
id <- c(0,0,2,100,200)
counselee.id <- data.frame(fam.id, id)
lfsChompret2009(fam.data, cancer.data, counselee.id)
}

