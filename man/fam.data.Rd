\name{fam.data}
\alias{fam.data}
\docType{data}
\title{
Family information data
}
\description{
Family information data used in function LFSpro.
}
\usage{fam.data}
\format{
  A data frame containing 9 variables.
  \describe{
    \item{\code{fam.id}}{family id}
    \item{\code{id}}{individual id}
    \item{\code{fid}}{father id}
    \item{\code{mid}}{mother id}
    \item{\code{gender}}{gender. 0: female, 1: male}
    \item{\code{age}}{If the individual was alive at the last contact date, use the age at that time. otherwise, use the age at death.}
  }
}
\details{
Family id could include characters, but id, father id and mother id must be interger larger than 0.
}
\seealso{
\code{\link{cancer.data}}
}
\examples{
fam.data
}
\keyword{datasets}
