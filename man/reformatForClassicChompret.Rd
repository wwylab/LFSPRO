\name{reformatForClassicChompret}
\alias{reformatForClassicChompret}

\title{
Reformat input data of LFSpro for evaluation using the Classic or the Chompret criteria
}

\description{
Reformat input data of LFSpro for evaluation using the Classic or the Chompret criteria
}

\usage{
reformatForClassicChompret(fam.cancer.data)
}

\arguments{
  \item{fam.cancer.data}{
Combined family and cancer information data for \emph{ONE FAMILY ONLY}. See \link{fam.cancer.data} for details.
}
}

\value{
Data with format used in LFSClassic and LFSChompret. See \code{\link{LFSClassic}} and \code{\link{LFSChompret}} for details.
}

\author{
Gang Peng
}
\note{
One family a time.
}

\seealso{
\code{\link{LFSClassic}}, \code{\link{LFSChompret}}, \code{\link{lfsproC}}
}
\examples{
fam.cancer.data <- CombineData(fam.data, cancer.data)
reformatForClassicChompret(fam.cancer.data[[1]])
}

\keyword{manip}
