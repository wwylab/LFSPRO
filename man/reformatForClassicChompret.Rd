\name{reformatForClassicChompret}
\alias{reformatForClassicChompret}

\title{
Reformat input data of LFSPRO for evaluation using the Classic or the Chompret criteria
}

\description{
Reformat input data of LFSPRO for evaluation using the Classic or the Chompret criteria
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
Data with format used in lfsClassic and lfsChompret. See \code{\link{lfsClassic}}, \code{\link{lfsChompret2015}} and \code{\link{lfsChompret2009}} for details.
}

\author{
Gang Peng
}
\note{
One family a time.
}

\seealso{
\code{\link{lfsClassic}}, \code{\link{lfsChompret2009}}, \code{\link{lfsChompret2015}},  \code{\link{lfsproC}}
}
\examples{
fam.cancer.data <- combinedata(fam.data, cancer.data)
reformatForClassicChompret(fam.cancer.data[[1]])
}

