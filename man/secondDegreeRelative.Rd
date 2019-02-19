\name{secondDegreeRelative}
\alias{secondDegreeRelative}

\title{
Second degree relative
}
\description{
Find out the second degree relatives of a given individuals in a family.
}
\usage{
secondDegreeRelative(pedigree, ii)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pedigree}{
A data frame including ID, fID and mID. These IDs should be positive integers.
}
  \item{ii}{
Individual \emph{index} (not ID here). If you want to get the second degree relatives of the second individual in the family, set ii as 2.
}
}

\value{
The \emph{indexes} of individuals who are the second degree relatives to the given individual.
}

\author{
Gang Peng
}

\seealso{
\code{\link{firstDegreeRelative}}, \code{\link{lfsClassic}}, \code{\link{lfsChompret2009}}, \code{\link{lfsChompret2015}}
}
\examples{
fam.cancer.data <- combinedata(fam.data, cancer.data)
lfsData <- reformatForClassicChompret(fam.cancer.data[[1]])
pedigree <- data.frame(lfsData$ID,lfsData$fID,lfsData$mID,stringsAsFactors=FALSE)
names(pedigree) <- c("ID","fID","mID")
secondDegreeRelative(pedigree,2)
}

