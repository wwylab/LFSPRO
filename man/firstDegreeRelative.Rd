\name{firstDegreeRelative}
\alias{firstDegreeRelative}
\title{
First degree relatives
}
\description{
Find out the first degree relatives of a given individuals in a family.
}
\usage{
firstDegreeRelative(pedigree, ii)
}

\arguments{
  \item{pedigree}{
A data frame including ID, fID and mID. These IDs should be a positive integer.
}
  \item{ii}{
Individual \emph{index} (not ID here). If you want to get the first degree relatives of the second individual in the family, set ii as 2.
}
}

\value{
The \emph{index} of individuals who are the first degree relatives to the given individual. 
}

\author{
Gang Peng
}

\seealso{
\code{\link{secondDegreeRelative}}, \code{\link{lfsClassic}}, \code{\link{lfsChompret2015}}
}
\examples{
fam.cancer.data <- combinedata(fam.data, cancer.data)
lfsData <- reformatForClassicChompret(fam.cancer.data[[1]])
pedigree <- data.frame(lfsData$ID,lfsData$fID,lfsData$mID,stringsAsFactors=FALSE)
names(pedigree) <- c("ID","fID","mID")
firstDegreeRelative(pedigree,2)
}

