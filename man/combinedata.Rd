\name{combinedata}
\alias{combinedata}
\title{
Combine the family information data and cancer information data
}
\description{
Combine the family information data and cancer information data into a single combined data object. The data object is organized family by family.
}
\usage{
combinedata(fam.data, cancer.data)
}

\arguments{
  \item{fam.data}{
  A data frame storing family information. See \link{fam.data} for details.
}
  \item{cancer.data}{
  A data frame storing cancer information. See \link{cancer.data} for details.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list. Each component in the list stores the family and cancer information for a family. See \link{fam.cancer.data} for details.
}

\author{
Gang Peng
}

\seealso{
\code{\link{lfspro}}, \code{\link{lfsproC}}
}
\examples{
cancer.data$cancer.type <- c(2,4,3,7,2,2,6,50,50,2,4,50,50,4,50) #cancer types
combinedata(fam.data, cancer.data)
}
