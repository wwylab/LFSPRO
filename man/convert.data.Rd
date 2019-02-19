\name{convert.data}
\alias{convert.data}
\title{
  Convert the combined family information data and cancer information data into separate data objects. 
}
\description{
  Separate the family information data and cancer information data into two data objects used in multiple primary cancer predictions. 
}
\usage{
  convert.data(fam.cancer.data)
}

\arguments{
  \item{fam.cancer.data}{
    Data frame storing family (ID, sex, father ID and mother ID) and cancer information (ID, time sex, test, death status (Yes/No) and counts of multiple primary cancers). 
  }
}
\value{
  A list. Each component in the list stores the family and cancer information for a family. See \link{fam.cancer.data} for details.
}

\author{
  Seung Jun Shin
}
\examples{
  fam.cancer.data <- combinedata(fam.data, cancer.data)
  data.obj <- convert.data(fam.cancer.data)
}







