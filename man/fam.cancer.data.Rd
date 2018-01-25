\name{fam.cancer.data}
\alias{fam.cancer.data}
\docType{data}
\title{
Combined data with family and cancer information
}
\description{
Combined data with family and cancer information. It is used in function LFSproC.
}
\usage{fam.cancer.data}
\format{
List. Each component of the list stores the family and cancer information of one family. Each component is also a list with the following variables:
  \describe{
    \item{\code{fam.id}}{family id}
    \item{\code{id}}{individual id}
    \item{\code{fid}}{father id}
    \item{\code{mid}}{mother id}
    \item{\code{gender}}{gender, 0: female, 1:male}
    \item{\code{age}}{If the individual is dead, it is death age, otherwise it is current age. We use last contact date to calculate current age.}
    \item{\code{cancer.info}}{Cancer informatice. List. Each component of the list stores cancer information for each individual in the family. For each component, if there there is no cancer for the individual, it is a 0 by 0 data frame. Otherwise, it is a data frame with cancer.type (coded number) and diag.age(the age at which the cancer was diagnosed).}
  }
}
\examples{
fam.cancer.data
}
\keyword{datasets}
