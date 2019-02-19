\name{fam.cancer.data}
\alias{fam.cancer.data}
\docType{data}
\title{
Data object with Family and Cancer Information Example
}
\description{
Data object example with integrated family and cancer information. After combining family and cancer information with \link{combinedata}, the data object should be similar to this example. 
}
\usage{fam.cancer.data}
\format{
List. Each component of the list stores the family and cancer information of one family. Each component is also a list with the following variables:
  \describe{
    \item{\code{fam.id}}{family id}
    \item{\code{id}}{individual id}
    \item{\code{fid}}{father id}
    \item{\code{mid}}{mother id}
    \item{\code{sex}}{sex, 0: female, 1:male}
    \item{\code{age}}{If the individual is dead, it is age at death, otherwise it is the age at last contact.}
    \item{\code{cancer.info}}{Cancer informatice. List. Each component of the list stores cancer information for each individual in the family. For each component, if there is no cancer for the individual, it is a 0 by 0 data frame. Otherwise, it is a data frame with two columns, cancer.type (coded number) and diag.age (the age at which the cancer was diagnosed).}
  }
}
\examples{
fam.cancer.data
}
