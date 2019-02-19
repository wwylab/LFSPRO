\name{fam.data}
\alias{fam.data}
\docType{data}
\title{
Built-in Family Information Data Example
}
\description{
This built-in family information dataset is used to demonstrate the examples of the pacakge. Input family data should maintain the same format, column names and data types of the data example to avoid warning or error messages. 
}
\usage{fam.data}
\format{
  A data frame containing 9 variables.
  \describe{
    \item{\code{fam.id}}{family id}
    \item{\code{id}}{individual id}
    \item{\code{fid}}{father id}
    \item{\code{mid}}{mother id}
    \item{\code{gender}}{sex. 0: female, 1: male}
    \item{\code{age}}{If the individual was alive at the last contact date, use the age at that time. otherwise, use the age at death. Note: do not record age at diagnosis, which shoule be recorded in cancer.data}
    \item{\code{vital}}{Vital sign. D: death, A: alive}
  }
}
\details{
Family id could include characters, but id, father id and mother id must be intergers larger than 0.
}
\seealso{
\code{\link{cancer.data}}
}
\examples{
fam.data
}
