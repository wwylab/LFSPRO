\name{cancer.data}
\alias{cancer.data}
\docType{data}
\title{
Built-in Cancer Information Data example
}
\description{
This built-in cancer information dataset is used to demonstrate the examples of the pacakge. Input cancer data should maintain the same format, column names and data types of the data example to avoid warning or error messages. 
}
\usage{cancer.data}
\format{
  A data frame with 4 variables.
  \describe{
    \item{\code{fam.id}}{Family id}
    \item{\code{id}}{Individual id}
    \item{\code{cancer.type}}{Cancer type. See \link{lfspro.cancer.type} for details}
    \item{\code{diag.age}}{The age when the individual was diagnosed with the corresponding cancer}
  }
}
\seealso{
\code{\link{lfspro.cancer.type}}, \code{\link{cancer.type.all}}, \code{\link{lfs.cut}}, \code{\link{invasive.cut}}
}
\examples{
cancer.data
}
