\name{lkNoneAffect}
\alias{lkNoneAffect}
\title{
Likelihood for non-affected individuals
}
\description{
Calculate the likelihood for non-affected individuals.
}
\usage{
lkNoneAffect(penetrance, age)
}
\arguments{
  \item{penetrance}{
The penetrance matrix for only male or female. See \link{LFSpenet.2010} for details.
}
  \item{age}{
The age when the individual is still healthy.
}
}

\value{
The likelihood (Pr(D|G)) for the none affected individual at 'age' years old. G can be TP530, TP531 and PT532.
}

\author{
Gang Peng
}

\seealso{
\code{\link{calLK}}
}
\examples{
lkNoneAffect(LFSpenet.2010$fFX, 50)
}

\keyword{model}
