\name{LFSPRO-package}
\alias{LFSPRO-package}
\docType{package}
\title{
TP53 mutation carrier estimation
}
\description{
TP53 mutation is a main cause of Li-Fraumeni syndrome. This package is designed to estimate the probability that the counselee is a TP53 mutation carrier on the basis of his/her family cancer history. We also included LFS classic and Chompret criteria in the package.
}
\details{
\tabular{ll}{
Package: \tab LFSPRO\cr
Type: \tab Package\cr
Version: \tab 1.0.5\cr
Date: \tab 2016-12-14\cr
License: \tab GPL-3 \cr
}
}
\author{
Gang Peng, Wenyi Wang

Maintainer: Gang Peng <michael.gang.peng@gmail.com>
}
\references{
Peng, G., Bojadzieva, J., Ballinger, M., Thomas, D., Strong, L. and Wang, W. Estimating TP53 Mutation Carrier Probability in Families with Li-Fraumeni Syndrome Using LFSpro. \bold{Submitted}.

Chen, S., Wang, W., Broman, K. and Parmigiani, G. (2004) BayesMendel: An R Environment for Mendelian Risk Prediction. \emph{Statistical Application in Genetics and Molecular Biology}, \bold{3(1):} Article 21.
}

\note{
New features in 1.0.5: (1) add function LFSChompret2015, (2) update the the default MAF and de novo mutation rate. 
}

\keyword{ package }
\seealso{
\code{\link{lfspro}}, \code{\link{LFSClassic}}, \code{\link{LFSChompret}}
}
\examples{
options(stringsAsFactors = FALSE)
fam.id <- c("fam1","fam2","fam2","fam2","fam2")
id <- c(0,0,2,100,200)
counselee.id <- data.frame(fam.id, id)
LFSClassic(fam.data, cancer.data, counselee.id)
LFSChompret(fam.data, cancer.data, counselee.id)

allef.g <- list(c(0.9997,0.0003))
mRate.g <- 6e-05

lfspro(fam.data, cancer.data, LFSpenet.2010, counselee.id, allef.g, 1,mRate.g)
}
