\name{lfspro}
\alias{lfspro}
\title{
Estimating TP53 mutation probability for families with Li-Fraumeni Syndrome
}
\description{
We use Mendelian risk prediction model to estimate the probability for the counselee as a TP53 mutation carrier on the basis of his/her family cancer history.
}
\usage{
lfspro(fam.data, cancer.data, penetrance.all, counselee.id, allef, nloci, mRate)
}
\arguments{
  \item{fam.data}{
Family information data. See \link{fam.data} for details.
}
  \item{cancer.data}{
Cancer information data. See \link{cancer.data} for details.
}
  \item{penetrance.all}{
Penetrance data. See \link{LFSpenet.2010} for details.
}
  \item{counselee.id}{
Data frame including two variables: fam.id (family id of counselees) and id (individual id of counselees).
}
  \item{allef}{
List. Allele frequency for each locus/gene. If there is only one gene and two alleles in the gene (allele frequency is 0.1 and 0.9), allef = list(c(0.1,0.9)), If there are two genes,two alleles (allele frequency is 0.1 and 0.9) for gene 1  and three alleles (allele frequncy is 0.2, 0.2 and 0.6) for gene 2, allef = list(c(0.1,0.9),c(0.2,0.2,0.6)). We suggeste to use list(c(1-maf,maf)). maf = 0.0001. 
}
  \item{nloci}{
Number of loci/genes in the model. It is set to be 1 here.
}
  \item{mRate}{
Mutation rate. We suggest to use 5e-4.
}
}

\value{
The TP53 mutation carrier probability for each counselee. It is a data frame with 3 variables: fam.id (family id), id (individual id) and pp(posterior probability that the counselee is a TP53 mutation carrier).
}

\references{
Peng, G., Bojadzieva, J., Ballinger, M., Thomas, D., Strong, L. and Wang, W. Estimating TP53 Mutation Carrier Probability in Families with Li-Fraumeni Syndrome Using LFSPRO. \bold{Submitted}.

Chen, S., Wang, W., Broman, K. and Parmigiani, G. (2004) BayesMendel: An R Environment for Mendelian Risk Prediction. \emph{Statistical Application in Genetics and Molecular Biology}, \bold{3(1):} Article 21.
}

\author{
Gang Peng, Wenyi Wang
}

\seealso{
\code{\link{LFSClassic}}, \code{\link{LFSChompret}}, \code{\link{lfsproC}}
}
\examples{
allef.g <- list(c(0.9997,0.0003))
mRate.g <- 6e-05
options(stringsAsFactors = FALSE)
fam.id <- c("fam1","fam2","fam2","fam2","fam2")
id <- c(0,0,2,100,200)
counselee.id <- data.frame(fam.id, id)
lfspro(fam.data, cancer.data, LFSpenet.2010, counselee.id, allef.g, 1,mRate.g)
}

\keyword{model}
