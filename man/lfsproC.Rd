\name{lfsproC}
\alias{lfsproC}
\title{
Calculate the posterior probability of p53 mutations on the basis of family history
}
\description{
lfsproC is used to calculate posterior probability of carrying p53 mutation on the baisis of counselees' family cancer history. 
}
\usage{
lfsproC(fam.cancer.data, penetrance.all, counselee.id, allef, nloci, mRate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fam.cancer.data}{
Combined family and cancer information data for \emph{ONE FAMILY ONLY}. See \link{fam.cancer.data} for details.
}
  \item{penetrance.all}{
Penetrance data. See \link{lfspenet.2010} for details.
}
  \item{counselee.id}{
Individual id for the counselee. If you want to estimate multiple samples at the same time, just set counselee.id as a vector of IDs for all.
}
  \item{allef}{
List. Allele frequency for each locus/gene. If there is only one gene and two alleles in the gene (allele frequency is 0.1 and 0.9), allef = list(c(0.1,0.9)), If there are two genes,two alleles (allele frequency is 0.1 and 0.9) for gene 1  and three alleles (allele frequncy is 0.2, 0.2 and 0.6) for gene 2, allef = list(c(0.1,0.9),c(0.2,0.2,0.6)). We set it as allef=list(c(0.9994,0.0006)) by default. 
}
  \item{nloci}{
Number of loci/genes in the model. It's set to be 1 here.
}
  \item{mRate}{
Mutation rate. We set it as 0.00012 by default.
}
}
\value{
A data frame of posterior probabilities of having wild type, one allele mutated and two allels mutated in TP53 for each counselee.
}
\references{
Peng, G., Bojadzieva, J., Ballinger, M. L., Li, J., Blackford, A. L., Mai, P. L., Savage, S. A., Thomas, D. M., Strong, L. C., … Wang, W. (2017). Estimating TP53 Mutation Carrier Probability in Families with Li-Fraumeni Syndrome Using LFSPRO. \emph{Cancer Epidemiology, Biomarker and Prevention}, \bold{26(6)}, 837-844. DOI: 10.1158/1055-9965.EPI-16-0695

Chen, S., Wang, W., Broman, K. and Parmigiani, G. (2004) BayesMendel: An R Environment for Mendelian Risk Prediction. \emph{Statistical Application in Genetics and Molecular Biology}, \bold{3(1):} Article 21. DOI: 10.2202/1544-6115.1063

}
\author{
Gang Peng, Wenyi Wang
}

\seealso{
\code{\link{lfsClassic}}, \code{\link{lfsChompret2009}}, \code{\link{lfsChompret2015}}, \code{\link{lfspro.mode}}, \code{\link{peelingRC}}
}
\examples{
# convert cancer type to specific number and check the cancer type
cancer.data$cancer.type <- c(2,4,3,7,2,2,6,50,50,2,4,50,50,4,50) #cancer types
allef <- allef.g <- list(c(0.9997,0.0003))
mRate.g <- 6e-05
fam.cancer.data <- combinedata(fam.data, cancer.data)
lfsproC(fam.cancer.data[[1]], lfspenet.2010, 0, allef.g, 1,mRate.g)
}

