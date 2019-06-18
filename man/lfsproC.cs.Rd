\name{lfsproC.cs}
\alias{lfsproC.cs}
\title{
Calculate the posterior probability of p53 mutations for cancer-specific risk prediction model on the basis of family history
}
\description{
lfsproC.cs is used to calculate posterior probability of p53 mutation based on cancer-specific risk prediction model specifically integrated counselees' family cancer history. 
}
\usage{
  lfsproC.cs(fam.cancer.data, penetrance.all, counselee.id, allef, nloci, mRate) 
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
Shin, S. J., Yuan, Y., Strong, L. C., Bojadzieva, J., and Wang, W. (2018). Bayesian Semiparametric Estimation of Cancer-Specific Age-at-Onset Penetrance With Application to Li-Fraumeni Syndrome. \emph{Journal of the American Statistical Association}, 1–12. https://doi.org/10.1080/01621459.2018.1482749
}
\author{
  Seung Jun Shin, Wenyi Wang
}

\seealso{
  \code{\link{lfsClassic}}, \code{\link{lfsChompret2009}}, \code{\link{lfsChompret2015}}, \code{\link{lfspro}}
}
\examples{
famdata <- fam.data[fam.data$fam.id == "fam2",] 
cancerdata <- cancer.data[cancer.data$fam.id == "fam2",] 
cancer.data$cancer.type <- c(2,4,3,7,2,2,6,50,50,2,4,50,50,4,50) #cancer types
nloci <- 1 
allef <- list(c(1 - 1.0e-4, 1.0e-4))
mRate <- 1.0e-6
cid <- c(1:4)
fam.cancer.data <- combinedata(fam.data, cancer.data) 
lfsproC.cs(fam.cancer.data[[1]], lfspenet.cs, 
           cid, allef, nloci, mRate)
}
