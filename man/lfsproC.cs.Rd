\name{lfsproC.cs}
\alias{lfsproC.cs}
\title{
Calculate the posterior probability of p53 mutations for cancer-specific risk prediction model on the basis of family history
}
\description{
lfsproC.cs is used to calculate posterior probability of p53 mutation based on cancer-specific risk prediction model specifically integrated counselees' family cancer history. 
}
\usage{
  lfsproC.cs(data1, data2, counselee.id, penetrance.all, allef, nloci, mRate, mut.info = TRUE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data1}{
    Family data that describe how the individuals are related, for \emph{ONE FAMILY ONLY}.
  }
  \item{data2}{
    A comprehensive dataset that describes each individual in details, including age, gender (1 for male, 2 for female), mutation status (0 for wildtype, 1 for mutation, NA if unknown), vital status (A if alive, D if dead), and cancer history (cancer type and age at diagnosis). Both data1 and data2 can be generated from family and cancer datasets using \link{combinedata} and \link{convert.data}. See the example below.
  }
  \item{counselee.id}{
    Individual id for the counselee. If you want to estimate multiple samples at the same time, just set counselee.id as a vector of IDs for all.
  }
  \item{penetrance.all}{
    Penetrance data. See \link{lfspenet.2010} for details.
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
  \item{mut.info}{
    Set to TRUE (default) if you want to use known genetic testing results of family members in the calculation of mutation probabilities. Set to FALSE otherwise. 
  }
}
\value{
A data frame of posterior probabilities of having wild type, one allele mutated and two allels mutated in TP53 for each counselee.
}
\references{
Shin, S. J., Yuan, Y., Strong, L. C., Bojadzieva, J., and Wang, W. (2018). Bayesian Semiparametric Estimation of Cancer-Specific Age-at-Onset Penetrance With Application to Li-Fraumeni Syndrome. \emph{Journal of the American Statistical Association}, 1–12. https://doi.org/10.1080/01621459.2018.1482749
}
\author{
  Seung Jun Shin, Nam Nguyen, Wenyi Wang
}

\seealso{
  \code{\link{lfsClassic}}, \code{\link{lfsChompret2009}}, \code{\link{lfsChompret2015}}, \code{\link{lfspro}}
}
\examples{
fam.cancer.data <- combinedata(fam.data, cancer.data)
data.obj <- convert.data(fam.cancer.data)
data.obj1 <- data.obj[[1]]
data.obj2 <- data.obj[[2]]

allef <- list(c(0.9997,0.0003))
nloci <-1
mRate <- 6e-05
cid <- c(5,6,9)
lfsproC.cs(data.obj1[[1]], data.obj2[[1]], cid, lfspenet.cs, allef, nloci, mRate)
}
