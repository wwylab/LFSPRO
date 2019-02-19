\name{lfsproC.mpc}
\alias{lfsproC.mpc}
\title{
Calculate the posterior probability of p53 mutations for multiple primary cancer on the basis of family history
}
\description{
lfsproC.mpc is used to calculate posterior probability of p53 mutation based on prediciting multiple primary cancer model specifically on the baisis of counselees' family cancer history. 
}
\usage{
  lfsproC.mpc(fam.cancer.data, parameter, data1, data2, counselee.id, allef, nloci, mRate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fam.cancer.data}{
    Combined family and cancer information data for \emph{ONE FAMILY ONLY}. See \link{fam.cancer.data} for details.
  }
  \item{parameter}{
    See \link{parameter.mpc} for details.
  }
  \item{data1}{sex, father ID and mother ID}
  \item{data2}{age of diagnosis, sex and multiple primary cancers index}
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
Shin, S. J., Ning, J., Bojadzieva, J., Strong, L. C., and Wang, W. (2018). Bayesian estimation of a semiparametric recurrent event model with applications to the penetrance estimation of multiple primary cancers in Li-Fraumeni syndrome. \emph{Biostatistics}, 00, 1–16. https://doi.org/10.1093/biostatistics/kxy066
}
\author{
  Seung Jun Shin, Wenyi Wang
}

\examples{

fam.cancer.data <- combinedata(fam.data, cancer.data)
data.obj <- convert.data(fam.cancer.data)
data.obj1 <- data.obj[[1]]
data.obj2 <- data.obj[[2]]

allef.g <- list(c(0.9997,0.0003))
nloci.g <-1
mRate.g <- 6e-05
cid <- c(0,2,100,200)
lfsproC.mpc(fam.cancer.data[[1]], parameter.mpc, 
                        data.obj1[[1]], data.obj2[[1]], cid, allef.g, nloci.g, mRate.g) 
}

