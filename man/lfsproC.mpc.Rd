\name{lfsproC.mpc}
\alias{lfsproC.mpc}
\title{
Calculate the posterior probability of p53 mutations for multiple primary cancer on the basis of family history
}
\description{
lfsproC.mpc is used to calculate posterior probability of p53 mutation based on prediciting multiple primary cancer model specifically on the baisis of counselees' family cancer history. 
}
\usage{
  lfsproC.mpc(data1, data2, counselee.id, parameter, allef, nloci, mRate, mut.info = TRUE)
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
  \item{parameter}{
    See \link{parameter.mpc} for details.
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
Shin, S. J., Ning, J., Bojadzieva, J., Strong, L. C., and Wang, W. (2018). Bayesian estimation of a semiparametric recurrent event model with applications to the penetrance estimation of multiple primary cancers in Li-Fraumeni syndrome. \emph{Biostatistics}, 00, 1–16. https://doi.org/10.1093/biostatistics/kxy066
}
\author{
  Seung Jun Shin, Nam Nguyen, Wenyi Wang
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
lfsproC.mpc(data.obj1[[1]], data.obj2[[1]], cid, parameter.mpc, allef, nloci, mRate)
}
