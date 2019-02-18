\name{lfsproC}
\alias{lfsproC}
\title{
Estimating TP53 mutation probability for families with Li-Fraumeni Syndrome for combined family and cancer data
}
\description{
We use Mendelian risk prediction model to estimate the probability for the counselee as a TP53 mutation carrier on the basis of his/her family cancer history.
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
Penetrance data. See \link{LFSpenet.2010} for details.
}
  \item{counselee.id}{
Individual id for the counselee. If you want to estimate multiple samples at the same time, just set counselee.id as a vector of IDs for all the samples you want to estimate.
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
The TP53 mutation carrier probability for each counselee.
}
\references{
Peng, G., Bojadzieva, J., Ballinger, M., Thomas, D., Strong, L. and Wang, W. Estimating TP53 Mutation Carrier Probability in Families with Li-Fraumeni Syndrome Using LFSPRO. \bold{Submitted}.

Chen, S., Wang, W., Broman, K. and Parmigiani, G. (2004) BayesMendel: An R Environment for Mendelian Risk Prediction. \emph{Statistical Application in Genetics and Molecular Biology}, \bold{3(1):} Article 21.
}
\author{
Gang Peng, Wenyi Wang
}

\seealso{
\code{\link{LFSClassic}}, \code{\link{LFSChompret}}, \code{\link{lfspro}}, \code{\link{peelingRC}}
}
\examples{
# convert cancer type to specific number and check the cancer type
num.cancer <- nrow(cancer.data)
cancer.type.num <- rep(-1, num.cancer)
for(i in 1:num.cancer){
  tmp <- LFSpro.cancer.type[cancer.data$cancer.type[i]]
  if(is.na(tmp)){
    print(paste("Cannot find cancer ", cancer.data$cancer.type[i], 
                  " in the LFSpro predefined cancer type", sep = ""))
    print("LFSpro predefined cancer types are: ")
    print(cancer.type.all)
    print("Please check the input cancer information data.")
      
    num.counselee <- nrow(counselee.id)
    pp <- rep(-1, num.counselee)
      
    rlt <- data.frame(cbind(counselee.id, pp),check.names = FALSE)
    colnames(rlt) <- c("fam.id", "id", "pp")
    return(rlt)
  }
  cancer.type.num[i] <- tmp
}
  
cancer.data$cancer.type <- cancer.type.num
  
allef <- allef.g <- list(c(0.9997,0.0003))
mRate.g <- 6e-05
fam.cancer.data <- CombineData(fam.data, cancer.data)
lfsproC(fam.cancer.data[[1]], LFSpenet.2010, 0, allef.g, 1,mRate.g)
}

\keyword{model}
