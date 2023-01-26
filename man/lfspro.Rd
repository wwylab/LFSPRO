\name{lfspro}
\alias{lfspro}
\title{
Estimate TP53 mutation probability and predict cancer risk for families with Li-Fraumeni Syndrome
}
\description{
We use a Mendelian risk prediction model to estimate the probability of the counselee being a TP53 germline mutation carrier on the basis of his/her family cancer history. lfspro estimates TP53 carrier probability for every counselee but predicts cancer risks based on cancer history of counselee. If the counselee never had primary cancer, lfspro estimates cancer-specific risks. If the counselee had one primary cancer, lfspro estimates multiple primary cancer risks. If the counselee had multiple primary cancer, lfspro cannot predict cancer risks. 
}
\usage{
  lfspro(fam.data, cancer.data, counselee.id, method, allef, nloci, mRate, mut.info)
}
\arguments{
  \item{fam.data}{
    Family information data. See \link{fam.data} for details.
  }
  \item{cancer.data}{
    Cancer information data. See \link{cancer.data} for details.
  }
    \item{counselee.id}{
    Data frame including two variables: fam.id (family id of counselees) and id (individual id of counselees).
  }
  \item{method}{
    The model used to calculate TP53 mutation probabilities. The default is "MPC", but can be switched to "CS".
  } 
  \item{allef}{
    List. Allele frequency for each locus/gene. If there is only one gene and two alleles in the gene (allele frequency is 0.1 and 0.9
    ), allef = list(c(0.1,0.9)), If there are two genes,two alleles (allele frequency is 0.1 and 0.9) for gene 1  and three alleles 
    (allele frequncy is 0.2, 0.2 and 0.6) for gene 2, allef = list(c(0.1,0.9),c(0.2,0.2,0.6)). We suggeste to use list(c(1-maf,maf)). We set it as allef = list(c(0.9994,0.0006)) by default. 
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
The output consists of 3 lists, TP53 carrier probability, cancer-specific risk and multiple primary cancer risk prediction.
TP53 carrier probability is a data frame with 3 variables: fam.id (family id), id (individual id) and pp (posterior probability that the counselee is a TP53 mutation carrier).
Cancer-specific risk prediction is a data frame of future risks of having breast cancer, sarcoma and other cancers for each counselee.
Multiple primary cancer risk prediction is a data frame of future cancer risk prediction of developing multiple primary cancer for each counselee. 
}

\references{
Peng, G., Bojadzieva, J., Ballinger, M. L., Li, J., Blackford, A. L., Mai, P. L., Savage, S. A., Thomas, D. M., Strong, L. C., … Wang, W. (2017). Estimating TP53 Mutation Carrier Probability in Families with Li-Fraumeni Syndrome Using LFSPRO. \emph{Cancer Epidemiology, Biomarker and Prevention}, \bold{26(6)}, 837-844. DOI: 10.1158/1055-9965.EPI-16-0695

Shin, S. J., Ning, J., Bojadzieva, J., Strong, L. C., and Wang, W. (2018). Bayesian estimation of a semiparametric recurrent event model with applications to the penetrance estimation of multiple primary cancers in Li-Fraumeni syndrome. \emph{Biostatistics}, 00, 1–16. https://doi.org/10.1093/biostatistics/kxy066

Shin, S. J., Yuan, Y., Strong, L. C., Bojadzieva, J., and Wang, W. (2018). Bayesian Semiparametric Estimation of Cancer-Specific Age-at-Onset Penetrance With Application to Li-Fraumeni Syndrome. \emph{Journal of the American Statistical Association}, 1–12. https://doi.org/10.1080/01621459.2018.1482749

Chen, S., Wang, W., Broman, K. and Parmigiani, G. (2004) BayesMendel: An R Environment for Mendelian Risk Prediction. \emph{Statistical Application in Genetics and Molecular Biology}, \bold{3(1):} Article 21. DOI: 10.2202/1544-6115.1063

}

\author{
  Gang Peng, Seung Jun Shin, Jingxiao Chen, Nam Nguyen, Wenyi Wang
}


\examples{
fam.id <- c("fam1","fam2","fam2","fam2","fam2")
id <- c(9,0,2,100,200)
counselee.id <- data.frame(fam.id, id)
lfspro(fam.data, cancer.data, counselee.id)
}
