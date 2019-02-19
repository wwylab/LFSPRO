\name{LFSPRO-package}
\alias{LFSPRO-package}
\docType{package}
\title{
TP53 mutation carrier estimation and cancer risk predictions
}
\description{
TP53 germline mutations are the main cause of Li-Fraumeni syndrome. This package is designed to estimate probabilities that:  1) the counselee is a TP53 germline mutation carrier, 2) the counselee develops any cancer in future, 3) the counselee develops breast cancer, sarcoma or any other cancers in future, 4) the counselee develops a first or second primary cancer in future, on the basis of his/her family cancer history. The package also provides functions for using the LFS classic and Chompret criteria.
}
\details{
\tabular{ll}{
Package: \tab LFSPRO\cr
Type: \tab Package\cr
Version: \tab 2.0.0\cr
Date: \tab 2018-12-04\cr
License: \tab GPL-3 \cr
}
}
\author{
Gang Peng, SeungJun Shin, Jasper Chen, Wenyi Wang

Maintainer: Wenyi Wang <wwang7@mdanderson.org>
}
\references{
Peng, G., Bojadzieva, J., Ballinger, M. L., Li, J., Blackford, A. L., Mai, P. L., Savage, S. A., Thomas, D. M., Strong, L. C., … Wang, W. (2017). Estimating TP53 Mutation Carrier Probability in Families with Li-Fraumeni Syndrome Using LFSPRO. \emph{Cancer Epidemiology, Biomarker and Prevention}, \bold{26(6)}, 837-844. DOI: 10.1158/1055-9965.EPI-16-0695

Shin, S. J., Ning, J., Bojadzieva, J., Strong, L. C., and Wang, W. (2018). Bayesian estimation of a semiparametric recurrent event model with applications to the penetrance estimation of multiple primary cancers in Li-Fraumeni syndrome. \emph{Biostatistics}, 00, 1–16. https://doi.org/10.1093/biostatistics/kxy066

Shin, S. J., Yuan, Y., Strong, L. C., Bojadzieva, J., and Wang, W. (2018). Bayesian Semiparametric Estimation of Cancer-Specific Age-at-Onset Penetrance With Application to Li-Fraumeni Syndrome. \emph{Journal of the American Statistical Association}, 1–12. https://doi.org/10.1080/01621459.2018.1482749

Chen, S., Wang, W., Broman, K. and Parmigiani, G. (2004) BayesMendel: An R Environment for Mendelian Risk Prediction. \emph{Statistical Application in Genetics and Molecular Biology}, \bold{3(1):} Article 21. DOI: 10.2202/1544-6115.1063

}

\note{
New features in 2.0.0: lfspro.mode is added to replace lfspro to account for additional utilities in cancer-specific and multiple primary cancer risk predictions. Simply specify parameter "mode" to call desired functions. 

New features in 1.0.5: (1) add function lfsChompret2015, (2) update the the default MAF and de novo mutation rate. 
}

\seealso{
\code{\link{lfspro.mode}}, \code{\link{lfsClassic}}, \code{\link{lfsChompret2009}}
}
\examples{
fam.id <- c("fam1","fam2","fam2","fam2","fam2")
id <- c(1,1,2,100,200)
counselee.id <- data.frame(fam.id, id)

# LFS classic criteria
lfsClassic(fam.data, cancer.data, counselee.id)
# Chompret criteria
lfsChompret2015(fam.data, cancer.data, counselee.id)

# "1st.all" predict the probability of carrying TP53 mutations
lfspro.mode(fam.data, cancer.data, counselee.id, "1st.all")
# "mpc" predict future risks of developing multiple primary cancers
lfspro.mode(fam.data, cancer.data, counselee.id, "mpc")
# "1st.cs" predict future risks of having breast cancer, sarcoma, other cancers, and death
lfspro.mode(fam.data, cancer.data, counselee.id, "1st.cs")
}
