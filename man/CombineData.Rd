\name{CombineData}
\alias{CombineData}
\title{
Combine the family information data and cancer information data
}
\description{
Combine the family information data and cancer information data into a single combined data. The data is organized family by family.
}
\usage{
CombineData(fam.data, cancer.data)
}

\arguments{
  \item{fam.data}{
  Data frame storing family information. See \link{fam.data} for details.
}
  \item{cancer.data}{
  Data frame storing cancer information. See \link{cancer.data} for details.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list. Each component in the list stores the family and cancer information for a family. See \link{fam.cancer.data} for details.
}

\author{
Gang Peng
}

\seealso{
\code{\link{lfspro}},\code{\link{lfsproC}}
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

CombineData(fam.data, cancer.data)
}

\keyword{manip}
