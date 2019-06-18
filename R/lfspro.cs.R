.lfspro.cs <-
function (fam.data, cancer.data, penetrance.all, counselee.id, 
          allef, nloci, mRate) 
{
      
  # add cancer type for death as 0
  
  num.cancer <- nrow(cancer.data) # number of cancers observed
  cancer.type.num <- rep(-1, num.cancer)
  colnames(counselee.id) <- c("fam.id", "id")
  # cancer.type.code -> number
  for (i in 1:num.cancer) {
    tmp <- lfspro.cancer.type[cancer.data$cancer.type[i]]
    
    if (is.na(tmp)) {
      print(paste("Cannot find cancer ", cancer.data$cancer.type[i], 
                  " in the LFSpro predefined cancer type", sep = ""))
      print("LFSpro predefined cancer types are: ")
      print(names(new.lfspro.cancer.type))
      print("Please check the input cancer information data.")
      num.counselee <- nrow(counselee.id)
      pp <- rep(-1, num.counselee)
      rlt <- data.frame(cbind(counselee.id, pp), check.names = FALSE)
      colnames(rlt) <- c("fam.id", "id", "pp")
      return(rlt)
    }
    cancer.type.num[i] <- tmp
  }
  cancer.data$cancer.type <- cancer.type.num
  fam.cancer.data <- combinedata(fam.data, cancer.data) # okay!
  
  num.fam <- length(fam.cancer.data)
  pp.all <- NULL
  risk.cs <- NULL
  for (i in 1:num.fam) {
    cid <- counselee.id$id[counselee.id$fam.id == fam.cancer.data[[i]]$fam.id[1]]
    if (length(cid) < 1) {
      print(paste("Cannot find any counselee id in family ", 
                  fam.cancer.data[[i]]$fam.id[1], sep = ""))
      next
    }
    pp <- lfsproC.cs(fam.cancer.data[[i]], penetrance.all, 
                     cid, allef, nloci, mRate)
    pp.tmp <- 1 - pp[,1]
    pp.all <- c(pp.all, pp.tmp)
    risk.temp <- risk.cs(fam.cancer.data[[i]], lfspenet.cs, cid, pp)
    if (is.null(risk.cs)) {
      risk.cs <- c(risk.cs, risk.temp)
    } else {
      risk.cs <- Map(list,risk.cs,risk.temp)
    }
  }
  pp <- data.frame(cbind(counselee.id, pp.all), check.names = F)
  names(pp)[3] <- "pp"
  output <- list(pp, risk.cs)
  return(output)
}
