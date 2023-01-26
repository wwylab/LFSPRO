.lfspro.cs <- function (fam.data, cancer.data, counselee.id, 
                        penetrance.all, allef, nloci, mRate, mut.info = TRUE) {
  num.cancer <- nrow(cancer.data) # number of cancers observed
  colnames(counselee.id) <- c("fam.id", "id")
  
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
  }
  
  fam.cancer.data <- combinedata(fam.data, cancer.data) 
  data.obj <- convert.data(fam.cancer.data)
  data.obj1 <- data.obj[[1]]
  data.obj2 <- data.obj[[2]]
  
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
    pp <- lfsproC.cs(data.obj1[[i]], data.obj2[[i]], cid, 
                     penetrance.all, allef, nloci, mRate, mut.info)
    pp.tmp <- 1 - pp[,1]
    pp.all <- c(pp.all, pp.tmp)
    risk.temp <- risk.cs(fam.cancer.data[[i]], penetrance.all, cid, pp)
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
