.lfspro.mpc <- function (fam.data, cancer.data, counselee.id, 
                         parameter, allef, nloci, mRate, mut.info = TRUE) {
  num.cancer <- nrow(cancer.data)
  colnames(counselee.id) <- c("fam.id", "id")
  
  for (i in 1:num.cancer) {
    tmp <- lfspro.cancer.type[cancer.data$cancer.type[i]]
    if (is.na(tmp)) {
      print(paste("Cannot find cancer ", cancer.data$cancer.type[i], 
                  " in the LFSPRO predefined cancer type", sep = ""))
      print("LFSPRO predefined cancer types are: ")
      print(cancer.type.all)
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
  risk_mpc <- NULL
  
  for (i in 1:num.fam) {
    cid <- counselee.id$id[counselee.id$fam.id == fam.cancer.data[[i]]$fam.id[1]]
    if (length(cid) < 1) {
      print(paste("Cannot find any counselee id in family ", 
                  fam.cancer.data[[i]]$fam.id[1], sep = ""))
      next
    }
    pp.tmp <- lfsproC.mpc(data.obj1[[i]], data.obj2[[i]], cid, 
                          parameter, allef, nloci, mRate, mut.info) 
    risk <- risk.mpc(fam.cancer.data[[i]], cancer.data, cid, data.obj2[[i]], parameter)
    risk_mpc <- rbind(risk_mpc, risk)
    pp.all <- rbind(pp.all, pp.tmp)
  }
  
  pp <- 1 - pp.all[, 1]
  rlt <- data.frame(cbind(counselee.id, pp), check.names = FALSE)
  risk_mpc <- data.frame(risk_mpc)
  colnames(rlt) <- c("fam.id", "id", "pp")
  colnames(risk_mpc) <- c("fam.id", "ID", "age","5 years(wildtype)", "5 years(mutation)", 
                     "10 years(wildtype)", "10 years(mutation)", "15 years(wildtype)", 
                     "15 years(mutation)")
  risk.all <- combined.risk.mpc(pp.all, risk_mpc, counselee.id)
  output <- list(rlt, risk.all)
  return(output)
}
