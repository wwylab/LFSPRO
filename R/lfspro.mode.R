lfspro.mode <- function(fam.data, cancer.data, counselee.id, mode, penetrance.all=NULL,
                        allef=list(c(0.9994,0.0006)), nloci=1, mRate=0.00012){
  # Aim: calculate the posterior probability of TP53 mutaitons on the basis of
  #    family history
  # Input:
  #    fam.data: data frame. family information including fam.id(family id), id,
  #      fid(father id), mid(mother id), gender(0: female, 1:male), 
  #      vital(vital status, A: alive, D: dead), dob (date of birth)
  #      dod (data of death), lcd (last contact date)
  #    cancer.data: dara frame. cancer information including fam.id(family id),
  #      id, diag.cancer(diagnosed cancer type), diag.date (diagnosed date)
  #    penetrance.all: penetrance table
  #    counselee.id: data frame. the family id and id for the counselees.
  #    allef: alle frequency
  #    nloci: number of loci
  #    mRate: mutation rate
  #    mode: type of LFS prediction model, "1st.all", "mpc" or "1st.cs"
  #    Output: the posterior probability as a TP53 mutation carrier for each counselees
  
  # convert cancer type to specific number and check the cancer type
  # if (mode == "1st.all"| mode == "mpc"){
  #   num.cancer <- nrow(cancer.data)
  #   cancer.type.num <- rep(-1, num.cancer)
  # } else if ((mode == "1st.cs")){
  #   num.cancer <- nrow(cancer.data) 
  #   cancer.type.num <- rep(-1, num.cancer)
  # } else stop(print("Please specify the correct LFS mode"))
  # 
  num.cancer <- nrow(cancer.data) 
  cancer.type.num <- rep(-1, num.cancer)
  colnames(counselee.id) <- c("fam.id", "id")
  for(i in 1:num.cancer){
    tmp <- lfspro.cancer.type[cancer.data$cancer.type[i]]
    if(is.na(tmp)){
      print(paste("Cannot find cancer ", cancer.data$cancer.type[i], 
                  " in the LFSPRO predefined cancer type", sep = ""))
      print("LFSPRO predefined cancer types are: ")
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
  fam.cancer.data <- combinedata(fam.data, cancer.data)
  
  num.fam <- length(fam.cancer.data)
  risk_mpc <- NULL
  risk.cs <- NULL
  pp.all <- NULL
  for(i in 1:num.fam){
    cid <- counselee.id$id[counselee.id$fam.id == fam.cancer.data[[i]]$fam.id[1]]
    if(length(cid)<1){
      print(paste("Cannot find any counselee id in family ", 
                  fam.cancer.data[[i]]$fam.id[1], sep=""))
      next
    }
    if (mode == "1st.all") {
      if (is.null(penetrance.all)){
        penetrance.all = lfspenet.2010
      }
      pp.tmp <- lfsproC(fam.cancer.data[[i]], penetrance.all, 
                        cid, allef, nloci,mRate)
      pp.all <- c(pp.all, pp.tmp)
      output <- data.frame(cbind(counselee.id,pp.all), check.names = FALSE)
      colnames(output) <- c("fam.id", "id", "pp")
    } else if (mode == "mpc") {
      if (is.null(penetrance.all)){
        penetrance.all = parameter.mpc
      }
      data.obj <- convert.data(fam.cancer.data)
      data.obj1 <- data.obj[[1]]
      data.obj2 <- data.obj[[2]]
      pp.tmp <- lfsproC.mpc(fam.cancer.data[[i]], penetrance.all, data.obj1[[i]],
                            data.obj2[[i]], cid, allef, nloci, mRate) 
      risk <- risk.mpc(fam.cancer.data[[i]],cancer.data, cid, data.obj2[[i]], penetrance.all)
      risk_mpc <- rbind(risk_mpc, risk)
      pp.all <- rbind(pp.all, pp.tmp)
    } else if (mode == "1st.cs") {
      if (is.null(penetrance.all)){
        penetrance.all = lfspenet.cs.nodeath
      }
      pp <- lfsproC.cs(fam.cancer.data[[i]], penetrance.all, 
                       cid, allef, nloci, mRate)
      pp.tmp <- 1 - pp[,1]
      pp.all <- c(pp.all, pp.tmp)
      risk.temp <- risk.cs(fam.cancer.data[[i]], lfspenet.cs.death, cid, pp)
      if (is.null(risk.cs)) {
        risk.cs <- c(risk.cs, risk.temp)
      } else {
        risk.cs <- Map(list,risk.cs,risk.temp)
      }
    } else stop(print("Please specify the correct LFS mode"))
  }
  
  if (mode == "mpc"){
    pp <- 1 - pp.all[, 1]
    rlt <- data.frame(cbind(counselee.id, pp), check.names = FALSE)
    risk_mpc <- data.frame(risk_mpc)
    colnames(rlt) <- c("fam.id", "id", "pp")
    colnames(risk_mpc) <- c("fam.id", "ID", "age","5 years(wildtype)", "5 years(mutation)", 
                            "10 years(wildtype)", "10 years(mutation)", "15 years(wildtype)", 
                            "15 years(mutation)")
    risk.all <- combined.risk.mpc(pp.all, risk_mpc, counselee.id)
    output <- list(rlt, risk.all)
  } else if (mode == "1st.cs"){
    pp <- data.frame(cbind(counselee.id, pp.all), check.names = FALSE)
    names(pp)[3] <- "pp"
    output <- list(pp, risk.cs)
  }
  return(output)
}
