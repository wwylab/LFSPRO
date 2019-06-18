lfspro <- function(fam.data, cancer.data, counselee.id, penetrance.all=NULL,
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
  risk.mpc.output <- NULL
  risk.cs.output <- NULL
  risk.mpc.final <- data.frame()
  pp.all <- NULL
  for(i in 1:num.fam){
    cid <- counselee.id$id[counselee.id$fam.id == fam.cancer.data[[i]]$fam.id[1]]
    if(length(cid)<1){
      print(paste("Cannot find any counselee id in family ", 
                  fam.cancer.data[[i]]$fam.id[1], sep=""))
      next
    }
    
    ## Carrier probability calculation with MPC
    if (is.null(penetrance.all)){
      penetrance.all = parameter.mpc
    }
    data.obj <- convert.data(fam.cancer.data)
    data.obj1 <- data.obj[[1]]
    data.obj2 <- data.obj[[2]]
    pp.tmp <- lfsproC.mpc(fam.cancer.data[[i]], penetrance.all, data.obj1[[i]],
                          data.obj2[[i]], cid, allef, nloci, mRate) 
    pp.all <- rbind(pp.all, pp.tmp)

    ## risk prediction
    cid_num.cancer <- fam.cancer.data[[i]]$num.cancer[which(fam.cancer.data[[i]]$id %in% cid)] 
    cid.na <- cid[which(cid_num.cancer==0)] #counselee without previous cancers
    cid.1 <- cid[which(cid_num.cancer>=1)] #counselee with previous primary cancer
    
    pp.na <- pp.tmp[which(cid_num.cancer==0),]
    dim(pp.na) <- c(sum(cid_num.cancer==0), 3)
    pp.1 <- pp.tmp[which(cid_num.cancer>=1),]
    dim(pp.1) <- c(sum(cid_num.cancer>=1), 3)
    
    if (length(cid.na)>0){
      risk.cs.temp <- risk.cs(fam.cancer.data[[i]], lfspenet.cs, cid.na, pp.na)
    } else {
      risk.cs.temp <- NULL
    }
    if (is.null(risk.cs.output)) {
      risk.cs.output <- c(risk.cs.output, risk.cs.temp)
    } else {
      risk.cs.output <- Map(list,risk.cs.output,risk.cs.temp)
    }
    
    if (length(cid.1)>0){
      risk.mpc.temp <- risk.mpc(fam.cancer.data[[i]],cancer.data, cid.1, data.obj2[[i]], penetrance.all)
      risk.mpc.output <- data.frame(risk.mpc.temp)
      colnames(risk.mpc.output) <- c("fam.id", "ID", "age","5 years(wildtype)", "5 years(mutation)", 
                              "10 years(wildtype)", "10 years(mutation)", "15 years(wildtype)", 
                              "15 years(mutation)")
      counselee.id[which(cid_num.cancer>=1),]
      counselee.id.1 <- data.frame(fam.id=fam.cancer.data[[i]]$fam.id[1],id=cid.1)
      risk.all <- combined.risk.mpc(pp.1, risk.mpc.output, counselee.id.1)
      risk.mpc.final <- rbind(risk.mpc.final, risk.all)
    }
  }
  pp <- 1 - pp.all[, 1]
  rlt <- data.frame(cbind(counselee.id, pp), check.names = FALSE)
  colnames(rlt) <- c("fam.id", "id", "mutation_probability")
  output <- list(rlt, risk.cs.output, na.omit(risk.mpc.final))
  names(output) <- c("Mutation_probability", "Cancer_specific_risks",
                     "Multiple_primary_cancer_risks")
  return(output)
}
