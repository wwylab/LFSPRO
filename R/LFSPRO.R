lfspro <- function(fam.data, cancer.data, penetrance.all, counselee.id,
                   allef, nloci, mRate){
  # Aim: calculate the posterior probability of TP53 mutaitons on the basis of
  #    family history
  # Author: Gang Peng
  #
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
  # 
  # Output: the posterior probability as a TP53 mutation carrier for each counselees
  
  # convert cancer type to specific number and check the cancer type
  num.cancer <- nrow(cancer.data)
  cancer.type.num <- rep(-1, num.cancer)
  for(i in 1:num.cancer){
    tmp <- LFSpro.cancer.type[cancer.data$cancer.type[i]]
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
  
  
  fam.cancer.data <- CombineData(fam.data, cancer.data)
  num.fam <- length(fam.cancer.data)
  
  pp.all <- NULL
  for(i in 1:num.fam){
    cid <- counselee.id$id[counselee.id$fam.id == fam.cancer.data[[i]]$fam.id[1]]
    if(length(cid)<1){
      print(paste("Cannot find any counselee id in family ", 
                  fam.cancer.data[[i]]$fam.id[1], sep=""))
      next
    }
    
    pp.tmp <- lfsproC(fam.cancer.data[[i]], penetrance.all, cid, allef, nloci,mRate)
    pp.all <- c(pp.all, pp.tmp)
  }

  rlt <- data.frame(cbind(counselee.id,pp.all),check.names = FALSE)
  colnames(rlt) <- c("fam.id", "id", "pp")
  return(rlt)
}

lfsproC <- function(fam.cancer.data, penetrance.all, counselee.id, allef, nloci,mRate){
  # Aim: calculate the posterior probability of p53 mutations on the basis of
  # family history for combined family and cancer data
  # Author: Gang Peng
  # 
  # Input:
  # fam.cancer.data: family and cancer data including at least id, gender, 
  # fid (father id), mid (mother id), age (current age or death age), 
  # num.cancer (number of invasive cancers), cancer.info (list, for each individual
  # store diag.cancer and diag.age)
  # Output:
  #
  # the posterior probability for p53 variant
  lik <- calLK(fam.cancer.data, penetrance.all)
  
  ##################################
  # convert data
  ##################################
  id <- as.integer(fam.cancer.data$id)
  fid <- as.integer(fam.cancer.data$fid)
  mid <- as.integer(fam.cancer.data$mid)
  counselee.id <- as.integer(counselee.id)
  
  if(min(id)==0){
    id <- id+1
    fid <- fid+1
    mid <- mid+1
    counselee.id <- counselee.id+1
  }
  fid[is.na(fid)] <- 0
  mid[is.na(mid)] <- 0
  
  ped <- data.frame(ID=id, Gender = fam.cancer.data$gender,
                    FatherID = fid, 
                    MotherID = mid, 
                    stringsAsFactors=FALSE)
  
#   print(allef)
#   print(lik)
#   print(ped)
#   print(counselee.id)
#   print(nloci)
#   print(mRate)
  pp <- peelingRC(allef, lik, ped,
                 counselee.id, nloci, mRate)
  return(1-pp[,1])
}

calLK <- function(fam.cancer.data, penetrance.all){
  # Aim: calculate the likelihood, Pr(D|G) for each individual
  # Author: Gang Peng
  #
  # input: 
  # fam.cancer.data: family data
  # penetrance.all: penetrance in the population for male and female
  #
  # output:
  # likelihood
  
  num.individual <- length(fam.cancer.data$id)
  lik <- matrix(NA, nrow=num.individual, ncol=3)
  
  # penetrance length, oldest age in the penetrance table
  length.pene <- nrow(penetrance.all$fMX)
  
  for(i in 1:num.individual){
    if(fam.cancer.data$num.cancer[i]==0){
      if(is.na(fam.cancer.data$gender[i])){
        pene.mix <- (penetrance.all$fMX + penetrance.all$fFX)/2
        lik[i,] <- lkNoneAffect(pene.mix, fam.cancer.data$age[i])
      }
      else{        
        if(fam.cancer.data$gender[i]==1){
          lik[i,] <- lkNoneAffect(penetrance.all$fMX, fam.cancer.data$age[i])
        }
        else{
          lik[i,] <- lkNoneAffect(penetrance.all$fFX, fam.cancer.data$age[i])
        }
      }
    }
    else{
      flag.benign <- TRUE
      for(j in 1:fam.cancer.data$num.cancer[i]){
        if(fam.cancer.data$cancer.info[[i]]$cancer.type[j] < invasive.cut){
          flag.benign <- FALSE
          break
        }
      }
      
      if(flag.benign){
        if(is.na(fam.cancer.data$gender[i])){
          pene.mix <- (penetrance.all$fMX + penetrance.all$fFX)/2
          lik[i,] <- lkNoneAffect(pene.mix, fam.cancer.data$age[i])
        }
        else{        
          if(fam.cancer.data$gender[i]==1){
            lik[i,] <- lkNoneAffect(penetrance.all$fMX, fam.cancer.data$age[i])
          }
          else{
            lik[i,] <- lkNoneAffect(penetrance.all$fFX, fam.cancer.data$age[i])
          }
        }
      }
      else
      {
        age.tmp <- min(fam.cancer.data$cancer.info[[i]]$diag.age)
        if(age.tmp < 1){
          age.tmp = 1
        }
        if(age.tmp > length.pene){
          age.tmp = length.pene
        }
        if(is.na(fam.cancer.data$gender[i])){
          lik[i,] <- (penetrance.all$fMX[age.tmp,]+penetrance.all$fFX[age.tmp,])/2
        }
        else{
          if(fam.cancer.data$gender[i]==1){
            lik[i,] <- penetrance.all$fMX[age.tmp,]
          }
          else{
            lik[i,] <- penetrance.all$fFX[age.tmp,]
          }
        }
      }
      
      
    }
  }
  
  return(lik)
}


lkNoneAffect <- function(penetrance, age){
  # calculate the likelihood for unaffected sample
  # Args:
  #   age: sample's age when he/she still doesn't get disease
  #
  # Returns:
  #   likelihood
  
  length.pene <- nrow(penetrance)
  if(age==0)
  {
    age <- 1
  }
  if(age > length.pene)
  {
    age <- length.pene
  }
  rlt <- rep(1,3)
  rlt[1] <- 1-sum(penetrance[1:age,1])
  rlt[2] <- 1-sum(penetrance[1:age,2])
  rlt[3] <- 1-sum(penetrance[1:age,3])
  #   for(i in 1:age){
  #     rlt[1] <- rlt[1]*(1-penetrance[i,1])
  #     rlt[2] <- rlt[2]*(1-penetrance[i,2])
  #     rlt[3] <- rlt[3]*(1-penetrance[i,3])
  #   }
  rlt
}