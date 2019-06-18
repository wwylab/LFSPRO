utils::globalVariables(c("lfspenet.cs.nodeath", "new.lfspro.cancer.type"))
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

