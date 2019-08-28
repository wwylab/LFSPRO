reformatForClassicChompret <- function(fam.cancer.data){
  # Aim: reformat new data type to the format suit for clinical criteria
  # Author: Gang Peng
  #
  # Input:
  # fam.cancer.data: A list, formatted for input in LFSpro.
  #
  # Output:
  # the data with format for clinical criteria
  # ID, fID, mID, primaryInfo
  ID <- fam.cancer.data$id
  fID <- fam.cancer.data$fid
  mID <- fam.cancer.data$mid
  
  primaryInfo <- list()
  for(i in 1:length(ID)){
    diagAge <- NULL
    cancerType <- NULL
    if(fam.cancer.data$num.cancer[i]>0){
      diagAge <- fam.cancer.data$cancer.info[[i]]$diag.age
      cancerType <- fam.cancer.data$cancer.info[[i]]$cancer.type
    }
    else{
      diagAge <- NA
      cancerType <- NA
    }
    primaryInfo[[i]] <- list(diagAge=diagAge, cancerType=cancerType)
  }
  
  return(list(ID=ID, fID=fID, mID=mID,primaryInfo=primaryInfo))
}

combinedata <- function(fam.data, cancer.data){
  # Aim: combine family information data and cancer information data into 
  #    one unified data type
  # Author: Gang Peng 
  # Input:
  #    fam.info: data frame containing family information
  #    cancer.info: data frame containing cancer information
  # Output:
  #    Combined family data including at least id, gender, fid (father id), 
  #    mid (mother id), age (current age or death age), 
  #    num.cancer (number of invasive cancers), 
  #    cancer.info (list, for each individual store diag.cancer and diag.age)
  #    This data is used in function LFSpro
  
  fam.id.all <- unique(fam.data$fam.id)
  num.fam <- length(fam.id.all)
  
  fam.cancer.data <- list(num.fam)
  
  for(i in 1:num.fam){
    fam.data.tmp <- fam.data[fam.data$fam.id == fam.id.all[i],]
    cancer.data.tmp <- cancer.data[cancer.data$fam.id == fam.id.all[i],]
    
    num.ind.tmp <- length(fam.data.tmp$fam.id)
    
#     age <- 1:num.ind.tmp
#     for(j in 1:num.ind.tmp){
#       if(fam.data.tmp$vital[j]=="A"){
#         age[j] = calAge(fam.data.tmp$lcd[j],fam.data.tmp$dob[j])
#       }
#       else{
#         age[j] = calAge(fam.data.tmp$dod[j],fam.data.tmp$dob[j])
#       }
#     }
#     
#     if(sum(is.na(age))>0){
#       print(paste("Errors in age for sample ", fam.data.tmp$fam.id[is.na(age)], 
#                   fam.data.tmp$id[is.na(age)],sep=" "))
#     }
    
    cancer.info <- list(num.ind.tmp)
    uniq.id.fam <- paste(fam.data.tmp$fam.id,fam.data.tmp$id,sep="")
    uniq.id.cancer <- paste(cancer.data.tmp$fam.id,cancer.data.tmp$id,sep="")
    num.cancer <- rep(0, num.ind.tmp)
    for(j in 1:num.ind.tmp){
      index <- uniq.id.cancer %in% uniq.id.fam[j]
      num.cancer.tmp <- sum(index)
      num.cancer[j] <- num.cancer.tmp
      if(num.cancer.tmp > 0){
        cancer.type <- cancer.data.tmp$cancer.type[index]
        #diag.date <- cancer.data.tmp$diag.date[index]
        diag.age <- cancer.data.tmp$diag.age[index]
#         diag.age <- rep(0,num.cancer.tmp)
#         for(k in 1:num.cancer.tmp){
#           diag.age[k] <- calAge(diag.date[k], fam.data.tmp$dob[j])
#         }
        
        cancer.info[[j]] <- data.frame(cancer.type,diag.age)
      }
      else
      {
        cancer.info[[j]] <- data.frame()
      }
    }
    
    fam.cancer.data[[i]] <- list(fam.id = fam.data.tmp$fam.id, 
                                 id = fam.data.tmp$id,
                                 fid = fam.data.tmp$fid, 
                                 mid = fam.data.tmp$mid,
                                 gender = fam.data.tmp$gender, 
                                 vital = fam.data.tmp$vital,
                                 age = fam.data.tmp$age, 
                                 num.cancer = num.cancer,
                                 cancer.info = cancer.info)
  }
  
  return(fam.cancer.data)
}
