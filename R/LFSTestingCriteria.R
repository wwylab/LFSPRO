# first degree relative
# input: 
# pedigree: data frame, ID, fID, mID (ID>=0)
# ii: sample ii (index)


# Avoid the binding error. 
 utils::globalVariables(c("lfspro.cancer.type", "lfs.cut", "cancer.type.all", "parameter.mpc", "lfspenet.2010",
                          "lfspenet.cs.death", "invasive.cut"))

# if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))


firstDegreeRelative=function(pedigree,ii)
{
  rlt=NULL
  numSample=dim(pedigree)[1]
  for(i in 1:numSample)
  {
    if(is.na(pedigree$fID[i]) && !is.na(pedigree$mID[i]))
    {
      pedigree$fID[i]=-999
    }
    if(!is.na(pedigree$fID[i]) && is.na(pedigree$mID[i]))
    {
      pedigree$mID[i]=-9999
    }
  }
  fID=pedigree$fID[ii]
  mID=pedigree$mID[ii]
  ID=pedigree$ID[ii]
  for(i in 1:numSample)
  {
    if(i==ii)
    {
      next
    }
    #sibling
    if((pedigree$fID[i]==fID && pedigree$mID[i]==mID) && !is.na((pedigree$fID[i]==fID && pedigree$mID[i]==mID)))
    {
      rltTmp=data.frame(matrix(data=NA,nrow=1,ncol=2))
      names(rltTmp)=c("index","lineage")
      rltTmp[1]=i
      rltTmp[2]="S"
      rlt=rbind(rlt,rltTmp)
      next
    }
    #father
    if(pedigree$ID[i]==fID && !is.na(pedigree$ID[i]==fID))
    {
      rltTmp=data.frame(matrix(data=NA,nrow=1,ncol=2))
      names(rltTmp)=c("index","lineage")
      rltTmp[1]=i
      rltTmp[2]="F"
      rlt=rbind(rlt,rltTmp)
      next
    }
    #mother
    if(pedigree$ID[i]==mID && !is.na(pedigree$ID[i]==mID))
    {
      rltTmp=data.frame(matrix(data=NA,nrow=1,ncol=2))
      names(rltTmp)=c("index","lineage")
      rltTmp[1]=i
      rltTmp[2]="M"
      rlt=rbind(rlt,rltTmp)
      next
    }
    #child
    if((pedigree$fID[i]==ID || pedigree$mID[i]==ID) && !is.na(pedigree$fID[i]==ID || pedigree$mID[i]==ID))
    {
      rltTmp=data.frame(matrix(data=NA,nrow=1,ncol=2))
      names(rltTmp)=c("index","lineage")
      rltTmp[1]=i
      rltTmp[2]="C"
      rlt=rbind(rlt,rltTmp)
      next
    }
  }
  rlt
}

# second degree relative
# input: 
# pedigree: data frame, ID, fID, mID (ID>=0)
# ii: sample ii (index)
secondDegreeRelative=function(pedigree,ii)
{
  numSample=dim(pedigree)[1]
  for(i in 1:numSample)
  {
    if(is.na(pedigree$fID[i]) && !is.na(pedigree$mID[i]))
    {
      pedigree$fID[i]=-999
    }
    if(!is.na(pedigree$fID[i]) && is.na(pedigree$mID[i]))
    {
      pedigree$mID[i]=-9999
    }
  }
  fdRelative=firstDegreeRelative(pedigree,ii)
  
  rlt=NULL
  if(is.null(fdRelative))
  {
    return (rlt)
  }
  #fID=pedigree$fID[ii]
  #mID=pedigree$mID[ii]
  numFDR=dim(fdRelative)[1]
  for(i in 1:numFDR)
  {
    if(fdRelative[i,2]=="F")
    {
      fID=pedigree$fID[fdRelative[i,1]]
      mID=pedigree$mID[fdRelative[i,1]]
      ID=pedigree$ID[fdRelative[i,1]]
      for(j in 1:numSample)
      {
        if(j==fdRelative[i,1] || j==ii)
        {
          next
        }
        #Fahter's parents
        if((pedigree$ID[j]==fID || pedigree$ID[j]==mID) && !is.na(pedigree$ID[j]==fID || pedigree$ID[j]==mID))
        {
          rltTmp=data.frame(matrix(data=NA,nrow=1,ncol=2))
          names(rltTmp)=c("index","lineage")
          rltTmp[1]=j
          rltTmp[2]="F"
          rlt=rbind(rlt,rltTmp)
          next
        }
        #Father's siblings
        if((pedigree$fID[j]==fID && pedigree$mID[j]==mID) && !is.na(pedigree$fID[j]==fID && pedigree$mID[j]==mID))
        {
          rltTmp=data.frame(matrix(data=NA,nrow=1,ncol=2))
          names(rltTmp)=c("index","lineage")
          rltTmp[1]=j
          rltTmp[2]="F"
          rlt=rbind(rlt,rltTmp)
          next
        }
        #Father's child (half siblings)
        if(pedigree$fID[j]==ID && !is.na(pedigree$fID[j]==ID))
        {
          if(!(j %in% fdRelative[,1]))
          {
            rltTmp=data.frame(matrix(data=NA,nrow=1,ncol=2))
            names(rltTmp)=c("index","lineage")
            rltTmp[1]=j
            rltTmp[2]="F"
            rlt=rbind(rlt,rltTmp)
            next
          }
        }
      }
      next
    }
    if(fdRelative[i,2]=="M")
    {
      fID=pedigree$fID[fdRelative[i,1]]
      mID=pedigree$mID[fdRelative[i,1]]
      ID=pedigree$ID[fdRelative[i,1]]
      for(j in 1:numSample)
      {
        if(j==fdRelative[i,1] || j==ii)
        {
          next
        }
        #Mohter's parents
        if((pedigree$ID[j]==fID || pedigree$ID[j]==mID) && !is.na(pedigree$ID[j]==fID || pedigree$ID[j]==mID))
        {
          rltTmp=data.frame(matrix(data=NA,nrow=1,ncol=2))
          names(rltTmp)=c("index","lineage")
          rltTmp[1]=j
          rltTmp[2]="M"
          rlt=rbind(rlt,rltTmp)
          next
        }
        #Mother's siblings
        if((pedigree$fID[j]==fID && pedigree$mID[j]==mID) && !is.na(pedigree$fID[j]==fID && pedigree$mID[j]==mID))
        {
          rltTmp=data.frame(matrix(data=NA,nrow=1,ncol=2))
          names(rltTmp)=c("index","lineage")
          rltTmp[1]=j
          rltTmp[2]="M"
          rlt=rbind(rlt,rltTmp)
          next
        }
        #Mother's child (half siblings)
        if(pedigree$mID[j]==ID && !is.na(pedigree$mID[j]==ID))
        {
          if(!(j %in% fdRelative[,1]))
          {
            rltTmp=data.frame(matrix(data=NA,nrow=1,ncol=2))
            names(rltTmp)=c("index","lineage")
            rltTmp[1]=j
            rltTmp[2]="M"
            rlt=rbind(rlt,rltTmp)
            next
          }
        }
      }
      next
    }
    if(fdRelative[i,2]=="S")
    {
      ID=pedigree$ID[fdRelative[i,1]]
      for(j in 1:numSample)
      {
        if(j==fdRelative[i,1] || j==ii)
        {
          next
        }
        #sibling's child
        if((pedigree$fID[j]==ID || pedigree$mID[j]==ID) && !is.na(pedigree$fID[j]==ID || pedigree$mID[j]==ID))
        {
          rltTmp=data.frame(matrix(data=NA,nrow=1,ncol=2))
          names(rltTmp)=c("index","lineage")
          rltTmp[1]=j
          rltTmp[2]="S"
          rlt=rbind(rlt,rltTmp)
          next
        }
      }
      next
    }
    if(fdRelative[i,2]=="C")
    {
      ID=pedigree$ID[fdRelative[i,1]]
      for(j in 1:numSample)
      {
        if(j==fdRelative[i,1] || j==ii)
        {
          next
        }
        #child's child
        if((pedigree$fID[j]==ID || pedigree$mID[j]==ID) && !is.na(pedigree$fID[j]==ID || pedigree$mID[j]==ID))
        {
          rltTmp=data.frame(matrix(data=NA,nrow=1,ncol=2))
          names(rltTmp)=c("index","lineage")
          rltTmp[1]=j
          rltTmp[2]="C"
          rlt=rbind(rlt,rltTmp)
          next
        }
      }
      next
    }
  }
  rlt
}

lfsClassic=function(fam.data, cancer.data, counselee.id){
  # Classic Li-Fraumeni syndrome criteria
  # Input:
  # fam.data: family information data
  # cancer.data: cancer information data
  # counselee.id: data frame. the family id and id for the counselees
  
  # convert cancer type to specific number and check the cancer type
  num.cancer <- nrow(cancer.data)
  cancer.type.num <- rep(-1, num.cancer)
  for(i in 1:num.cancer){
    tmp <- lfspro.cancer.type[cancer.data$cancer.type[i]]
    if(is.na(tmp)){
      print(paste("Cannot find cancer ", cancer.data$cancer.type[i], 
                  " in the LFSpro predefined cancer type", sep = ""))
      print("LFSpro predefined cancer types are: ")
      print(cancer.type.all)
      print("Please check the input cancer information data.")
      
      num.counselee <- nrow(counselee.id)
      result <- rep(NA, num.counselee)
      
      rlt <- data.frame(cbind(counselee.id, result),check.names = FALSE)
      colnames(rlt) <- c("fam.id", "id", "result")
      return(rlt)
    }
    cancer.type.num[i] <- tmp
  }
  
  cancer.data$cancer.type <- cancer.type.num
  
  fam.cancer.data <- combinedata(fam.data, cancer.data)
  num.fam <- length(fam.cancer.data)
  
  rlt.all <- NULL
  for(ii in 1:num.fam){
    
    cid <- counselee.id$id[counselee.id$fam.id == fam.cancer.data[[ii]]$fam.id[1]]
    if(length(cid)<1){
      print(paste("Cannot find any counselee id in family ", 
                  fam.cancer.data[[ii]]$fam.id[1], sep=""))
      next
    }
    
    lfsData <- reformatForClassicChompret(fam.cancer.data[[ii]])
    numSample=length(lfsData$ID)
    pedigree=data.frame(lfsData$ID,lfsData$fID,lfsData$mID,stringsAsFactors=FALSE)
    names(pedigree)=c("ID","fID","mID")
    rlt=rep(FALSE,numSample)
    c45=rep(FALSE,numSample)
    sacoma=rep(FALSE,numSample)
    for(i in 1:numSample)
    {
      for(j in 1:length(lfsData$primaryInfo[[i]]$cancerType))
      {
        if(!is.na(lfsData$primaryInfo[[i]]$cancerType[j]))
        {
          if(!is.na(lfsData$primaryInfo[[i]]$diagAge[j]) && lfsData$primaryInfo[[i]]$diagAge[j]<45
             && lfsData$primaryInfo[[i]]$cancerType[j] < invasive.cut)
          {
            c45[i]=TRUE
          }
          if(lfsData$primaryInfo[[i]]$cancerType[j] ==  lfspro.cancer.type["ost"] || 
               lfsData$primaryInfo[[i]]$cancerType[j] == lfspro.cancer.type["sts"])
          {
            sacoma[i]=TRUE
          }
        }
      }
    }
    for(i in 1:numSample)
    {
      flagSarcoma=FALSE
      if(sacoma[i] && c45[i])
      {
        fdRelative=firstDegreeRelative(pedigree,i)
        numFDR=dim(fdRelative)[1]
        cancerInfo1=NULL
        if(!is.null(numFDR) && numFDR>0){
          for(j in 1:numFDR)
          {
            if(c45[fdRelative[j,1]] || sacoma[fdRelative[j,1]])
            {
              cancerInfoTmp=as.data.frame(matrix(data=NA,nrow=1,ncol=3))
              names(cancerInfoTmp)=c("c45","sacoma","lineage")
              cancerInfoTmp$c45=c45[fdRelative[j,1]]
              cancerInfoTmp$sacoma=sacoma[fdRelative[j,1]]
              cancerInfoTmp$lineage=fdRelative[j,2]
              cancerInfo1=rbind(cancerInfo1,cancerInfoTmp)
            }
          }
        }
        sdRelative=secondDegreeRelative(pedigree,i)
        numSDR=dim(sdRelative)[1]
        cancerInfo2=NULL
        if(!is.null(numSDR) && numSDR>0){
          for(j in 1:numSDR)
          {
            cancerInfoTmp=as.data.frame(matrix(data=NA,nrow=1,ncol=3))
            names(cancerInfoTmp)=c("c45","sacoma","lineage")
            cancerInfoTmp$c45=c45[sdRelative[j,1]]
            cancerInfoTmp$sacoma=sacoma[sdRelative[j,1]]
            cancerInfoTmp$lineage=sdRelative[j,2]
            cancerInfo2=rbind(cancerInfo2,cancerInfoTmp)
          }
        }
        
        if(is.null(cancerInfo1))
        {
          next
        }
        for(j in 1:dim(cancerInfo1)[1])
        {
          if(cancerInfo1$c45[j] && cancerInfo1$sacoma[j])
          {
            lineage=cancerInfo1$lineage[j]
            for(k in 1:dim(cancerInfo1)[1])
            {
              if(k==j)
              {
                next
              }
              if(cancerInfo1$c45[k] || cancerInfo1$sacoma[k])
              {
                if(lineage=="S" || lineage=="C")
                {
                  rlt[i]=TRUE
                  break
                }
                else if(lineage=="F")
                {
                  if(!(cancerInfo1$lineage[k]=="M"))
                  {
                    rlt[i]==TRUE
                    break
                  }
                }
                else if(lineage=="M")
                {
                  if(!(cancerInfo1$lineage[k]=="F"))
                  {
                    rlt[i]==TRUE
                    break
                  }
                }
                else
                {
                  print("ERROR: unrecognized lineage type")
                }
              }
            }
            if(rlt[i])
            {
              break
            }
            
            if(is.null(cancerInfo2))
            {
              next
            }
            for(k in 1:dim(cancerInfo2)[1])
            {
              if(cancerInfo2$c45[k] || cancerInfo2$sacoma[k])
              {
                if(lineage=="S" || lineage=="C")
                {
                  rlt[i]=TRUE
                  break
                }
                else if(lineage=="F")
                {
                  if(!(cancerInfo2$lineage[k]=="M"))
                  {
                    rlt[i]==TRUE
                    break
                  }
                }
                else if(lineage=="M")
                {
                  if(!(cancerInfo2$lineage[k]=="F"))
                  {
                    rlt[i]==TRUE
                    break
                  }
                }
                else
                {
                  print("ERROR: unrecognized lineage type")
                }
              }
            }
            if(rlt[i])
            {
              break
            }
            
          }
        }
        
      }
    }
    
    rlt.all <- c(rlt.all, rlt[fam.cancer.data[[ii]]$id %in% cid])
  }
  
  rlt <- data.frame(cbind(counselee.id, rlt.all), check.names = FALSE)
  colnames(rlt) <- c("fam.id", "id", "result")
  return(rlt)
}


lfsChompret2009=function(fam.data, cancer.data, counselee.id)
{
  # Chompret criteria
  # Input:
  # fam.data: family information data
  # cancer.data: cancer information data
  # counselee.id: data frame. the family id and id for the counselees
  
  # convert cancer type to specific number and check the cancer type
  
  num.cancer <- nrow(cancer.data)
  cancer.type.num <- rep(-1, num.cancer)
  for(i in 1:num.cancer){
    tmp <- lfspro.cancer.type[cancer.data$cancer.type[i]]
    if(is.na(tmp)){
      print(paste("Cannot find cancer ", cancer.data$cancer.type[i], 
                  " in the LFSpro predefined cancer type", sep = ""))
      print("LFSpro predefined cancer types are: ")
      print(cancer.type.all)
      print("Please check the input cancer information data.")
      
      num.counselee <- nrow(counselee.id)
      result <- rep(NA, num.counselee)
      
      rlt <- data.frame(cbind(counselee.id, result),check.names = FALSE)
      colnames(rlt) <- c("fam.id", "id", "result")
      return(rlt)
    }
    cancer.type.num[i] <- tmp
  }
  
  cancer.data$cancer.type <- cancer.type.num
  
  fam.cancer.data <- combinedata(fam.data, cancer.data)
  num.fam <- length(fam.cancer.data)
  
  rlt.all <- NULL
  for(ii in 1:num.fam){
    
    cid <- counselee.id$id[counselee.id$fam.id == fam.cancer.data[[ii]]$fam.id[1]]
    if(length(cid)<1){
      print(paste("Cannot find any counselee id in family ", 
                  fam.cancer.data[[ii]]$fam.id[1], sep=""))
      next
    }
    
    lfsData <- reformatForClassicChompret(fam.cancer.data[[ii]])
    
    numSample=length(lfsData$ID)
    rlt=rep(FALSE,numSample)
    
    CancerInfo=NULL
    for(i in 1:numSample)
    {
      CancerInfoTmp=data.frame(matrix(data=NA,nrow=1,ncol=4))
      names(CancerInfoTmp)=c("bAge","oAge","multi","nonlfs")
      numCancer=0;
      for(j in 1:length(lfsData$primaryInfo[[i]]$cancerType))
      {
        cancerType <- lfsData$primaryInfo[[i]]$cancerType[j]
        if(!is.na(cancerType))
        {
          ageTmp=lfsData$primaryInfo[[i]]$diagAge[j]
          if(cancerType == lfspro.cancer.type["breast"])
          {
            if(is.na(CancerInfoTmp$bAge))
            {
              CancerInfoTmp$bAge=ageTmp
            }
            else
            {
              if(CancerInfoTmp$bAge>ageTmp)
              {
                CancerInfoTmp$bAge=ageTmp
              }
            }
            numCancer=numCancer+1
            next
          }
          
          if(cancerType < lfs.cut)
          {
            if(is.na(CancerInfoTmp$oAge))
            {
              CancerInfoTmp$oAge=ageTmp
            }
            else
            {
              if(CancerInfoTmp$oAge>ageTmp)
              {
                CancerInfoTmp$oAge=ageTmp
              }
            }
            numCancer=numCancer+1
            next
          }
        }
      }
      
      CancerInfoTmp$multi=numCancer
      CancerInfo=rbind(CancerInfo, CancerInfoTmp)
    }
    
    for(i in 1:numSample)
    {
      #term 3
      for(j in 1:length(lfsData$primaryInfo[[i]]$cancerType))
      {
        if(!is.na(lfsData$primaryInfo[[i]]$cancerType[j]))
        {
          if(lfsData$primaryInfo[[i]]$cancerType[j] == lfspro.cancer.type["acc"] || 
               lfsData$primaryInfo[[i]]$cancerType[j] == lfspro.cancer.type["choroid"])
          {
            rlt[i]=TRUE
            break
          }
        }
      }
      if(rlt[i])
      {
        next
      }
      
      #term 2
      if(length(lfsData$primaryInfo[[i]]$cancerType)>1)
      {
        lfs.cancer <- lfsData$primaryInfo[[i]]$cancerType[lfsData$primaryInfo[[i]]$cancerType<lfs.cut]
        num.breast.cancer <- sum(lfs.cancer == lfspro.cancer.type["breast"])
        if(num.breast.cancer > 0)
        {
          num.lfs.cancer= length(lfs.cancer) - num.breast.cancer + 1
        }
        else
        {
          num.lfs.cancer= length(lfs.cancer)
        }
        
        if(num.lfs.cancer > 1)
        {
          if(min(lfsData$primaryInfo[[i]]$diagAge[lfsData$primaryInfo[[i]]$cancerType<lfs.cut]) < 46)
          {
            rlt[i] <- TRUE
          }
        }
      }
      
      if(rlt[i])
      {
        next
      }
      
      #term 1
      
      if((CancerInfo$bAge[i]<46 || CancerInfo$oAge[i]<46) && !is.na(CancerInfo$bAge[i]<46 || CancerInfo$oAge[i]<46))
      {
        pedigree=data.frame(lfsData$ID,lfsData$fID,lfsData$mID,stringsAsFactors=FALSE)
        names(pedigree)=c("ID","fID","mID")
        fdRelative=firstDegreeRelative(pedigree,i)
        if(is.null(fdRelative))
        {
          next
        }
        sdRelative=secondDegreeRelative(pedigree,i)
        if(!is.null(fdRelative))
        {
          for(j in 1:dim(fdRelative)[1])
          {
            if(CancerInfo$multi[fdRelative[j,1]]>1)
            {
              rlt[i]=TRUE
              break
            }
            else
            {
              if(is.na(CancerInfo$bAge[i]))
              {
                if(!(is.na(CancerInfo$bAge[fdRelative[j,1]]<56 || CancerInfo$oAge[fdRelative[j,1]]<56)) && 
                  (CancerInfo$bAge[fdRelative[j,1]]<56 || CancerInfo$oAge[fdRelative[j,1]]<56))
                {
                  rlt[i]=TRUE
                  break
                }
              }
              else
              {
                if(!is.na(CancerInfo$oAge[fdRelative[j,1]]) && CancerInfo$oAge[fdRelative[j,1]]<56)
                {
                  rlt[i]=TRUE
                  break
                }
              }
            }
          }
          
          if(rlt[i])
          {
            next
          }
        }
        
        if(!is.null(sdRelative))
        {
          for(j in 1:dim(sdRelative)[1])
          {
            if(CancerInfo$multi[sdRelative[j,1]]>1)
            {
              rlt[i]=TRUE
              break
            }
            else
            {
              if(is.na(CancerInfo$bAge[i]))
              {
                if(!(is.na(CancerInfo$bAge[sdRelative[j,1]]<56 || CancerInfo$oAge[sdRelative[j,1]]<56)) && 
                  (CancerInfo$bAge[sdRelative[j,1]]<56 || CancerInfo$oAge[sdRelative[j,1]]<56))
                {
                  rlt[i]=TRUE
                  break
                }
              }
              else
              {
                if(!is.na(CancerInfo$oAge[sdRelative[j,1]]) && CancerInfo$oAge[sdRelative[j,1]]<56)
                {
                  rlt[i]=TRUE
                  break
                }
              }
            }
          }
        }
      }
      
    }
    rlt.all <- c(rlt.all, rlt[fam.cancer.data[[ii]]$id %in% cid])
  }
  
  rlt <- data.frame(cbind(counselee.id, rlt.all), check.names = FALSE)
  colnames(rlt) <- c("fam.id", "id", "result")
  return(rlt)
}


lfsChompret2015=function(fam.data, cancer.data, counselee.id)
{
  # Chompret criteria
  # Input:
  # fam.data: family information data
  # cancer.data: cancer information data
  # counselee.id: data frame. the family id and id for the counselees
  
  # convert cancer type to specific number and check the cancer type
  num.cancer <- nrow(cancer.data)
  cancer.type.num <- rep(-1, num.cancer)
  for(i in 1:num.cancer){
    tmp <- lfspro.cancer.type[cancer.data$cancer.type[i]]
    if(is.na(tmp)){
      print(paste("Cannot find cancer ", cancer.data$cancer.type[i], 
                  " in the LFSpro predefined cancer type", sep = ""))
      print("LFSpro predefined cancer types are: ")
      print(cancer.type.all)
      print("Please check the input cancer information data.")
      
      num.counselee <- nrow(counselee.id)
      result <- rep(NA, num.counselee)
      
      rlt <- data.frame(cbind(counselee.id, result),check.names = FALSE)
      colnames(rlt) <- c("fam.id", "id", "result")
      return(rlt)
    }
    cancer.type.num[i] <- tmp
  }
  
  cancer.data$cancer.type <- cancer.type.num
  
  fam.cancer.data <- combinedata(fam.data, cancer.data)
  num.fam <- length(fam.cancer.data)
  
  rlt.all <- NULL
  for(ii in 1:num.fam){
    
    cid <- counselee.id$id[counselee.id$fam.id == fam.cancer.data[[ii]]$fam.id[1]]
    if(length(cid)<1){
      print(paste("Cannot find any counselee id in family ", 
                  fam.cancer.data[[ii]]$fam.id[1], sep=""))
      next
    }
    
    lfsData <- reformatForClassicChompret(fam.cancer.data[[ii]])
    
    numSample=length(lfsData$ID)
    rlt=rep(FALSE,numSample)
    
    CancerInfo=NULL
    for(i in 1:numSample)
    {
      CancerInfoTmp=data.frame(matrix(data=NA,nrow=1,ncol=4))
      names(CancerInfoTmp)=c("bAge","oAge","multi","nonlfs")
      numCancer=0;
      for(j in 1:length(lfsData$primaryInfo[[i]]$cancerType))
      {
        cancerType <- lfsData$primaryInfo[[i]]$cancerType[j]
        if(!is.na(cancerType))
        {
          ageTmp=lfsData$primaryInfo[[i]]$diagAge[j]
          if(cancerType == lfspro.cancer.type["breast"])
          {
            if(is.na(CancerInfoTmp$bAge))
            {
              CancerInfoTmp$bAge=ageTmp
            }
            else
            {
              if(CancerInfoTmp$bAge>ageTmp)
              {
                CancerInfoTmp$bAge=ageTmp
              }
            }
            numCancer=numCancer+1
            next
          }
          
          if(cancerType < lfs.cut && cancerType != lfspro.cancer.type["leukemia"] && cancerType != lfspro.cancer.type["lung"])
          {
            if(is.na(CancerInfoTmp$oAge))
            {
              CancerInfoTmp$oAge=ageTmp
            }
            else
            {
              if(CancerInfoTmp$oAge>ageTmp)
              {
                CancerInfoTmp$oAge=ageTmp
              }
            }
            numCancer=numCancer+1
            next
          }
        }
      }
      
      CancerInfoTmp$multi=numCancer
      CancerInfo=rbind(CancerInfo, CancerInfoTmp)
    }
    
    for(i in 1:numSample)
    {
      #term breast
      for(j in 1:length(lfsData$primaryInfo[[i]]$cancerType))
      {
        if(!is.na(lfsData$primaryInfo[[i]]$cancerType[j]))
        {
          if(lfsData$primaryInfo[[i]]$cancerType[j] == lfspro.cancer.type["breast"] && 
             lfsData$primaryInfo[[i]]$diagAge[j] < 31)
          {
            rlt[i]=TRUE
            break
          }
        }
      }
      if(rlt[i])
      {
        next
      }
      #term 3
      for(j in 1:length(lfsData$primaryInfo[[i]]$cancerType))
      {
        if(!is.na(lfsData$primaryInfo[[i]]$cancerType[j]))
        {
          if(lfsData$primaryInfo[[i]]$cancerType[j] == lfspro.cancer.type["acc"] || 
             lfsData$primaryInfo[[i]]$cancerType[j] == lfspro.cancer.type["choroid"])
          {
            rlt[i]=TRUE
            break
          }
        }
      }
      if(rlt[i])
      {
        next
      }
      
      #term 2
      if(length(lfsData$primaryInfo[[i]]$cancerType)>1)
      {
        lfs.cancer <- lfsData$primaryInfo[[i]]$cancerType[lfsData$primaryInfo[[i]]$cancerType<lfs.cut & lfsData$primaryInfo[[i]]$cancerType != lfspro.cancer.type["leukemia"] & lfsData$primaryInfo[[i]]$cancerType != lfspro.cancer.type["lung"]]
        num.breast.cancer <- sum(lfs.cancer == lfspro.cancer.type["breast"])
        if(num.breast.cancer > 0)
        {
          num.lfs.cancer= length(lfs.cancer) - num.breast.cancer + 1
        }
        else
        {
          num.lfs.cancer= length(lfs.cancer)
        }
        
        if(num.lfs.cancer > 1)
        {
          if(min(lfsData$primaryInfo[[i]]$diagAge[lfsData$primaryInfo[[i]]$cancerType<lfs.cut & lfsData$primaryInfo[[i]]$cancerType != lfspro.cancer.type["leukemia"] & lfsData$primaryInfo[[i]]$cancerType != lfspro.cancer.type["lung"]]) < 46)
          {
            rlt[i] <- TRUE
          }
        }
      }
      
      if(rlt[i])
      {
        next
      }
      
      #term 1
      
      if((CancerInfo$bAge[i]<46 || CancerInfo$oAge[i]<46) && !is.na(CancerInfo$bAge[i]<46 || CancerInfo$oAge[i]<46))
      {
        pedigree=data.frame(lfsData$ID,lfsData$fID,lfsData$mID,stringsAsFactors=FALSE)
        names(pedigree)=c("ID","fID","mID")
        fdRelative=firstDegreeRelative(pedigree,i)
        if(is.null(fdRelative))
        {
          next
        }
        sdRelative=secondDegreeRelative(pedigree,i)
        if(!is.null(fdRelative))
        {
          for(j in 1:dim(fdRelative)[1])
          {
            if(CancerInfo$multi[fdRelative[j,1]]>1)
            {
              rlt[i]=TRUE
              break
            }
            else
            {
              if(is.na(CancerInfo$bAge[i]))
              {
                if(!(is.na(CancerInfo$bAge[fdRelative[j,1]]<56 || CancerInfo$oAge[fdRelative[j,1]]<56)) && 
                   (CancerInfo$bAge[fdRelative[j,1]]<56 || CancerInfo$oAge[fdRelative[j,1]]<56))
                {
                  rlt[i]=TRUE
                  break
                }
              }
              else
              {
                if(!is.na(CancerInfo$oAge[fdRelative[j,1]]) && CancerInfo$oAge[fdRelative[j,1]]<56)
                {
                  rlt[i]=TRUE
                  break
                }
              }
            }
          }
          
          if(rlt[i])
          {
            next
          }
        }
        
        if(!is.null(sdRelative))
        {
          for(j in 1:dim(sdRelative)[1])
          {
            if(CancerInfo$multi[sdRelative[j,1]]>1)
            {
              rlt[i]=TRUE
              break
            }
            else
            {
              if(is.na(CancerInfo$bAge[i]))
              {
                if(!(is.na(CancerInfo$bAge[sdRelative[j,1]]<56 || CancerInfo$oAge[sdRelative[j,1]]<56)) && 
                   (CancerInfo$bAge[sdRelative[j,1]]<56 || CancerInfo$oAge[sdRelative[j,1]]<56))
                {
                  rlt[i]=TRUE
                  break
                }
              }
              else
              {
                if(!is.na(CancerInfo$oAge[sdRelative[j,1]]) && CancerInfo$oAge[sdRelative[j,1]]<56)
                {
                  rlt[i]=TRUE
                  break
                }
              }
            }
          }
        }
      }
      
    }
    rlt.all <- c(rlt.all, rlt[fam.cancer.data[[ii]]$id %in% cid])
  }
  
  rlt <- data.frame(cbind(counselee.id, rlt.all), check.names = FALSE)
  colnames(rlt) <- c("fam.id", "id", "result")
  return(rlt)
}
