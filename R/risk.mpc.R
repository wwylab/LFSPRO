risk.mpc <- function(fam.cancer.data, cancer.data, cid, data2, parameter){
  # counselee information 
  data.c <- as.data.frame(matrix(unlist(data2[data2[,1] %in% cid,]), ncol=6))
  colnames(data.c) <- c("ID", "time", "gender", "test", "D", "Dp")
  famid <- unique(fam.cancer.data$fam.id)
  cancer.data.new <- cancer.data[cancer.data$fam.id == famid,]
  ## exclude counselees with multiple cancers 
  data.agg <- aggregate(Dp~ID+gender, sum, data=data.c) # find our who had MPC
  if (nrow(data.agg)==1) {
    if (data.agg$Dp > 1) {
      data.d <- data.c[nrow(data.c),] #last row
    } else {
      data.d <- data.c
    }
  } else if (all(data.agg$Dp==0)) {
    data.d <- data.c[order(data.c$ID),]
  } else {
    id.dup <- data.agg[which(data.agg$Dp > 1),]$ID 
    data.k <- NULL
    for (i in id.dup){
      id.kp <- which.max(data.c$time[data.c$ID == i])
      data.k <- rbind(data.k, data.c[data.c$ID == i,][id.kp,])
    }
    data.k$gender <- data.k$test <- data.k$Dp <- NA
    id.agg <- data.agg[-which(data.agg$Dp > 1),]$ID 
    data.d <- rbind(data.c[data.c$ID %in% id.agg,], data.k)
    data.d <- data.d[order(data.d$ID),]
  }
  data.e <- data.d[!duplicated(data.d$ID),]
  for (i in 1:nrow(data.e)) { # output NA for counselee older than 80 or MPC already
    if (data.e$time[i] > 80|data.e$D[i] == 1) {
      warning(paste0(famid," ID=",data.e$ID[i],": No multiple primary cancer risk predictions available for counselees with age over 80, never had or had multiple primary cancers!"))
      data.e[i,][3] <- NA
      data.e[i,][4] <- NA
      data.e[i,][5] <- NA
      data.e[i,][6] <- NA
    }
  }
  age.c <- fam.cancer.data$age[fam.cancer.data$id %in% data.e$ID] 
  if (all(is.na(data.d))){
    warning(paste0(famid," ID=",unique(data.c$ID),": No multiple primary cancer risk predictions available!"))
    lik <- as.matrix(cbind(famid, unique(data.c$ID), age.c, NA, NA, NA, 
                           NA, NA, NA, NA))
  } else {
    
  ## predict penetrance for the next 5, 10, 15 years
  ## add rows for counselee at different ages with the same information else
    
    ages.n <- c()
    data.f <- NULL
    
    if (any(duplicated(data.d$ID))) { #any counselee had one primary cancer
      for (i in 1:length(age.c)) {
        id <- data.e$ID[i]
        n.rep <- age.c[i]+15
        age.i <- seq(1, n.rep, 1)  #from age 1 to age in 15 years
        ages.n <- c(ages.n, age.i)
        data.g <- matrix(rep(as.matrix(data.e, ncol=6)[i,], each = n.rep), nrow = n.rep)
        
      if (fam.cancer.data$num.cancer[which(fam.cancer.data$id == id)]==1) #Dp marked as 1 for age had cancer
        {data.g[age.c[i], 6] <- 1}
        data.f <- rbind(data.f, data.g)
      }
        } else { # all counselee are healthy
      for (i in 1:nrow(data.d)){
      n.rep <- data.d[i,2]+15
      age.i <- seq(1, n.rep, 1)  #from age 1 to age in 15 years
      ages.n <- c(ages.n, age.i)
      data.g <- matrix(rep(as.matrix(data.e, ncol=6)[i,], each = n.rep), nrow = n.rep)
      data.f <- rbind(data.f, data.g)
      }
        }
    data.f[,2] <- ages.n
  beta <- parameter$beta
  gamma <- parameter$gamma
  M <- length(gamma)

  ###################
  ## Create the baseline hazard
  ###################
  

  time <- data.f[,2]
  time[time > 80 - 1.0e-8] <- 80
  id <- data.f[,1]
  time[time == 0] <- 1.0e-12
  tilde.t <- time/80
  gender <- data.f[,3]
  dp     <- data.f[,6]
  
  n2 <- length(time)
  
  # Bernstein Basis
  w <- unlist(lapply(1:M, function(k) dbeta(tilde.t, k, M-k+1)))    
  
  ft <- matrix(w, ncol = M)

  lambda0 <- (ft %*% gamma)/80      # baseline
  
  unique.id <- unique(id)
  
  # calculate likelihood
  test0 <- rep(0, n2) # wildtype
  test1 <- rep(1, n2) # mutation
  
  # X(t)
  xp.test0.m <- cbind(test0, rep(1, n2), rep(1, n2) * test0, dp, dp * test0)
  xp.test1.m <- cbind(test1, rep(1, n2), rep(1, n2) * test1, dp, dp * test1)
  
  xp.test0.f <- cbind(test0, rep(0, n2), rep(0, n2) * test0, dp, dp * test0)
  xp.test1.f <- cbind(test1, rep(0, n2), rep(0, n2) * test1, dp, dp * test1)
  
  xpbeta.test0.m <- xp.test0.m %*% beta
  xpbeta.test1.m <- xp.test1.m %*% beta
  
  xpbeta.test0.f <- xp.test0.f %*% beta
  xpbeta.test1.f <- xp.test1.f %*% beta
  
  
  # proposed semi-parametric model for the intensity function
  # hazard probability
  log.L0.m <- lambda0 * exp(xpbeta.test0.m)
  log.L1.m <- lambda0 * exp(xpbeta.test1.m)
  
  log.L0.f <- lambda0 * exp(xpbeta.test0.f)
  log.L1.f <- lambda0 * exp(xpbeta.test1.f)
  
  like0.5 <- like1.5 <- NULL
  like0.10 <- like1.10 <- NULL
  like0.15 <- like1.15 <- NULL
  
  lik <- matrix(NA, ncol = 9)
  for (ii in unique(id)) {
    temp.id <- which(id == ii)
    age.index <- which(fam.cancer.data$id == ii)
    age <- fam.cancer.data$age[age.index]
    age.diag <- cancer.data.new[cancer.data.new$id %in% ii,]$diag.age
    
    
    if (length(age.diag)==0){
      age.diag <- age
    } else if (length(age.diag) > 0){
      age.diag <- age.diag[1]
    }
      
    
    if (any(is.na(gender[temp.id]))) { # take average of the penetrance when gender missing
      temp.llike0.m.5 <- sum(log.L0.m[temp.id[1:(age+5)]])
      temp.llike1.m.5 <- sum(log.L1.m[temp.id[1:(age+5)]]) 
      
      temp.llike0.m.10 <- sum(log.L0.m[temp.id[1:(age+10)]]) 
      temp.llike1.m.10 <- sum(log.L1.m[temp.id[1:(age+10)]])
      
      temp.llike0.m.15 <- sum(log.L0.m[temp.id[1:(age+15)]]) 
      temp.llike1.m.15 <- sum(log.L1.m[temp.id[1:(age+15)]]) 
      
      temp.llike0.f.5 <- sum(log.L0.f[temp.id[1:(age+5)]])
      temp.llike1.f.5 <- sum(log.L1.f[temp.id[1:(age+5)]]) 
      
      temp.llike0.f.10 <- sum(log.L0.f[temp.id[1:(age+10)]]) 
      temp.llike1.f.10 <- sum(log.L1.f[temp.id[1:(age+10)]]) 
      
      temp.llike0.f.15 <- sum(log.L0.f[temp.id[1:(age+15)]]) 
      temp.llike1.f.15 <- sum(log.L1.f[temp.id[1:(age+15)]]) 
      
      temp.llike1.m.c <- sum(log.L1.m[temp.id[1:age]]) # penetrance at current age
      temp.llike0.m.c <- sum(log.L0.m[temp.id[1:age]]) # penetrance at current age
      temp.llike1.f.c <- sum(log.L1.f[temp.id[1:age]]) # penetrance at current age
      temp.llike0.f.c <- sum(log.L0.f[temp.id[1:age]]) # penetrance at current age
      
      temp.llike1.m.d <- sum(log.L1.m[temp.id[1:age.diag]]) # penetrance at diagosis age
      temp.llike0.m.d <- sum(log.L0.m[temp.id[1:age.diag]]) # penetrance at diagosis age
      temp.llike1.f.d <- sum(log.L1.f[temp.id[1:age.diag]]) # penetrance at diagosis age
      temp.llike0.f.d <- sum(log.L0.f[temp.id[1:age.diag]]) # penetrance at diagosis age
      
      temp.llike0.5 <- (temp.llike0.m.5 + temp.llike0.f.5)/2
      temp.llike1.5 <- (temp.llike1.m.5 + temp.llike1.f.5)/2
      
      temp.llike0.10 <- (temp.llike0.m.10 + temp.llike0.f.10)/2
      temp.llike1.10 <- (temp.llike1.m.10 + temp.llike1.f.10)/2
      
      temp.llike0.15 <- (temp.llike0.m.15 + temp.llike0.f.15)/2
      temp.llike1.15 <- (temp.llike1.m.15 + temp.llike1.f.15)/2
      
      temp.llike0.c <- (temp.llike0.m.c + temp.llike0.f.c)/2
      temp.llike1.c <- (temp.llike1.m.c + temp.llike1.f.c)/2
      
      temp.llike0.d <- (temp.llike0.m.d + temp.llike0.f.d)/2
      temp.llike1.d <- (temp.llike1.m.d + temp.llike1.f.d)/2
      
    } else if (all(gender[temp.id] == 1))   {
      
      temp.llike0.5 <- sum(log.L0.m[temp.id[1:(age+5)]])
      temp.llike1.5 <- sum(log.L1.m[temp.id[1:(age+5)]])
      
      temp.llike0.10 <- sum(log.L0.m[temp.id[1:(age+10)]]) 
      temp.llike1.10 <- sum(log.L1.m[temp.id[1:(age+10)]]) 
      
      temp.llike0.15 <- sum(log.L0.m[temp.id[1:(age+15)]])
      temp.llike1.15 <- sum(log.L1.m[temp.id[1:(age+15)]])
      
      temp.llike0.c <- sum(log.L0.m[temp.id[1:age]])
      temp.llike1.c <- sum(log.L1.m[temp.id[1:age]])
      
      temp.llike0.d <- sum(log.L0.m[temp.id[1:age.diag]])
      temp.llike1.d <- sum(log.L1.m[temp.id[1:age.diag]])
      
    } else if (all(gender[temp.id] == 2)) {

      temp.llike0.5 <- sum(log.L0.f[temp.id[1:(age+5)]]) 
      temp.llike1.5 <- sum(log.L1.f[temp.id[1:(age+5)]]) 
      
      temp.llike0.10 <- sum(log.L0.f[temp.id[1:(age+10)]]) 
      temp.llike1.10 <- sum(log.L1.f[temp.id[1:(age+10)]]) 
      
      temp.llike0.15 <- sum(log.L0.f[temp.id[1:(age+15)]]) 
      temp.llike1.15 <- sum(log.L1.f[temp.id[1:(age+15)]])
      
      temp.llike0.c <- sum(log.L0.f[temp.id[1:age]])
      temp.llike1.c <- sum(log.L1.f[temp.id[1:age]])
  
      temp.llike0.d <- sum(log.L0.f[temp.id[1:age.diag]])
      temp.llike1.d <- sum(log.L1.f[temp.id[1:age.diag]])
      
    } else stop("check gender code: 1 for female 2 for male!")
    
    if (fam.cancer.data$num.cancer[which(fam.cancer.data$id == ii)]==1) { # if had one primary cancer
      temp.like0.c <- 1 - exp(-temp.llike0.c)
      temp.like1.c <- 1 - exp(-temp.llike1.c)
       
      temp.like0.d <- 1 - exp(-temp.llike0.d)
      temp.like1.d <- 1 - exp(-temp.llike1.d)
      
      temp.like0.5 <- (1 - exp(-temp.llike0.5) - temp.like0.c)/(1 - temp.like0.c + temp.like0.d)
      temp.like1.5 <- (1 - exp(-temp.llike1.5) - temp.like1.c)/(1 - temp.like1.c + temp.like1.d)
      
      temp.like0.10 <- (1 - exp(-temp.llike0.10) - temp.like0.c)/(1 - temp.like0.c + temp.like0.d)
      temp.like1.10 <- (1 - exp(-temp.llike1.10) - temp.like1.c)/(1 - temp.like1.c + temp.like1.d)
      
      temp.like0.15 <- (1 - exp(-temp.llike0.15) - temp.like0.c)/(1 - temp.like0.c + temp.like0.d)
      temp.like1.15 <- (1 - exp(-temp.llike1.15) - temp.like1.c)/(1 - temp.like1.c + temp.like1.d)
    } else { # if counselee is healthy
      temp.like0.c <- 1 - exp(-temp.llike0.c)
      temp.like1.c <- 1 - exp(-temp.llike1.c)
      
      temp.like0.5 <- 1 - exp(-temp.llike0.5) - temp.like0.c
      temp.like1.5 <- 1 - exp(-temp.llike1.5) - temp.like1.c
      
      temp.like0.10 <- 1 - exp(-temp.llike0.10) - temp.like0.c
      temp.like1.10 <- 1 - exp(-temp.llike1.10) - temp.like1.c
      
      temp.like0.15 <- 1 - exp(-temp.llike0.15) - temp.like0.c
      temp.like1.15 <- 1 - exp(-temp.llike1.15) - temp.like1.c
    }
    
    like0.5 <- c(like0.5, temp.like0.5)
    like1.5 <- c(like1.5, temp.like1.5)
    
    like0.10 <- c(like0.10, temp.like0.10)
    like1.10 <- c(like1.10, temp.like1.10)
    
    like0.15 <- c(like0.15, temp.like0.15)
    like1.15 <- c(like1.15, temp.like1.15)
  }
  lik <- as.matrix(cbind(famid, unique(id), age.c, like0.5, like1.5, 
                         like0.10, like1.10, like0.15, like1.15))
  }
  return(lik)
}
