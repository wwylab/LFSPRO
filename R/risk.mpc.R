risk.mpc <- function(fam.cancer.data, counselee.id, parameter){
  beta <- parameter$beta
  gamma <- parameter$gamma
  M <- length(gamma)
  
  penet <- data.frame(matrix(NA, nrow = length(counselee.id), ncol = 11))
  colnames(penet) <- c('fam id', 'id', 'age', 'gender', 'first cancer', 
                       'risk0.5', 'risk0.10', 'risk0.15', 
                       'risk1.5', 'risk1.10', 'risk1.15')
  penet[, 1] <- fam.cancer.data$fam.id[1]
  penet[, 2] <- counselee.id
  for (i in 1:length(counselee.id)) {
    continue <- TRUE
    
    id <- counselee.id[i]
    index <- which(fam.cancer.data$id == id)
    n.cancer <- fam.cancer.data$num.cancer[index]
    age <- fam.cancer.data$age[index]
    penet[i, 3] <- age
    gender <- fam.cancer.data$gender[index]
    penet[i, 4] <- gender
    
    age.pred <- c(age + 5, age + 10, age + 15)
    
    if (n.cancer > 1) {
      penet[i, 5] <- fam.cancer.data$cancer.info[[index]]$diag.age[1]
      continue <- FALSE
      warning(paste0("ID = ", id, ": No MPC risk predictions available for counselees with MPC already"))
    }
    if (age > 80) {
      continue <- FALSE
      warning(paste0("ID = ", id, ": No MPC risk predictions available for counselees over 80 years old"))
    }
    
    if (continue ==  TRUE) {
      time <- 1:(age + 15)
      time[time > 80 - 1.0e-8] <- 80
      time[time == 0] <- 1.0e-12
      tilde.t <- time / 80
      n2 <- length(time)
      
      dp <- rep(0, n2)
      if (n.cancer == 1) {
        age.diag <- fam.cancer.data$cancer.info[[index]]$diag.age
        penet[i, 5] <- age.diag
        dp[time > age.diag] <- 1
      }
      
      # Bernstein Basis
      w <- unlist(lapply(1:M, function(k){dbeta(tilde.t, k, M - k + 1)}))    
      ft <- matrix(w, ncol = M)
      lambda0 <- as.numeric((ft %*% gamma) / 80) # baseline
      
      # X(t)
      test0 <- rep(0, n2) # wildtype
      test1 <- rep(1, n2) # mutation
      
      xp.test0.m <- cbind(test0, rep(1, n2), rep(1, n2) * test0, dp, dp * test0)
      xp.test1.m <- cbind(test1, rep(1, n2), rep(1, n2) * test1, dp, dp * test1)
      
      xp.test0.f <- cbind(test0, rep(0, n2), rep(0, n2) * test0, dp, dp * test0)
      xp.test1.f <- cbind(test1, rep(0, n2), rep(0, n2) * test1, dp, dp * test1)
      
      xpbeta.test.m <- cbind(xp.test0.m %*% beta, xp.test1.m %*% beta)
      xpbeta.test.f <- cbind(xp.test0.f %*% beta, xp.test1.f %*% beta)
      
      # Cox proportional hazards model
      L.m <- lambda0 * exp(xpbeta.test.m)
      L.f <- lambda0 * exp(xpbeta.test.f)
      
      temp.like.m <- temp.like.m.c <- temp.like.m.d <- matrix(0, nrow = 3, ncol = 2)
      temp.like.f <- temp.like.f.c <- temp.like.f.d <- matrix(0, nrow = 3, ncol = 2)
      for (j in 1:3) {
        temp.like.m[j,] <- apply(L.m[1:age.pred[j],], 2, sum) 
        temp.like.f[j,] <- apply(L.f[1:age.pred[j],], 2, sum)
      }
      temp.like.m.c <- matrix(rep(apply(L.m[1:age,], 2, sum), 3), nrow = 3, byrow = TRUE) # penetrance at the current age
      temp.like.f.c <- matrix(rep(apply(L.f[1:age,], 2, sum), 3), nrow = 3, byrow = TRUE) # penetrance at the current age
      if (n.cancer == 1){
        temp.like.m.d <- matrix(rep(apply(L.m[1:age.diag,], 2, sum), 3), nrow = 3, byrow = TRUE) # penetrance at the age of last diagnosis
        temp.like.f.d <- matrix(rep(apply(L.f[1:age.diag,], 2, sum), 3), nrow = 3, byrow = TRUE) # penetrance at the age of last diagnosis
      }
      
      penet.m1 <- 1 - exp(-(temp.like.m - temp.like.m.d)) 
      penet.m2 <- 1 - exp(-(temp.like.m.c - temp.like.m.d))
      penet.m <- (penet.m1 - penet.m2) / (1 - penet.m2)
      
      penet.f1 <- 1 - exp(-(temp.like.f - temp.like.f.d)) 
      penet.f2 <- 1 - exp(-(temp.like.f.c - temp.like.f.d))
      penet.f <- (penet.f1 - penet.f2) / (1 - penet.f2)
      
      if (is.na(gender)) {
        penet[i, 6:11] <- (as.vector(penet.m) + as.vector(penet.f)) / 2
      } else if (gender == 0) {
        penet[i, 6:11] <- as.vector(penet.f)
      } else {
        penet[i, 6:11] <- as.vector(penet.m)
      }
      penet[i, 6:11][rep(age.pred, 2) > 80] <- NA
    }
  }
  
  return(penet)
}
