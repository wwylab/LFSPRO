calLK.mpc <- function(data1, data2, parameter)
{
  beta <- parameter$beta
  gamma <- parameter$gamma
  
  M <- length(gamma)
  
  time <- data2[,2]
  time[time > 80 - 1.0e-8] <- 80
  id <- data2[,1]
  time[time == 0] <- 1.0e-12
  
  tilde.t <- time/80
  gender <- data2[,3]
  d      <- data2[,5]
  dp     <- data2[,6]
  
  n2 <- length(time)
  
  # Bernstein Basis
  W <- unlist(lapply(1:M, function(k) pbeta(tilde.t, k, M-k+1)))
  w <- unlist(lapply(1:M, function(k) dbeta(tilde.t, k, M-k+1)))    
  
  Ft <- matrix(W, ncol = M)
  ft <- matrix(w, ncol = M)
  
  # calculate likelihood
  test0 <- rep(0, n2)
  test1 <- rep(1, n2)
  
  
  xp.test0.m <- cbind(test0, rep(1, n2), rep(1, n2) * test0, dp, dp * test0)
  xp.test1.m <- cbind(test1, rep(1, n2), rep(1, n2) * test1, dp, dp * test1)
  
  xp.test0.f <- cbind(test0, rep(0, n2), rep(0, n2) * test0, dp, dp * test0)
  xp.test1.f <- cbind(test1, rep(0, n2), rep(0, n2) * test1, dp, dp * test1)
  
  xpbeta.test0.m <- xp.test0.m %*% beta
  xpbeta.test1.m <- xp.test1.m %*% beta
  
  xpbeta.test0.f <- xp.test0.f %*% beta
  xpbeta.test1.f <- xp.test1.f %*% beta
  
  Lambda <- (Ft %*% gamma)         # cumulative baseline
  lambda <- (ft %*% gamma)/80# baseline
  
  unique.id <- unique(id)
  diff.Lambda <- NULL
  for (ii in unique.id) {
    diff.Lambda <- c(diff.Lambda, diff(c(0, Lambda[id == ii])))
  }
  
  log.l0.m <- log(lambda) + xpbeta.test0.m
  log.l1.m <- log(lambda) + xpbeta.test1.m
  
  log.l0.f <- log(lambda) + xpbeta.test0.f
  log.l1.f <- log(lambda) + xpbeta.test1.f
  
  log.L0.m <- diff.Lambda * exp(xpbeta.test0.m)
  log.L1.m <- diff.Lambda * exp(xpbeta.test1.m)
  
  log.L0.f <- diff.Lambda * exp(xpbeta.test0.f)
  log.L1.f <- diff.Lambda * exp(xpbeta.test1.f)
  
  
  llike0 <- llike1 <- NULL
  for (ii in unique.id) {
    temp.id <- which(id == ii)
    
    if (any(is.na(gender[temp.id]))) {
      temp.llike0.m <- sum(log.l0.m[temp.id]) - log.l0.m[temp.id[length(temp.id)]] - sum(log.L0.m[temp.id])
      temp.llike1.m <- sum(log.l1.m[temp.id]) - log.l1.m[temp.id[length(temp.id)]] - sum(log.L1.m[temp.id])
      
      temp.llike0.f <- sum(log.l0.f[temp.id]) - log.l0.f[temp.id[length(temp.id)]] - sum(log.L0.f[temp.id])
      temp.llike1.f <- sum(log.l1.f[temp.id]) - log.l1.f[temp.id[length(temp.id)]] - sum(log.L1.f[temp.id])
      
      temp.llike0 <- (temp.llike0.m + temp.llike0.f)/2
      temp.llike1 <- (temp.llike1.m + temp.llike1.f)/2
    } else if (all(gender[temp.id] == 1))   {
      temp.llike0 <- sum(log.l0.m[temp.id]) - log.l0.m[temp.id[length(temp.id)]] - sum(log.L0.m[temp.id])
      temp.llike1 <- sum(log.l1.m[temp.id]) - log.l1.m[temp.id[length(temp.id)]] - sum(log.L1.m[temp.id])  
    } else if (all(gender[temp.id] == 2)) {
      temp.llike0 <- sum(log.l0.f[temp.id]) - log.l0.f[temp.id[length(temp.id)]] - sum(log.L0.f[temp.id])
      temp.llike1 <- sum(log.l1.f[temp.id]) - log.l1.f[temp.id[length(temp.id)]] - sum(log.L1.f[temp.id])
    } else stop("check gender code: 1 for female 2 for male!")
    
    llike0 <- c(llike0, temp.llike0)
    llike1 <- c(llike1, temp.llike1)
  }
  
  lik <- cbind(exp(llike0), exp(llike1), exp(llike1))
  
  return(lik)
}
