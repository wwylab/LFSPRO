calLK.cs <- function (ind.data, penetrance.all, mut.info = TRUE) {
  invasive.cut <- 99
  id <- unique(ind.data$id)
  num.individual <- length(id)
  lik <- matrix(NA, nrow = num.individual, ncol = 3)
  length.pene <- nrow(penetrance.all$fMX)
  for (i in 1:num.individual) {
    ind.data.i <- ind.data[ind.data$id == id[i],]
    num.cancer <- nrow(ind.data.i)-1
    gender <- ind.data.i$gender[1]
    age <- ind.data.i$time[num.cancer+1]
    if (num.cancer == 0) {
      if (is.na(gender)) {
        pene.mix <- (penetrance.all$fMX + penetrance.all$fFX)/2
        lik[i, ] <- lkNoneAffect.cs(pene.mix, age)
      } else {
        if (gender == 1) {
          lik[i, ] <- lkNoneAffect.cs(penetrance.all$fMX, age)
        } else {
          lik[i, ] <- lkNoneAffect.cs(penetrance.all$fFX, age)
        }
      }
    } else { 
      flag.benign <- TRUE
      for (j in 1:num.cancer) {
        cancer.type <- lfspro.cancer.type[ind.data.i$cancer.type[j]]
        if (cancer.type < invasive.cut) {
          flag.benign <- FALSE
          break
        }
      }
      if (flag.benign) {
        if (is.na(gender)) { 
          pene.mix <- (penetrance.all$fMX + penetrance.all$fFX)/2
          lik[i, ] <- lkNoneAffect.cs(pene.mix, age)
        } else { 
          if (gender == 1) {
            lik[i, ] <- lkNoneAffect.cs(penetrance.all$fMX, age)
          } else {
            lik[i, ] <- lkNoneAffect.cs(penetrance.all$fFX, age)
          }
        }
      } else {
        temp.id <- which.min(ind.data.i$time[1:num.cancer])
        age.tmp <- ind.data.i$time[temp.id]
        if (age.tmp < 1) {
          age.tmp <- 1
        }
        if (age.tmp > length.pene) {
          age.tmp <- length.pene
        }
        temp.diag.cancer <- lfspro.cancer.type[ind.data.i$cancer.type[temp.id]]
        {
          d <- if (temp.diag.cancer == 4) 1
          else if (temp.diag.cancer == 1) 2 
          else if (temp.diag.cancer == 2) 2
          else 3
        }
        if (is.na(gender)) { 
          lik[i, ] <- (penetrance.all$fMX[age.tmp, , d] + penetrance.all$fFX[age.tmp, , d])/2
        } else { 
          if (gender == 1) {
            lik[i, ] <- penetrance.all$fMX[age.tmp, , d]
          } else {
            lik[i, ] <- penetrance.all$fFX[age.tmp, , d]
          }
        }
      }
    }
    if (mut.info == TRUE) {
      mut <- ind.data.i$test[1]
      if (!is.na(mut)) {
        if (mut == 0) {
          lik[i,2:3] <- 0 
        } else if (mut == 1) {
          lik[i,1] <- 0
        }
      }
    }
  }
  return(lik)
}
