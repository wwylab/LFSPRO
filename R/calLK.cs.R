calLK.cs <-
function (fam.cancer.data, penetrance.all) 
{
  invasive.cut <- 99
  num.individual <- length(fam.cancer.data$id)
  lik <- matrix(NA, nrow = num.individual, ncol = 3)
  length.pene <- nrow(penetrance.all$fMX)
  for (i in 1:num.individual) {
    if (fam.cancer.data$num.cancer[i] == 0) # no.cancer
    {
      if (is.na(fam.cancer.data$gender[i])) # gender mssing
      {
        pene.mix <- (penetrance.all$fMX + penetrance.all$fFX)/2
        lik[i, ] <- lkNoneAffect.cs(pene.mix, fam.cancer.data$age[i])
      }
      else { # gender not missing
        if (fam.cancer.data$gender[i] == 1) # male
        {
          lik[i, ] <- lkNoneAffect.cs(penetrance.all$fMX, 
                                      fam.cancer.data$age[i])
        }
        else { # female
          lik[i, ] <- lkNoneAffect.cs(penetrance.all$fFX, 
                                      fam.cancer.data$age[i])
        }
      }
    }
    else { # if cancer diagnosed
      flag.benign <- TRUE
      for (j in 1:fam.cancer.data$num.cancer[i]) {
        if (fam.cancer.data$cancer.info[[i]]$cancer.type[j] < invasive.cut) {
          flag.benign <- FALSE
          break
        }
      }
      if (flag.benign) 
      { # if benign
        if (is.na(fam.cancer.data$gender[i])) { # gender missing
          pene.mix <- (penetrance.all$fMX + penetrance.all$fFX)/2
          lik[i, ] <- lkNoneAffect.cs(pene.mix, fam.cancer.data$age[i])
        }
        else { # gender not missing
          if (fam.cancer.data$gender[i] == 1) { # male
            lik[i, ] <- lkNoneAffect.cs(penetrance.all$fMX, 
                                        fam.cancer.data$age[i])
          }
          else { # female
            lik[i, ] <- lkNoneAffect.cs(penetrance.all$fFX, 
                                        fam.cancer.data$age[i])
          }
        }
      }
      else { # not benign
        # age at-diagnosed
        temp.id <- which.min(fam.cancer.data$cancer.info[[i]]$diag.age)
        age.tmp <- fam.cancer.data$cancer.info[[i]]$diag.age[temp.id]
        if (age.tmp < 1) {
          age.tmp = 1
        }
        if (age.tmp > length.pene) {
          age.tmp = length.pene
        }
        temp.diag.cancer <- fam.cancer.data$cancer.info[[i]]$cancer.type[temp.id]
        
        {
          d <- if (temp.diag.cancer == 4) 1
          else if (temp.diag.cancer == 1) 2 
          else if (temp.diag.cancer == 2) 2
          else 3
        }
        
        if (is.na(fam.cancer.data$gender[i])) 
        { # gender.missing
          lik[i, ] <- (penetrance.all$fMX[age.tmp, , d] + penetrance.all$fFX[age.tmp, , d])/2
        }
        else { # gender.not.missing
          if (fam.cancer.data$gender[i] == 1) { # male
            lik[i, ] <- penetrance.all$fMX[age.tmp, , d]
          }
          else { # female 
            lik[i, ] <- penetrance.all$fFX[age.tmp, , d]
          }
        }
      }
    }
  }
  return(lik)
}
