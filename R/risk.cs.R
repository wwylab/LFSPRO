risk.cs <- function (fam.cancer.data, penetrance.all, counselee.id, posterior){
    # cumulative sum of penetrance table of female
    penetrance.sum.F <- array(NA, dim = c(110, 3,4))
    penetrance.sum.F[,,1] <- matrix(cumsum(penetrance.all$fFX[,,1]), ncol=3)
    penetrance.sum.F[,,2] <- matrix(cumsum(penetrance.all$fFX[,,2]), ncol=3)
    penetrance.sum.F[,,3] <- matrix(cumsum(penetrance.all$fFX[,,3]), ncol=3)
    penetrance.sum.F[,,4] <- matrix(cumsum(penetrance.all$fFX[,,4]), ncol=3)
    # cumulative sum of penetrance table of Male
    penetrance.sum.M <- array(NA, dim = c(110, 3,4))
    penetrance.sum.M[,,1] <- matrix(cumsum(penetrance.all$fMX[,,1]), ncol=3)
    penetrance.sum.M[,,2] <- matrix(cumsum(penetrance.all$fMX[,,2]), ncol=3)
    penetrance.sum.M[,,3] <- matrix(cumsum(penetrance.all$fMX[,,3]), ncol=3)
    penetrance.sum.M[,,4] <- matrix(cumsum(penetrance.all$fMX[,,4]), ncol=3)
    
    num.counselee <- length(counselee.id)
    ages <- fam.cancer.data$age[fam.cancer.data$id %in% counselee.id]
    fam.id <- fam.cancer.data$fam.id
    br.risk.m <- br.risk.f <- br.risk <- matrix(NA, nrow = num.counselee, ncol = 4)    # breast cancer risk prediction table
    sar.risk.m <- sar.risk.f <- sar.risk <- matrix(NA, nrow = num.counselee, ncol = 4)   # sarcoma cancer risk prediction table
    other.risk.m <- other.risk.f <- other.risk <- matrix(NA, nrow = num.counselee, ncol = 4) # other cancers risk prediction table
    d.risk.m <- d.risk.f <- d.risk <- matrix(NA, nrow = num.counselee, ncol = 4)     # death risk prediction table
    c.risk.m <- c.risk.f <- c.risk <- matrix(NA, nrow = num.counselee, ncol = 4)     # risks at current age
    
    colnames(br.risk) <- c("breast in 5 yrs", "breast in 10 yrs", "breast in 15 yrs", "breast in 20 yrs")
    colnames(sar.risk) <- c("sarcoma in 5 yrs", "sarcoma in 10 yrs", "sarcoma in 15 yrs", "sarcoma in 20 yrs")
    colnames(other.risk) <- c("others in 5 yrs", "others in 10 yrs", "others in 15 yrs", "others in 20 yrs")
    colnames(d.risk) <- c("death in 5 yrs", "death in 10 yrs", "death in 15 yrs", "death in 20 yrs")
    colnames(c.risk) <- c("breast", "sarcoma", "others", "death")
    # calculate the risks of current age
    # current risk = 
    # cdf[age][1]*carrier probability of genotype1 + cdf[age][2]*carrier probability of genotype2 + cdf[age][2]*carrier probability of genotype3
    for (i in 1:num.counselee){
      id <- counselee.id[i]
      index <- which(fam.cancer.data$id == id)
      age <- fam.cancer.data$age[index]
      
      if (is.na(fam.cancer.data$gender[i])){ # gender is missing
        
        c.risk.m[i,1] <- penetrance.sum.M[,,1][age,][1]*posterior[i,][1] +
          penetrance.sum.M[,,1][age,][2]*posterior[i,][2] +
          penetrance.sum.M[,,1][age,][3]*posterior[i,][3]   
        c.risk.m[i,2] <- penetrance.sum.M[,,2][age,][1]*posterior[i,][1] +
          penetrance.sum.M[,,2][age,][2]*posterior[i,][2] +
          penetrance.sum.M[,,2][age,][3]*posterior[i,][3] 
        c.risk.m[i,3] <- penetrance.sum.M[,,3][age,][1]*posterior[i,][1] +
          penetrance.sum.M[,,3][age,][2]*posterior[i,][2] +
          penetrance.sum.M[,,3][age,][3]*posterior[i,][3] 
        c.risk.m[i,4] <- penetrance.sum.M[,,4][age,][1]*posterior[i,][1] +
          penetrance.sum.M[,,4][age,][2]*posterior[i,][2] +
          penetrance.sum.M[,,4][age,][3]*posterior[i,][3] 
        
        c.risk.f[i,1] <- penetrance.sum.F[,,1][age,][1]*posterior[i,][1] +
          penetrance.sum.F[,,1][age,][2]*posterior[i,][2] +
          penetrance.sum.F[,,1][age,][3]*posterior[i,][3] 
        c.risk.f[i,2] <- penetrance.sum.F[,,2][age,][1]*posterior[i,][1] +
          penetrance.sum.F[,,2][age,][2]*posterior[i,][2] +
          penetrance.sum.F[,,2][age,][3]*posterior[i,][3] 
        c.risk.f[i,3] <- penetrance.sum.F[,,3][age,][1]*posterior[i,][1] +
          penetrance.sum.F[,,3][age,][2]*posterior[i,][2] +
          penetrance.sum.F[,,3][age,][3]*posterior[i,][3] 
        c.risk.f[i,4] <- penetrance.sum.F[,,4][age,][1]*posterior[i,][1] +
          penetrance.sum.F[,,4][age,][2]*posterior[i,][2] +
          penetrance.sum.F[,,4][age,][3]*posterior[i,][3] 
        
        c.risk[i,1] <- (c.risk.m[i,1] + c.risk.f[i,1])/2
        c.risk[i,2] <- (c.risk.m[i,2] + c.risk.f[i,2])/2
        c.risk[i,3] <- (c.risk.m[i,3] + c.risk.f[i,3])/2
        c.risk[i,4] <- (c.risk.m[i,4] + c.risk.f[i,4])/2
       
      }
        else if (fam.cancer.data$gender[i] == 0) { #female
        c.risk[i,1] <- penetrance.sum.F[,,1][age,][1]*posterior[i,][1] +
          penetrance.sum.F[,,1][age,][2]*posterior[i,][2] +
          penetrance.sum.F[,,1][age,][3]*posterior[i,][3] 
        c.risk[i,2] <- penetrance.sum.F[,,2][age,][1]*posterior[i,][1] +
          penetrance.sum.F[,,2][age,][2]*posterior[i,][2] +
          penetrance.sum.F[,,2][age,][3]*posterior[i,][3] 
        c.risk[i,3] <- penetrance.sum.F[,,3][age,][1]*posterior[i,][1] +
          penetrance.sum.F[,,3][age,][2]*posterior[i,][2] +
          penetrance.sum.F[,,3][age,][3]*posterior[i,][3] 
        c.risk[i,4] <- penetrance.sum.F[,,4][age,][1]*posterior[i,][1] +
          penetrance.sum.F[,,4][age,][2]*posterior[i,][2] +
          penetrance.sum.F[,,4][age,][3]*posterior[i,][3] 
        }
        else { # male
          c.risk[i,1] <- penetrance.sum.M[,,1][age,][1]*posterior[i,][1] +
            penetrance.sum.M[,,1][age,][2]*posterior[i,][2] +
            penetrance.sum.M[,,1][age,][3]*posterior[i,][3]   
          c.risk[i,2] <- penetrance.sum.M[,,2][age,][1]*posterior[i,][1] +
            penetrance.sum.M[,,2][age,][2]*posterior[i,][2] +
            penetrance.sum.M[,,2][age,][3]*posterior[i,][3] 
          c.risk[i,3] <- penetrance.sum.M[,,3][age,][1]*posterior[i,][1] +
            penetrance.sum.M[,,3][age,][2]*posterior[i,][2] +
            penetrance.sum.M[,,3][age,][3]*posterior[i,][3] 
          c.risk[i,4] <- penetrance.sum.M[,,4][age,][1]*posterior[i,][1] +
            penetrance.sum.M[,,4][age,][2]*posterior[i,][2] +
            penetrance.sum.M[,,4][age,][3]*posterior[i,][3] 
        }
      }
      
      # calculate the risks of each competing risk
      for (i in 1:num.counselee){
        id <- counselee.id[i]
        index <- which(fam.cancer.data$id == id)
        age <- fam.cancer.data$age[index]
        if (age>80) age <- 80
        age.pred <- c(age+5, age+10, age+15, age+20)
        
        if (fam.cancer.data$num.cancer[index]>0) {
          warning("No risk predictions available for second primary cancer!")
          ### simply make age NA so the risk output is NA
          age.pred <- NA
        }
        # Example of breast cancer risk prediction: 
        # risk of breast cancer in n years = 
        #   (cdf.breast[age+n]-cdf.breast[age])*(1-cdf.sarcoma[age])*(1-cdf.others[age])*(1-cdf.death[age]) 
        
        if (is.na(fam.cancer.data$gender[i])){ # gender is missing
          for (j in 1:4){
            br.risk.m[i,j] <- ((penetrance.sum.M[,,1][age.pred[j],][1]*posterior[i,][1] + 
                                  penetrance.sum.M[,,1][age.pred[j],][2]*posterior[i,][2] +
                                  penetrance.sum.M[,,1][age.pred[j],][3]*posterior[i,][3]) -
                                 c.risk[i,1])*(1-c.risk[i,2])*(1-c.risk[i,3])*(1-c.risk[i,4])
            sar.risk.m[i,j] <- ((penetrance.sum.M[,,2][age.pred[j],][1]*posterior[i,][1] +
                                   penetrance.sum.M[,,2][age.pred[j],][2]*posterior[i,][2] +
                                   penetrance.sum.M[,,2][age.pred[j],][3]*posterior[i,][3]) -
                                  c.risk[i,2])*(1-c.risk[i,1])*(1-c.risk[i,3])*(1-c.risk[i,4])  
            other.risk.m[i,j] <- ((penetrance.sum.M[,,3][age.pred[j],][1]*posterior[i,][1] +
                                     penetrance.sum.M[,,3][age.pred[j],][2]*posterior[i,][2] +
                                     penetrance.sum.M[,,3][age.pred[j],][3]*posterior[i,][3])-
                                    c.risk[i,3])*(1-c.risk[i,1])*(1-c.risk[i,2])*(1-c.risk[i,4]) 
            d.risk.m[i,j] <- ((penetrance.sum.M[,,4][age.pred[j],][1]*posterior[i,][1] +
                                 penetrance.sum.M[,,4][age.pred[j],][2]*posterior[i,][2] +
                                 penetrance.sum.M[,,4][age.pred[j],][3]*posterior[i,][3]) -
                                c.risk[i,4])*(1-c.risk[i,1])*(1-c.risk[i,2])*(1-c.risk[i,3])
            
            br.risk.f[i,j] <- ((penetrance.sum.F[,,1][age.pred[j],][1]*posterior[i,][1] +
                                  penetrance.sum.F[,,1][age.pred[j],][2]*posterior[i,][2] +
                                  penetrance.sum.F[,,1][age.pred[j],][3]*posterior[i,][3])-
                                 c.risk[i,1])*(1-c.risk[i,2])*(1-c.risk[i,3])*(1-c.risk[i,4])
            sar.risk.f[i,j] <- ((penetrance.sum.F[,,2][age.pred[j],][1]*posterior[i,][1] +
                                   penetrance.sum.F[,,2][age.pred[j],][2]*posterior[i,][2] +
                                   penetrance.sum.F[,,2][age.pred[j],][3]*posterior[i,][3]) -
                                  c.risk[i,2])*(1-c.risk[i,1])*(1-c.risk[i,3])*(1-c.risk[i,4])  
            other.risk.f[i,j] <- ((penetrance.sum.F[,,3][age.pred[j],][1]*posterior[i,][1] +
                                     penetrance.sum.F[,,3][age.pred[j],][2]*posterior[i,][2] +
                                     penetrance.sum.F[,,3][age.pred[j],][3]*posterior[i,][3]) -
                                    c.risk[i,3])*(1-c.risk[i,1])*(1-c.risk[i,2])*(1-c.risk[i,4]) 
            d.risk.f[i,j] <- ((penetrance.sum.F[,,4][age.pred[j],][1]*posterior[i,][1] +
                                 penetrance.sum.F[,,4][age.pred[j],][2]*posterior[i,][2] +
                                 penetrance.sum.F[,,4][age.pred[j],][3]*posterior[i,][3]) -
                                c.risk[i,4])*(1-c.risk[i,1])*(1-c.risk[i,2])*(1-c.risk[i,3]) 
            
            br.risk[i,j] <- (br.risk.m[i,j] + br.risk.f[i,j])/2
            sar.risk[i,j] <- (sar.risk.m[i,j] + sar.risk.f[i,j])/2
            other.risk[i,j] <- (other.risk.m[i,j] + other.risk.f[i,j])/2
            d.risk[i,j] <- (d.risk.m[i,j] + d.risk.f[i,j])/2
          }
        }
        else if (fam.cancer.data$gender[i] == 0){ #female
          for (j in 1:4){
            br.risk[i,j] <- ((penetrance.sum.F[,,1][age.pred[j],][1]*posterior[i,][1] +
                                penetrance.sum.F[,,1][age.pred[j],][2]*posterior[i,][2] +
                                penetrance.sum.F[,,1][age.pred[j],][3]*posterior[i,][3])-
                               c.risk[i,1])*(1-c.risk[i,2])*(1-c.risk[i,3])*(1-c.risk[i,4])
            sar.risk[i,j] <- ((penetrance.sum.F[,,2][age.pred[j],][1]*posterior[i,][1] +
                                penetrance.sum.F[,,2][age.pred[j],][2]*posterior[i,][2] +
                                penetrance.sum.F[,,2][age.pred[j],][3]*posterior[i,][3]) -
                                c.risk[i,2])*(1-c.risk[i,1])*(1-c.risk[i,3])*(1-c.risk[i,4])  
            other.risk[i,j] <- ((penetrance.sum.F[,,3][age.pred[j],][1]*posterior[i,][1] +
                                penetrance.sum.F[,,3][age.pred[j],][2]*posterior[i,][2] +
                                penetrance.sum.F[,,3][age.pred[j],][3]*posterior[i,][3]) -
                                c.risk[i,3])*(1-c.risk[i,1])*(1-c.risk[i,2])*(1-c.risk[i,4]) 
            d.risk[i,j] <- ((penetrance.sum.F[,,4][age.pred[j],][1]*posterior[i,][1] +
                                penetrance.sum.F[,,4][age.pred[j],][2]*posterior[i,][2] +
                                penetrance.sum.F[,,4][age.pred[j],][3]*posterior[i,][3]) -
                                c.risk[i,4])*(1-c.risk[i,1])*(1-c.risk[i,2])*(1-c.risk[i,3])  
          }
        } else { # male
          for (j in 1:4){
            br.risk[i,j] <- ((penetrance.sum.M[,,1][age.pred[j],][1]*posterior[i,][1] + 
                                penetrance.sum.M[,,1][age.pred[j],][2]*posterior[i,][2] +
                                penetrance.sum.M[,,1][age.pred[j],][3]*posterior[i,][3]) -
                               c.risk[i,1])*(1-c.risk[i,2])*(1-c.risk[i,3])*(1-c.risk[i,4])
            sar.risk[i,j] <- ((penetrance.sum.M[,,2][age.pred[j],][1]*posterior[i,][1] +
                                 penetrance.sum.M[,,2][age.pred[j],][2]*posterior[i,][2] +
                                 penetrance.sum.M[,,2][age.pred[j],][3]*posterior[i,][3]) -
                                c.risk[i,2])*(1-c.risk[i,1])*(1-c.risk[i,3])*(1-c.risk[i,4])  
            other.risk[i,j] <- ((penetrance.sum.M[,,3][age.pred[j],][1]*posterior[i,][1] +
                                   penetrance.sum.M[,,3][age.pred[j],][2]*posterior[i,][2] +
                                   penetrance.sum.M[,,3][age.pred[j],][3]*posterior[i,][3])-
                                  c.risk[i,3])*(1-c.risk[i,1])*(1-c.risk[i,2])*(1-c.risk[i,4]) 
            d.risk[i,j] <- ((penetrance.sum.M[,,4][age.pred[j],][1]*posterior[i,][1] +
                               penetrance.sum.M[,,4][age.pred[j],][2]*posterior[i,][2] +
                               penetrance.sum.M[,,4][age.pred[j],][3]*posterior[i,][3]) -
                              c.risk[i,4])*(1-c.risk[i,1])*(1-c.risk[i,2])*(1-c.risk[i,3])
          }
        }
      }
    ## Invalid to predict cancer risks for patients with predicted ages beyond 80
    pred_age_mat <- matrix(rep(ages,4),ncol = 4) + t(matrix(rep(seq(5,20,5),length(ages)), nrow = 4))
    br.risk[which(pred_age_mat > 80)] <- NA
    sar.risk[which(pred_age_mat > 80)] <- NA
    other.risk[which(pred_age_mat > 80)] <- NA
    #breast <- data.frame(cbind(fam.id=unique(fam.id), counselee.id, ages, br.risk), check.names = FALSE, stringsAsFactors = F)
    breast <- cbind.data.frame(fam.id=unique(fam.id), counselee.id, ages, br.risk)
    #sarcoma <- data.frame(cbind(fam.id=unique(fam.id), counselee.id, ages, sar.risk), check.names = FALSE, stringsAsFactors = F)
    sarcoma <- cbind.data.frame(fam.id=unique(fam.id), counselee.id, ages, sar.risk)
    #others <- data.frame(cbind(fam.id=unique(fam.id), counselee.id, ages, other.risk), check.names = FALSE, stringsAsFactors = F)
    others <- cbind.data.frame(fam.id=unique(fam.id), counselee.id, ages, other.risk)
    output <- list(Breast_risks=breast, 
                   Sarcoma_risks=sarcoma, 
                   Other_cancers=others)
    return(output)
  }
    
    
