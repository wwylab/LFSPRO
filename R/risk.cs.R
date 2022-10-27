risk.cs <- function (fam.cancer.data, penetrance.all, counselee.id, posterior){
  
  # cumulative sum of penetrance table of female
  penetrance.sum.F <- array(NA, dim = c(110, 3, 4))
  penetrance.sum.F[,,1] <- matrix(cumsum(penetrance.all$fFX[,,1]), ncol=3)
  penetrance.sum.F[,,2] <- matrix(cumsum(penetrance.all$fFX[,,2]), ncol=3)
  penetrance.sum.F[,,3] <- matrix(cumsum(penetrance.all$fFX[,,3]), ncol=3)
  penetrance.sum.F[,,4] <- matrix(cumsum(penetrance.all$fFX[,,4]), ncol=3)
  
  # cumulative sum of penetrance table of Male
  penetrance.sum.M <- array(NA, dim = c(110, 3, 4))
  penetrance.sum.M[,,1] <- matrix(cumsum(penetrance.all$fMX[,,1]), ncol=3)
  penetrance.sum.M[,,2] <- matrix(cumsum(penetrance.all$fMX[,,2]), ncol=3)
  penetrance.sum.M[,,3] <- matrix(cumsum(penetrance.all$fMX[,,3]), ncol=3)
  penetrance.sum.M[,,4] <- matrix(cumsum(penetrance.all$fMX[,,4]), ncol=3)
  
  num.counselee <- length(counselee.id)
  ages <- fam.cancer.data$age[fam.cancer.data$id %in% counselee.id]
  fam.id <- fam.cancer.data$fam.id
  br.risk.m <- br.risk.f <- br.risk <- matrix(NA, nrow = num.counselee, ncol = 4)          # breast cancer risk prediction table
  sar.risk.m <- sar.risk.f <- sar.risk <- matrix(NA, nrow = num.counselee, ncol = 4)       # sarcoma cancer risk prediction table
  other.risk.m <- other.risk.f <- other.risk <- matrix(NA, nrow = num.counselee, ncol = 4) # other cancers risk prediction table
  d.risk.m <- d.risk.f <- d.risk <- matrix(NA, nrow = num.counselee, ncol = 4)             # death risk prediction table
  c.risk.m <- c.risk.f <- matrix(NA, nrow = 5, ncol = 4)                                   # risks at (current age, age + 5, age + 10, age + 15, age + 20)
  
  colnames(br.risk) <- c("breast in 5 yrs", "breast in 10 yrs", "breast in 15 yrs", "breast in 20 yrs")
  colnames(sar.risk) <- c("sarcoma in 5 yrs", "sarcoma in 10 yrs", "sarcoma in 15 yrs", "sarcoma in 20 yrs")
  colnames(other.risk) <- c("others in 5 yrs", "others in 10 yrs", "others in 15 yrs", "others in 20 yrs")
  colnames(d.risk) <- c("death in 5 yrs", "death in 10 yrs", "death in 15 yrs", "death in 20 yrs")
  
  for (i in 1:num.counselee){
    
    id <- counselee.id[i]
    index <- which(fam.cancer.data$id == id)
    age <- fam.cancer.data$age[index]
    if (age > 80) age <- 80
    age.pred <- c(age, age + 5, age + 10, age + 15, age + 20)
    
    if (fam.cancer.data$num.cancer[index] > 0) {
      
      warning(paste0("ID = ", id, ": No risk predictions available for second primary cancer!"))
      ### simply make age NA so the risk output is NA
      age.pred <- NA
      
    } else {
      
      # calculate the risk at a given age
      # risk = cdf[age][1] * carrier probability of genotype1 + 
      #        cdf[age][2] * carrier probability of genotype2 + 
      #        cdf[age][3] * carrier probability of genotype3
      
      for (c in 1:4) {
        
        c.risk.m[,c] <- penetrance.sum.M[,,c][age.pred,][,1] * posterior[i,][1] +
          penetrance.sum.M[,,c][age.pred,][,2] * posterior[i,][2] +
          penetrance.sum.M[,,c][age.pred,][,3] * posterior[i,][3]
        c.risk.f[,c] <- penetrance.sum.F[,,c][age.pred,][,1] * posterior[i,][1] +
          penetrance.sum.F[,,c][age.pred,][,2] * posterior[i,][2] +
          penetrance.sum.F[,,c][age.pred,][,3] * posterior[i,][3]
        
      }
      
      # Calculate the risks of each competing risk
      # Example of breast cancer risk prediction: 
      
      # A = {no cancer, no death by age}
      # Pr(A) = (1 - Pr(br < age)) * (1 - Pr(sar < age)) * (1 - Pr(others < age)) * (1 - Pr(d < age))
      
      # Risk of breast cancer in n years and no other cancers or death = 
      #   Pr(br =< age + n, sar > age + n, others > age + n, d > age + n | A)
      #   Pr(age < br =< age + n, sar > age + n, others > age + n, d > age + n) / P(A)
      #   Pr(age < br =< age + n, br is the first cancer) / P(A)
      #   (Pr(br =< age + n, br is the first cancer) - Pr(br < age, br is the first cancer)) / P(A)
      #   (br_penet(age + n) - br_penet(age)) / P(A)
      
      no.cancer.m <- (1 - c.risk.m[1,1]) * (1 - c.risk.m[1,2]) * (1 - c.risk.m[1,3]) * (1 - c.risk.m[1,4])
      no.cancer.f <- (1 - c.risk.f[1,1]) * (1 - c.risk.f[1,2]) * (1 - c.risk.f[1,3]) * (1 - c.risk.f[1,4])
      
      br.risk.m[i,] <- (c.risk.m[2:5,1] - c.risk.m[1,1]) / no.cancer.m
      sar.risk.m[i,] <- (c.risk.m[2:5,2] - c.risk.m[1,2]) / no.cancer.m
      other.risk.m[i,] <- (c.risk.m[2:5,3] - c.risk.m[1,3]) / no.cancer.m
      d.risk.m[i,] <- (c.risk.m[2:5,4] - c.risk.m[1,4]) / no.cancer.m
      
      br.risk.f[i,] <- (c.risk.f[2:5,1] - c.risk.f[1,1]) / no.cancer.f
      sar.risk.f[i,] <- (c.risk.f[2:5,2] - c.risk.f[1,2]) / no.cancer.f
      other.risk.f[i,] <- (c.risk.f[2:5,3] - c.risk.f[1,3]) / no.cancer.f
      d.risk.f[i,] <- (c.risk.f[2:5,4] - c.risk.f[1,4]) / no.cancer.f
      
      if (is.na(fam.cancer.data$gender[index])) { # gender is missing
        
        br.risk[i,] <- (br.risk.m[i,] + br.risk.f[i,]) / 2
        sar.risk[i,] <- (sar.risk.m[i,] + sar.risk.f[i,]) / 2
        other.risk[i,] <- (other.risk.m[i,] + other.risk.f[i,]) / 2
        d.risk[i,] <- (d.risk.m[i,] + d.risk.f[i,]) / 2
        
      } else if (fam.cancer.data$gender[index] == 0) { #female
        
        br.risk[i,] <- br.risk.f[i,]
        sar.risk[i,] <- sar.risk.f[i,]
        other.risk[i,] <- other.risk.f[i,]
        d.risk[i,] <- d.risk.f[i,]
        
      } else { # male
        
        br.risk[i,] <- br.risk.m[i,]
        sar.risk[i,] <- sar.risk.m[i,]
        other.risk[i,] <- other.risk.m[i,]
        d.risk[i,] <- d.risk.m[i,]
        
      }
      
    }
    
    ## Invalid to predict cancer risks for patients with predicted ages beyond 80
    pred_age_mat <- matrix(rep(ages, 4), ncol = 4) + t(matrix(rep(seq(5, 20, 5), length(ages)), nrow = 4))
    br.risk[which(pred_age_mat > 80)] <- NA
    sar.risk[which(pred_age_mat > 80)] <- NA
    other.risk[which(pred_age_mat > 80)] <- NA
    d.risk[which(pred_age_mat > 80)] <- NA
    
  }
  
  breast <- cbind.data.frame(fam.id = unique(fam.id), counselee.id, ages, br.risk)
  sarcoma <- cbind.data.frame(fam.id = unique(fam.id), counselee.id, ages, sar.risk)
  others <- cbind.data.frame(fam.id = unique(fam.id), counselee.id, ages, other.risk)
  death <- cbind.data.frame(fam.id = unique(fam.id), counselee.id, ages, d.risk)
  
  output <- list(Breast_risks = breast, 
                 Sarcoma_risks = sarcoma, 
                 Other_cancers = others, 
                 Death = death)
  
  return(output)
  
}


