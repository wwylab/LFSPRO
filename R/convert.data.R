convert.data <- function(fam.cancer.data)
{
  n.family <- length(fam.cancer.data)
  data.obj1 <- data.obj2 <- as.list(1:n.family)
  for (i in 1:n.family)
  {
    temp <- fam.cancer.data[[i]]
    ID <- temp$id
    Gender <- 2 - temp$gender
    FatherID <- temp$fid
    MotherID <- temp$mid
    data.obj1[[i]] <- data.frame(ID = ID,
                                 Gender = Gender,
                                 FatherID = FatherID, 
                                 MotherID = MotherID)

    id <- time <- gender <- D <- Dp <- NULL
    for (j in 1:length(ID))
    {
      cancer.temp <- temp$cancer.info[[j]]
      if (nrow(cancer.temp) == 0)
      {
        temp.id   <- temp$id[j]
        temp.time <- temp$age[j]
        temp.gender <- temp$gender[j]
        temp.D <- temp.Dp <- 0
      } else {
        temp.id     <- rep(temp$id[j], nrow(cancer.temp) + 1)
        temp.time   <- c(cancer.temp$diag.age, temp$age[j])
        temp.gender <- rep(temp$gender[j], nrow(cancer.temp) + 1)
        temp.D <- c(0, 0, rep(1, nrow(cancer.temp) - 1))
        temp.Dp <- c(0, rep(1, nrow(cancer.temp)))
      }
      id   <- c(id, temp.id)
      time <- c(time, temp.time)
      gender <- c(gender, temp.gender)
      D      <- c(D, temp.D)
      Dp     <- c(Dp, temp.Dp)
    }
    gender <- 2 - gender
    test <- rep(NA, length(gender))
    
    temp.result <- cbind(id, time, gender, test, D, Dp)
    colnames(temp.result) <- c("ID", "time", "gender", "test", "D", "Dp")
    data.obj2[[i]] <- temp.result
    }
  
  obj <- list(data.obj1 = data.obj1, data.obj2 = data.obj2)
  
}
