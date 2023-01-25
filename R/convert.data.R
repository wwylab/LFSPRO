convert.data <- function(fam.cancer.data) {
  n.family <- length(fam.cancer.data)
  data.obj1 <- data.obj2 <- as.list(1:n.family)
  for (i in 1:n.family) {
    temp <- fam.cancer.data[[i]]
    ID <- temp$id
    Gender <- 2-temp$gender
    Fid <- temp$fid
    Mid <- temp$mid
    data.obj1[[i]] <- data.frame(id = ID, gender = Gender, fid = Fid, mid = Mid)
    
    id <- time <- cancer.type <- gender <- test <- D <- Dp <- NULL
    for (j in 1:length(ID)) {
      cancer.temp <- temp$cancer.info[[j]]
      if (nrow(cancer.temp) == 0) {
        temp.id <- temp$id[j]
        temp.time <- temp$age[j]
        temp.cancer.type <- temp$vital[j]
        temp.gender <- 2-temp$gender[j]
        temp.test <- temp$test[j]
        temp.D <- temp.Dp <- 0
      } else {
        temp.id <- rep(temp$id[j], nrow(cancer.temp) + 1)
        temp.time <- c(cancer.temp$diag.age, temp$age[j])
        temp.cancer.type <- c(cancer.temp$cancer.type, temp$vital[j])
        temp.gender <- rep(2-temp$gender[j], nrow(cancer.temp) + 1)
        temp.test <- rep(temp$test[j], nrow(cancer.temp) + 1)
        temp.D <- c(0, 0, rep(1, nrow(cancer.temp) - 1))
        temp.Dp <- c(0, rep(1, nrow(cancer.temp)))
      }
      id  <- c(id, temp.id)
      time <- c(time, temp.time)
      cancer.type <- c(cancer.type, temp.cancer.type)
      gender <- c(gender, temp.gender)
      test <- c(test, temp.test)
      D <- c(D, temp.D)
      Dp <- c(Dp, temp.Dp)
    }
    data.obj2[[i]] <- data.frame(id, time, cancer.type, gender, test, D, Dp)
  }
  
  obj <- list(data.obj1 = data.obj1, data.obj2 = data.obj2)
}
