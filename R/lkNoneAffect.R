utils::globalVariables(c("lfspenet.cs.nodeath", "new.lfspro.cancer.type"))
lkNoneAffect <- function(penetrance, age){
  # calculate the likelihood for unaffected sample
  # Args:
  #   age: sample's age when he/she still doesn't get disease
  #
  # Returns:
  #   likelihood
  
  length.pene <- nrow(penetrance)
  if(age==0)
  {
    age <- 1
  }
  if(age > length.pene)
  {
    age <- length.pene
  }
  rlt <- rep(1,3)
  rlt[1] <- 1-sum(penetrance[1:age,1])
  rlt[2] <- 1-sum(penetrance[1:age,2])
  rlt[3] <- 1-sum(penetrance[1:age,3])
  #   for(i in 1:age){
  #     rlt[1] <- rlt[1]*(1-penetrance[i,1])
  #     rlt[2] <- rlt[2]*(1-penetrance[i,2])
  #     rlt[3] <- rlt[3]*(1-penetrance[i,3])
  #   }
  rlt
}
