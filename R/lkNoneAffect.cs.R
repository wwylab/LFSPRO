lkNoneAffect.cs <- 
function (penetrance, age) 
{
  length.pene <- nrow(penetrance)
  if (age == 0) {
    age <- 1
  }
  if (age > length.pene) {
    age <- length.pene
  }
  collapsed.penetrance <- apply(penetrance, c(1,2), sum)
  rlt <- rep(1, 3)
  rlt[1] <- 1 - sum(collapsed.penetrance[1:age, 1])
  rlt[2] <- 1 - sum(collapsed.penetrance[1:age, 2])
  rlt[3] <- 1 - sum(collapsed.penetrance[1:age, 3])
  rlt
}
