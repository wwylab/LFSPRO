combined.risk.mpc <- function(posterior, risk, counselee.id){
  # risk in 5 years
  risk.5 <- apply(as.matrix(cbind(risk[4], risk[5], risk[5])), 2, as.numeric)
  overall.5 <- risk.5 %*% t(posterior)
  
  # risk in 10 years
  risk.10 <- apply(as.matrix(cbind(risk[6], risk[7], risk[7])), 2, as.numeric)
  overall.10 <- risk.10 %*% t(posterior)
  
  # risk in 15 years
  risk.15 <- apply(as.matrix(cbind(risk[8], risk[9], risk[9])), 2, as.numeric)
  overall.15 <- risk.15 %*% t(posterior)
  
  risk <- data.frame(counselee.id, diag(overall.5), diag(overall.10), diag(overall.15))
  colnames(risk) <- c("fam.id", "id","5 years", "10 years", "15 years")
  return(risk)
}
