lfsproC.cs <- function (ped.data, ind.data, counselee.id, penetrance.all, 
                        allef, nloci, mRate, mut.info = TRUE) {
  lik <- calLK.cs(ind.data, penetrance.all, mut.info)
  id <- as.integer(ped.data$id)
  fid <- as.integer(ped.data$fid)
  mid <- as.integer(ped.data$mid)
  counselee.id <- as.integer(counselee.id)
  if (min(id) == 0) {
    id <- id + 1
    fid <- fid + 1
    mid <- mid + 1
    counselee.id <- counselee.id + 1
  }
  fid[is.na(fid)] <- 0
  mid[is.na(mid)] <- 0
  ped <- data.frame(ID = id, Gender = 2-ped.data$gender, 
                    FatherID = fid, MotherID = mid, stringsAsFactors = FALSE)
  pp <- peelingRC(allef, lik, ped, counselee.id, nloci, mRate)
  return(pp)
}
