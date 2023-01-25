lfsproC.cs <- function (data1, data2, counselee.id, penetrance.all, 
                        allef, nloci, mRate, mut.info = TRUE) {
  lik <- calLK.cs(data1, data2, penetrance.all, mut.info)
  id <- as.integer(data1$id)
  fid <- as.integer(data1$fid)
  mid <- as.integer(data1$mid)
  counselee.id <- as.integer(counselee.id)
  if (min(id) == 0) {
    id <- id + 1
    fid <- fid + 1
    mid <- mid + 1
    counselee.id <- counselee.id + 1
  }
  fid[is.na(fid)] <- 0
  mid[is.na(mid)] <- 0
  ped <- data.frame(ID = id, Gender = 2-data1$gender, 
                    FatherID = fid, MotherID = mid, stringsAsFactors = FALSE)
  pp <- peelingRC(allef, lik, ped, counselee.id, nloci, mRate)
  return(pp)
}
