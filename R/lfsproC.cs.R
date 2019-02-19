#fam.cancer.data <- fam.cancer.data[[1]]
#counselee.id <- cid
lfsproC.cs <- 
function (fam.cancer.data, penetrance.all, counselee.id, allef, 
          nloci, mRate) 
{
  lik <- calLK.cs(fam.cancer.data, penetrance.all)
  id <- as.integer(fam.cancer.data$id)
  fid <- as.integer(fam.cancer.data$fid)
  mid <- as.integer(fam.cancer.data$mid)
  counselee.id <- as.integer(counselee.id)
  if (min(id) == 0) {
    id <- id + 1
    fid <- fid + 1
    mid <- mid + 1
    counselee.id <- counselee.id + 1
  }
  fid[is.na(fid)] <- 0
  mid[is.na(mid)] <- 0
  ped <- data.frame(ID = id, Gender = fam.cancer.data$gender, 
                    FatherID = fid, MotherID = mid, stringsAsFactors = FALSE)
  pp <- peelingRC(allef, lik, ped, counselee.id, nloci, mRate)
  return(pp)
}
