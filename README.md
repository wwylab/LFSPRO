# LFSPRO Package

TP53 germline mutations are the main cause of Li-Fraumeni syndrome. This package is designed to estimate probabilities that:  1) the counselee is a TP53 germline mutation carrier, 2) the counselee develop any cancer in future, 3) the counselee develops breast cancer, sarcoma or any other cancers in future, 4) the counselee develops a second primary cancer in future, on the basis of his/her family cancer history. The package also provides functions for using the LFS classic and Chompret criteria.

## Installation

library(devtools)

install_github("wwylab/LFSPRO")

## Sample codes

The following code briefly illustrates how to use the major functions in the package. Check out the links below for details.

[LFSPRO Github](https://github.com/wwylab/LFSPRO)

[MD Anderson LFSPRO Website](https://bioinformatics.mdanderson.org/public-software/lfspro)

```
fam.id <- c("fam1","fam2","fam2","fam2","fam2")
id <- c(9,1,2,100,200)
counselee.id <- data.frame(fam.id, id)

# LFS classic criteria
lfsClassic(fam.data, cancer.data, counselee.id)
# Chompret criteria
lfsChompret2015(fam.data, cancer.data, counselee.id)

lfspro(fam.data, cancer.data, counselee.id)
```


