
devtools::install_github('DenisRustand/INLAjoint')



library(INLAjoint)
library(JM) # This package contains the dataset

data(pbc2) # dataset
pbc2 <- pbc2[which(pbc2$id %in% c(1:20)),]
# extract some variable of interest without missing values
Longi <- na.omit(pbc2[, c("id", "years", "status","drug","age", 
                          "sex","year","serBilir","SGOT", "albumin", "edema",
                          "platelets", "alkaline","spiders", "ascites")])
Surv <- Longi[c(which(diff(as.numeric(Longi[,which(colnames(Longi)=="id")]))==1),
                length(Longi[,which(colnames(Longi)=="id")])),-c(7:10, 12:16)]
Surv$death <- ifelse(Surv$status=="dead",1,0) # competing event 1
Surv$trans <- ifelse(Surv$status=="transplanted",1,0) # competing event 2

summary(Longi)
summary(Surv)







