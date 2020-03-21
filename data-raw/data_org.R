library(data.table)

dir <- "data-raw/data_org/"
TrainSurv <- fread(paste0(dir, "TrainSurv.csv"))
ValidSurv <- fread(paste0(dir, "ValidSurv.csv"))
TrainCode <- fread(paste0(dir, "TrainCode.csv"))
ValidCode <- fread(paste0(dir, "ValidCode.csv"))
TrainN <- fread(paste0(dir, "TrainN.csv"))
ValidN <- fread(paste0(dir, "ValidN.csv"))
list(TrainSurv, ValidSurv)
list(TrainCode, ValidCode)
list(TrainN, ValidN)

colnames(TrainSurv) <- colnames(ValidSurv) <- c(
  "case", "delta", "sx", "sc", "base_pred1", "base_pred2", "base_pred3")
colnames(TrainCode)[1:2] <- colnames(ValidCode)[1:2] <- c("case", "analysisfu")
colnames(TrainN)[1] <- colnames(ValidN)[1] <- "case"
data_org <- list(
  TrainSurv = TrainSurv,
  ValidSurv = ValidSurv,
  TrainCode = TrainCode,
  ValidCode = ValidCode,
  TrainN = TrainN,
  ValidN = ValidN)
str(data_org)

usethis::use_data(data_org, overwrite = TRUE)
