library(data.table)

data_org <- list(
  TrainSurv = fread("data-raw/data_org/TrainSurv.csv"),
  ValidSurv = fread("data-raw/data_org/ValidSurv.csv"),
  TrainCode = fread("data-raw/data_org/TrainCode.csv"),
  ValidCode = fread("data-raw/data_org/ValidCode.csv"),
  TrainN = fread("data-raw/data_org/TrainN.csv"),
  ValidN = fread("data-raw/data_org/ValidN.csv"))

usethis::use_data(data_org, overwrite = TRUE)
