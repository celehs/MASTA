library(data.table)

dir <- "data-raw/data_org/"
data_org <- list(
  TrainSurv = fread(paste0(dir, "TrainSurv.csv")),
  ValidSurv = fread(paste0(dir, "ValidSurv.csv")),
  TrainCode = fread(paste0(dir, "TrainCode.csv")),
  ValidCode = fread(paste0(dir, "ValidCode.csv")))

usethis::use_data(data_org, overwrite = TRUE)
