library(data.table)

dir <- "data-raw/data_org/"

data_org <- list(
  TrainSurv = fread(paste0(dir, "TrainSurv.csv")),
  ValidSurv = fread(paste0(dir, "ValidSurv.csv")),
  TrainCode = fread(paste0(dir, "TrainCode.csv")),
  ValidCode = fread(paste0(dir, "ValidCode.csv")))

# colnames(data_org$TrainSurv) <- colnames(data_org$ValidSurv) <- c(
#   "case", "delta", "sx", "sc", "base_pred1", "base_pred2", "base_pred3")
# colnames(data_org$TrainCode)[1:2] <- colnames(data_org$ValidCode)[1:2] <- c("case", "analysisfu")

usethis::use_data(data_org, overwrite = TRUE)
