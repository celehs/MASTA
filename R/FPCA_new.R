#' @title Functional Principal Component Analysis (FPCA)
#' @description Performs FPCA to extract features from longitudinal encounter data.
#' @param time longitudinal encounter times. They should be greater than or equal to 1.
#' @param fu_train follow-up time (training)
#' @param fu_valid follow-up time (validation)
#' @param K.select characters indicating which method to choose the number of principal components K.
#'  Default is K.select="PropVar", and K.select="PPIC" is also available.
#' @param n.grid an integer value for grid points used in estimating covariance function g. Default is \code{401}.
#' @param propvar a proportion of variation used to select number of FPCs. Default is \code{0.85}.
#' @export
fpca.new <- function(time, fu_train, fu_valid,
                     K.select = "PropVar", n.grid = 401, propvar = 0.85) {
  
  #--- check input ---
  fpca.check(time, fu_train, fu_valid)
  
  #--- data preparation ---
  data <- fpca.pre(time, fu_train, fu_valid, Tend=1, as_int = TRUE)
  h1 <- bw.nrd(data$time_std)
  h2 <- bw.nrd(data$time_std)^(5 / 6)
  
  #--- Fit FPCA ---
  tmp <- PP_FPCA_new(
    data$time_std, h1 = h1, h2 = h2,
    data$count, bw = "nrd", ngrid = n.grid, Tend = 1,
    K.select = K.select, propvar = propvar,
    #density.method = "local linear",
    derivatives = TRUE
  )
  
  #--- generate summary statistic ---
  out <- fpca.summary(data, tmp, fu_train, fu_valid)
  
  return(out)
}