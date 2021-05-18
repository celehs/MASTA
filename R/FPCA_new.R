
fpca.new <- function(time, fu_train, fu_valid,
                     K.select = "PropVar", Kmax = 5, n.grid = 401, propvar = 0.85) {
  
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
    K.select = K.select, Kmax = Kmax, propvar = propvar,
    #density.method = "local linear",
    derivatives = TRUE
  )
  
  #--- generate summary statistic ---
  out <- fpca.summary(data, tmp, fu_train, fu_valid)
  
  return(out)
}