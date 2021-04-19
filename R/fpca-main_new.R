######################################################################
### Second version of the functional principal component analysis on rare events (Wu et al.,	2013).
######################################################################
## FPCA approach by Wu et al (2013)
## n = num of patients
## t: observed event times of all the individuals, can have duplicates
## h: bandwidth
## N: vector, contains num of observed event of each patient
## index: index of the patient; default is NULL, i.e., patient are labelled from 1 to n.
## ngrid: number of grid points in estimating covariance function g
## subdivisions: number of subdivisions used in GCV
## propvar: proportion of variation used to select number of FPCs
## PPIC: if TRUE, the propvar will be ignored and PPIC will be used to select number of FPCS
## polybinx: if use the same partition (x) for the polynomial regression
########################################################################
#' @title Functional Principal Component Analysis (FPCA)
#' @description Performs FPCA to estimate the density function for each subject.
#' @param t named vector; standardized longitudinal encounter times for a code with the corresponding patient name. They should be less than or equal to 1.
#' @param h1,h2 integers; bandwidth used to estimate the mean intensity function and the covariance function. Default=null.
#' @param N named vector; the number of observed event with the corresponding patient name.
#' @param bw a character string; bandwidth estimating method when h1 and h2 are null. Default="ucv", but can also be "nrd0", "nrd", "bcv","SJ-dpi" and "SJ-ste".
#' @param Tend numeric; the upper bound of the encounter time for the estimated density function. Default=1.
#' @param N named vector; the number of observed event with the corresponding patient name.
#' @param n.grid an integer value for grid points used in estimating covariance function g. Default is \code{101}.
#' @param K.select characters indicating which method to choose the number of principal components K.
#'  Default is K.select="PropVar", and K.select="PPIC" is also available.
#' @param propvar a proportion of variation used to select number of FPCs. Default is \code{0.85}.
#' @param density.method a character string; the method of estimating density function when \code{K.select}="PPIC". Default is "kernal", but can also be "local linear".
#' @param polybinx logical; if use the same partition (x) for the polynomial regression when \code{density.method}="local linear". Default is FALSE.
#' @param derivatives logical; whether to estimate the first derivatives of the density function. Default is TRUE.
#' @export 
PP_FPCA_new <- function(t, 
                        h1 = NULL, h2 = NULL, N = NULL, bw = "ucv", Tend = 1, 
                        ngrid = 101, K.select = c("PropVar", "PPIC"), Kmax = 10,
                        propvar = 0.85, 
                        density.method = c("kernel", "local linear"), ## PPIC
                        polybinx = FALSE, derivatives = TRUE) {
  ## generate appropriate parameters
  Tend = 1
  par = gen.par(N,h1,h2,bw,ngrid,Tend)
  N = par$N; h1 = par$h1; h2 = par$h2; x = par$x

  ## get f_i and G
  fG.est = fpca.est.fG(x, t, N, h1, h2, Kmax)
  
  ## select number of FPCs
  K = fpca.est.K(K.select, Kmax, propvar,polybinx, density.method, x, t, N, 
                 fG.est$delta, fG.est$f_mu, fG.est$G.eigen, fG.est$xi)
  
  ## get density functions and scores using est K
  fi_Kadj <- fpca.Kadj.f(K, Kmax, ngrid, derivatives, x, N,
                             fG.est$delta, fG.est$f_mu, 
                             fG.est$G.eigen, fG.est$xi)
  
  result = list(
    mean = data.frame(x = x, f_mu = fG.est$f_mu),
    cov = fG.est$G,
    baseline = fG.est$baseline,
    K = K,
    densities = fi_Kadj$fi,
    derivatives = fi_Kadj$fi_d,
    scores = fi_Kadj$scores,
    basis = fi_Kadj$basis
  )
  return(result)
}
