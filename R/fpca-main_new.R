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
PP_FPCA_new <- function(t, 
                        h1 = NULL, h2 = NULL, N = NULL, bw = "ucv", Tend = 1, # assume it is transformed to [0,1]
                        ngrid = 101, K.select = c("PropVar", "PPIC"), Kmax = 10,
                        propvar = 0.9, ## For K.select == PropVar
                        density.method = c("kernel", "local linear"), ## PPIC
                        polybinx = FALSE, derivatives = FALSE) {
  ## generate appropriate parameters
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
