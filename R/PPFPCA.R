VTM <- function(vc, dm) {
  matrix(vc, ncol = length(vc), nrow = dm, byrow = TRUE)
}

FPC_Kern_S <- function(x, t, N, h1, h2) {
  grp <- rep(seq(N), N)
  M <- outer(t, x, "-")
  D1 <- dnorm(M, 0, h1)
  D2 <- if (h1 == h2) D1 else dnorm(M, 0, h2) 
  v <- rowsum(D2, grp)
  list(f_mu = colSums(D1), 
       Cov_G = crossprod(v) - crossprod(D2))
}

FPC.Kern.S <- function(x, t, N, h1 = NULL, h2 = NULL, bw = "ucv", nsubs = 10, n_core = NULL) {
  h <- switch(bw, 
              "nrd0" = bw.nrd0(t),
              "nrd" = bw.nrd(t),
              "ucv" = bw.ucv(t), # leave one out cross validation
              "bcv" = bw.bcv(t), # biased cross validation            
              "SJ-ste" = bw.SJ(t, method = "ste"), 
              "SJ-dpi" = bw.SJ(t, method = "dpi"))
  if (is.null(h1)) h1 <- h
  if (is.null(h2)) h2 <- h
  n <- length(N)
  a <- n %% nsubs
  b <- floor(n / nsubs)
  subsize <- rep((b + 1):b, c(a, nsubs - a))
  cumsumsub <- cumsum(c(0, subsize))
  grp <- rep(seq(subsize), subsize)
  cumsumN <- c(0, cumsum(rowsum(N, grp))) 
  if (is.null(n_core)) n_core <- parallel::detectCores()
  registerDoParallel(cores = n_core)  
  tmplist <- foreach(i = 1:nsubs) %dopar% {
    tsub <- t[(cumsumN[i] + 1):cumsumN[i + 1]]
    Nsub <- N[(cumsumsub[i] + 1):cumsumsub[i + 1]]
    FPC_Kern_S(x, tsub, Nsub, h1, h2)
  }
  L1 <- lapply(tmplist, `[[`, 1)
  L2 <- lapply(tmplist, `[[`, 2)
  f_mu <- Reduce("+", L1) / sum(N)
  G <- Reduce("+", L2) / sum(N * (N - 1)) - outer(f_mu, f_mu)
  list(f_mu = f_mu, G = G, cumsumsub = cumsumsub, cumsumN = cumsumN)
}

######################################################################
### Second version of the functional principal component analysis 
### on rare events (Wu et al.,	2013). 

########################################################################

########################################################################
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
############                    Parallel                    ############     
########################################################################
PP_FPCA_Parallel <- function(t, h1 = NULL, h2 = NULL, N, bw = "ucv", Tend = 1, # assume it is transformed to [0,1]
                             ngrid = 101, K.select = c('PropVar', 'PPIC'), Kmax = 10,
                             propvar = 0.9, ## For K.select == PropVar
                             density.method = c("kernel", "local linear"), ## PPIC
                             polybinx = FALSE, derivatives = FALSE,
                             nsubs = 10, subsize = NULL, PPIC.sub = TRUE) {
  ## eleminate patients with 0 observed event
  if (sum(N == 0) > 0) {
    NN <- N
    N.index <- N != 0
    N <- N[N.index]
    # cat("Note: patients with zero observed event have been eliminated from analysis!","\n")
  }
  n <- length(N) # number of patients with at least one observed event
  ## get a fine grid (x,y)
  # tmp  = range(t)
  x <- y <- seq(0, Tend, length.out = ngrid)
  # print(system.time(tmplist <- FPC.Kern.S(t, N, h1, h2, bw, Tend, ngrid))) 
  # f_mu <- tmplist$f_mu
  # G <- tmplist$G

  print(system.time(tmp <- FPC.Kern.S(x, t, N, h1, h2, bw = bw, nsubs = nsubs)))
  f_mu <- tmp$f_mu
  G <- tmp$G  
  cumsumsub <- tmp$cumsumsub
  cumsumN <- tmp$cumsumN

  G.eigen <- svd(G)
  delta <- sqrt(x[2] - x[1])
  baseline <- as.numeric(t(G.eigen$u[, 1:Kmax]) %*% f_mu * delta) # second term in xi
  cumsumN2 <- c(0, cumsum(N))
  # interpolate the eigenvectors to get eigenfunctions and then get the xi's
  xi <- foreach(i = 1:nsubs, .combine = cbind) %dopar% {
    tmp2 <- apply(G.eigen$u[, 1:Kmax], 2, function(s) approx(x = x, y = s, xout = t[{cumsumN[i] + 1}:cumsumN[i + 1]])$y) / delta 
    indexi <- rep({1 + cumsumsub[i]}:cumsumsub[i + 1], N[{cumsumsub[i] + 1}:cumsumsub[i + 1]])
    -baseline + t(apply(tmp2, 2, FUN = function(x) tapply(x, indexi, mean)))
  }
  attr(xi, "dimnames") <- NULL
  ## select number of FPCs
  K.select <- match.arg(K.select)
  if (K.select == "PropVar") {
    # method 1: proportion of variation >= 90%
    K <- cumsum(G.eigen$d) / sum(G.eigen$d)
    K <- min(which(K >= propvar))
    if (K > Kmax) K <- Kmax
  }
  ## get density functions
  if (K == 1) {
    tmp <- foreach(i = 1:nsubs, .combine = cbind) %dopar% {
      rawden <- f_mu + outer(G.eigen$u[, 1], xi[1, {cumsumsub[i] + 1}:cumsumsub[i + 1]]) / delta # density functions
      ## make adjustments to get valid density functions (nonnegative+integrate to 1)
      apply(rawden, 2, function(x) {
        x[x < 0] <- 0     # non-negative
        x <- x / sum(x) # integrate to delta^2
        return(x)
      }) / delta^2
    }
    scores <- data.frame(xi[1:max(K, Kmax), ])
  } else {
    tmp <- foreach(i = 1:nsubs, .combine = cbind) %dopar% {
      rawden <- f_mu + { G.eigen$u[, 1:K] / delta} %*% xi[1:K, {cumsumsub[i] + 1}:cumsumsub[i + 1]] # density functions
      # make adjustments to get valid density functions (nonnegative+integrate to 1)
      apply(rawden, 2, function(x) {
        x[x < 0] <- 0 # non-negative
        x <- x / sum(x) # integrate to delta^2
        return(x)
      }) / delta^2
    }
    scores <- data.frame(t(xi[1:max(K, Kmax), ]))
  }
  names(scores) <- paste0("score_", seq(1:max(K, Kmax))) # FPC scores
  basis <- data.frame(cbind(x, G.eigen$u[, 1:max(K, Kmax)] / delta))
  names(basis) <- c("x", paste0("basis_", seq(1:max(K, Kmax)))) # FPC eigenfunctions
  ## name the density functions 
  if (exists("N.index")) {
    colnames(tmp) <- paste0("f", which(N.index))
    scores <- data.frame(id = as.character(which(N.index)), 
                         scores, stringsAsFactors = FALSE) # add id to scores
  } else{
    colnames(tmp) <- paste0("f",seq(1,n))
    scores <- data.frame(id = as.character(seq(1, n)), scores, stringsAsFactors = FALSE)
  }
  tmp <- as.data.frame(tmp)
  tmp <- data.frame(x = x, tmp)
  # ## return K and prop of var by now
  # tmp     = data.frame(tmp,info=c(K,sum(G.eigen$d[1:K])/sum(G.eigen$d),rep(0,ngrid-2)))
  ## get derivatives if derivatives = TRUE
  if (derivatives) {
    tmp2 <- foreach(i = 1:nsubs, .combine = cbind) %dopar% {
      { tmp[-c(1:2), { cumsumsub[i] + 1}:cumsumsub[i + 1] + 1] - tmp[-c(1:2 + ngrid - 2), {cumsumsub[i] + 1}:cumsumsub[i + 1] + 1]} / diff(x, lag = 2) 
    }
    tmp2 <- as.data.frame(tmp2)
    tmp2 <- data.frame(x = x[-c(1, ngrid)], tmp2)
    colnames(tmp2) <- paste0("d", colnames(tmp))
    return(list(scores = scores, 
                densities = tmp, 
                derivatives = tmp2, 
                mean = data.frame(x = x, f_mu = f_mu),
                cov = G,
                basis = basis, 
                baseline = baseline,
                K = K))
  }
  else{
    return(list(scores = scores, 
                densities = tmp, 
                mean = data.frame(x = x, f_mu = f_mu), 
                cov  = G,
                basis = basis, 
                baseline = baseline,
                K = K))
  }
}


########################################################################
############        Predict Scores for New Subjects         ############     
########################################################################
PP_FPCA_Pred <- function(t, N, mean.fun, eigen.fun, K, parallel = FALSE) {
  delta <- mean.fun[2, 1] - mean.fun[1, 1]
  tmp <- as.numeric(t(eigen.fun[, -1]) %*% mean.fun[, -1] * delta) # second term in xi
  tmp2 <- apply(eigen.fun[, -1], 2, function(s) approx(x = mean.fun$x, y = s, xout = t)$y) ### check whether we need parallel this
  indx <- rep(1:length(N), N)
  xi <- -tmp + t(apply(tmp2, 2, FUN = function(x) tapply(x, indx, mean))) # FPC scores, ith column for ith patient
  if (K == 1) {
    tmp <- mean.fun[, -1] + outer(eigen.fun[, 2], xi[1, ]) # density functions
  } else {
    tmp <- as.numeric(mean.fun[, -1]) + as.matrix(eigen.fun[, 1:K + 1]) %*% xi[1:K, ] # density functions
  }
  tmp <- apply(tmp, 2, function(x) {
    x[x < 0] <- 0     # non-negative
    x <- x / sum(x) # integrate to delta^2
    return(x)
  })
  tmp <- tmp / delta
  tmp2 <- {tmp[-c(1:2), ] - tmp[-c(1:2 + length(mean.fun[, 1]) - 2), ]} / diff(mean.fun[, 1], lag = 2) 
  list(scores = t(xi), 
       densities = cbind(mean.fun[, 1], tmp), 
       derivatives = cbind({mean.fun[-c(1:2), 1] + mean.fun[-c(1:2 + length(mean.fun[, 1]) - 2), 1]} / 2, tmp2))
}

### make prediction for patients with or without codes
PP_FPCA_Pred2 <- function(t, N, mean.fun, eigen.fun, K, parallel = FALSE){
  NZ <- N == 0
  delta <- mean.fun[2, 1] - mean.fun[1, 1]
  baseline <- as.numeric(t(eigen.fun[, -1]) %*% mean.fun[, -1] * delta) # second term in xi
  tmp2 <- apply(eigen.fun[, -1], 2, function(s) approx(x = mean.fun$x, y = s, xout = t)$y) ### check whether we need parallel this
  indx <- rep(1:length(N[!NZ]), N[!NZ])
  xi <- -baseline + t(apply(tmp2, 2, FUN = function(x) tapply(x, indx, mean))) # FPC scores, ith column for ith patient
  if (K == 1) {
    tmp <- mean.fun[, -1] + outer(eigen.fun[, 2], xi[1, ]) # density functions
  } else {
    tmp <- as.numeric(mean.fun[, -1]) + as.matrix(eigen.fun[, 1:K+1]) %*% xi[1:K, ] # density functions
  }
  tmp <- apply(tmp, 2, function(x) {
    x[x < 0] <- 0 # non-negative
    x <- x / sum(x) # integrate to delta^2
    return(x)
  })
  tmp <- tmp / delta
  tmp2 <- {tmp[-c(1:2), ] -tmp[-c(1:2 + length(mean.fun[, 1]) - 2), ]} / diff(mean.fun[, 1], lag = 2) 
  list(scores = t(xi), 
       densities = cbind(mean.fun[, 1], tmp), 
       derivatives = cbind({mean.fun[-c(1:2), 1] + mean.fun[-c(1:2 + length(mean.fun[, 1]) - 2), 1]} / 2, tmp2), 
       baseline = baseline)
}
