gen.par <- function(t, N, h1, h2, bw, ngrid, Tend){
  ## eliminate patients with 0 observed event
  if (is.null(N)){
    N <- tapply(t, names(t), function(x){length(x)})
  }else{
    if (sum(N == 0) > 0) {
      NN <- N
      N.index <- N != 0
      N <- N[N.index]
    }
  }
  
  n <- length(N) # number of patients with at least one observed event
  ## if h is null then set it to be the optimal bandwidth under GCV
  if (is.null(h1) & is.null(h2)) {
    if (bw == "nrd0") {
      h1 <- h2 <- bw.nrd0(t)
    } else if (bw == "nrd") {
      h1 <- h2 <- bw.nrd(t)
    } else if (bw == "ucv") { # leave one out cross validation
      h1 <- h2 <- bw.ucv(t)
    } else if (bw == "bcv") { # biased cross validation
      h1 <- h2 <- bw.bcv(t)
    } else if (bw == "SJ-ste") {
      h1 <- h2 <- bw.SJ(t, method = "ste")
    } else if (bw == "SJ-dpi") {
      h1 <- h2 <- bw.SJ(t, method = "dpi")
    } else {
      h1 <- h2 <- bw.ucv(t)
    }
  }
  
  ## get a fine grid (x,y)
  x <- seq(0, Tend, length.out = ngrid)
  return(list(`N`=N,`h1`=h1,`h2`=h2,`x`=x))
}

fpca.est.fG <- function(x, t, N, h1, h2, Kmax){
  ## estimate the mean density f_mu and the intermediate value g that is used to
  tmp <- FPC_Kern_S(x, t, N, h1, h2)
  f_mu <- as.numeric(tmp$f_mu) / sum(N)
  G <- tmp$Cov_G / sum(N * (N - 1)) - outer(f_mu, f_mu)
  G.eigen <- svd(G)
  delta <- sqrt(x[2] - x[1])
  baseline <- as.numeric(t(G.eigen$u[, 1:Kmax]) %*% f_mu * delta) # second term in xi
  # interpolate the eigenvectors to get eigenfunctions
  tmp2 <- apply(G.eigen$u[, 1:Kmax], 2, function(s) approx(x = x, y = s, xout = t)$y) / delta 
  n <- length(N)
  xi <- -baseline + t(apply(tmp2, 2, FUN = function(x) tapply(x, rep(1:n, N), mean))) # FPC scores, ith column for ith patient
  return(list(`delta`=delta, `f_mu`=f_mu,`G`=G, `G.eigen`=G.eigen, 
              `xi` = xi, `baseline`=baseline))
}

fpca.est.K <- function(K.select = c("PropVar", "PPIC"), Kmax, 
                       propvar, polybinx, density.method = c("kernel", "local linear"), 
                       x, t, N, delta, f_mu, G.eigen, xi){
  K.select <- match.arg(K.select)
  if (K.select == "PropVar") {
    # method 1: proportion of variation >= 90%
    K <- cumsum(G.eigen$d) / sum(G.eigen$d)
    # cat("Propvar=", K, "\n")
    K <- min(which(K >= propvar))
    if (K > Kmax) K <- Kmax
  } else if (K.select == "PPIC") {
    K <- Kmax
    if (missing(density.method)) density.method <- "kernel" else density.method <- match.arg(density.method)
    if (density.method == "local linear") {
      if (polybinx) {
        f_locpoly <- den_locpoly(partition = x, t = t, N = N)
      } else {
        f_locpoly <- den_locpoly2(t = t, N = N)
      }
    } else {
      cumsumN2 <- c(0, cumsum(N))
      n <- length(N)
      f_locpoly <- sapply(1:n, function(j) {
        tmp <- density(t[{
          cumsumN2[j] + 1
        }:cumsumN2[j + 1]], bw = "nrd0") # nrd or nrd0?
        tmp$y[sapply(t[{
          cumsumN2[j] + 1
        }:cumsumN2[j + 1]], function(s) which.min(abs(s - tmp$x)))]
      })
    }
    K <- sapply(1:K, function(k) {
      PPIC(
        K = k, f_locpoly = f_locpoly, t = t, N = N,
        f_mu = f_mu, G.eigen_v = G.eigen$u, xi = xi, xgrid = x, delta = delta
      )
    })
    K <- which.min(K)
  }
  return(K)
}

fpca.Kadj.f <- function(K, Kmax, ngrid, derivatives,
                        x, N, delta, f_mu, G.eigen, xi){
  if (K == 1) {
    fi <- f_mu + outer(G.eigen$u[, 1], xi[1, ]) / delta # density functions
    scores <- data.frame(xi[1:min(K, Kmax), ])
  } else {
    fi <- f_mu + (G.eigen$u[, 1:K] / delta) %*% xi[1:K, ] # density functions
    scores <- data.frame(t(xi[1:min(K, Kmax), ]))
  }
  
  ## make adjustments to get valid density functions (nonnegative+integrate to 1)
  fi <- apply(fi, 2, function(x) {
    x[x < 0] <- 0 # non-negative
    x <- x / sum(x) # integrate to delta^2
    return(x)
  })
  fi <- fi / delta^2
  
  names(scores) <- paste0("score_", seq(1:min(K, Kmax))) # FPC scores
  basis <- data.frame(cbind(x, G.eigen$u[, 1:min(K, Kmax)] / delta))
  names(basis) <- c("x", paste0("basis_", seq(1:min(K, Kmax)))) # FPC eigenfunctions
  
  ## name the density functions
  n = length(N)
  colnames(fi) <- paste0("f", seq(1, n))
  scores <- data.frame(
    id = as.character(names(N)),
    scores, stringsAsFactors = FALSE
  )
  
  fi <- as.data.frame(fi)
  fi <- data.frame(x = x, fi)
  
  ## get derivatives if requested
  if (derivatives) {
    fi_d <- {
      fi[-c(1:2), -1] - fi[-c(1:2 + ngrid - 2), -1]
    } / diff(x, lag = 2)
    fi_d <- as.data.frame(fi_d)
    fi_d <- data.frame(x = x[-c(1, ngrid)], fi_d)
    colnames(fi_d) <- paste0("d", colnames(fi))
  }else{
    fi_d = NULL
  }
  return(list(`fi`=fi, `fi_d`=fi_d, `scores`=scores, `basis`=basis))
}