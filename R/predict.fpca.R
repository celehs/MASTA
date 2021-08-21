predict_fpca <- function(object, newdata) {

  FPCA <- object$FPCA

  fu <- newdata$fu
  time <- newdata$time
  names_time <- names(time)
  time_std <- time / fu[names_time]
  #- time_std <- newdata$time_std

  count <- tapply(time, names_time, length) 
  NN <- rep(0, length(fu))
  N <- data.frame(id = as.character(names(fu)), pred_total = NN[names(fu)])
  PKTS2 <- GetPK(
    id = names_time, ### INT/CHAR ###
    t = time,
    tseq = 0:floor(max(fu)),
    fu = fu
  )
  time_std <- time_std[names(time_std) %in% names(fu)]
  time1 <- tapply(time, names_time, min)
  tmp <- PP.FPCA.Pred(time_std, count, FPCA$mean, FPCA$basis, FPCA$K)
  FPCA$ValidPred <- tmp
  ft.e2 <- cbind(
    matrix(fu, nrow = length(fu), ncol = 3),
    -tmp$baseline[1], log(1 + count)
  )
  pos <- count > 0
  locm <- unlist(apply(tmp$densities[, 1:sum(pos) + 1], 2, which.max))
  #-- just patch ( will come back to PP.FPCA.Pred later
  locm <- pmin(locm, nrow(tmp$derivatives))
  ft.e2[pos, 2] <- ft.e2[pos, 2] * tmp$densities[locm, 1]
  ft.e2[pos, 3] <- ft.e2[pos, 3] * tmp$derivatives[
    sapply(1:sum(pos), function(i) {
      which.max(tmp$derivatives[1:locm[i], i + 1])
    }), 1
    ]
  ft.e2[pos, 4] <- tmp$scores[, 2]
  ft.e.S2 <- cbind(
    ft.e2[, 1], VTM(-tmp$baseline[1:4], length(fu)), log(1 + count)
  )
  ft.e.S2[pos, 2:5] <- as.matrix(tmp$scores[, 2:5])
  colnames(ft.e2) <- c("1stCode", "Pk", "ChP", "1stScore", "logN")
  colnames(ft.e.S2) <- c("1stCode", "1stScore", "2ndScore", "3rdScore", "4thScore", "logN")
  rownames(ft.e.S2) <- rownames(ft.e2) <- names(PKTS2) <- names(fu)
  list(
    N = N,
    Ft = ft.e2,
    Sc = ft.e.S2,
    PK = PKTS2
  )
}
