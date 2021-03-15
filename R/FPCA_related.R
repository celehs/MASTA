#' @importFrom utils head
#' @importFrom graphics hist
#' @importFrom stats IQR approx bw.SJ bw.bcv bw.nrd bw.nrd0 bw.ucv density dnorm optimize
NULL

VTM <- function(vc, dm) matrix(vc, ncol = length(vc), nrow = dm, byrow = TRUE)

GetPK <- function(id, t, tseq, fu) {
  id <- as.character(id)
  PK <- rep(NA, length(fu))
  names(PK) <- names(fu)
  e <- diff(tseq)
  aa <- tapply(t, id, FUN = function(t) {
    x <- hist(t, breaks = tseq, plot = FALSE)$counts
    avg_sp <- cumsum(x) / tseq[-1]
    avg_sp2 <- c(0, head(avg_sp, -1))
    which.max((avg_sp - avg_sp2) / (avg_sp2 + e))
  })
  PK[unique(id)] <- tseq[aa + 1]
  PK
}





fpca.check <- function(time, fu_train, fu_valid){
  #-1-- minimum in "time" should be greater than or equal to 1
  if (!is.vector(time)) stop("Data Entry Issue: 'time' should be a vector")
  chk <- min(time)
  if (chk < 1) {
    stop("Data Entry Issue: The minimum of 'time' should be 1. Please do not code Month 0 for Days 1 to 30 but Month 1.")
  }
  
  #-2-- fu_train and fu_valid should be all vecotrs and one-record per subject and mutually exclusive
  if (!is.vector(fu_train)) stop("Data Entry Issue: 'fu_train' should be a vector")
  if (!is.vector(fu_valid)) stop("Data Entry Issue: 'fu_valid' should be a vector")
  chk1 <- unique(names(fu_train))
  chk2 <- unique(names(fu_valid))
  if (length(fu_train) != length(chk1)) stop("Data Entry Issue: More than one entry from one subject in 'fu_train'")
  if (length(fu_valid) != length(chk2)) stop("Data Entry Issue: More than one entry from one subject in 'fu_valid'")
  chk3 <- c(names(fu_train), names(fu_valid))
  if (length(chk3) != length(unique(chk3))) stop("Data Entry Issue: There are subjects who are in both 'fu_train' and 'fu_valid'")
  
  #-3-- subjects in time should be embedded by fu_train or fu_valid
  chk <- match(names(time), c(names(fu_train), names(fu_valid)))
  if (sum(is.na(chk) != 0)) stop("Data Entry Issue: Some subjects in 'time' do not have follow-up time information in 'fu_train' or 'fu_valid'")
}





fpca.pre <- function(time, fu_train, fu_valid, Tend, as_int, least_c=1, least_uc=1){
  fu <- c(fu_train, fu_valid)
  uni.count <- tapply(time, names(time), function(x){length(unique(x))})
  count <- tapply(time, names(time), function(x){length(x)})
  names_least <- names(uni.count[uni.count>=least_uc & count>=least_c])
  names_fu = names(fu)
  time <- time[(names(time) %in% names_fu) 
               & (names(time) %in%  names_least)]    #only keep data with follow-up info
  names_time <- names(time)
  count <- tapply(time, names_time, length)
  time_std <- time / fu[names_time]
  if (as_int) {
    names_time <- as.integer(names_time)
  }
  count_all <- 0 * fu
  count_all[names(count)] = count
  out = list(`fu`=fu,
             `time`=time,`time_std`=time_std,`names_time`=names_time,
             `count`=count,`count_all`=count_all)
  return(out)
}






fpca.summary <- function(data, tmp, fu_train, fu_valid){
  
  PKTS <- GetPK(
    id = data$names_time, ### INT/CHAR ###
    t = data$time,
    tseq = 0:floor(max(data$fu)),
    fu = data$fu
  )
  
  ft.e <- cbind(
    matrix(data$fu, nrow = length(data$fu), ncol = 3),
    -tmp$baseline[1], log(1 + data$count_all)
  )                                             # if longitudinal data is missing - use baseline 
  rownames(ft.e) = names(data$count_all)
  
  pos <- names(data$count)
  ft.e[pos, 1] <- tapply(data$time, data$names_time, min)
  locm <- unlist(apply(tmp$densities[, 1:length(pos) + 1], 2, which.max))  # peak
  ft.e[pos, 2] <- ft.e[pos, 2] * tmp$densities[locm, 1]  #density[,1] is the x-coordinate
  ft.e[pos, 3] <- ft.e[pos, 3] * tmp$derivatives[
    sapply(1:length(pos), function(i) {
      #which.max(tmp$derivatives[1:locm[i], i + 1])  #--??
      which.max(tmp$derivatives[, i + 1])
    }), 1]
  ft.e[pos, 4] <- tmp$scores[, 2]
  colnames(ft.e) <- c("1stCode", "Pk", "ChP", "1stScore", "logN")
  TrainN = data.frame("id" = names(fu_train),
                      "pred_total" = data$count_all[names(fu_train)])
  ValidN = data.frame("id" = names(fu_valid),
                      "pred_total" = data$count_all[names(fu_valid)])
  
  out <- list(
    TrainN = TrainN,
    ValidN = ValidN,
    TrainFt = ft.e[names(fu_train),],
    ValidFt = ft.e[names(fu_valid),],
    TrainPK = PKTS[names(fu_train)],
    ValidPK = PKTS[names(fu_valid)]
  )
  class(out) <- "fpca"
  return(out)
}




#' @title Functional Principal Component Analysis (FPCA)
#' @description Performs FPCA to extract features from longitudinal encounter data.
#' @param time_code longitudinal encounter times. They should be greater than or equal to 1.
#' @param follow_up_time follow-up time (training)
#' @param K.select characters indicating which method to choose the number of principal components K.
#'  Default is K.select="PropVar", and K.select="PPIC" is also available.
#' @param n.grid an integer value for grid points used in estimating covariance function g. Default is \code{401}.
#' @param propvar a proportion of variation used to select number of FPCs. Default is \code{0.85}.
#' @export
fpca.combine <- function(time_code, follow_up_time, 
                         K.select = "PropVar", n.grid = 401, propvar = 0.85){
  longitudinal$code = as.numeric(longitudinal$code)
  codes = unique(longitudinal$code)
  if(min(codes)!=1) stop("Codes should start from 1.")
  
  diff.index = codes[2:length(codes)] - codes[1:(length(codes)-1)]
  if(sum(diff.index != 1) > 0) stop("Codes must be incremental with no skipping number.")
 
  D = sapply(codes, function(i){
    x = longitudinal$time[longitudinal$code == i]
    names(x) = longitudinal$id[longitudinal$code == i]
    return(x)
  })
  names(D) = codes
  follow_up_train <- follow_up_time[follow_up_time$train_valid==1, ]
  follow_up_valid <- follow_up_time[follow_up_time$train_valid==2, ]
  
  fu_train <- follow_up_train$fu_time
  fu_valid <- follow_up_valid$fu_time
  names(fu_train) <- follow_up_train$id
  names(fu_valid) <- follow_up_valid$id

  TrainFt <- c()
  ValidFt <- c()
  TrainPK <- c()
  ValidPK <- c()
  TrainN <- c()
  ValidN <- c()

  for (i in 1:length(codes)) {
    time <- D[[i]]
    ans <- fpca.new(time, fu_train, fu_valid, K.select = "PropVar") 
    #--- create an object for fitting --
    Ft_name <- colnames(ans$TrainFt)
    Ft_name <- paste0(Ft_name, codes[i])
    colnames(ans$TrainFt) = colnames(ans$ValidFt) = Ft_name
    TrainFt <- cbind(TrainFt, ans$TrainFt)
    ValidFt <- cbind(ValidFt, ans$ValidFt)
    
    TrainPK <- cbind(TrainPK, ans$TrainPK)
    ValidPK <- cbind(ValidPK, ans$ValidPK)
    if (i == 1) {
      TrainN <- ans$TrainN
      ValidN <- ans$ValidN
    }else{
      TrainN <- cbind(TrainN, ans$TrainN[, 2])
      ValidN <- cbind(ValidN, ans$ValidN[, 2])
    }
  }
  
  colnames(TrainPK) <- colnames(ValidPK) <- paste0("pred", codes)
  colnames(TrainN) <- colnames(ValidN) <- c("id", paste0("pred", codes, "_total"))
  
  Z <- list()
  Z$TrainFt <- TrainFt
  Z$ValidFt <- ValidFt
  Z$TrainPK <- TrainPK
  Z$ValidPK <- ValidPK
  Z$TrainN <- TrainN
  Z$ValidN <- ValidN
  return(Z)
}