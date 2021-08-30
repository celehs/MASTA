#' @title Function for validating the MASTA algorithm with new data
#' @description This function builds an algorithm to identify the occurrence of event outcome from trajectories of several predictors.
#' @param object The object returned by the \code{masta.fit} function
#' @param new_longitudinal longitudinal encounter times for the validation. See an example as \code{longitudinal} data.
#' @param new_follow_up_time the follow-up data for the validation
#' @param new_survival the labeled data for the validation. The columns should be 1) id, 2) event indicator, 3) event time, followed by baseline predictors.
#' @return A list with components:
#' @return \item{result_valid}{Performance of the derived algorithm. C-statistics, etc.}
#' @return \item{risk_score_valid}{Predcited Score}
#' @return \item{pred_surv_valid}{Predicted incidence curves}
#' @export
masta_validation <- function(object, new_longitudinal, new_follow_up_time, new_survival) {

#--- data derivation -- 
   org_longitudinal <- object$fpca_obj$data$longitudinal
   org_follow_up_time <- object$fpca_obj$data$follow_up_time
   org_survival <- object$data_survival

   # training data
   train_id <- org_follow_up_time$id[org_follow_up_time$train_valid == 1]

   # replace the validation data with the new data
   new_follow_up_time$train_valid <- 2
   follow_up_time_val <- rbind(org_follow_up_time[org_follow_up_time$id %in% train_id, ], new_follow_up_time)
   longitudinal_val <- rbind(org_longitudinal[org_longitudinal$id %in% train_id, ], new_longitudinal)
   survival_val <- rbind(org_survival[org_survival$id %in% train_id, ], new_survival)

  #--- parameters used for FPCA and Fit
  #- list(K.select = K.select, Kmax = Kmax , n.grid = n.grid, propvar = propvar)
  #- list(Tend=Tend, cov_group = cov_group, thresh = thresh, PCAthresh = PCAthresh, seed = seed)
  
  fp <- fpca.combine(
    longitudinal_val,
    follow_up_time_val,
    K.select = object$fpca_obj$parm$K.select,
    Kmax = object$fpca_obj$parm$Kmax,
    n.grid = object$fpca_obj$parm$n.grid,
    propvar = object$fpca_obj$parm$propvar
  )
  
  b <- masta.fit(
    fp,
    survival_val,
    follow_up_time_val,
    Tend = object$parm$Tend,
    cov_group = object$parm$cov_group,
    thresh = object$parm$thresh,
    PCAthresh = object$parm$PCAthresh,
    seed = object$parm$seed
  )
  # ================
  # Output
  # ================
  out <- list()

  #--- result with the validation data 
  out$result_valid <- b$result_valid

  #---- predicted score (validation data)
  out$risk_score_valid <- b$risk_score_valid
  out$pred_surv_valid <- b$pred_surv_valid
 return(out)
}
