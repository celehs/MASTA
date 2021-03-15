#' @name MASTA-package
#' @aliases MASTA
#' @docType package
#' @title MASTA Package
#' @description
#' Implements an algorithm to identify the occurrence of event outcome from trajectories of several predictors.
#' @details Visit the documentation website for more details. https://celehs.github.io/MASTA/
#' @import survival doParallel foreach data.table survC1 rootSolve splines glmnet gglasso rpart rpart.utils
#' @importFrom graphics hist
#' @importFrom utils head
#' @importFrom stats IQR aggregate approx as.formula binomial bw.SJ bw.bcv bw.nrd bw.nrd0 bw.ucv
#' coef cor density dnorm glm knots median optim optimize prcomp predict quantile runif sd stepfun uniroot var
NULL


#' @title Sample Data of the MASTA Pacakge
#' @description A list of training and validation data for illustrating the use of the MASTA algorithm
#' @format A list with 4 data tables:
#' \describe{
#'   \item{TrainSurv}{baseline survival data for training (labeled); 1st colum: patient id, 2nd colum: event indicator (1=event, 0=censoring), 3rd colum: event time, 4th colum: follow-up time, 5th colum--: covariates.}
#'   \item{ValidSurv}{baseline survival data for validation; baseline survival data for training (labeled); 1st colum: patient id, 2nd colum: event indicator (1=event, 0=censoring), 3rd colum: event time, 4th colum: follow-up time, 5th colum--: covariates.}
#'   \item{TrainCode}{1st colum: patient id, 2nd colum: follow-up time, 3rd colum: time label (month), 4th colum-- : predictors.}
#'   \item{ValidCode}{1st colum: patient id, 2nd colum: follow-up time, 3rd colum: time label (month), 4th colum-- : predictors.}
#' }
"longitudinal"

#' @title Sample Data of the MASTA Pacakge
#' @description A list of training and validation data for illustrating the use of the MASTA algorithm
#' @format A list with 4 data tables:
#' \describe{
#'   \item{TrainSurv}{baseline survival data for training (labeled); 1st colum: patient id, 2nd colum: event indicator (1=event, 0=censoring), 3rd colum: event time, 4th colum: follow-up time, 5th colum--: covariates.}
#'   \item{ValidSurv}{baseline survival data for validation; baseline survival data for training (labeled); 1st colum: patient id, 2nd colum: event indicator (1=event, 0=censoring), 3rd colum: event time, 4th colum: follow-up time, 5th colum--: covariates.}
#'   \item{TrainCode}{1st colum: patient id, 2nd colum: follow-up time, 3rd colum: time label (month), 4th colum-- : predictors.}
#'   \item{ValidCode}{1st colum: patient id, 2nd colum: follow-up time, 3rd colum: time label (month), 4th colum-- : predictors.}
#' }
"survival"

#' @title Sample Data of the MASTA Pacakge
#' @description A list of training and validation data for illustrating the use of the MASTA algorithm
#' @format A list with 4 data tables:
#' \describe{
#'   \item{TrainSurv}{baseline survival data for training (labeled); 1st colum: patient id, 2nd colum: event indicator (1=event, 0=censoring), 3rd colum: event time, 4th colum: follow-up time, 5th colum--: covariates.}
#'   \item{ValidSurv}{baseline survival data for validation; baseline survival data for training (labeled); 1st colum: patient id, 2nd colum: event indicator (1=event, 0=censoring), 3rd colum: event time, 4th colum: follow-up time, 5th colum--: covariates.}
#'   \item{TrainCode}{1st colum: patient id, 2nd colum: follow-up time, 3rd colum: time label (month), 4th colum-- : predictors.}
#'   \item{ValidCode}{1st colum: patient id, 2nd colum: follow-up time, 3rd colum: time label (month), 4th colum-- : predictors.}
#' }
"follow_up_time"

VTM <- function(vc, dm) matrix(vc, ncol = length(vc), nrow = dm, byrow = TRUE)
