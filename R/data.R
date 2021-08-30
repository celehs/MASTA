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
#' @description A long form longitudinal data for illustrating the use of the MASTA algorithm
#' @format A data.frame with 3 columns:
#' \describe{
#'   \item{code}{integer indicator for the type of codes. It must be incremental and start from 1, which contains no skipping number.}
#'   \item{id}{patient id.}
#'   \item{time}{encounter time with the code. One patient can encounter with different types of codes multiple times.}
#' }
"longitudinal"

#' @title Sample Data of the MASTA Pacakge
#' @description The follow up time of training and validation data for illustrating the use of the MASTA algorithm
#' @format A data.frame with 3 columns:
#' \describe{
#'   \item{id}{patient id.}
#'   \item{fu_time}{the follow up time for the patient.}
#'   \item{train_valid}{indicator for patients in the training data set; 1=training, 0=test}
#' }
"follow_up_time"

#' @title Sample Data of the MASTA Pacakge
#' @description Survival data for illustrating the use of the MASTA algorithm
#' @format A data.frame with 6 columns:
#' \describe{
#'   \item{id}{patient id.}
#'   \item{event_ind}{indicator for events; 1=event, 0=cencoring}
#'   \item{event_time}{the event time if \code{event_ind}=1.}
#'   \item{cov_1,cov_2,cov_3}{3 covariates used in a survival regression model.}
#' }
"survival"


#' @title Sample Data of the MASTA Pacakge
#' @description A long form longitudinal data for illustrating how to validate the MASTA algorithm with a new data
#' @format A data.frame with 3 columns:
#' \describe{
#'   \item{code}{integer indicator for the type of codes. It must be incremental and start from 1, which contains no skipping number.}
#'   \item{id}{patient id.}
#'   \item{time}{encounter time with the code. One patient can encounter with different types of codes multiple times.}
#' }
"new_longitudinal"

#' @title Sample Data of the MASTA Pacakge
#' @description The follow up time of training and validation data for illustrating how to validate the MASTA algorithm with a new data
#' @format A data.frame with 2 columns:
#' \describe{
#'   \item{id}{patient id.}
#'   \item{fu_time}{the follow up time for the patient.}
#' }
"new_follow_up_time"

#' @title Sample Data of the MASTA Pacakge
#' @description Survival data for illustrating how to validate the MASTA algorithm with a new data
#' @format A data.frame with 6 columns:
#' \describe{
#'   \item{id}{patient id.}
#'   \item{event_ind}{indicator for events; 1=event, 0=cencoring}
#'   \item{event_time}{the event time if \code{event_ind}=1.}
#'   \item{cov_1,cov_2,cov_3}{3 covariates used in a survival regression model.}
#' }
"new_survival"
