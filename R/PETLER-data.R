#' @title Sample Data of the PETLER Pacakge
#' @description A list of training and validation data for illustrating the use of the PETLER algorithm
#' @format A list with 4 data tables:
#' \describe{
#'   \item{TrainSurv}{baseline survival data for training (labeled); 1st colum: patient id, 2nd colum: event indicator (1=event, 0=censoring), 3rd colum: event time, 4th colum: follow-up time, 5th colum--: covariates.}
#'   \item{ValidSurv}{baseline survival data for validation; baseline survival data for training (labeled); 1st colum: patient id, 2nd colum: event indicator (1=event, 0=censoring), 3rd colum: event time, 4th colum: follow-up time, 5th colum--: covariates.}
#'   \item{TrainCode}{1st colum: patient id, 2nd colum: follow-up time, 3rd colum: time label (month), 4th colum-- : predictors.}
#'   \item{ValidCode}{1st colum: patient id, 2nd colum: follow-up time, 3rd colum: time label (month), 4th colum-- : predictors.}
#' }
"data_org"


#' @title OLD Sample Data of the PETLER Pacakge
#' @description A list of training and validation data for illustrating the use of the PETLER algorithm
#' @format A list with 4 data tables:
#' \describe{
#'   \item{TrainSurv}{baseline survival data for training (labeled); 1st colum: patient id, 2nd colum: event indicator (1=event, 0=censoring), 3rd colum: event time, 4th colum: follow-up time, 5th colum--: covariates.}
#'   \item{ValidSurv}{baseline survival data for validation; baseline survival data for training (labeled); 1st colum: patient id, 2nd colum: event indicator (1=event, 0=censoring), 3rd colum: event time, 4th colum: follow-up time, 5th colum--: covariates.}
#'   \item{TrainCode}{1st colum: patient id, 2nd colum: follow-up time, 3rd colum: time label (month), 4th colum-- : predictors.}
#'   \item{ValidCode}{1st colum: patient id, 2nd colum: follow-up time, 3rd colum: time label (month), 4th colum-- : predictors.}
#' }
"data_org_old"
