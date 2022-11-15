#############################################
#' Argo float dataset
#' @description Argo floats measurements.
#' @docType data
#' @usage data("argo_float")
#' @format Data frame containing 274 observations on 4 variables.
#' \describe{
#'   \item{lat}{Latitude.}
#'   \item{lon}{Longitude.}
#'   \item{sal}{Salinity.}
#'   \item{temp}{Temperature.}
#' }
#' 
#' @details 
#' The floats have a pressure case made of aluminium that is about 1.3m long and about 20cm diameter. They weigh about 40kg. On the top is an antenna to communicate with the satellites that fix the float's position and receive the data. Also on the top are the temperature, salinity and pressure sensors.
#' 
#' @source Data can be obtained from \href{https://argo.ucsd.edu/}{Argo float website}.
"argo_float"