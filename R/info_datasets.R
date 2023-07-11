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


#' The swamp of Cienaga Grande in Santa Marta, Colombia
#'
#'There is a total of 114 locations where some properties of the
#' swamp were measured. Those measurements were taken twice, however there is no information available about their chronological order so this data cannot
#' be treated as spatiotemporal, despite that, the multiple measurements can be treated as replicates.
#'
#' @format A data frame with 218 rows and 6 columns.
#' \describe{
#'   \item{East, North}{location}
#'   \item{depth}{depth of the swamp}
#'   \item{temp}{temperature}
#'   \item{oxyg}{oxygen}
#'   \item{measurement}{1 means the first measurement, 2 the second}
#' }
#' @source ..
"cienaga"

#' The x y location of the border of the swamp of Cienaga Grande in Santa Marta, Colombia
#'
#' The data is of dimension 472 * 2. It contains the x and y coordinates of the border of the swamp.
#'
#' @format A data frame with 472 rows and 2 columns.
#' \describe{
#'   \item{East, North}{location}
#' }
#' @source ..
"cienaga.border"