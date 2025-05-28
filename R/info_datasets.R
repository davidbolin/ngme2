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


#'srft data
#'
#' @description srft data
#' A data-set that contains for estimated glomerular filtration rate and
#' a number of explanatory variables for 22,910 primary care patients
#'
#' @format A data frame with 392870 rows and 6 variables:
#' \describe{
#' \item{id}{identification number}
#' \item{egfr}{estimated glomerular filtration rate,
#' in mL/min per 1.73\eqn{m^2} of body surface area)}
#' \item{sex}{sex of the patient; 0: male, 1: female}
#' \item{bage}{baseline age, in yars}
#' \item{fu}{follow-up time; in years}
#' \item{pwl}{a piece-wise linear term defined as max(age - 56.5, 0),
#' where age is age at measurement}
#' }
#'@source \url{https://doi.org/10.1093/biostatistics/kxu053}
"srft_data"
