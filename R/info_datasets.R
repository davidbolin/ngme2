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


#' Quarterly PM2.5 monitor averages over CONUS in 2022
#'
#' A processed example dataset containing quarterly fine particulate matter
#' (PM2.5) averages for monitoring sites in the contiguous United States.
#' The data are derived from the US EPA Air Quality System and AirData
#' pre-generated monitor files and aggregated to site-quarter means.
#'
#' @docType data
#' @usage data("pm25_quarterly_2022")
#' @format A data frame with 2640 rows and 9 columns.
#' \describe{
#'   \item{site_id}{EPA monitor site identifier.}
#'   \item{state}{State name.}
#'   \item{county}{County name.}
#'   \item{lat}{Monitor latitude in decimal degrees.}
#'   \item{lon}{Monitor longitude in decimal degrees.}
#'   \item{quarter}{Calendar quarter, coded as \code{Q1}, \code{Q2}, \code{Q3}, or \code{Q4}.}
#'   \item{pm25_q}{Quarterly average PM2.5 concentration in micrograms per cubic meter.}
#'   \item{n_days}{Number of daily monitor summaries contributing to the quarterly average.}
#'   \item{replicate}{Integer replicate index corresponding to quarter.}
#' }
#'
#' @source US EPA Air Quality System Data Mart and EPA AirData pre-generated
#' files for PM2.5 monitor summaries.
#' \url{https://www.epa.gov/aqs}
#' \url{https://aqs.epa.gov/aqsweb/airdata/download_files.html}
"pm25_quarterly_2022"


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
