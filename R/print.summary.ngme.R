#' print ngme summary
#'
#' @param object summary.ngme object
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 


print.summary.ngme <- function(object, ...) {

#   cat("Call:\n")
#   print(object$call)

#   if(object$estimate_fisher){#fisher matrix is estimated
#     cat("\n"); cat("\n")
#     cat("Fixed effects:\n")
#     cat("\n")
#     printCoefmat(object$fixed_results, P.valueS = T, has.Pvalue = T)

#     cat("\n"); cat("\n")
#     cat("Random effects:\n")
#     cat("\n")
#     cat("Sigma matrix:\n")
#     print(object$random_results$Sigma_est)
#     cat("\n")
#     cat("Standard errors for Sigma matrix:\n")
#     print(object$random_results$Sigma_est_se)

#     if(object$random_distr %in% c("NIG", "tdist")){
#       mu_results <-
#         matrix(
#         c(object$random_results$mu_est,
#           object$random_results$mu_est_se),
#         ncol = 2, byrow = F,
#         dimnames = list(paste0("mu_", 1:length(object$random_results$mu_est)), c("Estimate", "SE"))
#       )

#       nu_results <-
#         matrix(c(object$random_results$nu_est,
#                  object$random_results$nu_est_se), nrow = 1,
#                dimnames = list(c("nu"), c("Estimate", "SE")))

#       cat("\n")
#       print(mu_results)

#       cat("\n")
#       print(nu_results)
#     }

#     if(object$use_process){
#       cat("\n"); cat("\n")
#       cat("Operator:\n")
#       cat("\n")
#       print(object$operator_results)

#       if(object$process_distr %in% c("NIG", "GAL")){
#         cat("\n"); cat("\n")
#         cat("Process:\n")
#         cat("\n")
#         print(object$process_results)
#       }
#     }

#     cat("\n"); cat("\n")
#     cat("Measurement error:\n")
#     cat("\n")
#     print(object$meas_error_results)

#   }else{#fisher is not estimated

#     cat("\n"); cat("\n")
#     cat("Fixed effects:\n")
#     cat("\n")
#     print(object$fixed_results)

#     cat("\n"); cat("\n")
#     cat("Random effects:\n")
#     cat("\n")
#     cat("Sigma matrix:\n")
#     print(object$random_results$Sigma_est)

#     if(object$random_distr %in% c("NIG", "tdist")){
#       mu_results <-
#         matrix(
#           object$random_results$mu_est,
#           ncol = 1,
#           dimnames = list(paste0("mu_", 1:length(object$random_results$mu_est)), c("Estimate"))
#         )

#       nu_results <-
#         matrix(object$random_results$nu_est,
#                nrow = 1,
#                dimnames = list(c("nu"), c("Estimate")))

#       cat("\n")
#       print(mu_results)

#       cat("\n")
#       print(nu_results)
#     }

#     if(object$use_process){

#       cat("\n"); cat("\n")
#       cat("Operator:\n")
#       cat("\n")
#       print(object$operator_results)

#       if(object$process_distr %in% c("NIG", "GAL")){
#         cat("\n"); cat("\n")
#         cat("Process:\n")
#         cat("\n")
#         print(object$process_results)
#       }
#     }

#     cat("\n"); cat("\n")
#     cat("Measurement error:\n")
#     cat("\n")
#     print(object$meas_error_results)

#   }
}