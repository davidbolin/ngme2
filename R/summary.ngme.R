# #' ngme summary
# #'
# #' @param ngme.object
# #'
# #' @return list of objects
# #' @export
# #'
# #' @examples
# #'
# summary.ngme <- function(object) {
#     ret <- object

#     class(ret) <- "summary.ngme"
#     return (ret)
# }

# print.summary.ngme <- function(object) {
#     cat("************ NGME fit ************ \n")

#     m_err = object$estimates$m_err
#     l_est = object$estimates$latents_est
#     f_eff = object$estimates$f_eff

#     # f_eff
#     cat("\nFIXED EFFECTS:\n")
#     cat(f_eff); cat("\n")

#     # m_err
#     cat("\nMEASUREMENT ERROR:\n")
#     cat("Type: "); cat(object$family); cat("\n")
#     cat(m_err); cat("\n")

#     # process
#     cat("\nLATENT PROCESS:\n")
#     proc.names = paste(ngme_out$model.types, ngme_out$var.types)
#     for (i in seq_along(proc.names)) {
#         cat(paste0("Prcess ", i, ": ", proc.names[i])); cat("\n")

#         cat("Estimates:  \n")
#         est = unlist(l_est[[i]])
#         names = format(names(est), width = 7, justify = "right")
#         est   = format(est, width = 7, justify = "right", digits=4)

#         cat(names); cat("\n")
#         cat(est); cat("\n")
#     }

#     ret = object$estimates
# }


