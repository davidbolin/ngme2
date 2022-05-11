#' ngme summary
#'
#' @param object 
#'
#' @return list of objects
#' @export
#'
#' @examples
#' 
#' 

summary.ngme <- function(object, ...) {
    ret <- list()
    
    # ret["m_err"] = object$estimates$m_err
    ret = object$estimates
    
    class(ret) <- "summary.ngme"
    return (ret)
}
