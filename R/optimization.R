# different methods for SGD
# vanilla
# momentum
# adagrad
# rmsprop
# adam
# adamW

#' Vanilla SGD optimization
#'
#' @param stepsize stepsize for SGD
#'
#' @return a list of control variables for optimization 
#' (used in \code{control_opt} function)
#' @export
vanilla <- function(
  stepsize = 0.5
) {

  ret <- list(
    method         = "vanilla",
    stepsize       = stepsize,
    sgd_parameters = NULL
  )
  class(ret) <- "ngme_optimization"
  ret
}

#' Momentum SGD optimization
#'
#' @details
#' The update rule for momentum is:
#'  v_t = beta1 * v_{t-1} + beta2 * g_t
#'  x_{t+1} = x_t - stepsize * v_t
#'
#' @param stepsize stepsize for SGD
#' @param beta1 beta1 for momentum
#' @param beta2 beta2 for momentum
#'
#' @return a list of control variables for optimization 
#' (used in \code{control_opt} function)
#' @export
momentum <- function(
  stepsize = 0.05,
  beta1 = 0.9,
  beta2 = 1 - beta1
) {
  ret <- list(
    method         = "momentum",
    stepsize       = stepsize,
    sgd_parameters = c(beta1, beta2)
  )
  class(ret) <- "ngme_optimization"
  ret
}


#' AdaGrad SGD optimization
#'
#' @details
#' The update rule for AdaGrad is:
#'  v_t = v_{t-1} + g_t^2
#'  x_{t+1} = x_t - stepsize * g_t / (sqrt(v_t) + epsilon)
#'
#' @param stepsize stepsize for SGD
#' @param epsilon epsilon for numerical stability
#'
#' @return a list of control variables for optimization
#' (used in \code{control_opt} function)
#' @export
adagrad <- function(
  stepsize = 0.05,
  epsilon = 1e-8
) {
  ret <- list(
    method         = "adagrad",
    stepsize       = stepsize,
    sgd_parameters = epsilon
  )
  class(ret) <- "ngme_optimization"
  ret
}


#' Root Mean Square Propagation (RMSProp) SGD optimization
#'
#' @details
#' The update rule for RMSProp is:
#'  v_t = beta1 * v_{t-1} + (1-beta1) * g_t^2
#'  x_{t+1} = x_t - stepsize * g_t / (sqrt(v_t) + epsilon)
#'
#' @param stepsize stepsize for SGD
#' @param beta1 beta1 for momentum
#' @param epsilon epsilon for numerical stability
#'
#' @return a list of control variables for optimization
#' (used in \code{control_opt} function)
#' @export
rmsprop <- function(
  stepsize = 0.05,
  beta1 = 0.9,
  epsilon = 1e-8
) {
  ret <- list(
    method         = "rmsprop",
    stepsize       = stepsize,
    sgd_parameters = c(beta1, epsilon)
  )
  class(ret) <- "ngme_optimization"
  ret
}

#' Adam SGD optimization
#'
#' @details
#' The update rule for Adam is:
#'  m_t = beta1 * m_{t-1} + (1 - beta1) * g_t
#'  v_t = beta2 * v_{t-1} + (1 - beta2) * g_t^2
#'  m_t_hat = m_t / (1 - beta1^t)
#'  v_t_hat = v_t / (1 - beta2^t)
#'  x_{t+1} = x_t - stepsize * m_t_hat / (sqrt(v_t_hat) + epsilon)
#' 
#' @param stepsize stepsize for SGD
#' @param beta1 beta1 for Adam
#' @param beta2 beta2 for Adam
#' @param epsilon epsilon for numerical stability
#'
#' @return a list of control variables for optimization 
#' (used in \code{control_opt} function)
#' @export
adam <- function(
  stepsize = 0.05,
  beta1 = 0.9,
  beta2 = 0.999,
  epsilon = 1e-8
) {
  ret <- list(
    method         = "adam",
    stepsize       = stepsize,
    sgd_parameters = c(beta1, beta2, epsilon)
  )
  class(ret) <- "ngme_optimization"
  ret
}

#' AdamW SGD optimization
#'
#' @details
#' The update rule for AdamW is:
#'  m_t = beta1 * m_{t-1} + (1 - beta1) * g_t
#'  v_t = beta2 * v_{t-1} + (1 - beta2) * g_t^2
#'  m_t_hat = m_t / (1 - beta1^t)
#'  v_t_hat = v_t / (1 - beta2^t)
#'  x_{t+1} = x_t - stepsize * (lambda * x_t + m_t_hat / (sqrt(v_t_hat) + epsilon))
#'
#' @param stepsize stepsize for SGD
#' @param beta1 beta1 for AdamW
#' @param beta2 beta2 for AdamW
#' @param lambda lambda (weight decay) for AdamW
#' @param epsilon epsilon for numerical stability
#'
#' @return a list of control variables for optimization
#' (used in \code{control_opt} function)
#' @export
adamW <- function(
  stepsize = 0.05,
  beta1 = 0.9,
  beta2 = 0.999,
  lambda = 0.01,
  epsilon = 1e-8
) {
  ret <- list(
    method         = "adamW",
    stepsize       = stepsize,
    sgd_parameters = c(beta1, beta2, lambda, epsilon)
  )
  class(ret) <- "ngme_optimization"
  ret
}