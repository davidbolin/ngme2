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
  class(ret) <- "ngme_optimizer"
  ret
}

#' Momentum SGD optimization
#'
#' @details
#' The update rule for momentum is:
#' \deqn{v_t = \beta_1 v_{t-1} + \beta_2 g_t}
#' \deqn{x_{t+1} = x_t - \text{stepsize} * v_t}
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
  class(ret) <- "ngme_optimizer"
  ret
}


#' AdaGrad SGD optimization
#'
#' @details
#' The update rule for AdaGrad is:
#' \deqn{v_t = v_{t-1} + g_t^2}
#' \deqn{x_{t+1} = x_t - \text{stepsize} * \frac{g_t}{\sqrt{v_t} + \epsilon}}
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
  class(ret) <- "ngme_optimizer"
  ret
}


#' Root Mean Square Propagation (RMSProp) SGD optimization
#'
#' @details
#' The update rule for RMSProp is:
#' \deqn{v_t = \beta_1 v_{t-1} + (1 - \beta_1) g_t^2}
#' \deqn{x_{t+1} = x_t - \text{stepsize} * \frac{g_t}{\sqrt{v_t} + \epsilon}}
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
  class(ret) <- "ngme_optimizer"
  ret
}

#' Adam SGD optimization
#'
#' @details
#' The update rule for Adam is:
#' \deqn{m_t = \beta_1 m_{t-1} + (1 - \beta_1) g_t}
#' \deqn{v_t = \beta_2 v_{t-1} + (1 - \beta_2) g_t^2}
#' \deqn{\hat{m_t} = m_t / (1 - \beta_1^t)}
#' \deqn{\hat{v_t} = v_t / (1 - \beta_2^t)}
#' \deqn{x_{t+1} = x_t - \text{stepsize} * \frac{\hat{m_t}}{\sqrt{\hat{v_t}} + \epsilon}}
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
  class(ret) <- "ngme_optimizer"
  ret
}

#' AdamW SGD optimization
#'
#' @details
#' The update rule for AdamW is:
#' \deqn{m_t = \beta_1 m_{t-1} + (1 - \beta_1) g_t}
#' \deqn{v_t = \beta_2 v_{t-1} + (1 - \beta_2) g_t^2}
#' \deqn{\hat{m_t} = m_t / (1 - \beta_1^t)}
#' \deqn{\hat{v_t} = v_t / (1 - \beta_2^t)}
#' \deqn{x_{t+1} = x_t - \text{stepsize} * \left( \lambda x_t + \frac{\hat{m_t}}{\sqrt{\hat{v_t}} + \epsilon} \right)}
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
  class(ret) <- "ngme_optimizer"
  ret
}