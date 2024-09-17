# different methods for SGD
# preconditioner_sgd
# momentum
# adagrad
# rmsprop
# adam
# adamW

#' Preconditioner SGD optimization
#'
#' @param stepsize stepsize for SGD
#' @param preconditioner  preconditioner, can be c("none", "fast", "full")
#' "none" means no preconditioner, i.e., vanilla SGD,
#' "full" means compute numerical hessian as preconditioner
#' "fast" means compute numerical hessian except for the parameter of K matrix (for speed reason) 
#' @param numerical_eps   numerical, the gap used for estimate preconditioner, default is 1e-5
#' @param precond_by_diff_chain logical, if TRUE, use different chains to estimate preconditioner (only computed at check points), if FALSE, use the same chain to estimate preconditioner (computed at each iteration)
#' @param compute_precond_each_iter logical, if TRUE, compute preconditioner at each iteration, if FALSE, only compute preconditioner at check points (if has only 1 chain running, it will be set TRUE)
#'
#' @return a list of control variables for optimization 
#' (used in \code{control_opt} function)
#' @export
precond_sgd <- function(
  stepsize = 0.5,
  # preconditioner related
  preconditioner    = "fast",
  numerical_eps       = 1e-5,
  precond_by_diff_chain = FALSE,
  compute_precond_each_iter = FALSE
) {

  ret <- list(
    # sgd related
    method         = "precond_sgd",
    stepsize       = stepsize,
    sgd_parameters = NULL,
    # preconditioner related
    preconditioner  = preconditioner,
    numerical_eps     = numerical_eps,
    precond_by_diff_chain = precond_by_diff_chain,
    compute_precond_each_iter = compute_precond_each_iter
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


#' BFGS optimization
#'
#' @param line_search line search method, can be c("backtracking", "wofle")
#' "bfgs" means use BFGS preconditioner
# ' The update rule for BFGS is:
# ' \deqn{H_{t+1} = H_t + \frac{y_t y_t^T}{y_t^T s_t} - \frac{H_t s_t s_t^T H_t}{s_t^T H_t s_t}}
# ' \deqn{x_{t+1} = x_t - H_{t+1} g_t}
#' @return a list of control variables for optimization
#' (used in \code{control_opt} function)
#' @export
bfgs <- function(
  line_search = "wolfe"
) {
  stopifnot(line_search %in% c("backtracking", "wolfe"))

  ret <- list(
    method         = "bfgs",
    stepsize       = 1,
    line_search    = line_search,
    sgd_parameters = NULL
  )
  class(ret) <- "ngme_optimizer"
  ret
}

#' Adaptive gradient descent
#' From the paper: https://arxiv.org/pdf/1910.09529
#' The update rule for adaptive gradient descent is:
#' \deqn{\lambda_k = \min(\sqrt{1 + \theta_{k-1}} \lambda_{k-1}, \frac{||x_k - x_{k-1}||}{2 ||\nabla f(x_k) - \nabla f(x_{k-1})||} )}
#' \deqn{x_{k+1} = x_k - \lambda_k \nabla f(x_k)}
#' \deqn{\theta_k = \lambda_k / \lambda_{k-1}}
#'  
#' @param stepsize initial stepsize for SGD
#'
#' @return a list of control variables for optimization
#' (used in \code{control_opt} function)
#' @export
adaptive_gd <- function(
  stepsize = 0.01
) {
  ret <- list(
    method         = "adaptive_gd",
    stepsize       = stepsize,
    sgd_parameters = double(0)
  )
  class(ret) <- "ngme_optimizer"
  ret
}