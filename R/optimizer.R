# different methods for SGD
# sgd (vanilla)
# preconditioner_sgd
# momentum
# adagrad
# rmsprop
# adam
# adamW

#' Preconditioner SGD optimization
#'
#' @details
#' Preconditioner SGD is using fisher information matrix as preconditioner (natural gradient descent).
#'
#' @param stepsize stepsize for SGD
#' @param numerical_eps   numerical, the gap used for estimate preconditioner, default is 1e-5
#'
#' @return a list of control variables for optimization
#' (used in \code{control_opt} function)
#' @export
precond_sgd <- function(
    stepsize = 1,
    numerical_eps = 1e-5
    # precond_by_diff_chain = FALSE
    ) {
  ret <- list(
    # sgd related
    method = "precond_sgd",
    stepsize = stepsize,
    sgd_parameters = NULL,
    numerical_eps = numerical_eps
  )
  class(ret) <- "ngme_optimizer"
  ret
}

#' Vanilla SGD optimization
#'
#' Simple stochastic gradient descent without momentum or preconditioning.
#'
#' @param stepsize stepsize for SGD
#'
#' @return a list of control variables for optimization
#' (used in \code{control_opt} function)
#' @export
sgd <- function(
    stepsize = 0.001) {
  ret <- list(
    method         = "sgd",
    stepsize       = stepsize,
    sgd_parameters = double(0)
  )
  class(ret) <- "ngme_optimizer"
  ret
}

#' Stochastic Gradient Langevin Dynamics (SGLD) optimization
#'
#' @details
#' SGLD adds Gaussian noise to vanilla SGD updates:
#' \deqn{x_{t+1} = x_t - \eta_t \nabla U(x_t) + \sqrt{2 T \eta_t}\,\xi_t, \quad \xi_t \sim \mathcal{N}(0, I)}
#' where \eqn{T} is \code{temperature}. The implementation applies this
#' component-wise using the current effective stepsize.
#'
#' @param stepsize base stepsize
#' @param temperature non-negative Langevin temperature
#'
#' @return a list of control variables for optimization
#' (used in \code{control_opt} function)
#' @export
sgld <- function(
    stepsize = 0.001,
    temperature = 1) {
  stopifnot(
    is.numeric(stepsize),
    length(stepsize) == 1,
    is.finite(stepsize),
    stepsize > 0,
    is.numeric(temperature),
    length(temperature) == 1,
    is.finite(temperature),
    temperature >= 0
  )

  ret <- list(
    method         = "sgld",
    stepsize       = stepsize,
    sgd_parameters = temperature
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
    beta2 = 1 - beta1) {
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
    epsilon = 1e-8) {
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
    epsilon = 1e-8) {
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
    epsilon = 1e-8) {
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
    epsilon = 1e-8) {
  ret <- list(
    method         = "adamW",
    stepsize       = stepsize,
    sgd_parameters = c(beta1, beta2, lambda, epsilon)
  )
  class(ret) <- "ngme_optimizer"
  ret
}

#' Stepsize decay schedule
#'
#' @details
#' Define a stepsize decay strategy for use in \code{control_opt()}.
#' Currently supports the plateau-based schedule using grad.norm().
#'
#' @param method decay strategy. "none" disables decay; "grad_norm_plateau"
#'   decays when mean grad.norm() across chains fails to decrease for a number
#'   of epochs.
#' @param patience number of consecutive epochs without grad.norm() improvement
#'   before decaying stepsize.
#' @param gamma decay factor applied when triggered (0 < gamma < 1).
#' @param min_delta minimum required decrease in grad.norm() to be counted as improvement.
#' @param warmup number of initial epochs to skip decay checks.
#' @param min_stepsize lower bound for stepsize after decay (absolute value).
#'
#' @return a list of control variables for stepsize decay (used in \code{control_opt}).
#' @export
stepsize_decay <- function(
    method = c("none", "grad_norm_plateau"),
    patience = 3,
    gamma = 0.5,
    min_delta = 0,
    warmup = 0,
    min_stepsize = 0) {
  method <- match.arg(method)
  ret <- list(
    method = method,
    patience = patience,
    gamma = gamma,
    min_delta = min_delta,
    warmup = warmup,
    min_stepsize = min_stepsize
  )
  class(ret) <- "ngme_stepsize_decay"
  ret
}

#' Stepsize schedule
#'
#' @details
#' Define an iteration-dependent stepsize schedule for use in \code{control_opt()}.
#' This schedule is independent of \code{stepsize_decay()} and applies a direct
#' multiplicative scaling by iteration index:
#' \deqn{\eta_t = \eta_0 (t + t_0)^{-\alpha}}.
#'
#' @param method schedule type. "constant" keeps \eqn{\eta_t = \eta_0} and
#'   "poly" uses polynomial decay.
#' @param alpha polynomial exponent used when \code{method = "poly"}.
#'   Must satisfy \eqn{1/2 < \alpha < 1}.
#' @param t0 non-negative offset in iteration index.
#' @param burnin_iter non-negative integer. During these initial iterations,
#'   schedule scaling is fixed to 1. Afterward, polynomial decay is applied
#'   with reset local time index.
#'
#' @return a list of control variables for stepsize schedule
#' (used in \code{control_opt}).
#' @export
stepsize_schedule <- function(
    method = c("constant", "poly"),
    alpha = 0.501,
    t0 = 1,
    burnin_iter = 0) {
  method <- match.arg(method)

  stopifnot(
    is.numeric(alpha),
    length(alpha) == 1,
    is.numeric(t0),
    length(t0) == 1,
    t0 >= 0,
    is.numeric(burnin_iter),
    length(burnin_iter) == 1,
    is.finite(burnin_iter),
    burnin_iter >= 0,
    abs(burnin_iter - round(burnin_iter)) < sqrt(.Machine$double.eps)
  )
  if (method == "poly") {
    stopifnot(alpha > 0.5, alpha < 1)
  }

  ret <- list(
    method = method,
    alpha = alpha,
    t0 = t0,
    burnin_iter = as.integer(round(burnin_iter))
  )
  class(ret) <- "ngme_stepsize_schedule"
  ret
}

.default_stepsize_control <- function() {
  ret <- list(
    schedule = stepsize_schedule(method = "constant"),
    decay = stepsize_decay(method = "none")
  )
  class(ret) <- "ngme_stepsize_control"
  ret
}

#' Unified stepsize control
#'
#' @details
#' Bundle an iteration schedule (\code{stepsize_schedule}) and an optional
#' checkpoint-based decay rule (\code{stepsize_decay}) into one object for
#' \code{control_opt()}.
#'
#' @param schedule schedule component. Either an object from
#'   \code{stepsize_schedule()} or a method string accepted by
#'   \code{stepsize_schedule()}.
#' @param decay decay component. Either an object from
#'   \code{stepsize_decay()} or a method string accepted by
#'   \code{stepsize_decay()}.
#'
#' @return a bundled stepsize-control object for \code{control_opt()}.
#' @export
stepsize_control <- function(
    schedule = stepsize_schedule(method = "constant"),
    decay = stepsize_decay(method = "none")) {
  if (is.character(schedule)) {
    stopifnot(length(schedule) == 1)
    schedule <- stepsize_schedule(method = schedule)
  }
  if (is.character(decay)) {
    stopifnot(length(decay) == 1)
    decay <- stepsize_decay(method = decay)
  }

  stopifnot(
    inherits(schedule, "ngme_stepsize_schedule"),
    inherits(decay, "ngme_stepsize_decay")
  )

  ret <- list(
    schedule = schedule,
    decay = decay
  )
  class(ret) <- "ngme_stepsize_control"
  ret
}

#' Polynomial schedule helper
#'
#' @details
#' Convenience helper that enables polynomial schedule and disables
#' checkpoint-based decay.
#'
#' @param alpha polynomial exponent in \eqn{(1/2,1)}.
#' @param t0 non-negative schedule offset.
#' @param burnin_iter non-negative integer. Initial iterations without
#'   polynomial schedule scaling.
#'
#' @return a \code{stepsize_control()} object.
#' @export
poly_decay <- function(alpha = 0.501, t0 = 1, burnin_iter = 0) {
  stepsize_control(
    schedule = stepsize_schedule(
      method = "poly",
      alpha = alpha,
      t0 = t0,
      burnin_iter = burnin_iter
    ),
    decay = stepsize_decay(method = "none")
  )
}

#' Batch/checkpoint decay helper
#'
#' @details
#' Convenience helper that enables grad-norm plateau decay and keeps
#' iteration schedule constant.
#'
#' @param patience number of consecutive checkpoints without improvement before decay.
#' @param gamma decay factor in \eqn{(0,1)}.
#' @param min_delta minimum required decrease to count as improvement.
#' @param warmup number of initial checkpoints to skip.
#' @param min_stepsize lower bound for absolute stepsize.
#'
#' @return a \code{stepsize_control()} object.
#' @export
batch_decay <- function(
    patience = 3,
    gamma = 0.5,
    min_delta = 0,
    warmup = 0,
    min_stepsize = 0) {
  stepsize_control(
    schedule = stepsize_schedule(method = "constant"),
    decay = stepsize_decay(
      method = "grad_norm_plateau",
      patience = patience,
      gamma = gamma,
      min_delta = min_delta,
      warmup = warmup,
      min_stepsize = min_stepsize
    )
  )
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
    line_search = "wolfe") {
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
#'
#' @description
#' Adaptive gradient descent optimizer.
#'
#' @details
#' Based on the method described in
#' \url{https://arxiv.org/pdf/1910.09529}.
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
    stepsize = 0.01) {
  ret <- list(
    method         = "adaptive_gd",
    stepsize       = stepsize,
    sgd_parameters = double(0)
  )
  class(ret) <- "ngme_optimizer"
  ret
}

#' List supported optimizers
#'
#' @description
#' This function returns a list of supported optimizers in the \code{ngme} package.
#' The optimizers are categorized into three groups:
#' \itemize{
#' \item \strong{Gradient descent}: \code{"sgd"}, \code{"sgld"}, \code{"momentum"}, \code{"adaptive_gd"}
#' \item \strong{Adaptive learning rate}: \code{"adagrad"}, \code{"rmsprop"}, \code{"adam"}, \code{"adamW"}
#' \item \strong{Preconditioner}: \code{"precond_sgd"}, \code{"bfgs"}
#' }
#'
#' @return a character vector of supported optimizers
#' @export
ngme_optimizers <- function() {
  c("sgd", "sgld", "precond_sgd", "momentum", "adagrad", "rmsprop", "adam", "adamW", "bfgs", "adaptive_gd")
}
