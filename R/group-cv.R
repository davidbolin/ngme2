#' Leave-group-out cross-validation
#'
#' Scores each fold against a posterior that never saw it, conditioning on the
#' fitted parameters as [cross_validation()] does, but reaching that posterior
#' by the cheapest route the model allows.
#'
#' It is meant as a drop-in for [cross_validation()]: it takes that function's
#' arguments under the same names, so an existing call can be switched over by
#' changing the function name alone. A few of them describe machinery this
#' function does not have and are either reported as ignored or refused outright.
#'
#' @details
#' `method = "auto"` picks between three strategies:
#'
#' * **Fully Gaussian fit** -- closed form. The leave-group-out posterior is
#'   available exactly from the factorization already computed, so no sampling
#'   happens at all and `n_gibbs_samples` is ignored.
#' * **Gaussian latent field, non-Gaussian measurement noise** -- one
#'   full-data Gibbs chain serves every fold. Conditional on the mixing
#'   variables the latent posterior is Gaussian with precision `QQ`, and
#'   dropping a group subtracts a rank-|I| term from it, so each fold costs a
#'   handful of triangular solves against the existing factor. Those draws
#'   target `p(. | y)` rather than `p(. | y_-I)`, and are corrected by
#'   Pareto-smoothed importance weights `1 / p(y_I | y_-I, V)`.
#' * **Anything else** (a non-Gaussian latent field) -- a real leave-group-out
#'   Gibbs chain per fold, run inside C++ and parallel across folds.
#'
#' Only the second route approximates anything, and its error is diagnosed per
#' group rather than assumed away.
#'
#' As |I| grows both the linear algebra advantage and the importance weights
#' degrade together, so `"auto"` hands groups larger than `max_lowrank_group`
#' to the per-fold chains instead.
#'
#' Like [cross_validation()], parameters are **not** re-estimated per fold:
#' everything conditions on the fitted values in `ngme`.
#'
#' @param ngme a fitted `ngme` object.
#' @param type how to form the folds, as in [cross_validation()]: `"loo"` (one
#'   group per observation), `"k-fold"`, `"lpo"` (leave-`percent`-out, repeated
#'   `times`), or `"custom"` to supply the design directly.
#' @param test_idx for `type = "custom"`, a list of integer vectors -- one per
#'   fold, holding the observations to be **scored**. As in
#'   [cross_validation()].
#' @param train_idx for `type = "custom"`, a list of the same length holding
#'   the observations each fold may **condition on**. Optional: left out, each
#'   fold trains on everything it does not score.
#' @param seed random seed. If `NULL`, drawn from the current R stream.
#' @param print logical; emit the advisory warnings (ignored arguments, a
#'   `min_ess` that will flag everything, a method that will degenerate).
#'   `FALSE` silences them, as it silences [cross_validation()]'s progress
#'   output. It does not affect what is returned.
#' @param N_sim accepted as a synonym for `n_chains`, which is the same
#'   quantity: `N_sim` independent runs of `n_gibbs_samples` draws each. An
#'   explicit `n_chains` takes precedence. It has an effect only where folds
#'   have their own chains -- see `n_chains`.
#' @param metric **not supported.** A custom scoring function; this function
#'   reports a fixed set of scores. Supplying it is an error rather than a
#'   warning, because quietly returning different scores than were asked for
#'   would be worse.
#' @param keep_pred logical; if `TRUE`, also return draws from each fold's
#'   predictive, as [cross_validation()] does. A shorthand for `n_draw_out`,
#'   which sets how many.
#' @param parallel,cores_layer1,cores_layer2,thining_gap,chain_combine
#'   **ignored.** Parallelism happens inside C++ across folds and replicates --
#'   see `num_threads` -- so there are no nested worker pools to configure, and
#'   there is one chain rather than several to thin or combine.
#' @param max_num_threads accepted as a synonym for `num_threads`.
#' @param merge_groups,merged_group_name,data **not supported.** Paired splits
#'   for grouped data; use [cross_validation()].
#' @param k number of folds for `type = "k-fold"`.
#' @param percent proportion held out per fold for `type = "lpo"`.
#' @param times number of repetitions for `type = "lpo"`.
#' @param transform optional function applied to the predictive and to the
#'   observation before scoring, so scores are reported on a different scale.
#'   It may take a second argument, the observation index, for a
#'   position-dependent transform. A named list of functions scores the same
#'   folds on each scale and returns one result per entry. `log_score` is only
#'   defined on the untransformed scale.
#' @param method force a strategy instead of dispatching on model class. See
#'   Details. `"lowrank"` on a non-Gaussian latent field warns, because that is
#'   where the weights degenerate.
#' @param max_lowrank_group largest group the low-rank route will accept under
#'   `method = "auto"` before deferring to per-fold chains.
#' @param n_score_draws draws used for CRPS and sCRPS when a `transform` is in
#'   force, where the scores have no closed form. Unused otherwise.
#' @param n_gibbs_samples number of posterior draws to retain. **Ignored for a
#'   fully Gaussian model**, where the answer is exact from the factorization
#'   and one draw suffices; honouring a large value there would cost orders of
#'   magnitude more for the same number. What governs the error elsewhere is
#'   not this count but the r_eff-corrected effective sample size it yields,
#'   which is a fraction of it because the draws come from a Markov chain.
#'   Check the `ess` column rather than assuming. The default is higher than
#'   [cross_validation()]'s, and `n_burnin` likewise: that function splits its
#'   budget over `N_sim` repeats of each fold, while here one chain carries
#'   all of it. This is the only default that differs between the two.
#' @param n_q number of prior draws of the held-out group's measurement mixing
#'   variable taken per retained posterior draw. The reported `mean` and `sd`
#'   do not depend on this where the measurement family has closed-form mixing
#'   moments (normal, NIG, GAL) -- there it is used only for the proper scores,
#'   which need the whole predictive mixture rather than its first two moments,
#'   and for the importance weights. Raise it when the measurement noise is
#'   strongly heavy tailed. Ignored for Gaussian measurement noise, where the
#'   mixing variable is the constant 1.
#' @param crps_max budget, in mixture-component pairs, for the `E|X - X'|`
#'   term behind CRPS and sCRPS. The exact double sum is taken whenever it fits
#'   inside the budget; past that the term is estimated from this many sampled
#'   pairs, which is unbiased and draws on every component. It is what limits
#'   the precision of those two scores; the other scores do not use it.
#' @param n_burnin number of burn-in sweeps.
#' @param num_threads number of OpenMP threads. How they are spent depends on
#'   the route: the low-rank one parallelises across replicates first and then
#'   across groups within a chain; the per-fold routes run replicates one after
#'   another and parallelise across **(fold, chain) pairs**.
#' @param chunk_cols bounds how many right-hand sides are handed to the solver
#'   at once, and so the peak memory of the dense solve. Larger is faster until
#'   the dense block stops fitting in cache.
#' @param n_chains run this many independent chains per fold, of
#'   `n_gibbs_samples` draws each, started from the stored per-chain states of
#'   the fit where the object kept them. Their between-chain spread is reported
#'   as `rhat`, and is a real convergence signal precisely because the starts
#'   are independent; folds that have not settled are warned about, never
#'   dropped.
#'
#'   This applies to the routes that give each fold its own chain. The
#'   low-rank route draws one full-data chain that every fold shares, so there
#'   is nothing per-fold to replicate: `n_chains` there is reported as ignored,
#'   and `n_gibbs_samples` is the knob that governs its error.
#' @param min_ess groups whose r_eff-corrected effective sample size falls below
#'   this are reported as unreliable. **This is the trigger that matters.**
#'   Pareto k diagnoses the tail of the weight distribution; what actually
#'   fails on the low-rank route is how slowly the chain mixes, and only the
#'   corrected ESS sees that. Groups where the estimator is badly biased can
#'   leave k untouched while their corrected ESS collapses.
#'
#'   A flagged group is a group to re-score exactly, not one to throw draws at.
#'   Where the weights are degenerate the effective sample size does not grow
#'   with `n_gibbs_samples`: more draws can shrink the spread while the
#'   estimate concentrates on the wrong value. Low variance and large bias
#'   together is the dangerous failure mode, and it looks like convergence.
#' @param khat_threshold groups whose Pareto k exceeds this are also reported as
#'   unreliable. 0.7 is the usual PSIS cutoff. Keep it, but do not rely on it
#'   alone -- see `min_ess`.
#' @param n_draw_out if positive, also return `n_draw_out` draws from each
#'   group's leave-group-out predictive, as a list of `|I| x n_draw_out` matrices
#'   in `attr(x, "draws")`. The predictive is a weighted mixture of Gaussians, so
#'   these are exact draws from it rather than a re-run of any sampler. Useful
#'   for plotting it against what [cross_validation()] simulates.
#' @param r_eff logical; correct the effective sample size for the
#'   autocorrelation of the chain. Leave this on: the weight-based ESS alone
#'   assumes independent draws and will call a group reliable when it is not.
#' @return An object of class `ngme_group_cv`: a data frame with one row per
#'   fold, holding the fold's identity (`group`, `size`, `n_held_out`), the
#'   predictive summary (`mean`, `sd`), the scores (`log_score`, `MAE`, `MSE`,
#'   `neg.CRPS`, `neg.sCRPS` -- all lower-is-better despite the `neg.`, which
#'   is inherited from [cross_validation()]), and the diagnostics (`khat`,
#'   `ess`, `r_eff`, `rhat`, `exact`, and the `low_ess` / `high_khat` /
#'   `unreliable` flags). A `summary` attribute carries the score means, and
#'   `print()` reports which folds are untrustworthy and why.
#'
#'   When `transform` is a named list, a list of such objects is returned, one
#'   per scale.
#' @seealso [cross_validation()], which re-scores by refitting each fold and
#'   covers models this does not (correlated measurement noise).
#' @export
group_cv <- function(
    ngme,
    type = c("k-fold", "loo", "lpo", "custom"),
    seed = NULL,
    print = TRUE,
    N_sim = 5,
    n_gibbs_samples = 1000,
    n_burnin = 200,
    k = 5,
    percent = 0.2,
    times = 10,
    metric = NULL,
    transform = identity,
    test_idx = NULL,
    train_idx = NULL,
    keep_pred = FALSE,
    parallel = FALSE,
    thining_gap = 1,
    max_num_threads = NULL,
    cores_layer1 = 1,
    cores_layer2 = 1,
    merge_groups = FALSE,
    merged_group_name = NULL,
    data = NULL,
    chain_combine = c("param_mean", "predictive_average"),
    # -- specific to this function ------------------------------------------
    method = c("auto", "lowrank", "exact"),
    max_lowrank_group = 30L,
    n_score_draws = 4000L,
    n_q = 8L,
    crps_max = 200000L,
    num_threads = 1,
    chunk_cols = 256,
    n_chains = 1L,
    min_ess = 100,
    khat_threshold = 0.7,
    n_draw_out = 0L,
    r_eff = TRUE) {
  type <- match.arg(type)
  method <- match.arg(method)
  if (is.null(seed)) seed <- ngme_random_seed()

  # Drop-in compatibility with cross_validation(): take its arguments under
  # its own names, honour the ones that have a counterpart here, and say
  # plainly which do not rather than ignoring them behind the user's back.
  if (isTRUE(merge_groups))
    stop("group_cv() does not support `merge_groups`; use cross_validation().",
         call. = FALSE)
  if (!is.null(metric))
    stop("group_cv() reports a fixed set of scores and does not take a custom ",
         "`metric`; use cross_validation().", call. = FALSE)
  if (!is.null(max_num_threads)) num_threads <- max_num_threads
  # N_sim is cross_validation()'s replication count: N_sim runs of
  # n_gibbs_samples draws each. `n_chains` is the same quantity here -- n_chains
  # independent chains of n_gibbs_samples each -- so honour it under that name
  # rather than dropping it. An explicit n_chains wins.
  if (!missing(N_sim) && missing(n_chains)) n_chains <- as.integer(N_sim)
  if (isTRUE(keep_pred) && n_draw_out < 1L) n_draw_out <- 1000L
  ignored <- c(

    if (isTRUE(parallel)) "parallel",
    if (!missing(thining_gap)) "thining_gap",
    if (!missing(cores_layer1)) "cores_layer1",
    if (!missing(cores_layer2)) "cores_layer2",
    if (!is.null(merged_group_name)) "merged_group_name",
    if (!is.null(data)) "data",
    if (!missing(chain_combine)) "chain_combine")
  if (length(ignored) && isTRUE(print))
    warning("group_cv() ignores ", paste(ignored, collapse = ", "),
            ": it draws once and parallelises inside C++ (see `num_threads`), ",
            "so there are no repeated simulations or nested worker pools to ",
            "configure.", call. = FALSE)

  # Several scales: score the same folds on each, as cross_validation() does.
  if (is.list(transform) && !is.null(names(transform)) &&
      all(nzchar(names(transform))) &&
      all(vapply(transform, is.function, logical(1)))) {
    # Recurse explicitly rather than re-evaluating match.call(): inside lapply,
    # parent.frame() is the lapply frame, so a call naming the fit by symbol
    # cannot be resolved there.
    return(structure(stats::setNames(lapply(transform, function(tf)
      group_cv(ngme = ngme, type = type, test_idx = test_idx,
        train_idx = train_idx, k = k, percent = percent, times = times,
        transform = tf, print = print, keep_pred = keep_pred,
        method = method, max_lowrank_group = max_lowrank_group,
        n_score_draws = n_score_draws, n_q = n_q, crps_max = crps_max,
        n_gibbs_samples = n_gibbs_samples,
        n_burnin = n_burnin, seed = seed, num_threads = num_threads,
        chunk_cols = chunk_cols, n_chains = n_chains, min_ess = min_ess,
        khat_threshold = khat_threshold, n_draw_out = n_draw_out,
        r_eff = r_eff)),
      names(transform)), class = "ngme_group_cv_list"))
  }

  # Several models: the SAME folds for each, or the comparison is meaningless.
  if (!inherits(ngme, "ngme")) {
    if (!is.list(ngme) || !all(vapply(ngme, inherits, logical(1), "ngme"))) {
      stop("`ngme` must be a fitted ngme object or a list of them.")
    }
    if (is.null(names(ngme)))
      names(ngme) <- paste0("model_", seq_along(ngme))
    n_shared <- attr(ngme[[1]], "fit")$n_data
    shared <- resolve_cv_type(type, test_idx, train_idx, n_shared,
                              k, percent, times, seed)
    # Pass the drawn folds on as an explicit custom design, so a randomised
    # type is not re-drawn per model.
    shared_train <- lapply(shared$hold, function(h) setdiff(seq_len(n_shared), h))
    return(structure(stats::setNames(lapply(ngme, function(m)
      group_cv(ngme = m, type = "custom", test_idx = shared$test,
               train_idx = shared_train,
               transform = transform, print = print, keep_pred = keep_pred,
               k = k, percent = percent, times = times,
               method = method, max_lowrank_group = max_lowrank_group,
               n_score_draws = n_score_draws, n_q = n_q, crps_max = crps_max,
               n_gibbs_samples = n_gibbs_samples, n_burnin = n_burnin,
               seed = seed, num_threads = num_threads,
               chunk_cols = chunk_cols, n_chains = n_chains, min_ess = min_ess,
               khat_threshold = khat_threshold, n_draw_out = n_draw_out,
               r_eff = r_eff)),
      names(ngme)), class = "ngme_group_cv_list"))
  }
  stopifnot("Not a ngme object." = inherits(ngme, "ngme"))
  seed_int <- as.integer(abs(seed) %% 2147483647)

  reps <- ngme$replicates

  # ---- pick the method from the model class ----
  # Three regimes, each with exactly one method that suits it:
  #   fully Gaussian            : QQ constant, weights uniform -> exact, 1 draw
  #   Gaussian latent + non-normal noise : weights well behaved (ESS 2700-3900
  #                               of 4000 measured) -> low-rank, no gate needed
  #   non-Gaussian latent field : weights degenerate at influential points
  #                               (ESS -> 0, and MORE draws make it worse), so
  #                               run the exact per-fold chains instead -- in
  #                               C++ and across threads, which is where that
  #                               estimator's speedup actually comes from.
  latent_gaussian <- length(reps) > 0 && all(vapply(reps, function(r) {
    length(r$models) == 0 ||
      all(vapply(r$models, function(m) all(m$noise$noise_type == "normal"),
                 logical(1)))
  }, logical(1)))

  # A fully Gaussian fit has no mixing variables: QQ is identical across draws,
  # so one draw is already the exact answer and every further draw is waste --
  # orders of magnitude of runtime for a result that agrees to rounding.
  all_gauss <- length(reps) > 0 &&
    all(vapply(reps, function(r) isTRUE(r$all_gaussian), logical(1))) &&
    all(vapply(reps, function(r) identical(r$noise$noise_type, "normal"), logical(1)))
  if (all_gauss && (n_gibbs_samples > 1L || n_burnin > 0L)) {
    n_gibbs_samples <- 1L
    n_burnin <- 0L
  }
  n_data <- attr(ngme, "fit")$n_data
  if (is.null(n_data)) stop("`ngme` does not carry fit information.")

  # Groups are given in the ORIGINAL data ordering; each replicate holds its own
  # slice of the data, so they have to be translated per replicate and groups
  # that straddle two replicates rejected rather than silently split.
  if (min_ess > n_gibbs_samples / 2 && !isTRUE(reps[[1]]$all_gaussian)) {
    if (isTRUE(print)) warning("min_ess (", min_ess, ") is more than half of n_gibbs_samples (",
            n_gibbs_samples, "). The effective sample size is bounded by the ",
            "draw count and is typically a small fraction of it, so most ",
            "groups will be flagged. Raise n_gibbs_samples instead.",
            call. = FALSE)
  }

  folds <- resolve_cv_type(type, test_idx, train_idx, n_data, k, percent,
                           times, seed)
  groups <- folds$test      # scored
  hold <- folds$hold        # conditioned away: the complement of train_idx
  per_rep <- map_groups_to_replicates(ngme, hold, groups)
  max_k <- max(lengths(hold))

  if (method == "auto") {
    # The low-rank path is exact for a fully Gaussian fit at any group size --
    # there are no weights to degrade. Otherwise its weights get worse as the
    # group grows, and past the measured break-even rank (~31 solve columns per
    # factorization) it is not even cheaper, so hand large groups to the exact
    # chains.
    method <- if (all_gauss) "lowrank"
      else if (latent_gaussian && max_k <= max_lowrank_group) "lowrank"
      else "exact"
  }
  # The low-rank route draws ONE full-data chain that every fold shares, so
  # there is nothing per-fold to replicate; saying so beats dropping it.
  if (method == "lowrank" && n_chains > 1L && isTRUE(print)) {
    warning("n_chains (or N_sim) is ignored by method = \"lowrank\": one ",
            "full-data chain serves every fold, so there are no per-fold ",
            "chains to replicate. Raise n_gibbs_samples, or use ",
            "method = \"exact\" to run independent chains per fold.",
            call. = FALSE)
  }
  if (method == "lowrank" && !latent_gaussian) {
    if (isTRUE(print)) warning("method = \"lowrank\" with a non-Gaussian latent field: the ",
            "importance weights degenerate at influential groups and more ",
            "draws do not help. Check the `ess` column.", call. = FALSE)
  }

  # Per-chain starting states from the fit, when it kept them: independent
  # starts are what make the between-chain spread a real convergence signal.
  chain_starts <- list()
  if (method == "exact" && n_chains > 1L) {
    cf <- try(resolve_ngme_chain_fits(ngme), silent = TRUE)
    if (!inherits(cf, "try-error") && length(cf) > 1L) {
      chain_starts <- lapply(seq_along(reps), function(i)
        lapply(cf, function(z) {
          w <- lapply(z$replicates[[i]]$models, function(mm) mm$W)
          if (any(vapply(w, is.null, TRUE))) numeric(0) else unlist(w)
        }))
    }
  }

  res <- if (method == "lowrank") {
    group_cv_cpp(
      ngme_replicates = lapply(seq_along(reps), function(i) reps[[i]]),
      groups_per_rep = lapply(per_rep, `[[`, "local"),
      n = as.integer(n_gibbs_samples),
      n_burnin = as.integer(n_burnin),
      seed = seed_int,
      num_threads = as.integer(num_threads),
      chunk_cols = as.integer(chunk_cols)
    )
  } else {
    group_cv_exact_cpp(
      ngme_replicates = lapply(seq_along(reps), function(i) reps[[i]]),
      groups_per_rep = lapply(per_rep, `[[`, "local"),
      n = as.integer(n_gibbs_samples),
      n_burnin = as.integer(n_burnin),
      seed = seed_int,
      num_threads = as.integer(num_threads),
      n_chains = as.integer(n_chains),
      chain_starts = chain_starts
    )
  }

  scored <- list()
  for (i in seq_along(reps)) {
    if (length(per_rep[[i]]$local) == 0) next
    scored[[length(scored) + 1]] <- if (method == "lowrank") {
      score_group_cv_replicate(
        rep_model = reps[[i]],
        draws = res[[i]],
        local_groups = per_rep[[i]]$local,
        global_ids = per_rep[[i]]$global_id,
        score_pos = per_rep[[i]]$score_pos,
        score_global = per_rep[[i]]$score_global,
        use_r_eff = r_eff,
        n_draw_out = as.integer(n_draw_out),
        transform = transform,
        n_score_draws = as.integer(n_score_draws),
        n_q = as.integer(n_q),
        crps_max = as.integer(crps_max)
      )
    } else {
      score_group_cv_exact(
        rep_model = reps[[i]],
        draws = res[[i]],
        local_groups = per_rep[[i]]$local,
        global_ids = per_rep[[i]]$global_id,
        score_pos = per_rep[[i]]$score_pos,
        score_global = per_rep[[i]]$score_global,
        n_draw_out = as.integer(n_draw_out),
        transform = transform,
        n_score_draws = as.integer(n_score_draws),
        n_q = as.integer(n_q),
        crps_max = as.integer(crps_max),
        n_chains = as.integer(n_chains)
      )
    }
  }
  if (length(scored) == 0) stop("No group fell inside a replicate.")

  draws_out <- NULL
  if (n_draw_out > 0L) {
    draws_out <- do.call(c, lapply(scored, attr, "draws"))
  }

  out <- do.call(rbind, scored)
  ord <- order(out$group)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  if (!is.null(draws_out)) draws_out <- draws_out[ord]

  # An exact (fully Gaussian) fit is never unreliable, whatever the draw count.
  out$low_ess <- !out$exact & out$ess < min_ess
  out$high_khat <- !out$exact & !is.na(out$khat) & out$khat > khat_threshold
  out$unreliable <- out$low_ess | out$high_khat

  structure(
    out,
    class = c("ngme_group_cv", "data.frame"),
    summary = colMeans(out[, c("MAE", "MSE", "neg.CRPS", "neg.sCRPS",
                               "log_score")], na.rm = TRUE),
    n_unreliable = sum(out$unreliable),
    n_low_ess = sum(out$low_ess),
    n_high_khat = sum(out$high_khat),
    min_ess = min_ess,
    khat_threshold = khat_threshold,
    n_draws = n_gibbs_samples,
    method = method,
    transform = !identical(transform, identity),
    draws = draws_out
  )
}


resolve_cv_groups <- function(groups, n_data) {
  if (identical(groups, "loo")) {
    return(lapply(seq_len(n_data), function(i) i))
  }
  if (!is.list(groups)) {
    stop("`groups` must be \"loo\" or a list of integer vectors.")
  }
  groups <- lapply(groups, function(g) {
    g <- as.integer(g)
    if (anyNA(g) || any(g < 1L) || any(g > n_data)) {
      stop("`groups` contains indices outside 1:", n_data, ".")
    }
    if (anyDuplicated(g)) stop("`groups` contains a repeated index within a group.")
    sort(g)
  })
  if (any(lengths(groups) == 0)) stop("`groups` contains an empty group.")
  groups
}


# Translate globally-indexed groups into each replicate's own row numbering.
# `hold` is what is removed from the likelihood; `score` is the subset of it
# that is actually scored. They differ whenever a buffer is used: INLA's group
# CV conditions on y_{-G_i} but scores only y_i, and a spatial leave-one-out is
# only meaningful that way, since otherwise the neighbours give the answer away.
map_groups_to_replicates <- function(ngme, hold, score) {
  lapply(seq_along(ngme$replicates), function(i) {
    data_idx <- ngme$replicates[[i]]$data_idx
    local <- list(); score_pos <- list(); score_global <- list()
    global_id <- integer()
    for (g in seq_along(hold)) {
      hit <- match(hold[[g]], data_idx)
      if (all(is.na(hit))) next
      if (anyNA(hit)) {
        stop("Group ", g, " spans more than one replicate. Leave-group-out ",
             "conditions on a single replicate's latent field, so split the ",
             "group or use cross_validation().")
      }
      sp <- match(score[[g]], hold[[g]])
      if (anyNA(sp)) {
        stop("Group ", g, ": every scored index must also be held out.")
      }
      local[[length(local) + 1]] <- as.integer(hit)
      score_pos[[length(score_pos) + 1]] <- as.integer(sp)
      score_global[[length(score_global) + 1]] <- as.integer(score[[g]])
      global_id <- c(global_id, g)
    }
    list(local = local, score_pos = score_pos, score_global = score_global,
         global_id = global_id)
  })
}


# Turn a cross_validation-style `type` into an explicit list of held-out groups.
# Folds in the same terms cross_validation() uses: `test_idx` is what gets
# scored, `train_idx` is what the posterior may condition on. Anything in
# neither is held out without being scored -- the exclusion buffer, which needs
# no argument of its own because the complement of `train_idx` already says it.
resolve_cv_type <- function(type, test_idx, train_idx, n_data, k, percent,
                            times, seed) {
  idx <- seq_len(n_data)
  if (identical(type, "custom") || !is.null(test_idx) || !is.null(train_idx)) {
    if (is.null(test_idx))
      stop("type = \"custom\" requires `test_idx`.", call. = FALSE)
    if (!is.list(test_idx))
      stop("`test_idx` must be a list.", call. = FALSE)
    # A design with no buffer trains on everything it does not score, so
    # train_idx carries no information and may be left out. Supplying it is
    # what creates a buffer, exactly as in cross_validation().
    if (is.null(train_idx))
      train_idx <- lapply(test_idx, function(te) setdiff(idx, as.integer(te)))
    if (!is.list(train_idx))
      stop("`train_idx` must be a list.", call. = FALSE)
    if (length(test_idx) != length(train_idx))
      stop("`test_idx` and `train_idx` must have the same length (",
           length(test_idx), " vs ", length(train_idx), ").", call. = FALSE)
    te <- resolve_cv_groups(test_idx, n_data)
    hold <- lapply(seq_along(te), function(j) {
      tr <- as.integer(train_idx[[j]])
      if (length(tr) && (anyNA(tr) || any(tr < 1L) || any(tr > n_data)))
        stop("`train_idx[[", j, "]]` contains indices outside 1:", n_data, ".",
             call. = FALSE)
      h <- sort(setdiff(idx, tr))
      if (!all(te[[j]] %in% h))
        stop("`test_idx[[", j, "]]` overlaps `train_idx[[", j, "]]`: a scored ",
             "observation cannot also be conditioned on.", call. = FALSE)
      h
    })
    return(list(test = te, hold = hold))
  }
  te <- switch(type,
    loo = lapply(idx, function(i) i),
    "k-fold" = {
      stopifnot("k must be at least 2" = k >= 2)
      withr_seed(seed, {
        fold <- cut(sample(idx), breaks = k, labels = FALSE)
        lapply(seq_len(k), function(j) idx[fold == j])
      })
    },
    lpo = {
      stopifnot("percent must be in (0, 1)" = percent > 0 && percent < 1)
      m <- max(1L, round(percent * n_data))
      withr_seed(seed, lapply(seq_len(times), function(j) sort(sample(idx, m))))
    })
  # The built-in designs train on everything they do not score, so there is no
  # buffer and the held-out set is the scored set.
  list(test = te, hold = te)
}

# Draw folds reproducibly without disturbing the caller's RNG stream.
withr_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  set.seed(seed)
  on.exit({
    if (is.null(old)) rm(".Random.seed", envir = .GlobalEnv)
    else assign(".Random.seed", old, envir = .GlobalEnv)
  }, add = TRUE)
  force(expr)
}


#' @export
print.ngme_group_cv_list <- function(x, ...) {
  cat("Leave-group-out cross-validation for", length(x), "models\n\n")
  tab <- do.call(rbind, lapply(x, function(z) attr(z, "summary")))
  rownames(tab) <- names(x)
  print(round(tab, 5))
  bad <- vapply(x, function(z) attr(z, "n_unreliable"), numeric(1))
  if (any(bad > 0))
    cat("\nunreliable groups per model:",
        paste(names(x), bad, sep = "=", collapse = "  "), "\n")
  invisible(x)
}
