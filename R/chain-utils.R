# Internal helpers for chain-wise post-processing

strip_ngme_chain_attributes <- function(object) {
  attr(object, "chain_outputs") <- NULL
  attr(object, "chain_fits") <- NULL
  object
}


build_ngme_from_chain_output <- function(object, chain_output) {
  stopifnot(
    inherits(object, "ngme"),
    is.list(chain_output)
  )

  chain_fit <- object
  stopifnot(length(chain_fit$replicates) == length(chain_output))

  for (i in seq_along(chain_fit$replicates)) {
    chain_fit$replicates[[i]] <- update_ngme_est(
      chain_fit$replicates[[i]],
      chain_output[[i]]
    )
  }

  strip_ngme_chain_attributes(chain_fit)
}


get_ngme_chain_fits <- function(object) {
  stopifnot(inherits(object, "ngme"))

  chain_fits <- attr(object, "chain_fits")
  if (is.list(chain_fits) && length(chain_fits) > 0) {
    return(chain_fits)
  }

  chain_outputs <- attr(object, "chain_outputs")
  if (!is.list(chain_outputs) || length(chain_outputs) == 0) {
    return(list(object))
  }

  fits <- lapply(chain_outputs, function(chain_output) {
    build_ngme_from_chain_output(object, chain_output)
  })

  fits
}
