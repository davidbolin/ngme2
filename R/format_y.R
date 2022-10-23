#' Format the input data
#'
#' Function to reformat the input data, Y and supply to ngme2.R and model.R
#' The function does not change any observations, but restrucutres data
#' to correctly handle replicates internally.
#'
#' @param Y1      a list of observations for the first field,
#'                where each list may have different length, meaning that
#'                replicates are not uniformly observed at same locations
#' @param Y2      a list of observations for the second field
#'
#' @return Y       a vector Y of observations stacked in the format specific
#' @return nrep    a vector specifying number of observations for each replicate
#' for the projector matrix, A specified in ngme2.R
#' @export
#'
#' @examples
#' Y <- format_y(Y1, Y2, nrep)
#'
format_y <- function(Y1, Y2) {
    nrep1 <- length(Y1)
    nrep2 <- length(Y2)
    # -------------  CHECK INPUT ---------------
    if (nrep1 != nrep2) {
        stop("Replicates for two fields differ\n")
    }
    nrep <- matrix(NA, nrep1)
    Y1_full <- NULL
    Y2_full <- NULL
    Y <- NULL

    for (i in 1:nrep1) {
        Y1_full <- c(Y1_full, Y1[[i]])
        Y2_full <- c(Y2_full, Y2[[i]])
        Y <- c(Y, Y1[[i]], Y2[[i]])
        nrep[i] <- length(Y1[[i]])
    }

#returns a stacked vector Y = [Y11^T, Y21^T, ..., Y1nrep^T, Y2nrep^T]
#where Y11 is a vector of observations for first replicate at locations of replicate 1 # nolint
# NOTE Y is not stacked by observations: Y = [Y11(s1),Y21(s1),Y11(s2),Y21(s2)...]^T # nolint
    return(list(Y, nrep))
}
