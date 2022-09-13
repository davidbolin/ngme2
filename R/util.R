#' Make observation matrix for time series
#'
#' @param loc   integers (after sorting, no gaps > 1)
#' @param replicates indicating replicate measure at same location
#' @param range range for the mesh
#'  by default range=(min(loc), max(loc))
#'
#' @return A matrix (length(loc) * length(unique(loc)))
#' @export
#'
#' @examples ngme.ts.make.A(c(1, 2, 2), end = 5, replicates = c(1, 1, 2))
#'
ngme.ts.make.A <- function(loc, replicates = NULL, range = c(1, max(loc))) {
  if (is.null(loc)) return (NULL)

  n_loc <- length(loc)
  nrep <- 1

  start <- range[1]; end <- range[2]
  n_range <- end - start + 1

  if (is.null(replicates)) {
    replicates <- rep(1, n_range)
  }

  unique_rep <- unique(replicates)
  nrep <- length(unique_rep)

  A <- matrix(0, nrow = n_loc, ncol = n_range * nrep)

  for (i in 1:n_loc) {
    ncol_rep <- which(unique_rep == replicates[i])
    A[i, (ncol_rep - 1) * n_range + loc[i] - start + 1] <- 1
  }
  as(A, "dgCMatrix")
}


# ngme.ts.make.A <- function(loc) {
#   n_loc = length(loc)
#   n_range = max(loc)-min(loc)+1
#   if (any((diff(sort(loc))) > 1)) stop("no gaps allowed")
#
#   A = matrix(0, nrow=n_loc, ncol=n_range)
#   for (i in 1:n_loc) {
#     A[i, loc[i]-min(loc)+1] = 1
#   }
# #Do an na.rm
#   as(A, "dgCMatrix")
# }

#' NGME starting point for block model
#'
#' @param sigma.eps mesurement noise
#' @param fixed.effects fixed effects
#' @param W initial W for the block model
#'
#' @return
#' @export
#'
#' @examples
ngme.start <- function(
  fixed.effects = NULL,
  sigma.eps     = NULL,
  W = NULL
) {

  start = list(
    fixed.effects     = fixed.effects,
    mesurement.noise  = sigma.eps,
    block.W           = W
  )

  return (start)
}

#' Convert sparse matrix into sparse dgCMatrix
#'
#' @return sparse dgCMatrix
#' @export
#'
#' @examples
ngme.as.sparse <- function(G) {
  tryCatch(
    expr={
      G <- as(G, "dgCMatrix")
    },
    error = function(e) {
      G <- as(G, "dgTMatrix")
      idx <- which(G@i <= G@j)
      G = Matrix::sparseMatrix(i=G@i[idx], j = G@j[idx], x= G@x[idx],
                               symmetric=FALSE, index1 = FALSE)
      G = as(G, "dgCMatrix")
    },
    finally = {
      G
    }
  )
}

#' Make index for the matern model
#'
#' @param name
#' @param n.spde
#' @param n.repl
#' @param mesh
#' @param dim
#'
#' @return
#' @export
#'
#' @examples
ngme.matern.make.index <- function(n.spde=NULL,
                                   n.repl = 1,
                                   mesh = NULL,
                                   dim = NULL){

  if(is.null(n.spde)&&is.null(mesh)){
    stop("You should provide either n.spde or mesh!")
  }

  if(!is.null(mesh)){
    n_mesh = mesh$n

    if(mesh$manifold == "R1"){
      dim = 1
    } else if(mesh$manifold == "R2"){
      dim = 2
    } else{
      stop("The domain must be flat manifolds of dimension 1 or 2, that is,
         the domain must be a line or a plane.")
    }

  } else{
    n_mesh <- n.spde
    if(is.null(dim)) {
      stop("You should provide the dimension d!")
    }
  }

  out <- list()
  out$index <- rep(1:n_mesh, times = n.repl)

  out$replicates <- rep(1:n.repl, each = n_mesh)
  return(out)
}


#' Split the data according to NA in response variable
#'
#' @param formula
#' @param data
#'
#' @return
#' @export
#'
#' @examples
parse_formula_NA <- function(formula, data) {
  stopifnot("Please provide a valid formula" = inherits(formula, "formula"))

  Y_full <- model.frame(formula, data, na.action = NULL)[[1]]

  mf <- model.frame(formula, data)
  Y_data <- mf[[1]]
  index_data  <- as.numeric(rownames(mf))
  X_data      <- model.matrix((terms(formula)), data)
  
  # watch out! X_full can be 0*0 
  X_full      <- model.matrix(delete.response(terms(formula)), data)
  contain_NA  <- !is.null(attr(mf, "na.action"))
  
  # if there is NA in response variable
  if (!is.null(attr(mf, "na.action"))) {
    index_NA     <- as.numeric(attr(mf, "na.action"))
    
    # X_NA         <- X_full[index_NA, , drop = FALSE]
    if (all(dim(X_full) == c(0, 0)))
      X_NA      <- X_full
    else 
      X_NA      <- X_full[index_NA, , drop = FALSE]
  } else {
    index_NA <- X_NA <- NULL
  }

  list(
    length      = length(Y_data),
    Y_data      = Y_data,
    X_data      = X_data,
    index_data  = index_data,
    contain_NA  = contain_NA,
    index_NA    = index_NA,
    X_NA        = X_NA

    # Y_full      = Y_full # only for f to use
  )
}

# compare numerical grad. and ana. grad.
# compare_grad_ngme <- function(
#   formula,
#   data,
#   controls      = ngme.control(),
#   debug         = ngme.debug(),
#   noise         = ngme.noise(),
#   last_fit      = NULL,
#   beta          = NULL,
#   seed          = NULL
# ) {
#   control_ana <- 
#   control_num <- 
# }

# draft
# acceptNA <- function(Y) {
#   Y_clean <- na.rm(Y)
  
#   # compute A matrix
#   A <- makeA(loc = Y, mesh)
#   A_clean <- makeA(loc = Y_clean, mesh)
  
#   ngme1 <- ngme(
#     formula = Y_clean ~ f(A = A_clean),
#     control = ngme.control(estimation = TRUE)
#   )

#   # Y has NA
#   ngme2 <- ngme(
#     formula = Y_clean ~ f(A = A),
#     control = ngme.control(estimation = FALSE),
#     last_fit = ngme1
#   )

#   ngme2$output$W
# }

# post.sampleW(ngme) {
#   W
# }

# plot_density <- function(ngme) {

# }