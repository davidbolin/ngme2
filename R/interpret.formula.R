#' Interpret the formula
#'
#' @param gf formula
#' @param data data.frame
#' @param debug
#'
#' @return 1. plain formula without f 2. latents_in - a list of lists
#' @export
#'
#' @examples
ngme.interpret.formula <- function(gf,
                                   data,
                                   debug=FALSE) {

  tf <- terms.formula(gf, specials = c("f"))
  terms <- attr(tf, "term.labels")
  intercept <- attr(tf, "intercept")
  nt <- length(terms)

  latents_in = list()
  # order of f terms in labels
  spec.order = attr(tf, "specials")$f - 1
  for (i in spec.order) {
    str = gsub("^f\\(", "ngme2::f(", terms[i])
    res = eval(parse(text = str), envir = data)
    latents_in[[length(latents_in) + 1]] = res
  }

  fixf = terms[-spec.order]


  # construct formula without f
    fm = as.character(attr(tf, "variables")[[2]])
    fm = paste(fm, "~", intercept, paste(c("", fixf), collapse = " + "))
# paste(c("", "x1", "x2"), collapse = "+")
  return(list(latents_in=latents_in,
              plain.fm = formula(fm)))
}

# ngme.interpret.formula2 <- function(gf,
#                                    debug = FALSE,
#                                    data = NULL,
#                                    data.model=NULL,
#                                    parent.frame = NULL
#                                   ) {
#   if (!is.null(data.model)) {
#     ## make a copy of the environment
#     p.env <- new.env(hash = TRUE, parent = environment(gf))
#     ## then add the entries in data.model to that environment
#     for (nm in names(data.model)) {
#       idx <- which(names(data.model) == nm)
#       assign(nm, data.model[[idx]], envir = p.env)
#     }
#   } else {
#     ## otherwise, just use (for read-only) this environment
#     if (missing(parent.frame)) {
#       p.env <- environment(gf)
#     } else {
#       p.env <- parent.frame
#     }
#   }

# ## check f func
#   tf <- terms.formula(gf, specials = c("f"), data = NULL)
#   terms <- attr(tf, "term.labels")
#   nt <- length(terms)

#   if (attr(tf, "response") > 0) {
#     ## fixf formula with ONLY fixed effects.  randf formula with
#     ## ONLY random effect.  weightf formula where are the names of
#     ## the (possible) weigths for the covariates
#     response <- as.character(attr(tf, "variables")[2])
#     fixf <- randf <- weightf <- paste(response, "~", sep = "")
#   } else {
#     stop("\n\tA response variable has to be present")
#   }

#   rt <- attr(tf, "specials")$f
#   vtab <- attr(tf, "factors")
#   if (length(rt) > 0) {
#     for (i in 1:length(rt)) {
#       ind <- (1:nt)[as.logical(vtab[rt[i], ])]
#       rt[i] <- ind
#     }
#   }

#   k <- ks <- kp <- 1
#   len.rt <- length(rt)
#   random.spec <- list()
#   if (nt > 0) {
#     for (i in 1:nt) {
#       if (k <= len.rt && ((ks <= len.rt && rt[ks] == i))) {
# # print(data)
#         st <- eval(parse(text = gsub("^f\\(", "ngme2::f(", terms[i])), envir = data, enclos = p.env)
#         # st <- eval(parse(text = gsub("^f\\(", "f(", terms[i])), envir = data, enclos = p.env)
#         random.spec[[k]] <- st
#         if (ks <= len.rt && rt[ks] == i) {
#           ks <- ks + 1
#         } else {
#           kt <- kt + 1
#         }
#         k <- k + 1
#       } else {
#         if (kp > 1) {
#           fixf <- paste(fixf, " + ", terms[i], sep = "")
#         } else {
#           fixf <- paste(fixf, terms[i], sep = "")
#         }
#         kp <- kp + 1
#       }
#     }
#   }

#   if (debug) {
#     print(fixf)
#     # print(random.spec)
#   }

# # build in_list
#   # create design matrix for fix effect
#   # fix is string "y~x"
#   X = NULL
#   tryCatch({
#     X = model.matrix(formula(fixf), data=data)
#   },
#   error=function(cond) {
#     if (attr(tf, "intercept")==1)
#       fixf <- paste(fixf, 1, sep = "")
#     else
#       fixf <- paste(fixf, -1, sep = "")
#     X <<- model.matrix(formula(fixf), data=data)
#   },
#   finally = {
#   })

#   n_regs = Reduce("+", lapply(random.spec, function(l) {l$n_reg}))
#   Y = eval(parse(text=response), env=data)
#   general_in=list(X=X, Y=Y, n=length(Y), n_regs=n_regs)

#   return(list(general_in=general_in,
#               latents_in=random.spec))
# }

# test
# a = 3;
# environment(y~x)
# evalq(1+1, list(`+`=`-`))
#
# eval({ xx <- pi; xx^2}) ; xx
# ngme.interpret.formula(y~x+f(x,model="ar1"), debug=T)

# evalq(x+x, envir=list(x=1))



