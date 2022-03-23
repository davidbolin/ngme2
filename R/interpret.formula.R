ngme.interpret.formula <- function(
                                  gf,
                                  debug = FALSE,
                                  data = NULL,
                                  data.model=NULL,
                                  parent.frame = NULL
                                  ) {
  if (!is.null(data.model)) {
    ## make a copy of the environment
    p.env <- new.env(hash = TRUE, parent = environment(gf))
    ## then add the entries in data.model to that environment
    for (nm in names(data.model)) {
      idx <- which(names(data.model) == nm)
      assign(nm, data.model[[idx]], envir = p.env)
    }
  } else {
    ## otherwise, just use (for read-only) this environment
    if (missing(parent.frame)) {
      p.env <- environment(gf)
    } else {
      p.env <- parent.frame
    }
  }

## check f func
  tf <- terms.formula(gf, specials = c("f"), data = NULL)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)

  rt <- attr(tf, "specials")$f
  vtab <- attr(tf, "factors")
  if (length(rt) > 0) {
    for (i in 1:length(rt)) {
      ind <- (1:nt)[as.logical(vtab[rt[i], ])]
      rt[i] <- ind
    }
  }

  k <- ks <- kp <- 1
  len.rt <- length(rt)
  random.spec <- list()
  if (nt > 0) {
    for (i in 1:nt) {
      if (k <= len.rt && ((ks <= len.rt && rt[ks] == i))) {
        st <- eval(parse(text = gsub("^f\\(", "ngme2::f(", terms[i])), envir = data, enclos = p.env)
        # st <- eval(parse(text = gsub("^f\\(", "f(", terms[i])), envir = data, enclos = p.env)
        random.spec[[k]] <- st
        if (ks <= len.rt && rt[ks] == i) {
          ks <- ks + 1
        } else {
          kt <- kt + 1
        }
        k <- k + 1
      } else {
        if (kp > 1) {
          fixf <- paste(fixf, " + ", terms[i], sep = "")
        } else {
          fixf <- paste(fixf, terms[i], sep = "")
        }
        kp <- kp + 1
      }
    }
  }

## check response
  if (attr(tf, "response") > 0) {
    ## fixf formula with ONLY fixed effects.  randf formula with
    ## ONLY random effect.  weightf formula where are the names of
    ## the (possible) weigths for the covariates
    response <- as.character(attr(tf, "variables")[2])
    fixf <- randf <- weightf <- paste(response, "~", sep = "")
  } else {
    stop("\n\tA response variable has to be present")
  }

  if (debug) {
    print(response)
    # print(random.spec)
  }

# build in_list
  n_regs = Reduce("+", lapply(random.spec, function(l) {l$n_reg}))
  Y = eval(parse(text=response))
  general_in=list(Y=Y, n=length(Y), n_regs=n_regs)

  return(list(general_in=general_in, latents_in=random.spec))
}

# terms.formula(y~x)

