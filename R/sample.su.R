#' @export
getdenspars <- function(sigma.sq.beta, kappa, fam = "power") {

  length.key <- 10000000
  if (fam == "power") {
    qs <- exp(seq(log(0.1), log(6), length.out = length.key))
    qstoks <- gamma(5/qs)*gamma(1/qs)/(gamma(3/qs)^2) - 3
    q <- qs[min(which(qstoks <= kappa))]
    lambda <- (sigma.sq.beta*gamma(1/q)/gamma(3/q))^(-q/2)
  } else if (fam == "bessel") {
    q <- (6/kappa - 1)/2
    lambda <- sqrt(sigma.sq.beta/(2*q + 1))
  }  else if (fam == "pearson") {
    q <- 5/2 + 3/kappa
    lambda <- sqrt(sigma.sq.beta*(2*q - 3))
  } else if (fam == "weibull") {
    # pdf is infinite at 0 when kappa > 3
    qs <- exp(seq(log(0.1), log(10), length.out = length.key))
    qstoks <- (gamma(4/qs + 1) - 3*gamma(2/qs + 1)^2)/gamma(2/qs + 1)^2
    q <- qs[min(which(qstoks <= kappa))]
    lambda <- sqrt(sigma.sq.beta/gamma(2/q + 1))
  } else if (fam == "dl") {
    qs <- exp(seq(log(0.000001), log(999999), length.out = length.key))
    qstoks <- 3*((qs^2 + 9*qs + 12)/(qs*(qs + 1)))
    q <- qs[min(which(qstoks <= kappa))]
    lambda <- sigma.sq.beta
  }
  return(list("q" = q, "lambda" = lambda))
}

#' @export
shrinkdens <- function(x, sigma.sq.beta, kappa, fam = "power", pars = NULL, log.nocons = FALSE) {

  if (is.null(pars)) {pars <- getdenspars(sigma.sq.beta = sigma.sq.beta, kappa = kappa, fam = fam)}

  q <- pars[["q"]]
  lambda <- pars[["lambda"]]

  if (fam == "power") {
    if (!log.nocons) {
      return(q*lambda^(1/q)*exp(-lambda*abs(x)^q)/(2*gamma(1/q)))
    } else {
      return(-lambda*abs(x)^q)
    }
  } else if (fam == "bessel") {
    if (!log.nocons) {
      return(abs(x)^q/(sqrt(pi)*2^q*lambda^(q + 1)*gamma(q + 1/2))*besselK(abs(x)/lambda, q))
    } else {
      return(q*log(abs(x)) + log(besselK(abs(x)/lambda, q)))
    }
  }  else if (fam == "pearson") {
    if (!log.nocons) {
      return((1 + (x/lambda)^2)^(-q)/(lambda*beta(q - 1/2, 1/2)))
    } else {
      return(-q*log(1 + (x/lambda)^2))
    }
  } else if (fam == "weibull") {
    # pdf is infinite at 0 when kappa > 3
    if (!log.nocons) {
      return((q/lambda)*(abs(x)/lambda)^(q - 1)*exp(-(abs(x)/lambda)^c)/2)
    } else {
      return((q - 1)*log(abs(x)) -(abs(x)/lambda)^q)
    }
  } else if (fam == "dl") {
    alpha <- sqrt(1/(2^3*q*(q + 1)))*sqrt(lambda)
    if (!log.nocons) {
      return((2^((1 + q)/2)*gamma(q))^(-1)*abs(x/alpha)^((q - 1)/2)*besselK(sqrt(2)*sqrt(abs(x/alpha)), 1 - q)/alpha)
    } else {
      return((q - 1)/2)*log(abs(x)) + log(besselK(sqrt(2)*sqrt(abs(x/alpha), 1 - q)))
    }
  }
}

sample.u <- function(XtX, Xty, s, u, sigma.sq.z) {

  A <- XtX*(s%*%t(s))
  B <- Xty*s

  for (i in 1:length(u)) {
    probs <- exp(-1/(2*sigma.sq.z)*((A[i, i] + 2*c(-1, 1)*sum(A[i, -i]*u[-i])) - 2*c(-1, 1)*B[i]))
    probs[is.infinite(probs)] <- 1
    u[i] <- sample(c(-1, 1), 1,
                   prob = probs)
  }
  return(u)

}

sample.s <- function(XtX, Xty, u, sigma.sq.z, sigma.sq.beta, kappa = 3, s.old,
                     sing.x = FALSE,
                     epsilon = 0, fam = "power", pars = NULL, proposal = "marginal", delta = 0, mean.adj = 0) {

  A <- XtX*(u%*%t(u))/sigma.sq.z

  b <- (Xty*u/sigma.sq.z) - mean.adj

  e.A <- eigen(A)

  V <- solve(A + epsilon*diag(ncol(A)))
  m <- V%*%b

  # e.V <- eigen(V)
  # rt.V <- e.V$vectors[, e.V$values > 0]%*%diag(sqrt(e.V$values[e.V$values > 0]))%*%t(e.V$vectors[, e.V$values > 0])
  #
  # s <- rep(-1, length(m))
  #
  # while(min(s) < 0) {
  #   s <- m + rt.V%*%rnorm(length(m))
  # }

  acc <- rep(1, p)
  s <- s.old


  if (proposal == "joint") {

    V.eig <- eigen(V)
    V.inv <- V.eig$vectors[, V.eig$values > 0]%*%diag(1/V.eig$values[V.eig$values > 0])%*%t(V.eig$vectors[, V.eig$values > 0])
    s <- rtmvnorm(1, mean = as.numeric(m), lower = rep(0, p), algorithm = "gibbs",
                  H = V.inv)[1, ]
    while(sum(is.infinite(s) | is.na(s)) > 0) {
      s <- rtmvnorm(1, mean = as.numeric(m), lower = rep(0, p), algorithm = "gibbs",
                    H = V.inv, start.value = abs(m)/2)[1, ]
    }


    lik.new <- -(1/2)*(crossprod(t(crossprod(s, A)), s) - 2*crossprod(s, b)) + sum(shrinkdens(s, sigma.sq.beta = sigma.sq.beta, kappa = kappa, fam = fam,
                                                                                              pars = pars, log.nocons = TRUE))
    lik.old <- -(1/2)*(crossprod(t(crossprod(s.old, A)), s.old) - 2*crossprod(s.old, b)) + sum(shrinkdens(s.old, sigma.sq.beta = sigma.sq.beta, kappa = kappa, fam = fam,
                                                                                                          pars = pars, log.nocons = TRUE))

    diff <- exp(lik.new - lik.old)[1, 1]

    if (diff < 1 & runif(1, 0, 1) > diff) {
      s <- s.old
      acc <- rep(0, p)
    }


  } else if (proposal == "marginal" | proposal == "conditional") {

    for (i in 1:p) {
      if (proposal == "marginal") {
        mm <- m[i]
        vv <- V[i, i]
      } else if (proposal == "conditional") {
        BB <- crossprod(solve(V[-i, -i]), V[i, -i])
        mm <- as.numeric(m[i] + crossprod(BB, s[-i] - m[-i]))
        vv <- as.numeric(V[i, i] - crossprod(BB, V[i, -i]))
      }

      s.new <- rtmvnorm(1, mean = mm, sigma = vv, lower = 0, algorithm = "gibbs",
                        burn.in.samples = 1000)
      # while(is.infinite(s.new)) {
      #   s.new <- rtmvnorm(1, mean = mm, sigma = vv, lower = 0, algorithm = "gibbs",
      #                     start.value =  0)
      # }
      # Get extra samples from prior for densities with infinite p(|x|) at x = 0
      if (is.infinite(s.new)) {
        s[i] <- s.old[i]
        acc[i] <- 0
      } else {
        extra.zero <- rbinom(1, 1, delta[i])
        if ((fam == "bessel" | fam == "dl") & extra.zero > 0) {
          p.samp <- -1
          while (p.samp < 0) {
            if (fam == "bessel") {
              p.samp <- sqrt(sigma.sq.beta)*sqrt(2/(pars[["q"]] + 1/2))*sqrt(rgamma(1, pars[["q"]] + 1/2, 2))*rnorm(1)
            } else if (fam == "dl") {
              p.samp <- sqrt(sigma.sq.beta)*(1/2)*rexp(1, 1)*(-1)^rbinom(1, 1, 0.5)*rgamma(1, shape = pars[["q"]], rate = 1/2)
            }
          }
        } else {
          p.samp <- 0
        }
        s[i] <- s.new*(1 - extra.zero) + (extra.zero)*p.samp

        lik.new <- -(1/2)*(crossprod(t(crossprod(s, A)), s) - 2*crossprod(s, b)) + sum(shrinkdens(s, sigma.sq.beta = sigma.sq.beta, kappa = kappa, fam = fam,
                                                                                                pars = pars, log.nocons = TRUE))
        lik.old <- -(1/2)*(crossprod(t(crossprod(s.old, A)), s.old) - 2*crossprod(s.old, b)) + sum(shrinkdens(s.old, sigma.sq.beta = sigma.sq.beta, kappa = kappa, fam = fam,
                                                                                                            pars = pars, log.nocons = TRUE))

        diff <- exp(lik.new - lik.old)[1, 1]
        if (diff < 1 & runif(1, 0, 1) > diff) {
          s[i] <- s.old[i]
          acc[i] <- 0
        }
      }
    }
  }

  return(list("s" = s, "acc" = acc))
}

#' @export
sample.su <- function(X, y, sigma.sq.z, sigma.sq.beta, kappa = 3,
                      epsilon = 0,
                      num.samp = 1000, print.iter = FALSE,
                      fam = "power", delta = 0, proposal = "marginal", kl.tweak = FALSE) {

  library(tmvtnorm)


  if (fam == "dl" & kappa <= 3) {
    cat("Values of excess kurtosis less than 3 cannot be represented by the DL prior.\n")
  }
  pars <- getdenspars(sigma.sq.beta = sigma.sq.beta, kappa = kappa, fam = fam)
  pars.pow <- getdenspars(sigma.sq.beta = sigma.sq.beta, kappa = kappa, fam = "power") # Can use to adjust sampler

  if (kl.tweak) {
    q.mommatch <- pars.pow[["q"]]
    mean.adj <- sqrt(gamma(1/q.mommatch)*gamma(3/q.mommatch)/(gamma(2/q.mommatch)^2))/sqrt(sigma.sq.beta)
    print(mean.adj)
  } else {mean.adj = 0}

  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  sing.x <- min(eigen(XtX)$values) <= 0

  ridge.est <- solve(XtX + 10^(-1)*diag(ncol(X)))%*%(Xty)
  u <- sign(ridge.est)
  s <- abs(ridge.est)

  samples.s <- samples.u <- samples.beta <- matrix(nrow = num.samp, ncol = ncol(X))
  samples.s[1, ] <- s
  samples.u[1, ] <- u
  samples.beta[1, ] <- s*u

  acc <- matrix(NA, nrow = num.samp, ncol = ncol(X))
  for (i in 2:num.samp) {
    if (print.iter) { cat("i = ", i, "\n") }

    s.s <- sample.s(XtX, Xty,  u = samples.u[i - 1, ], sigma.sq.z = sigma.sq.z,
                    sigma.sq.beta = sigma.sq.beta, kappa = kappa,
                    s.old = samples.s[i - 1, ],
                    sing.x = sing.x,
                    epsilon = epsilon,
                    fam = fam, pars = pars, proposal = proposal, delta = delta, 
                    mean.adj = mean.adj)
    samples.s[i, ] <- s.s$s
    acc[i, ] <- s.s$acc
    if (print.iter & max(s.s$acc) == 1) {cat("Accepted!\n")}

    samples.u[i, ] <- sample.u(XtX, Xty, samples.s[i, ], samples.u[i - 1, ], sigma.sq.z)
    samples.beta[i, ] <- samples.s[i, ]*samples.u[i, ]
  }
  return(list("s" = samples.s[-1, ], "u" = samples.u[-1, ],
              "beta" = samples.beta[-1, ], "acc" = acc[-1, ]))

}


