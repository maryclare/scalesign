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
    if (!log.nocons) {
      return((2^((1 + q)/2)*gamma(q))^(-1)*(sqrt(sigma.sq.beta))^((q - 1)/2)*abs(x/sqrt(sigma.sq.beta))^((q - 1)/2)*besselK(sqrt(2)*sqrt(sqrt(sigma.sq.beta))*sqrt(abs(x/sqrt(sigma.sq.beta))), 1 - q))
    } else {
      return(((q - 1)/2)*log(abs(x/sqrt(sigma.sq.beta))) + log(besselK(sqrt(2)*sqrt(sqrt(sigma.sq.beta))*sqrt(abs(x/sqrt(sigma.sq.beta))), 1 - q)))
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
                     epsilon = 0, fam = "power", pars = NULL, proposal = "marginal") {

  A <- XtX*(u%*%t(u))/sigma.sq.z

  b <- Xty*u/sigma.sq.z

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

  for (i in 1:p) {

    if (proposal == "marginal") {
      mm <- m[i]
      vv <- V[i, i]
    } else if (proposal == "conditional") {
      BB <- crossprod(solve(V[-i, -i]), V[i, -i])
      mm <- m[i] + crossprod(BB, s[-i] - m[-i])
      vv <- V[i, i] - crossprod(BB, V[i, -i])
    }

    # Compute mean and variance of this distribution
    alpha <- -mm/sqrt(vv)
    log.z <- pnorm(-alpha, log = TRUE)
    alpha.z <- exp(dnorm(alpha, log = TRUE)- log.z)
    mean <- (mm + alpha.z*sqrt(vv))
    var <- vv*(1 + alpha*alpha.z - (alpha.z)^2)

    k <- (var/mean)
    theta <- (mean/k)
    s.new.inv <- rinvgamma(1, shape = theta, scale = 1/k)
    s.new <- 1/s.new.inv

    s[i] <- s.new
    # prior.new <- ifelse(prior == "power", -(sqrt(sigma.sq.beta*gamma(1/q)/gamma(3/q)))^(-q)*sum(s^q),
    #                     ifelse(prior == "pearson", sum(-q*log(1 + (s/(sqrt(sigma.sq.beta)*sqrt(2*q - 3)))^2)),
    #                            ifelse(prior == "bessel",
    #                                   sum(q*log(s) + log(besselK(sqrt(1 + q)*s/sqrt(sigma.sq.beta), q))), NA)))
    # prior.old <- ifelse(prior == "power", -(sqrt(sigma.sq.beta*gamma(1/q)/gamma(3/q)))^(-q)*sum(s.old^q),
    #                     ifelse(prior == "pearson",
    #                            sum(-q*log(1 + (s.old/(sqrt(sigma.sq.beta)*sqrt(2*q - 3)))^2)),
    #                          ifelse(prior == "bessel",
    #                                 sum(q*log(s.old) + log(besselK(sqrt(1 + q)*s.old/sqrt(sigma.sq.beta), q))),
    #                                 NA)))
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

  return(list("s" = s, "acc" = acc))
}

#' @export
sample.su <- function(X, y, sigma.sq.z, sigma.sq.beta, kappa = 3,
                      epsilon = 0,
                      num.samp = 1000, print.iter = FALSE,
                      fam = "power", delta = 10^(-7), proposal = "marginal") {

  if (fam == "dl" & kappa <= 3) {
    cat("Values of excess kurtosis less than 3 cannot be represented by the DL prior.\n")
  }
  pars <- getdenspars(sigma.sq.beta = sigma.sq.beta, kappa = kappa, fam = fam)

  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  sing.x <- min(eigen(XtX)$values) <= 0

  ridge.est <- solve(XtX + delta*diag(ncol(X)))%*%(Xty)
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
                    fam = fam, pars = pars, proposal = proposal)
    samples.s[i, ] <- s.s$s
    acc[i, ] <- s.s$acc
    if (print.iter & max(s.s$acc) == 1) {cat("Accepted!\n")}

    samples.u[i, ] <- sample.u(XtX, Xty, samples.s[i, ], samples.u[i - 1, ], sigma.sq.z)
    samples.beta[i, ] <- samples.s[i, ]*samples.u[i, ]
  }
  return(list("s" = samples.s[-1, ], "u" = samples.u[-1, ],
              "beta" = samples.beta[-1, ], "acc" = acc[-1, ]))

}


