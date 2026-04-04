#####################################################
#############     Univariate GLMMs     ##############
#####################################################



# Univariate Weibull PH 
weibullPH_lognormal <- function() {
  L <- stats::make.link("identity")
  list(
    family  = "WeibullPH_lognormal",
    link = "identity",
    linkfun = L$linkfun,
    linkinv = L$linkinv,
    log_dens = function(y, eta, mu_fun, phis, eta_zi) {
      if (is.null(phis)) stop("'phis' must be supplied for WeibullPH_lognormal")
      y <- .coerce_y3(y, "WeibullPH_lognormal")
      start  <- as.numeric(y[, 1]); stop <- as.numeric(y[, 2]); status <- as.numeric(y[, 3])
      eps <- 1e-12; stop <- pmax(stop, eps); start <- pmax(start, 0)
      lambda <- exp(phis[1]); rho <- exp(phis[2])     # scale, shape
      H0_inc <- lambda * (stop^rho - pmax(start, eps)^rho)
      logh0  <- log(lambda) + log(rho) + (rho - 1) * log(stop)
      eta <- pmin(pmax(eta, -15), 15); haz <- exp(eta)
      if (is.matrix(eta)) {
        out <- (matrix(status, nrow(eta), ncol(eta)) *
                  (matrix(logh0, nrow(eta), ncol(eta)) + eta)) -
               matrix(H0_inc, nrow(eta), ncol(eta)) * haz
      } else {
        out <- status * (logh0 + eta) - H0_inc * haz
      }
      attr(out, "mu_y") <- haz
      out
    },
    score_eta_fun = function(y, mu, phis, eta_zi) {
      if (is.null(phis)) stop("'phis' must be supplied for WeibullPH_lognormal (score)")
      y <- .coerce_y3(y, "WeibullPH_lognormal")
      start  <- as.numeric(y[, 1]); stop <- as.numeric(y[, 2]); status <- as.numeric(y[, 3])
      eps <- 1e-12; stop <- pmax(stop, eps); start <- pmax(start, 0)
      lambda <- exp(phis[1]); rho <- exp(phis[2])
      H0_inc <- lambda * (stop^rho - pmax(start, eps)^rho)
      status - H0_inc * mu
    }
  )
}

# Univariate cumulative probit
cumulativeprobit <- function(K) {
  log_dens <- function(y, eta, mu_fun, phis, eta_zi) {
    n <- length(y)
    if (!is.matrix(eta)) eta <- matrix(eta, nrow = n)
    if (nrow(eta) != n) stop("nrow(eta) deve essere uguale a length(y)")
    if (any(y < 1L | y > K)) stop("y deve essere in 1..K")
    tau <- as.numeric(phis)
    if (length(tau) != K - 1L) stop("length(phis) deve essere K-1")
    tau_u <- rep( Inf, n); iu <- which(y <  K); tau_u[iu] <- tau[y[iu]]
    tau_l <- rep(-Inf, n); il <- which(y > 1L); tau_l[il] <- tau[y[il] - 1L]
    U <- matrix(tau_u, nrow = n, ncol = ncol(eta))
    L <- matrix(tau_l, nrow = n, ncol = ncol(eta))
    P <- pnorm(U - eta) - pnorm(L - eta)
    eps <- .Machine$double.eps
    P[P <  eps] <- eps; P[P > 1 - eps] <- 1 - eps
    out <- log(P)
    attr(out, "mu_y") <- P
    out
  }
  structure(list(
    family = "cumulative link",
    link= "probit",
    linkfun = function(mu) mu,
    linkinv = function(eta) eta,
    log_dens = log_dens
  ), class = "family")
}









.coerce_y3 <- function(y, who = "WeibullPH_lognormal") {
  if (is.matrix(y)) {
    if (ncol(y) < 3L) stop(sprintf("[%s] y must be cbind(start, stop, status).", who))
    y[, 1:3, drop = FALSE]
  } else {
    if (length(y) != 3L) stop(sprintf("[%s] y must have 3 elements.", who))
    matrix(y, nrow = 1L)
  }
}
