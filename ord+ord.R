#######################################################
########                                ###############
####### Code for a pair of ordinal items ##############
########                                ###############
#######################################################
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(stringr)
library(GLMMadaptive)
library(MASS)
library(mvtnorm)  
library(pbivnorm)



# ── family: cumulative probit (y stackata; K variabile per item) ────────────────
biv_cumulative_probit_glmmadaptive <- function(K_item_named, link = "identity") {
  if (length(K_item_named) != 2L) stop("K_item_named deve avere 2 elementi.")

  K1 <- as.integer(K_item_named[1])
  K2 <- as.integer(K_item_named[2])

  if (K1 < 2L || K2 < 2L) stop("each item must have at least 2 categories.")

  stats <- make.link(link)
  len_th1 <- K1 - 1L
  len_th2 <- K2 - 1L
  n_phis <- len_th1 + len_th2 + 1L
  eps <- .Machine$double.eps

  F_vec <- function(x, y, rho) {
  out <- numeric(length(x))

  m1 <- (is.infinite(x) & x < 0) | (is.infinite(y) & y < 0)
  out[m1] <- 0

  m2 <- (is.infinite(x) & x > 0) & (is.infinite(y) & y > 0)
  out[m2] <- 1

  m3 <- (is.infinite(x) & x > 0) & is.finite(y)
  out[m3] <- pnorm(y[m3])

  m4 <- (is.infinite(y) & y > 0) & is.finite(x)
  out[m4] <- pnorm(x[m4])

  mf <- is.finite(x) & is.finite(y)
  if (any(mf)) {
    xf <- x[mf]
    yf <- y[mf]
    tmp <- pbivnorm::pbivnorm(xf, yf, rho)
    bad <- !is.finite(tmp)
    if (any(bad)) {
      Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
      tmp[bad] <- vapply(
        seq_along(xf[bad]),
        function(i) as.numeric(mvtnorm::pmvnorm(
          lower = c(-Inf, -Inf),
          upper = c(xf[bad][i], yf[bad][i]),
          sigma = Sigma
        )),
        numeric(1)
      )
    }
    out[mf] <- tmp
  }

  out[out < 0] <- 0
  out[out > 1] <- 1
  out
}

  log_dens <- function(y, eta, mu_fun, phis, eta_zi) {
    y <- as.integer(y)
    n <- length(y)

    if (length(phis) != n_phis) stop(sprintf("phis deve avere lunghezza %d.", n_phis))

    th1 <- phis[seq_len(len_th1)]
    th2 <- phis[len_th1 + seq_len(len_th2)]
    rho <- tanh(phis[n_phis])

    alpha <- c(-Inf, th1, Inf)
    beta  <- c(-Inf, th2, Inf)

    mu1 <- mu_fun(eta)
    mu2 <- mu_fun(eta_zi)

    decode_y <- function(y) {
      idx1  <- which(y >= 1L & y <= K1)
      idx2  <- which(y > K1 & y <= K1 + K2)
      idx12 <- which(y > K1 + K2 & y <= K1 + K2 + K1 * K2)

      c1 <- if (length(idx1)) y[idx1] else integer(0)
      d2 <- if (length(idx2)) y[idx2] - K1 else integer(0)

      if (length(idx12)) {
        tmp <- y[idx12] - (K1 + K2)
        c12 <- ((tmp - 1L) %% K1) + 1L
        d12 <- ((tmp - 1L) %/% K1) + 1L
      } else {
        c12 <- integer(0)
        d12 <- integer(0)
      }

      list(idx1 = idx1, idx2 = idx2, idx12 = idx12, c1 = c1, d2 = d2, c12 = c12, d12 = d12)
    }

    one_eval <- function(m1, m2) {
      dec <- decode_y(y)
      P <- rep(1, n)

      if (length(dec$idx1)) {
        ii <- dec$idx1
        L1 <- alpha[dec$c1]
        U1 <- alpha[dec$c1 + 1L]
        P[ii] <- pnorm(U1 - m1[ii]) - pnorm(L1 - m1[ii])
      }

      if (length(dec$idx2)) {
        ii <- dec$idx2
        L2 <- beta[dec$d2]
        U2 <- beta[dec$d2 + 1L]
        P[ii] <- pnorm(U2 - m2[ii]) - pnorm(L2 - m2[ii])
      }

      if (length(dec$idx12)) {
        ii <- dec$idx12

        L1 <- alpha[dec$c12]
        U1 <- alpha[dec$c12 + 1L]
        L2 <- beta[dec$d12]
        U2 <- beta[dec$d12 + 1L]

        a1 <- U1 - m1[ii]
        a2 <- U2 - m2[ii]
        P12 <- F_vec(a1, a2, rho)

        ok1 <- is.finite(L1)
        if (any(ok1)) {
          b1 <- L1[ok1] - m1[ii][ok1]
          P12[ok1] <- P12[ok1] - F_vec(b1, a2[ok1], rho)
        }

        ok2 <- is.finite(L2)
        if (any(ok2)) {
          b2 <- L2[ok2] - m2[ii][ok2]
          P12[ok2] <- P12[ok2] - F_vec(a1[ok2], b2, rho)
        }

        ok12 <- ok1 & ok2
        if (any(ok12)) {
          b1 <- L1[ok12] - m1[ii][ok12]
          b2 <- L2[ok12] - m2[ii][ok12]
          P12[ok12] <- P12[ok12] + F_vec(b1, b2, rho)
        }

        P[ii] <- P12
      }

      P[P < eps] <- eps
      P[P > 1 - eps] <- 1 - eps
      P
    }

    if (!is.matrix(mu1)) {
      m1 <- as.numeric(mu1)
      m2 <- as.numeric(mu2)
      if (length(m1) != n) m1 <- rep(m1, length.out = n)
      if (length(m2) != n) m2 <- rep(m2, length.out = n)
      P <- one_eval(m1, m2)
      out <- log(P)
      attr(out, "mu_y") <- P
      return(out)
    }

    nGH <- ncol(mu1)
    if (!is.matrix(mu2)) mu2 <- matrix(mu2, nrow = n, ncol = nGH)

    out <- matrix(0, n, nGH)
    Pmat <- matrix(1, n, nGH)

    for (k in seq_len(nGH)) {
      Pk <- one_eval(mu1[, k], mu2[, k])
      out[, k] <- log(Pk)
      Pmat[, k] <- Pk
    }

    attr(out, "mu_y") <- Pmat
    out
  }

  structure(
    list(
      family = "bivariate cumulative probit",
      link = stats$name,
      linkfun = stats$linkfun,
      linkinv = stats$linkinv,
      log_dens = log_dens
    ),
    class = "family"
  )
}

get_K <- function(data, items) {
  setNames(
    vapply(items, function(nn) {
      x <- data[[nn]]
      if (is.ordered(x) || is.factor(x)) nlevels(x) else max(as.integer(x), na.rm = TRUE)
    }, integer(1)),
    items
  )
}


make_pair_data_biv_code <- function(dat2, items) {
  stopifnot(length(items) == 2L)

  item1 <- items[1]
  item2 <- items[2]
  K_item <- get_K(dat2, items)
  K1 <- K_item[1]
  K2 <- K_item[2]

  dp <- dat2 %>%
    dplyr::select(
      ID, visit_idx, I0, I6, I12, Ipost,
      waiting, prost_uni, loghosp_post, age, sex01, diag_oa,
      dplyr::all_of(items)
    ) %>%
    dplyr::mutate(
      y1_int = as.integer(.data[[item1]]),
      y2_int = as.integer(.data[[item2]]),
      y_code = dplyr::case_when(
        !is.na(y1_int) & is.na(y2_int) ~ y1_int,
        is.na(y1_int) & !is.na(y2_int) ~ K1 + y2_int,
        !is.na(y1_int) & !is.na(y2_int) ~ K1 + K2 + (y2_int - 1L) * K1 + y1_int,
        TRUE ~ NA_integer_
      ),
      wI6 = waiting * I6,
      wI12 = waiting * I12,
      pI6 = prost_uni * I6,
      pI12 = prost_uni * I12,

      z1 = as.integer(!is.na(y1_int)),
      z2 = as.integer(!is.na(y2_int)),

      z1I0 = z1 * I0,
      z1Ipost = z1 * Ipost,
      z2I0 = z2 * I0,
      z2Ipost = z2 * Ipost
    ) %>%
    dplyr::filter(!is.na(y_code))

  attr(dp, "K_item") <- K_item
  dp
}

.start_zeta_vec <- function(x, K) {
  y_fac <- ordered(as.integer(x), levels = seq_len(K))
  fit0 <- tryCatch(
    MASS::polr(y_fac ~ 1, method = "probit", Hess = FALSE),
    error = function(e) NULL,
    warning = function(w) invokeRestart("muffleWarning")
  )
  if (!is.null(fit0)) return(as.numeric(fit0$zeta))
  qnorm(seq(1 / K, (K - 1) / K, length.out = K - 1L))
}

start_zetas_pair_wide <- function(dp, items, K_item) {
  c(
    .start_zeta_vec(dp[[items[1]]], K_item[1]),
    .start_zeta_vec(dp[[items[2]]], K_item[2])
  )
}

extract_re_summary <- function(fit) {
  V <- tryCatch(VarCorr(fit)$ID, error = function(e) NULL)
  if (is.null(V)) V <- tryCatch(fit$D, error = function(e) NULL)
  if (is.null(V)) {
    return(list(
      sds = c(z1I0 = NA_real_, z1Ipost = NA_real_, z2I0 = NA_real_, z2Ipost = NA_real_),
      corr = c(item1_I0_vs_Ipost = NA_real_, item2_I0_vs_Ipost = NA_real_,
               cross_I0 = NA_real_, cross_Ipost = NA_real_)
    ))
  }

  V <- as.matrix(V)
  nm <- rownames(V)

  map <- c(
    z1I0 = "z1I0",
    z1Ipost = "z1Ipost",
    z2I0 = if ("zi_z2I0" %in% nm) "zi_z2I0" else "z2I0",
    z2Ipost = if ("zi_z2Ipost" %in% nm) "zi_z2Ipost" else "z2Ipost"
  )

  get_sd <- function(k) {
    kk <- map[[k]]
    if (!(kk %in% nm)) return(NA_real_)
    sqrt(V[kk, kk])
  }

  get_cor <- function(a, b) {
    aa <- map[[a]]
    bb <- map[[b]]
    if (!(aa %in% nm) || !(bb %in% nm)) return(NA_real_)
    V[aa, bb] / sqrt(V[aa, aa] * V[bb, bb])
  }

  list(
    sds = c(
      z1I0 = get_sd("z1I0"),
      z1Ipost = get_sd("z1Ipost"),
      z2I0 = get_sd("z2I0"),
      z2Ipost = get_sd("z2Ipost")
    ),
    corr = c(
      item1_I0_vs_Ipost = get_cor("z1I0", "z1Ipost"),
      item2_I0_vs_Ipost = get_cor("z2I0", "z2Ipost"),
      cross_I0 = get_cor("z1I0", "z2I0"),
      cross_Ipost = get_cor("z1Ipost", "z2Ipost")
    )
  )
}

make_initial_values_biv <- function(dp, items, K_item, start_from = NULL, initial_values = NULL) {
  p <- 10L
  n_phis <- sum(K_item - 1L) + 1L

  out <- list(
    phis = c(start_zetas_pair_wide(dp, items, K_item), 0)
  )

  if (!is.null(start_from)) {
    fit0 <- if (!is.null(start_from$fit)) start_from$fit else start_from

    if (!is.null(fit0$D) && is.matrix(fit0$D) && all(dim(fit0$D) == c(4L, 4L))) {
      out$D <- fit0$D
    }

    if (!is.null(fit0$phis)) {
      if (length(fit0$phis) == n_phis) out$phis <- unname(fit0$phis)
      if (length(fit0$phis) == n_phis - 1L) out$phis <- c(unname(fit0$phis), 0)
    }

    if (!is.null(fit0$coefficients)) {
      nm <- names(fit0$coefficients)

      old_covars <- c("I6","I12","wI6","wI12","pI6","pI12","age","sex01","diag_oa")
      new_covars <- c("I6","I12","wI6","wI12","pI6","pI12","loghosp_post","age","sex01","diag_oa")

      nm1_old <- paste0("item", items[1], ":", old_covars)
      nm2_old <- paste0("item", items[2], ":", old_covars)

      if (!is.null(nm) && all(nm1_old %in% nm) && all(nm2_old %in% nm)) {
        b1 <- setNames(rep(0, length(new_covars)), new_covars)
        b2 <- setNames(rep(0, length(new_covars)), new_covars)
        b1[old_covars] <- fit0$coefficients[nm1_old]
        b2[old_covars] <- fit0$coefficients[nm2_old]
        out$betas <- unname(b1)
        out$gammas <- unname(b2)
      } else {
        if (length(fit0$coefficients) == p) out$betas <- unname(fit0$coefficients)
        if (!is.null(fit0$gammas) && length(fit0$gammas) == p) out$gammas <- unname(fit0$gammas)
        if (is.null(out$gammas) && !is.null(fit0$coefficients_zi) && length(fit0$coefficients_zi) == p) {
          out$gammas <- unname(fit0$coefficients_zi)
        }
      }
    }
  }

  if (!is.null(initial_values)) {
    out <- modifyList(out, initial_values)
  }

  if (!is.null(out$betas) && length(out$betas) != p) stop("betas must be of length 10.")
  if (!is.null(out$gammas) && length(out$gammas) != p) stop("gammas  must be of length 10.")
  if (!is.null(out$D) && (!is.matrix(out$D) || any(dim(out$D) != c(4L, 4L)))) stop("D  must be 4x4.")
  if (length(out$phis) != n_phis) stop("phis has wrong length.")

  out
}

fit_one_pair_biv_cor <- function(dat2, items,
                                 ctrl = list(iter_EM = 0, nAGQ = 1, verbose = FALSE),
                                 start_from = NULL,
                                 initial_values = NULL) {
  stopifnot(length(items) == 2L)

  dp <- make_pair_data_biv_code(dat2, items)
  K_item <- attr(dp, "K_item")
  fam <- biv_cumulative_probit_glmmadaptive(K_item)

  fixed_fml <- y_code ~ 0 + I6 + I12 + wI6 + wI12 +
    pI6 + pI12 + loghosp_post + age + sex01 + diag_oa

  random_fml <- ~ 0 + z1I0 + z1Ipost | ID

  zi_fixed_fml <- ~ 0 + I6 + I12 + wI6 + wI12 +
    pI6 + pI12 + loghosp_post + age + sex01 + diag_oa

  zi_random_fml <- ~ 0 + z2I0 + z2Ipost | ID

  inits <- make_initial_values_biv(
    dp = dp,
    items = items,
    K_item = K_item,
    start_from = start_from,
    initial_values = initial_values
  )

  fit <- tryCatch(
    GLMMadaptive::mixed_model(
      fixed = fixed_fml,
      random = random_fml,
      data = dp,
      family = fam,
      zi_fixed = zi_fixed_fml,
      zi_random = zi_random_fml,
      n_phis = sum(K_item - 1L) + 1L,
      initial_values = inits,
      control = modifyList(list(iter_EM = 0, nAGQ = 11, verbose = FALSE), ctrl)
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(
      ok = FALSE,
      pair = tibble::tibble(
        item1 = items[1],
        item2 = items[2],
        logLik = NA_real_,
        AIC = NA_real_,
        sd1I0 = NA_real_,
        sd1Ipost = NA_real_,
        sd2I0 = NA_real_,
        sd2Ipost = NA_real_,
        corr_1_I0Ipost = NA_real_,
        corr_2_I0Ipost = NA_real_,
        corr_cross_I0 = NA_real_,
        corr_cross_Ipost = NA_real_,
        rho_resid = NA_real_,
        conv = FALSE,
        msg = conditionMessage(fit)
      ),
      fit = NULL
    ))
  }

  re <- extract_re_summary(fit)
  s <- re$sds
  crr <- re$corr
  rho_hat <- tanh(tail(fit$phis, 1))

  out_row <- tibble::tibble(
    item1 = items[1],
    item2 = items[2],
    logLik = as.numeric(logLik(fit)),
    AIC = AIC(fit),
    sd1I0 = unname(s["z1I0"]),
    sd1Ipost = unname(s["z1Ipost"]),
    sd2I0 = unname(s["z2I0"]),
    sd2Ipost = unname(s["z2Ipost"]),
    corr_1_I0Ipost = unname(crr["item1_I0_vs_Ipost"]),
    corr_2_I0Ipost = unname(crr["item2_I0_vs_Ipost"]),
    corr_cross_I0 = unname(crr["cross_I0"]),
    corr_cross_Ipost = unname(crr["cross_Ipost"]),
    rho_resid = rho_hat,
    conv = isTRUE(fit$converged),
    msg = NA_character_
  )

  list(ok = isTRUE(fit$converged), pair = out_row, fit = fit)
}

# UTILITIES
.safe_saveRDS <- function(object, file, compress = "xz") {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(file, ".tmp_", sprintf("%06d", sample.int(1e6, 1)))
  on.exit(try(unlink(tmp), silent = TRUE), add = TRUE)
  saveRDS(object, tmp, compress = compress)
  file.rename(tmp, file)
}


#  PARALLEL PER COPPIA con progress bar + NOMI + conv + log + resume
fit_all_pairs_parallel <- function(dat2, itemnames,
                                   ctrl = list(iter_EM = 5, nAGQ = 11),
                                   out_dir = "pairwise_fits",
                                   resume = TRUE,
                                   n_cores = max(1L, parallel::detectCores(logical = TRUE) - 1L),
                                   dat2_path = NULL,     # se dato, i worker leggono da file
                                   progress_log = file.path(out_dir, "_progress.log"),
                                   print_names = TRUE) {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(dat2_path)) dat2_path <- normalizePath(dat2_path, mustWork = TRUE)
  
  pairs <- combn(itemnames, 2, simplify = FALSE)
  N <- length(pairs)
  
  .log_line <- function(...) {
    line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", sprintf(...), "\n")
    cat(line, file = progress_log, append = TRUE)
    invisible(NULL)
  }
  
  stubs <- vapply(seq_len(N), function(i) {
    items <- pairs[[i]]
    sprintf("%03d_%s__%s.rds",
            i,
            gsub("[^A-Za-z0-9]+", "_", items[1]),
            gsub("[^A-Za-z0-9]+", "_", items[2]))
  }, character(1))
  fpaths <- file.path(out_dir, stubs)
  todo_idx <- which(!(resume & file.exists(fpaths)))
  skip_idx <- setdiff(seq_len(N), todo_idx)
  
  if (print_names && length(skip_idx)) {
    for (k in skip_idx) {
      items <- pairs[[k]]
      cat(sprintf("[SKIP] %03d/%03d :: %s vs %s (already done)\n", k, N, items[1], items[2]))
    }
    flush.console()
  }
  
  if (length(todo_idx) == 0L) {
    res_loaded <- lapply(fpaths, readRDS)
    summary_tbl <- dplyr::bind_rows(lapply(res_loaded, `[[`, "pair"))
    fits_list   <- setNames(lapply(res_loaded, `[[`, "fit"),
                            vapply(pairs, paste, character(1), collapse = " vs "))
    return(list(summary = summary_tbl, fits = fits_list, out_dir = out_dir,
                progress_log = progress_log))
  }
  
  has_future   <- requireNamespace("future.apply", quietly = TRUE)
  has_progress <- requireNamespace("progressr",    quietly = TRUE)
  
  if (has_future && has_progress) {
    future::plan(future::multisession, workers = min(n_cores, length(todo_idx)))
    on.exit({ try(future::plan(future::sequential), silent = TRUE) }, add = TRUE)
    
    if (requireNamespace("progress", quietly = TRUE)) {
      progressr::handlers("progress")
    } else {
      progressr::handlers("txtprogressbar")
    }
    progressr::handlers(global = TRUE)
    
    worker_fun <- function(i) {
      suppressWarnings(suppressPackageStartupMessages({
        library(GLMMadaptive); library(MASS)
        library(dplyr); library(tidyr); library(tibble); library(stringr)
        library(mvtnorm)  ## *** MOD ***
      }))
      items   <- pairs[[i]]
      fpath   <- fpaths[i]
      errfile <- sub("\\.rds$", ".err.txt", fpath)
      
      if (resume && file.exists(fpath)) {
        out <- tryCatch(readRDS(fpath), error = function(e) NULL)
        if (!is.null(out)) return(out)
      }
      
      .log_line("START %03d :: %s vs %s", i, items[1], items[2])
      
      out <- tryCatch({
        dat2_local <- if (!is.null(dat2_path)) readRDS(dat2_path) else dat2
        fit_one_pair_biv_cor(dat2_local, items, ctrl = ctrl) 
      }, error = function(e) {
        tib <- tibble::tibble(
          item1 = items[1], item2 = items[2],
          logLik = NA_real_, AIC = NA_real_,
          sd1I0 = NA_real_, sd1Ipost = NA_real_, sd2I0 = NA_real_, sd2Ipost = NA_real_,
          corr_1_I0Ipost = NA_real_, corr_2_I0Ipost = NA_real_,
          corr_cross_I0 = NA_real_, corr_cross_Ipost = NA_real_,
          conv = FALSE, msg = conditionMessage(e)
        )
        list(ok = FALSE, pair = tib, fit = NULL)
      })
      
      .safe_saveRDS(out, fpath, compress = "xz")
      if (!isTRUE(out$ok)) cat(paste0(out$pair$msg, "\n"), file = errfile, append = FALSE)
      
      .log_line("DONE  %03d :: %s vs %s (conv=%s)", i, items[1], items[2],
                ifelse(isTRUE(out$ok) && isTRUE(out$pair$conv), "TRUE", "FALSE"))
      out
    }
    
    res <- progressr::with_progress({
      p <- progressr::progressor(steps = length(todo_idx))
      future.apply::future_lapply(todo_idx, function(i) {
        items <- pairs[[i]]
        lbl <- sprintf("%03d/%03d :: %s vs %s", i, N, items[1], items[2])
        p(lbl)
        if (print_names) { cat("[START] ", lbl, "\n", sep = ""); flush.console() }
        out <- worker_fun(i)
        conv <- ifelse(isTRUE(out$ok) && isTRUE(out$pair$conv), "TRUE", "FALSE")
        aic  <- suppressWarnings(tryCatch(as.numeric(out$pair$AIC), error = function(e) NA_real_))
        if (print_names) {
          cat(sprintf("[DONE ] %03d/%03d :: %s vs %s  conv=%s  AIC=%s\n",
                      i, N, items[1], items[2], conv, ifelse(is.na(aic), "NA", sprintf("%.1f", aic))))
          flush.console()
        }
        out
      }, future.seed = TRUE)
    })
    
    all_res <- vector("list", N)
    if (length(skip_idx)) all_res[skip_idx] <- lapply(fpaths[skip_idx], readRDS)
    all_res[todo_idx] <- res
    
  } else {
    if (!requireNamespace("pbapply", quietly = TRUE))
      stop("Instal 'future.apply' + 'progressr' or at least 'pbapply' for progress bar.")
    
    cl <- parallel::makeCluster(min(n_cores, length(todo_idx)))
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    
    parallel::clusterEvalQ(cl, {
      library(dplyr)
      library(tidyr)
      library(GLMMadaptive)
      library(MASS)
      library(tibble)
      library(stringr)
      library(mvtnorm)
      library(pbivnorm)
      NULL
    })

    parallel::clusterExport(
      cl,
      varlist = c(
        "dat2","pairs","fpaths","ctrl","resume","dat2_path",".log_line",
        ".safe_saveRDS",
        "biv_cumulative_probit_glmmadaptive",
        "make_pair_data_biv_code",
        "get_K",
        "fit_one_pair_biv_cor",
        "start_zetas_pair_wide",
        ".start_zeta_vec",
        "extract_re_summary",
        "make_initial_values_biv"
      ),
      envir = environment()
    )
    
    if (print_names) {
      for (i in todo_idx) {
        items <- pairs[[i]]
        cat(sprintf("[QUEUE] %03d/%03d :: %s vs %s\n", i, N, items[1], items[2]))
      }
      flush.console()
    }
    
    res_todo <- pbapply::pblapply(todo_idx, cl = cl, FUN = function(i) {
      items   <- pairs[[i]]
      fpath   <- fpaths[i]
      errfile <- sub("\\.rds$", ".err.txt", fpath)
      
      if (resume && file.exists(fpath)) {
        out <- tryCatch(readRDS(fpath), error = function(e) NULL)
        if (!is.null(out)) return(out)
      }
      
      .log_line("START %03d :: %s vs %s", i, items[1], items[2])
      
      out <- tryCatch({
        dat2_local <- if (!is.null(dat2_path)) readRDS(dat2_path) else dat2
        fit_one_pair_biv_cor(dat2_local, items, ctrl = ctrl)
      }, error = function(e) {
        tib <- tibble::tibble(
          item1 = items[1], item2 = items[2],
          logLik = NA_real_, AIC = NA_real_,
          sd1I0 = NA_real_, sd1Ipost = NA_real_, sd2I0 = NA_real_, sd2Ipost = NA_real_,
          corr_1_I0Ipost = NA_real_, corr_2_I0Ipost = NA_real_,
          corr_cross_I0 = NA_real_, corr_cross_Ipost = NA_real_,
          conv = FALSE, msg = conditionMessage(e)
        )
        list(ok = FALSE, pair = tib, fit = NULL)
      })
      
      .safe_saveRDS(out, fpath, compress = "xz")
      if (!isTRUE(out$ok)) cat(paste0(out$pair$msg, "\n"), file = errfile, append = FALSE)
      
      .log_line("DONE  %03d :: %s vs %s (conv=%s)", i, items[1], items[2],
                ifelse(isTRUE(out$ok) && isTRUE(out$pair$conv), "TRUE", "FALSE"))
      out
    })
    
    all_res <- vector("list", N)
    if (length(skip_idx)) all_res[skip_idx] <- lapply(fpaths[skip_idx], readRDS)
    all_res[todo_idx] <- res_todo
  }
  
  summary_tbl <- dplyr::bind_rows(lapply(all_res, `[[`, "pair"))
  fits_list   <- setNames(lapply(all_res, `[[`, "fit"),
                          vapply(pairs, paste, character(1), collapse = " vs "))
  list(summary = summary_tbl, fits = fits_list, out_dir = out_dir, progress_log = progress_log)
}

# 
itemnamesEQ <- c('Mobility','SelfCare','UsualActivities','PainorDiscomfort','Anxietyordepression')
itemnamesKOOS <- c('Gettingoutofbed','Puttingonsocks','Gettingupfromsitting',
                   'Bendingtowardthefloor_pickinguponobjectfromtheground',
                   'Twisting_pivotingontheinjuredknee','Kneeling','Squatting')
itemnames <- c(itemnamesEQ, itemnamesKOOS)

res_pairs <- fit_all_pairs_parallel(
  dat2, itemnames,
  ctrl = list(iter_EM = 0, nAGQ = 11),
  out_dir = "pairwise_EQ_ER_11_bivar",
  resume = TRUE,
  n_cores = 66,
  print_names = TRUE
)




# model.start <- readRDS("pairwise_EQ_ER_11_for_failures/001_Mobility__SelfCare.rds")

# res_biv <- fit_one_pair_biv_cor(
#   dat2,
#   items = c("Mobility","SelfCare"),
#   ctrl = list(iter_EM = 0, nAGQ = 11, verbose = TRUE)
# )
