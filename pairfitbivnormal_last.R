library(dplyr)
library(MASS)
library(GLMMadaptive)

dataknee <- readRDS("dat2_with_hosp30g.rds")

itemnames <- c("Mobility","SelfCare","UsualActivities","PainorDiscomfort","Anxietyordepression",
               "Gettingoutofbed","Puttingonsocks","Gettingupfromsitting",
               "Bendingtowardthefloor_pickinguponobjectfromtheground",
               "Twisting_pivotingontheinjuredknee","Kneeling","Squatting")

dataknee[itemnames] <- data.frame(lapply(dataknee[itemnames], as.ordered))

dataknee$ID <- factor(dataknee$ID)
dataknee$sex01 <- as.integer(dataknee$sex == "M")
dataknee$prost_uni <- as.integer(dataknee$prosthesis_type == "unicompartimental")
dataknee$diag_oa <- as.integer(dataknee$diagnosis != "osteoarthritis_knee")


dat2 <- dat2[dat2$residence == "Emilia-Romagna", ]

dat2 <- dataknee %>%
  group_by(ID) %>%
  arrange(times, .by_group = TRUE) %>%
  mutate(
    t_within = times - first(times),
    visit_idx = row_number() - 1L,
    post = as.integer(visit_idx >= 1L),
    time_post = ifelse(post == 1L, t_within, 0),
    loghosp = log(1 + hospitalizationdays),
    waiting = as.numeric(difftime(date_intervention, date_baseline, units = "days"))/30.4375,
    time_discrete = factor(visit_idx, levels = c(0,1,2), labels = c("t0","t6","t12")),
    I0 = as.integer(time_discrete == "t0"),
    I6 = as.integer(time_discrete == "t6"),
    I12 = as.integer(time_discrete == "t12"),
    Ipost = I6 + I12,
    loghosp_post = loghosp * Ipost
  ) %>%
  ungroup()


age_by_id <- dat2 %>%
  distinct(ID, age) %>%
  mutate(age_sc = as.numeric(scale(age)))

dat2 <- dat2 %>%
  select(-age) %>%
  left_join(age_by_id, by = "ID") %>%
  rename(age = age_sc)
  
event_by_id <- dat2 %>%
  transmute(ID, t_event_days = as.numeric(difftime(date_revision, date_intervention, units="days"))) %>%
  mutate(t_event_days = ifelse(is.na(t_event_days) | t_event_days < 0, Inf, t_event_days)) %>%
  group_by(ID) %>%
  summarise(T_event_days = min(t_event_days), .groups="drop")

dat2 <- dat2 %>%
  dplyr::left_join(event_by_id, by="ID") %>%
  dplyr::filter(
    visit_idx == 0L |
      is.infinite(T_event_days) |
      (visit_idx == 1L & T_event_days > 6*30.4375) |
      (visit_idx == 2L & T_event_days > 12*30.4375)
  ) %>%
  dplyr::select(-T_event_days)


make_pair_data_ord_norm_code <- function(dat2, ord_item, cont_item, scale100 = TRUE, eps_y = 0) {
  if (!(is.factor(dat2[[ord_item]]) || is.ordered(dat2[[ord_item]]))) dat2[[ord_item]] <- as.ordered(dat2[[ord_item]])
  K_ord <- nlevels(dat2[[ord_item]])
  cont_shift <- K_ord + 1

  dp <- dat2 %>%
    dplyr::select(
      ID, visit_idx, I0, I6, I12, Ipost,
      waiting, prost_uni, loghosp_post, age, sex01, diag_oa,
      dplyr::all_of(c(ord_item, cont_item))
    ) %>%
    dplyr::mutate(
      y_ord  = as.integer(as.ordered(.data[[ord_item]])),
      y_cont = as.numeric(.data[[cont_item]]) / if (scale100) 100 else 1,
      y_cont = ifelse(is.na(y_cont), NA_real_, pmin(pmax(y_cont, 0), 1 - eps_y)),
      wI6  = waiting * I6,
      wI12 = waiting * I12,
      pI6  = prost_uni * I6,
      pI12 = prost_uni * I12,
      y_mix = dplyr::case_when(
        !is.na(y_ord) &  is.na(y_cont) ~ as.numeric(y_ord),
         is.na(y_ord) & !is.na(y_cont) ~ cont_shift + y_cont,
        !is.na(y_ord) & !is.na(y_cont) ~ cont_shift + y_ord + y_cont,
        TRUE                           ~ NA_real_
      ),
      z1 = as.integer(!is.na(y_ord)),
      z2 = as.integer(!is.na(y_cont)),
      z1I0    = z1 * I0,
      z1Ipost = z1 * Ipost,
      z2I0    = z2 * I0,
      z2Ipost = z2 * Ipost
    ) %>%
    dplyr::filter(!is.na(y_mix))

  attr(dp, "K_ord") <- K_ord
  attr(dp, "cont_shift") <- cont_shift
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

start_phis_ord_norm_rho <- function(dp, ord_item, K_ord, cont_shift) {
  zeta0 <- .start_zeta_vec(dp[[ord_item]], K_ord)

  tmp <- dp$y_mix[dp$z2 == 1] - cont_shift
  ycont <- tmp - floor(tmp)
  ycont <- ycont[is.finite(ycont)]
  s0 <- if (length(ycont) >= 5L) stats::sd(ycont) else 0.2
  s0 <- max(s0, 1e-3)
  logsigma0 <- log(s0)

  c(zeta0, logsigma0, 0)
}

biv_probit_plus_gaussian_rho <- function(K_ord, cont_shift = K_ord + 1, link = "identity", sigma_min = 1e-6) {
  stats <- make.link(link)
  eps <- .Machine$double.eps

  log_dens <- function(y, eta, mu_fun, phis, eta_zi) {
    y <- as.numeric(y)
    mu1 <- mu_fun(eta)
    mu2 <- mu_fun(eta_zi)

    n_th <- K_ord - 1L
    if (length(phis) != (n_th + 2L)) stop(sprintf("length(phis) deve essere %d.", n_th + 2L))

    zeta <- phis[seq_len(n_th)]
    logsigma <- phis[n_th + 1L]
    z_rho <- phis[n_th + 2L]

    sigma <- sigma_min + exp(logsigma)
    rho <- tanh(z_rho)
    s <- sqrt(1 - rho^2)

    ok <- !is.na(y)
    if (is.matrix(mu1)) {
      Q <- ncol(mu1)
      out <- matrix(0, nrow = length(y), ncol = Q)
    } else {
      out <- numeric(length(y))
    }
    if (!any(ok)) return(out)

    yy <- y[ok]
    idx_ok <- which(ok)

    is_ord  <- yy < cont_shift
    is_cont <- (yy >= cont_shift) & (yy < cont_shift + 1)
    is_both <- yy >= (cont_shift + 1)

    calc_vec <- function(et1, et2) {
      ll <- numeric(length(y))
      ll[!ok] <- 0

      if (any(is_ord)) {
        idx <- idx_ok[is_ord]
        ycat <- as.integer(round(yy[is_ord]))
        if (any(ycat < 1L | ycat > K_ord)) stop("ordinale fuori range.")

        theta_u <- ifelse(ycat < K_ord, zeta[ycat],  Inf)
        theta_l <- ifelse(ycat > 1L,    zeta[ycat - 1L], -Inf)

        P <- pnorm(theta_u - et1[idx]) - pnorm(theta_l - et1[idx])
        P[P < eps] <- eps
        ll[idx] <- log(P)
      }

      if (any(is_cont)) {
        idx <- idx_ok[is_cont]
        ycont <- yy[is_cont] - cont_shift
        dens <- dnorm(ycont, mean = et2[idx], sd = sigma, log = TRUE)
        dens[!is.finite(dens)] <- log(eps)
        ll[idx] <- dens
      }

      if (any(is_both)) {
        idx <- idx_ok[is_both]
        tmp <- yy[is_both] - cont_shift
        k <- floor(tmp)
        ycont <- tmp - k
        k <- as.integer(k)

        if (any(k < 1L | k > K_ord)) stop("k fuori range nel caso both.")
        ycont <- pmin(pmax(ycont, 0), 1 - 1e-12)

        v <- (ycont - et2[idx]) / sigma

        theta_u <- ifelse(k < K_ord, zeta[k],  Inf)
        theta_l <- ifelse(k > 1L,    zeta[k - 1L], -Inf)

        uu <- (theta_u - et1[idx] - rho * v) / s
        llow <- (theta_l - et1[idx] - rho * v) / s
        Pk <- pnorm(uu) - pnorm(llow)
        Pk[Pk < eps] <- eps

        ll_y <- dnorm(v, mean = 0, sd = 1, log = TRUE) - log(sigma)
        ll[idx] <- ll_y + log(Pk)
      }

      ll
    }

    if (!is.matrix(mu1)) {
      et1 <- as.numeric(mu1); if (length(et1) != length(y)) et1 <- rep(et1, length.out = length(y))
      et2 <- as.numeric(mu2); if (length(et2) != length(y)) et2 <- rep(et2, length.out = length(y))
      return(calc_vec(et1, et2))
    }

    Q <- ncol(mu1)
    if (!is.matrix(mu2)) mu2 <- matrix(mu2, nrow = length(y), ncol = Q)
    for (q in seq_len(Q)) out[, q] <- calc_vec(mu1[, q], mu2[, q])
    out
  }

  structure(list(
    family = "ordinal probit + gaussian with resid rho",
    link = stats$name,
    linkfun = stats$linkfun,
    linkinv = stats$linkinv,
    log_dens = log_dens
  ), class = "family")
}

.default_inits_ord_norm_rho <- function(dp, ord_item, fixed_fml, zi_fixed_fml, K_ord, cont_shift) {
  X <- model.matrix(delete.response(terms(fixed_fml)), dp)
  X_zi <- model.matrix(terms(zi_fixed_fml), dp)
  phis0 <- start_phis_ord_norm_rho(dp, ord_item, K_ord, cont_shift)
  list(
    betas = rep(0, ncol(X)),
    gammas = rep(0, ncol(X_zi)),
    D = diag(0.1, 4),
    phis = phis0
  )
}

fit_one_pair_ord_norm_cor <- function(dat2, ord_item, cont_item,
                                      ctrl = list(iter_EM = 0, nAGQ = 3, verbose = FALSE),
                                      scale100 = TRUE) {
  dp <- make_pair_data_ord_norm_code(dat2, ord_item, cont_item, scale100 = scale100)
  K_ord <- attr(dp, "K_ord")
  cont_shift <- attr(dp, "cont_shift")

  fam <- biv_probit_plus_gaussian_rho(K_ord = K_ord, cont_shift = cont_shift, link = "identity")
  n_phis <- (K_ord - 1L) + 2L

  fixed_fml  <- y_mix ~ 0 + I6 + I12 + wI6 + wI12 + pI6 + pI12 + loghosp_post + age + sex01 + diag_oa
  random_fml <- ~ 0 + z1I0 + z1Ipost | ID

  zi_fixed_fml  <- ~ 1 + I6 + I12 + wI6 + wI12 + pI6 + pI12 + loghosp_post + age + sex01 + diag_oa
  zi_random_fml <- ~ 0 + z2I0 + z2Ipost | ID

  inits <- .default_inits_ord_norm_rho(dp, ord_item, fixed_fml, zi_fixed_fml, K_ord, cont_shift)

  fit <- GLMMadaptive::mixed_model(
    fixed = fixed_fml,
    random = random_fml,
    data = dp,
    family = fam,
    zi_fixed = zi_fixed_fml,
    zi_random = zi_random_fml,
    n_phis = n_phis,
    initial_values = inits,
    control = ctrl
  )

  list(
    fit = fit,
    rho_resid = tanh(fit$phis[n_phis]),
    sigma = 1e-6 + exp(fit$phis[(K_ord - 1L) + 1L])
  )
}


.safe_saveRDS <- function(object, file, compress = "xz") {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(file, ".tmp_", sprintf("%06d", sample.int(1e6, 1)))
  on.exit(try(unlink(tmp), silent = TRUE), add = TRUE)
  saveRDS(object, tmp, compress = compress)
  file.rename(tmp, file)
}

#wrap:converte il tuo output (list(fit=...,rho_resid=...,sigma=...)) in list(ok,pair,fit)
fit_one_ord_norm_wrap <- function(dat2, ord_item, cont_item,
                                 ctrl = list(iter_EM = 0, nAGQ = 3, verbose = FALSE),
                                 scale100 = TRUE,
                                 fit_verbose_to_file = FALSE,
                                 fitlog_path = NULL) {
  out <- tryCatch({
    if (isTRUE(fit_verbose_to_file) && !is.null(fitlog_path)) {
      zz <- file(fitlog_path, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      on.exit({
        try(sink(type = "message"), silent = TRUE)
        try(sink(type = "output"),  silent = TRUE)
        try(close(zz), silent = TRUE)
      }, add = TRUE)
      ctrl$verbose <- TRUE
    }

    rr <- fit_one_pair_ord_norm_cor(
      dat2 = dat2,
      ord_item = ord_item,
      cont_item = cont_item,
      ctrl = ctrl,
      scale100 = scale100
    )

    fit <- rr$fit
    pair_row <- tibble::tibble(
      ord_item = ord_item,
      cont_item = cont_item,
      logLik = suppressWarnings(as.numeric(logLik(fit))),
      AIC = suppressWarnings(AIC(fit)),
      sigma = rr$sigma,
      rho_resid = rr$rho_resid,
      conv = isTRUE(fit$converged),
      msg = NA_character_
    )

    list(ok = isTRUE(fit$converged), pair = pair_row, fit = fit)
  }, error = function(e) {
    pair_row <- tibble::tibble(
      ord_item = ord_item,
      cont_item = cont_item,
      logLik = NA_real_,
      AIC = NA_real_,
      sigma = NA_real_,
      rho_resid = NA_real_,
      conv = FALSE,
      msg = conditionMessage(e)
    )
    list(ok = FALSE, pair = pair_row, fit = NULL)
  })

  out
}

fit_all_ord_norm_parallel <- function(dat2,
                                      ord_items,
                                      cont_items,
                                      ctrl = list(iter_EM = 0, nAGQ = 3, verbose = FALSE),
                                      out_dir = "ord_norm_fits",
                                      resume = TRUE,
                                      n_cores = max(1L, parallel::detectCores(logical = TRUE) - 1L),
                                      progress_log = file.path(out_dir, "_progress.log"),
                                      print_names = TRUE,
                                      fit_verbose_to_file = FALSE,
                                      scale100 = TRUE) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  grid <- expand.grid(ord_item = ord_items, cont_item = cont_items, stringsAsFactors = FALSE)
  N <- nrow(grid)

  .log_line <- function(...) {
    line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", sprintf(...), "\n")
    cat(line, file = progress_log, append = TRUE)
    invisible(NULL)
  }

  stubs <- vapply(seq_len(N), function(i) {
    o <- gsub("[^A-Za-z0-9]+", "_", grid$ord_item[i])
    c <- gsub("[^A-Za-z0-9]+", "_", grid$cont_item[i])
    sprintf("%03d_%s__%s.rds", i, o, c)
  }, character(1))

  fpaths <- file.path(out_dir, stubs)
    exists_vec <- file.exists(fpaths)
  if (isTRUE(resume)) {
    todo_idx <- which(!exists_vec)
    skip_idx <- which(exists_vec)
  } else {
    todo_idx <- seq_len(N)
    skip_idx <- integer(0)
  }

  if (print_names && length(skip_idx)) {
    for (k in skip_idx) {
      cat(sprintf("[SKIP] %03d/%03d :: %s + %s\n", k, N, grid$ord_item[k], grid$cont_item[k]))
    }
    flush.console()
  }

  if (length(todo_idx) == 0L) {
    res_loaded <- lapply(fpaths, readRDS)
    summary_tbl <- dplyr::bind_rows(lapply(res_loaded, `[[`, "pair"))
    fits_list <- setNames(lapply(res_loaded, `[[`, "fit"),
                          paste(grid$ord_item, grid$cont_item, sep = " + "))
    return(list(summary = summary_tbl, fits = fits_list, out_dir = out_dir, progress_log = progress_log))
  }

  has_future <- requireNamespace("future.apply", quietly = TRUE)
  has_progress <- requireNamespace("progressr", quietly = TRUE)

  if (!has_future || !has_progress) {
    stop("per avere progress bar in parallelo installa future.apply e progressr.")
  }

  future::plan(future::multisession, workers = min(n_cores, length(todo_idx)))
  on.exit({ try(future::plan(future::sequential), silent = TRUE) }, add = TRUE)
  progressr::handlers("txtprogressbar")
  progressr::handlers(global = TRUE)

 #dat2 catturato qui, prima dei worker
dat2_local <- dat2

worker_fun <- function(i) {
  suppressWarnings(suppressPackageStartupMessages({
    library(GLMMadaptive)
    library(dplyr)
    library(tibble)
  }))

  ord_item <- grid$ord_item[i]
  cont_item <- grid$cont_item[i]
  fpath <- fpaths[i]
  errfile <- sub("\\.rds$", ".err.txt", fpath)
  fitlog <- sub("\\.rds$", ".fitlog.txt", fpath)

  if (isTRUE(resume) && file.exists(fpath)) {
    out <- tryCatch(readRDS(fpath), error = function(e) NULL)
    if (!is.null(out)) return(out)
  }

  .log_line("START %03d :: %s + %s", i, ord_item, cont_item)

  out <- fit_one_ord_norm_wrap(
    dat2 = dat2_local,
    ord_item = ord_item,
    cont_item = cont_item,
    ctrl = ctrl,
    scale100 = scale100,
    fit_verbose_to_file = fit_verbose_to_file,
    fitlog_path = fitlog
  )

  .safe_saveRDS(out, fpath, compress = "xz")
  if (!isTRUE(out$ok)) cat(paste0(out$pair$msg, "\n"), file = errfile, append = FALSE)

  aic <- suppressWarnings(as.numeric(out$pair$AIC))
  .log_line("DONE  %03d :: %s + %s (conv=%s AIC=%s)", i, ord_item, cont_item,
            ifelse(isTRUE(out$pair$conv), "TRUE", "FALSE"),
            ifelse(is.na(aic), "NA", sprintf("%.1f", aic)))
  out
}

  res <- progressr::with_progress({
    p <- progressr::progressor(steps = length(todo_idx))
    future.apply::future_lapply(todo_idx, function(i) {
      lbl <- sprintf("%03d/%03d :: %s + %s", i, N, grid$ord_item[i], grid$cont_item[i])
      p(lbl)
      if (print_names) { cat("[START] ", lbl, "\n", sep = ""); flush.console() }
      out <- worker_fun(i)
      conv <- ifelse(isTRUE(out$pair$conv), "TRUE", "FALSE")
      aic <- suppressWarnings(as.numeric(out$pair$AIC))
      if (print_names) {
        cat(sprintf("[DONE ] %03d/%03d :: %s + %s  conv=%s  AIC=%s  rho=%s\n",
                    i, N, grid$ord_item[i], grid$cont_item[i], conv,
                    ifelse(is.na(aic), "NA", sprintf("%.1f", aic)),
                    ifelse(is.na(out$pair$rho_resid), "NA", sprintf("%.3f", out$pair$rho_resid))))
        flush.console()
      }
      out
    }, future.seed = TRUE)
  })

  all_res <- vector("list", N)
  if (length(skip_idx)) all_res[skip_idx] <- lapply(fpaths[skip_idx], readRDS)
  all_res[todo_idx] <- res

  summary_tbl <- dplyr::bind_rows(lapply(all_res, `[[`, "pair"))
  fits_list <- setNames(lapply(all_res, `[[`, "fit"),
                        paste(grid$ord_item, grid$cont_item, sep = " + "))
  list(summary = summary_tbl, fits = fits_list, out_dir = out_dir, progress_log = progress_log)
}


itemnamesEQ <- c('Mobility','SelfCare','UsualActivities','PainorDiscomfort','Anxietyordepression')
itemnamesKOOS <- c('Gettingoutofbed','Puttingonsocks','Gettingupfromsitting',
                   'Bendingtowardthefloor_pickinguponobjectfromtheground',
                   'Twisting_pivotingontheinjuredknee','Kneeling','Squatting')
itemnames <- c(itemnamesEQ, itemnamesKOOS)

# res <- fit_all_ord_norm_parallel(
#   dat2,
#   ord_items = itemnames,
#   cont_items = "Perceivedcurrenthealthstatus",
#   ctrl = list(iter_EM = 0, nAGQ = 11, verbose = FALSE),
#   out_dir = "pairwise_VAS_vs_ord_norm_resid",

#   resume = TRUE,
#   n_cores = 1,
#   print_names = TRUE,
#   fit_verbose_to_file = TRUE
# )


# #test singolo
res <- fit_one_pair_ord_norm_cor(
  dat2,
  ord_item  = "Twisting_pivotingontheinjuredknee",
  cont_item = "Perceivedcurrenthealthstatus",
  ctrl      = list(iter_EM = 0, nAGQ = 12, verbose = TRUE)
)
saveRDS(res,"/home/niccolo.cao/pairfit/pairwise_VAS_vs_ord_norm_resid/010_Twisting_pivotingontheinjuredknee__Perceivedcurrenthealthstatus.rds")
# res$rho_resid
# res$sigma
# summary(res$fit)
# res$rho_resid
# res$sigma
# summary(res$fit)

