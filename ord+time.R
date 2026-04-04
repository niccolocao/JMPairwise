#####################################################################
#################                              ######################
###############  Pair ordinal + time-to-event    ####################
#################                              ######################
#####################################################################
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(GLMMadaptive)
library(MASS)
library(survival)



cprobit_weibull_indexed <- function(K, type_vec, start_vec, stop_vec, status_vec, link="identity") {
  stats <- make.link(link); len_th <- K - 1L; eps <- 1e-12
  force(type_vec); force(start_vec); force(stop_vec); force(status_vec)
  
  log_dens <- function(Y, eta, mu_fun, phis, eta_zi) {
    Y <- as.matrix(Y)
    if (ncol(Y) != 2L) stop("Y deve essere n x 2: cbind(y_main, idx).")
    y_main <- as.numeric(Y[,1]); idx <- as.integer(Y[,2])
    if (length(phis) != len_th + 2L)
      stop(sprintf("phis attesi = %d (K-1 soglie + 2 Weibull).", len_th + 2L))
    
    mu  <- mu_fun(eta)
    n   <- length(y_main)
    out <- if (is.matrix(mu)) matrix(0, n, ncol(mu)) else numeric(n)
    
    type <- type_vec[idx]
    
    # ordinal
    io <- which(type == 0L)
    if (length(io)) {
      Ki <- K
      y_cat <- as.integer(y_main[io])
      bad <- which(!is.finite(y_cat) | y_cat < 1L | y_cat > K)
      if (length(bad))
        stop(sprintf("[ordinale] y fuori da 1..%d alle righe: %s", K, paste(io[bad], collapse=", ")))
      
      get_zeta <- function(c) phis[c]
      th_u <- rep( Inf, length(y_cat))
      iu <- which(y_cat < Ki); if (length(iu)) th_u[iu] <- vapply(y_cat[iu], get_zeta, numeric(1))
      th_l <- rep(-Inf, length(y_cat))
      il <- which(y_cat > 1L); if (length(il)) th_l[il] <- vapply(y_cat[il]-1L, get_zeta, numeric(1))
      
      if (is.matrix(mu)) {
        et <- mu[io,, drop=FALSE]
        P  <- pnorm(th_u - et) - pnorm(th_l - et)
        P[P < eps] <- eps; P[P > 1-eps] <- 1-eps
        out[io,] <- log(P)
      } else {
        et <- mu[io]
        P  <- pnorm(th_u - et) - pnorm(th_l - et)
        P[P < eps] <- eps; P[P > 1-eps] <- 1-eps
        out[io] <- log(P)
      }
    }
    
    #Weibull PH
    isv <- which(type == 1L)
    if (length(isv)) {
      start <- pmax(start_vec[idx[isv]], 0)
      stop  <- pmax(stop_vec[idx[isv]],  1e-12)
      stat  <- status_vec[idx[isv]]
      
      log_lambda <- phis[len_th + 1L]; log_rho <- phis[len_th + 2L]
      lambda <- exp(log_lambda); rho <- exp(log_rho)
      
      H0    <- lambda * ( exp(rho*log(stop)) - exp(rho*log(pmax(start, 1e-12))) )
      logh0 <- log_lambda + log_rho + (rho - 1)*log(stop)
      
      if (is.matrix(mu)) {
        et <- mu[isv,, drop=FALSE]
        STATUS <- matrix(stat,  length(isv), ncol(et))
        LOGH0  <- matrix(logh0, length(isv), ncol(et))
        H0M    <- matrix(H0,    length(isv), ncol(et))
        out[isv,] <- STATUS*(LOGH0 + et) - H0M*exp(et)
      } else {
        et <- mu[isv]
        out[isv] <- stat*(logh0 + et) - H0*exp(et)
      }
    }
    
    attr(out, "mu_y") <- mu
    out
  }
  
  structure(list(
    family   = "cprobit + WeibullPH (indexed)",
    link     = stats$name,
    linkfun  = stats$linkfun,
    linkinv  = stats$linkinv,
    log_dens = log_dens
  ), class="family")
}


.safe_saveRDS <- function(object, file, compress = "xz") {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(file, ".tmp_", sprintf("%06d", sample.int(1e6, 1)))
  on.exit(try(unlink(tmp), silent = TRUE), add = TRUE)
  saveRDS(object, tmp, compress = compress)
  file.rename(tmp, file)
}


.extract_re_summary_item_surv <- function(fit) {
  V <- tryCatch(VarCorr(fit)$ID, error = function(e) NULL)
  if (is.null(V)) V <- tryCatch(fit$D, error = function(e) NULL)
  if (is.null(V)) {
    return(list(
      sds  = c(z1I0 = NA_real_, z1Ipost = NA_real_, z2 = NA_real_),
      corr = c(item_I0_vs_Ipost = NA_real_, I0_vs_frailty = NA_real_, Ipost_vs_frailty = NA_real_)
    ))
  }
  V <- as.matrix(V); nm <- rownames(V)
  getc <- function(a,b) if (!(a %in% nm) || !(b %in% nm)) NA_real_ else V[a,b]/sqrt(V[a,a]*V[b,b])
  sds <- setNames(rep(NA_real_,3), c("z1I0","z1Ipost","z2"))
  d   <- diag(V); names(d) <- nm
  present <- intersect(names(sds), names(d)); sds[present] <- sqrt(d[present])

  corr <- c(
    item_I0_vs_Ipost = getc("z1I0","z1Ipost"),
    I0_vs_frailty    = getc("z1I0","z2"),
    Ipost_vs_frailty = getc("z1Ipost","z2")
  )
  list(sds = sds, corr = corr)
}


                                
fit_one_item_surv <- function(item_name,
                              ctrl_joint = list(iter_EM = 0, nAGQ = 7),
                              dat_path   = NULL,
                              defs_file  = NULL) {

  
  if (!is.null(defs_file)) {
    source(defs_file, local = TRUE)
  }

  
  if (!is.null(dat_path)) {
    dataknee_ER_first_hosp <- readRDS(dat_path)
  } else if (!exists("dataknee30g_ER_first_hosp", inherits = TRUE)) {
    stop("Oggetto 'dataknee30g_ER_first_hosp' non trovato e 'dat_path' non fornito.")
  }

  
  dataknee_ER_first_hosp$ID        <- factor(dataknee_ER_first_hosp$ID)
  dataknee_ER_first_hosp$sex01     <- as.integer(dataknee_ER_first_hosp$sex == "M")
  dataknee_ER_first_hosp$prost_uni <- as.integer(dataknee_ER_first_hosp$prosthesis_type == "unicompartimental")
  dataknee_ER_first_hosp$diag_oa   <- as.integer(dataknee_ER_first_hosp$diagnosis != "osteoarthritis_knee")
  dataknee_ER_first_hosp <- dataknee_ER_first_hosp %>% filter(!if_all(everything(), is.na))

  age_ref <- dataknee_ER_first_hosp %>%
    dplyr::distinct(ID, age) %>%
    dplyr::filter(!is.na(age))

  age_center <- mean(age_ref$age)
  age_scale  <- sd(age_ref$age)

  dataknee_ER_first_hosp <- dataknee_ER_first_hosp %>%
    dplyr::mutate(
      age = (as.numeric(age) - age_center) / age_scale
    )

 base_long <- dataknee_ER_first_hosp %>%
  group_by(ID) %>%
  arrange(times, .by_group = TRUE) %>%
  mutate(
    t_within   = times - first(times),
    visit_idx  = row_number() - 1L,
    post       = as.integer(visit_idx >= 1L),
    time_post  = ifelse(post == 1L, t_within, 0),
    loghosp    = log(1 + hospitalizationdays),
    waiting    = as.numeric(difftime(date_intervention, date_baseline, units = "days"))/30.4375,
    time_discrete = factor(visit_idx, levels = c(0,1,2), labels = c("t0","t6","t12")),
    I0    = as.integer(time_discrete == "t0"),
    I6    = as.integer(time_discrete == "t6"),
    I12   = as.integer(time_discrete == "t12"),
    Ipost = I6 + I12,
    loghosp_post = loghosp * Ipost
  ) %>%
  ungroup()

  


  # prosthesis failure
  admin_date <- as.POSIXct("2025-01-01", tz = attr(dataknee_ER_first_hosp$date_intervention, "tzone"))
  
  surv_core <- dataknee_ER_first_hosp %>%
    mutate(
      t_event_days = as.numeric(difftime(date_revision, date_intervention, units = "days"))/ 30.4375,
       t_admin_days = as.numeric(difftime(admin_date, date_intervention, units = "days"))/ 30.4375
    ) %>%
    group_by(ID) %>%
    summarise(
      t_event = {
        v <- t_event_days[!is.na(t_event_days) & t_event_days >= 0]
        if (length(v)) min(v) else NA_real_
      },
       t_admin = min(t_admin_days, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      status = as.integer(!is.na(t_event) & (t_event <= t_admin)),
      stop   = ifelse(status == 1L, t_event, t_admin),
      start  = 0
    ) %>%
    transmute(
      ID     = factor(ID),
      start  = start,
      stop   = pmax(as.numeric(stop), 1e-8),
      status = status
    )

dataknee_ER_first_hosp2 <- base_long %>%
  left_join(
    surv_core %>% mutate(t_event = ifelse(status == 1L, stop, NA_real_)),
    by = "ID"
  ) %>%
  relocate(start, stop, status, t_event, .after = ID) %>%
  tidyr::drop_na(start, stop, status, sex01, prost_uni, loghosp_post, age, diag_oa)


  
  items_ord <- c(item_name)

  dat0 <- dataknee_ER_first_hosp2 %>%
    mutate(
      ID        = as.factor(ID),
      visit_idx = coalesce(visit_idx, 0L),
      I0        = as.integer(I0),
      I6        = as.integer(I6),
      I12       = as.integer(I12),
      Ipost     = as.integer(Ipost),
      prost_uni = as.integer(prosthesis_type == "unicompartimental")
    )


ord <- dat0 %>%
  dplyr::select(
    ID, visit_idx, I0, I6, I12, Ipost,
    waiting, loghosp_post, age, sex01, diag_oa, prost_uni,
    start, stop, status,
    dplyr::all_of(items_ord)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(items_ord),
    names_to = "item",
    values_to = "y_raw"
  ) %>%
  dplyr::filter(!is.na(y_raw)) %>%
  dplyr::mutate(
    visit_m = dplyr::case_when(
      visit_idx == 0L ~ 0,
      visit_idx == 1L ~ 6,
      visit_idx == 2L ~ 12,
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::filter(is.na(stop) | visit_m <= stop + 1e-8) %>%
  dplyr::mutate(
    type = 0L,
    time = dplyr::case_when(
      visit_idx == 0 ~ 0L,
      visit_idx == 1 ~ 1L,
      visit_idx == 2 ~ 2L,
      TRUE ~ NA_integer_
    ),
    y_raw = as.integer(y_raw),
    y     = as.integer(y_raw - min(y_raw, na.rm = TRUE) + 1L),
    item  = factor(item, levels = items_ord),
    wI6   = waiting * I6,
    wI12  = waiting * I12,
    pI6   = prost_uni * I6,
    pI12  = prost_uni * I12,
    z1I0    = ifelse(I0    == 1L, 1, 0),
    z1Ipost = ifelse(Ipost == 1L, 1, 0)
  )

  surv <- dat0 %>%
  group_by(ID) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    ID,
    visit_idx = NA_integer_,
    I0, I6, I12, Ipost,
    waiting, wI6 = waiting * I6, wI12 = waiting * I12,
    pI6 = prost_uni * I6, pI12 = prost_uni * I12,
    loghosp_post, sex01, diag_oa,
    age, prost_uni,
    item = factor(items_ord[1], levels = items_ord),
    type = 1L, time = NA_integer_,
    y = 0L,
    start, stop, status
  )

  long <- ord %>%
    mutate(start = NA_real_, stop = NA_real_, status = NA_integer_) %>%
    bind_rows(surv) %>%
    arrange(ID, type, time, item)

  K <- max(long$y[long$type == 0L], na.rm = TRUE)

long <- long %>%
  mutate(
    I6    = I6,
    I12   = I12,
    wI6   = wI6,
    wI12  = wI12,
    pI6   = pI6,
    pI12  = pI12,
    loghosp_post = loghosp_post,
    sex01   = sex01,
    diag_oa = diag_oa,
    is_ord  = ifelse(type == 0L, 1, 0),
    is_surv = ifelse(type == 1L, 1, 0)
  ) %>%
  mutate(
    z1I0    = ifelse(type == 0L & I0    == 1L, 1, 0),
    z1Ipost = ifelse(type == 0L & Ipost == 1L, 1, 0),
    z2      = ifelse(type == 1L, 1, 0),
    age     = age
  )

  long$row_id <- seq_len(nrow(long))
  long$Y <- cbind(y = ifelse(long$type == 0L, long$y, 0), idx = long$row_id)

  
  fam <- cprobit_weibull_indexed(
    K          = K,
    type_vec   = long$type,
    start_vec  = long$start,
    stop_vec   = long$stop,
    status_vec = long$status
  )

  zeta_init <- tryCatch({
    fit0 <- suppressWarnings(polr(ordered(long$y[long$type==0L]) ~ 1, method = "probit", Hess = FALSE))
    as.numeric(fit0$zeta)
  }, error = function(e) qnorm(seq(1/K, (K-1)/K, length.out = K - 1L)))

  surv_dat <- long[long$type == 1L, ]
  fit_wb0  <- survreg(Surv(stop, status) ~ age, data = surv_dat, dist = "weibull")
  log_rho_init    <- log(1 / fit_wb0$scale)
  log_lambda_init <- log( exp(-fit_wb0$coefficients[1] * (1 / fit_wb0$scale)) )
  phis_init <- c(zeta_init, log_lambda_init, log_rho_init)
  stopifnot(length(phis_init) == (K - 1L + 2L))

  fixed_fml <- Y ~ 0 +
  is_ord:(I6 + I12 +
            wI6 + wI12 +
            pI6 + pI12 + diag_oa +
            sex01 + loghosp_post +
            age) +
  is_surv:(age)

  random_fml <- ~ 0 + z1I0 + z1Ipost + z2 | ID

                        
  univariateGLMM <- mixed_model(
  fixed  = y ~ 0 + I6 + I12 +
    wI6 + wI12 + pI6 + pI12 + loghosp_post + sex01 + age + diag_oa,
  random = ~ 0 + z1I0 + z1Ipost | ID,
  data   = ord,
  family = cumulativeprobit(K),
  n_phis = K - 1,
  initial_values = list(phis = zeta_init),
  control = list(iter_EM = 0, nAGQ = 11)
)

  
  wb0 <- readRDS("wb0fail.rds") # a PH weibull model with age as fixed effect
  bet0 <- stats::coef(wb0)[3:length(coef(wb0))]
  phi0 <- c(wb0$res.t["scale","est"], wb0$res.t["shape","est"])

  weibullPHfrailty <- GLMMadaptive::mixed_model(
    fixed      = cbind(start, stop, status) ~ 0 + age,
    random     = ~ 1 | ID,
    data       = surv_dat,
    family     = weibullPH_lognormal(),
    n_phis     = 2L,
    control    = list(iter_EM = 0, nAGQ = 20),
    initial_values = list(betas = bet0, phis = phi0, D = matrix(0.4^2, 1, 1))
  )

 phis_initsurv <- c(univariateGLMM$phis, weibullPHfrailty$phis)
    names(phis_initsurv) <- paste0("phi_", c(1:(K - 1 + 2)))

    X_joint <- model.matrix(fixed_fml, data = long)
    bn <- colnames(X_joint)
    betas_init <- setNames(rep(0, length(bn)), bn)

    b_ord <- univariateGLMM$coefficients
    names(b_ord) <- paste0("is_ord:", names(b_ord))
    betas_init[names(b_ord)] <- b_ord

    b_surv <- weibullPHfrailty$coefficients
    names(b_surv) <- paste0("is_surv:", names(b_surv))
    betas_init[names(b_surv)] <- b_surv

    D0s <- rbind(
    cbind(univariateGLMM$D, -0.1),
    matrix(c(rep(-0.1, 2), 1.30), nrow = 1)
    )
    colnames(D0s)[3] <- "z2"
    rownames(D0s)[3] <- "z2"

  fit_joint <- mixed_model(
  fixed  = fixed_fml,
  random = random_fml,
  data   = long,
  family = fam,
  n_phis = (K - 1) + 2,
  initial_values = list(
    betas = unname(betas_init[colnames(X_joint)]),
    phis  = phis_initsurv,
    D     = D0s
  ),
  control = ctrl_joint
)

                        
  re  <- .extract_re_summary_item_surv(fit_joint)
  s   <- re$sds; crr <- re$corr

  out_row <- tibble::tibble(
    item  = item_name,
    logLik = as.numeric(logLik(fit_joint)),
    AIC    = AIC(fit_joint),
    sd_z1I0    = unname(s["z1I0"]),
    sd_z1Ipost = unname(s["z1Ipost"]),
    sd_z2      = unname(s["z2"]),
    corr_item_I0_vs_Ipost = unname(crr["item_I0_vs_Ipost"]),
    corr_I0_vs_frailty    = unname(crr["I0_vs_frailty"]),
    corr_Ipost_vs_frailty = unname(crr["Ipost_vs_frailty"]),
    conv = isTRUE(fit_joint$converged),
    msg  = NA_character_
  )

  list(ok = isTRUE(fit_joint$converged), pair = out_row, fit = fit_joint)
}


                        
fit_all_items_vs_surv_parallel <- function(itemnames,
                                           out_dir = "pairwise_item_surv",
                                           resume  = TRUE,
                                           n_cores = max(1L, parallel::detectCores(logical = TRUE) - 1L),
                                           dat_path = NULL,
                                           ctrl_joint = list(iter_EM = 0, nAGQ = 7),
                                           progress_log = file.path(out_dir, "_progress.log"),
                                           print_names = TRUE,
                                           defs_file = NULL) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(dat_path)) dat_path <- normalizePath(dat_path, mustWork = TRUE)

  N <- length(itemnames)
  .log_line <- function(...) {
    line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", sprintf(...), "\n")
    cat(line, file = progress_log, append = TRUE)
    invisible(NULL)
  }

  stubs  <- vapply(seq_len(N), function(i) sprintf("%03d_%s__SURV.rds",
                                                   i, gsub("[^A-Za-z0-9]+", "_", itemnames[i])), character(1))
  fpaths <- file.path(out_dir, stubs)
  todo_idx <- which(!(resume & file.exists(fpaths)))
  skip_idx <- setdiff(seq_len(N), todo_idx)

  if (print_names && length(skip_idx)) {
    for (k in skip_idx) {
      cat(sprintf("[SKIP] %03d/%03d :: %s vs SURV (già presente)\n", k, N, itemnames[k]))
    }
    flush.console()
  }

  if (length(todo_idx) == 0L) {
    res_loaded <- lapply(fpaths, readRDS)
    summary_tbl <- dplyr::bind_rows(lapply(res_loaded, `[[`, "pair"))
    fits_list   <- setNames(lapply(res_loaded, `[[`, "fit"),
                            vapply(itemnames, function(x) paste(x, "vs SURV"), character(1)))
    return(list(summary = summary_tbl, fits = fits_list, out_dir = out_dir, progress_log = progress_log))
  }

  has_future   <- requireNamespace("future.apply", quietly = TRUE)
  has_progress <- requireNamespace("progressr",    quietly = TRUE)

  if (has_future && has_progress) {
    future::plan(future::multisession, workers = min(n_cores, length(todo_idx)))
    on.exit({ try(future::plan(future::sequential), silent = TRUE) }, add = TRUE)
    if (requireNamespace("progress", quietly = TRUE)) progressr::handlers("progress") else progressr::handlers("txtprogressbar")
    progressr::handlers(global = TRUE)

    worker_fun <- function(i) {
      suppressWarnings(suppressPackageStartupMessages({
        library(dplyr); library(tidyr); library(tibble); library(stringr)
        library(GLMMadaptive); library(MASS); library(survival)
      }))
      item   <- itemnames[i]
      fpath  <- fpaths[i]
      errfile <- sub("\\.rds$", ".err.txt", fpath)

      if (resume && file.exists(fpath)) {
        out <- tryCatch(readRDS(fpath), error = function(e) NULL)
        if (!is.null(out)) return(out)
      }

      .log_line("START %03d :: %s vs SURV", i, item)

                        
      if (!is.null(defs_file)) source(defs_file, local = TRUE)

      out <- tryCatch({
        fit_one_item_surv(item_name = item,
                          ctrl_joint = ctrl_joint,
                          dat_path   = dat_path,
                          defs_file  = defs_file)
      }, error = function(e) {
        tib <- tibble::tibble(
          item = item,
          logLik = NA_real_, AIC = NA_real_,
          sd_z1I0 = NA_real_, sd_z1Ipost = NA_real_, sd_z2 = NA_real_,
          corr_item_I0_vs_Ipost = NA_real_, corr_I0_vs_frailty = NA_real_, corr_Ipost_vs_frailty = NA_real_,
          conv = FALSE, msg = conditionMessage(e)
        )
        list(ok = FALSE, pair = tib, fit = NULL)
      })

      .safe_saveRDS(out, fpath, compress = "xz")
      if (!isTRUE(out$ok)) cat(paste0(out$pair$msg, "\n"), file = errfile, append = FALSE)
      .log_line("DONE  %03d :: %s vs SURV (conv=%s)", i, item,
                ifelse(isTRUE(out$ok) && isTRUE(out$pair$conv), "TRUE", "FALSE"))
      out
    }

    res <- progressr::with_progress({
      p <- progressr::progressor(steps = length(todo_idx))
      future.apply::future_lapply(todo_idx, function(i) {
        lbl <- sprintf("%03d/%03d :: %s vs SURV", i, N, itemnames[i])
        p(lbl)
        if (print_names) { cat("[START] ", lbl, "\n", sep = ""); flush.console() }
        out <- worker_fun(i)
        aic <- suppressWarnings(tryCatch(as.numeric(out$pair$AIC), error = function(e) NA_real_))
        if (print_names) {
          cat(sprintf("[DONE ] %03d/%03d :: %s vs SURV  conv=%s  AIC=%s\n",
                      i, N, itemnames[i],
                      ifelse(isTRUE(out$ok) && isTRUE(out$pair$conv), "TRUE", "FALSE"),
                      ifelse(is.na(aic), "NA", sprintf("%.1f", aic))))
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
      stop("Instal 'future.apply' + 'progressr' or at least 'pbapply' for the progress bar.")

    cl <- parallel::makeCluster(min(n_cores, length(todo_idx)))
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

    parallel::clusterEvalQ(cl, {
      library(dplyr); library(tidyr); library(tibble); library(stringr)
      library(GLMMadaptive); library(MASS); library(survival)
      NULL
    })
    parallel::clusterExport(
      cl,
      varlist = c("itemnames","fpaths","resume","dat_path",".log_line",
                  "fit_one_item_surv","ctrl_joint",".safe_saveRDS",
                  ".extract_re_summary_item_surv","defs_file"),
      envir = environment()
    )

    if (print_names) {
      for (i in todo_idx) cat(sprintf("[QUEUE] %03d/%03d :: %s vs SURV\n", i, N, itemnames[i]))
      flush.console()
    }

    res_todo <- pbapply::pblapply(todo_idx, cl = cl, FUN = function(i) {
      item   <- itemnames[i]
      fpath  <- fpaths[i]
      errfile <- sub("\\.rds$", ".err.txt", fpath)

      if (resume && file.exists(fpath)) {
        out <- tryCatch(readRDS(fpath), error = function(e) NULL)
        if (!is.null(out)) return(out)
      }

      .log_line("START %03d :: %s vs SURV", i, item)

      ## FIX 4: source defs anche in PSOCK
      if (!is.null(defs_file)) source(defs_file, local = TRUE)

      out <- tryCatch({
        fit_one_item_surv(item_name = item,
                          ctrl_joint = ctrl_joint,
                          dat_path   = dat_path,
                          defs_file  = defs_file)
      }, error = function(e) {
        tib <- tibble::tibble(
          item = item,
          logLik = NA_real_, AIC = NA_real_,
          sd_z1I0 = NA_real_, sd_z1Ipost = NA_real_, sd_z2 = NA_real_,
          corr_item_I0_vs_Ipost = NA_real_, corr_I0_vs_frailty = NA_real_, corr_Ipost_vs_frailty = NA_real_,
          conv = FALSE, msg = conditionMessage(e)
        )
        list(ok = FALSE, pair = tib, fit = NULL)
      })

      .safe_saveRDS(out, fpath, compress = "xz")
      if (!isTRUE(out$ok)) cat(paste0(out$pair$msg, "\n"), file = errfile, append = FALSE)
      .log_line("DONE  %03d :: %s vs SURV (conv=%s)", i, item,
                ifelse(isTRUE(out$ok) && isTRUE(out$pair$conv), "TRUE", "FALSE"))
      out
    })

    all_res <- vector("list", N)
    if (length(skip_idx)) all_res[skip_idx] <- lapply(fpaths[skip_idx], readRDS)
    all_res[todo_idx] <- res_todo
  }

  summary_tbl <- dplyr::bind_rows(lapply(all_res, `[[`, "pair"))
  fits_list   <- setNames(lapply(all_res, `[[`, "fit"),
                          vapply(itemnames, function(x) paste(x, "vs SURV"), character(1)))
  list(summary = summary_tbl, fits = fits_list, out_dir = out_dir, progress_log = progress_log)
}

                                 
itemnames_all <- c('Mobility','SelfCare','UsualActivities','PainorDiscomfort','Anxietyordepression',
                   'Gettingoutofbed','Puttingonsocks','Gettingupfromsitting',
                   'Bendingtowardthefloor_pickinguponobjectfromtheground',
                   'Twisting_pivotingontheinjuredknee','Kneeling','Squatting')

defs_file <- normalizePath("univariateFuns.R", mustWork = TRUE)
res_item_surv <- fit_all_items_vs_surv_parallel(
  itemnames = itemnames_all,
  out_dir   = "pairwise_failure_hospost",
  resume    = TRUE,
  n_cores   = 12,
  dat_path  = "dat2_with_hosp30g.rds", 
  ctrl_joint = list(iter_EM = 0, nAGQ = 25),
  defs_file = defs_file 
)

                                 


# items_test <- c("SelfCare")
# defs_file <- normalizePath("jm_defs.R", mustWork = TRUE)
# res = fit_one_item_surv(item_name = items_test,
#                        ctrl_joint = list(iter_EM = 0, nAGQ = 24, verbose = TRUE),
#                          dat_path   = "dat2_with_hosp30g.rds",
#                          defs_file  = defs_file)


                                 
