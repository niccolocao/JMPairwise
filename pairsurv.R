## ================================================================
##  PAIRWISE (ITEM + SURVIVAL)
## ================================================================



###################################################################
#                 ATTENZIONE PER 30 GIONRI DA INTERVENTO QUI      #
###################################################################



suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(stringr)
  library(GLMMadaptive); library(MASS); library(survival)
})

## -- safe save 
.safe_saveRDS <- function(object, file, compress = "xz") {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(file, ".tmp_", sprintf("%06d", sample.int(1e6, 1)))
  on.exit(try(unlink(tmp), silent = TRUE), add = TRUE)
  saveRDS(object, tmp, compress = compress)
  file.rename(tmp, file)
}

## -- estrazione RE summary per (z1I0, z1Ipost, z2)
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


## --- FAMILY INDICIZZATA (usa idx per rileggere type/start/stop/status) ---
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
    
    ## ------ ORDINALE ------
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
    
    ## ------ SURVIVAL (Weibull–PH) ------
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

## ----------------------------------------------------------------
## fit di UNA coppia (item + survival) 
## ----------------------------------------------------------------
fit_one_item_surv <- function(item_name,
                              ctrl_joint = list(iter_EM = 0, nAGQ = 7),
                              dat_path   = NULL,
                              defs_file  = NULL) {

  ## FIX 1: sorgenti funzioni custom SE presenti
  if (!is.null(defs_file)) {
    source(defs_file, local = TRUE)
  }

  ## FIX 2: leggi il dataset se passato, altrimenti usa quello già in RAM
  if (!is.null(dat_path)) {
    dataknee_ER_first_hosp <- readRDS(dat_path)
  } else if (!exists("dataknee30g_ER_first_hosp", inherits = TRUE)) {
    stop("Oggetto 'dataknee30g_ER_first_hosp' non trovato e 'dat_path' non fornito.")
  }

  ## ---------- COPIA 1: tue trasformazioni dati (invariato) ----------
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

  # surv_core <- dataknee_ER_first_hosp %>%
  #   mutate(t_event = as.numeric(difftime(date_hosp30g, date_intervention, units = "days"))) %>%
  #   group_by(ID) %>%
  #   summarise(
  #     t_first = { v <- t_event[!is.na(t_event) & t_event >= 0]; if (length(v)) min(v) else NA_real_ },
  #     .groups = "drop"
  #   ) %>%
  #   mutate(
  #   status = as.integer(!is.na(t_first) & t_first <= 30),
  #   stop   = ifelse(status == 1L, pmin(t_first, 30), 30),
  #   start  = 0
  # ) %>%
  #   transmute(
  #     ID    = factor(ID),
  #     start = start,
  #     stop  = pmax(as.numeric(stop), 1e-8),
  #     status = as.integer(status)
  #   )



  # 2) per FALLIMENTO PROTESI
  admin_date <- as.POSIXct("2025-01-01", tz = attr(dataknee_ER_first_hosp$date_intervention, "tzone"))
  
  surv_core <- dataknee_ER_first_hosp %>%
    ## tempo evento in giorni (come prima)
    mutate(
      t_event_days = as.numeric(difftime(date_revision, date_intervention, units = "days"))/ 30.4375,
      ## tempo fino alla censura amministrativa (sempre >=0 se l'intervento è prima del 2025-01-01)
      t_admin_days = as.numeric(difftime(admin_date, date_intervention, units = "days"))/ 30.4375
    ) %>%
    group_by(ID) %>%
    summarise(
      ## per ogni ID prendo il PRIMO intervento che hai registrato (come facevi prima)
      t_event = {
        v <- t_event_days[!is.na(t_event_days) & t_event_days >= 0]
        if (length(v)) min(v) else NA_real_
      },
      ## anche la censura la voglio per ID
      t_admin = min(t_admin_days, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      ## evento osservato se c'è un evento e avviene PRIMA della data amministrativa
      status = as.integer(!is.na(t_event) & (t_event <= t_admin)),
      ## stop = min(evento, censura amm.)
      stop   = ifelse(status == 1L, t_event, t_admin),
      start  = 0
    ) %>%
    transmute(
      ID     = factor(ID),
      start  = start,
      stop   = pmax(as.numeric(stop), 1e-8),
      status = status
    )

  ## 2) per EARLY READMISSION entro 30 giorni (date_hosp30g)
  # surv_core <- dataknee_ER_first_hosp %>%
  #   mutate(
  #     t_event = as.numeric(difftime(date_hosp30g, date_intervention, units = "days"))/ 30.4375
  #   ) %>%
  #   group_by(ID) %>%
  #   summarise(
  #     t_first = {
  #       v <- t_event[!is.na(t_event) & is.finite(t_event) & t_event >= 0]
  #       if (length(v)) min(v) else NA_real_
  #     },
  #     .groups = "drop"
  #   ) %>%
  #   mutate(
  #     status = as.integer(!is.na(t_first) & t_first <= 30),
  #     stop   = ifelse(status == 1L, pmin(t_first, 30), 30),
  #     start  = 0
  #   ) %>%
  #   transmute(
  #     ID     = factor(ID),
  #     start  = start,
  #     stop   = pmax(as.numeric(stop), 1e-8),
  #     status = as.integer(status)
  #   )


# dataknee_ER_first_hosp2 <- base_long %>%
#   left_join(surv_core, by = "ID") %>%
#   relocate(start, stop, status, .after = ID) %>%
#   tidyr::drop_na(start, stop, status, sex01, prost_uni, loghosp, age, diag_oa)

# PER Y23
# dataknee_ER_first_hosp2 <- base_long %>%
#   left_join(surv_core %>% mutate(t_event = ifelse(is.na(stop), NA_real_, ifelse(status == 1L, stop, NA_real_))),
#             by = "ID") %>%
#   # Qui sopra: t_event = stop solo se status==1, altrimenti NA (censurati)
#   relocate(start, stop, status, t_event, .after = ID) %>%
#   tidyr::drop_na(start, stop, status, sex01, prost_uni, loghosp, age, diag_oa)

dataknee_ER_first_hosp2 <- base_long %>%
  left_join(
    surv_core %>% mutate(t_event = ifelse(status == 1L, stop, NA_real_)),
    by = "ID"
  ) %>%
  relocate(start, stop, status, t_event, .after = ID) %>%
  tidyr::drop_na(start, stop, status, sex01, prost_uni, loghosp_post, age, diag_oa)


  ## ---------- COPIA 2: preparazione long item+surv (invariato) ----------
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

  # ord <- dat0 %>%
  #   dplyr::select(
  #     ID, visit_idx, I0, I6, I12, Ipost,
  #     waiting, loghosp, age, sex01, diag_oa, prost_uni,
  #     dplyr::all_of(items_ord)
  #   ) %>%
  #   tidyr::pivot_longer(
  #     cols = dplyr::all_of(items_ord),
  #     names_to = "item",
  #     values_to = "y_raw"
  #   ) %>%
  #   dplyr::filter(!is.na(y_raw)) %>%
  #   dplyr::mutate(
  #     type = 0L,
  #     time = dplyr::case_when(
  #       visit_idx == 0 ~ 0L,
  #       visit_idx == 1 ~ 1L,
  #       visit_idx == 2 ~ 2L,
  #       TRUE ~ NA_integer_
  #     ),
  #     y_raw = as.integer(y_raw),
  #     y     = as.integer(y_raw - min(y_raw, na.rm = TRUE) + 1L),
  #     item  = factor(item, levels = items_ord),
  #     wI6   = waiting * I6,
  #     wI12  = waiting * I12,
  #     pI6   = prost_uni * I6,
  #     pI12  = prost_uni * I12,
  #     z1I0    = ifelse(I0    == 1L, 1, 0),
  #     z1Ipost = ifelse(Ipost == 1L, 1, 0)
  #   )

# PER Y23
# ord <- dat0 %>%
#   dplyr::select(
#     ID, visit_idx, I0, I6, I12, Ipost,
#     waiting, loghosp, age, sex01, diag_oa, prost_uni,
#     start, stop, status, t_event,
#     dplyr::all_of(items_ord)
#   ) %>%
#   tidyr::pivot_longer(
#     cols = dplyr::all_of(items_ord),
#     names_to = "item",
#     values_to = "y_raw"
#   ) %>%
#   dplyr::filter(!is.na(y_raw)) %>%
#   dplyr::mutate(
#     visit_days = dplyr::case_when(
#       visit_idx == 0L ~ 0,
#       visit_idx == 1L ~ 6  * 30.4375,
#       visit_idx == 2L ~ 12 * 30.4375,
#       TRUE ~ NA_real_
#     ),
#     t_event_days = ifelse(is.na(t_event) | t_event < 0, Inf, t_event),

#     type = 0L,
#     time = dplyr::case_when(
#       visit_idx == 0 ~ 0L,
#       visit_idx == 1 ~ 1L,
#       visit_idx == 2 ~ 2L,
#       TRUE ~ NA_integer_
#     ),
#     y_raw = as.integer(y_raw),
#     y     = as.integer(y_raw - min(y_raw, na.rm = TRUE) + 1L),
#     item  = factor(item, levels = items_ord),
#     wI6   = waiting * I6,
#     wI12  = waiting * I12,
#     pI6   = prost_uni * I6,
#     pI12  = prost_uni * I12,
#     z1I0    = ifelse(I0    == 1L, 1, 0),
#     z1Ipost = ifelse(Ipost == 1L, 1, 0)
#   )

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
#PER Y23
# bad <- ord %>% filter(!is.na(stop), visit_days > stop + 1e-8)
# stopifnot(nrow(bad) == 0)

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

  ## ---------- COPIA 3: family + start values (invariato) ----------
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

  ## ---------- COPIA 4: formule (invariato) ----------
  # fixed_fml <- Y ~ 0 +
  #   is_ord:(I6 + I12 +
  #             wI6 + wI12 +
  #             pI6 + pI12 + diag_oa +
  #             sex01 + loghosp +
  #             age) +
  #   is_surv:(age) #+ prost_uni)
  fixed_fml <- Y ~ 0 +
  is_ord:(I6 + I12 +
            wI6 + wI12 +
            pI6 + pI12 + diag_oa +
            sex01 + loghosp_post +
            age) +
  is_surv:(age)

  random_fml <- ~ 0 + z1I0 + z1Ipost + z2 | ID

  ## ---------- COPIA 5: start separati (invariato) ----------
  # univariateGLMM <- mixed_model(
  #   fixed  = y ~ 0 + I6 + I12 + wI6 + wI12 + pI6 + pI12 + loghosp + sex01 + age + diag_oa,
  #   random = ~ 0 + z1I0 + z1Ipost | ID,
  #   data   = ord,
  #   family = cumulativeprobit(K),
  #   n_phi  = K - 1,
  #   initial_values = list(phis = zeta_init),
  #   control = list(iter_EM = 10, nAGQ = 11)
  # )

# PER Y23
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


  wb0 <- readRDS("wb0fail.rds")
  bet0 <- stats::coef(wb0)[3:length(coef(wb0))]
  phi0 <- c(wb0$res.t["scale","est"], wb0$res.t["shape","est"])

  weibullPHfrailty <- GLMMadaptive::mixed_model(
    fixed      = cbind(start, stop, status) ~ 0 + age,# + prost_uni,
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

  ## ---------- COPIA 6: fit joint (invariato) ----------
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

  ## ---------- riassunto ----------
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

  list(ok = TRUE, pair = out_row, fit = fit_joint)
}

## ----------------------------------------------------------------
## parallel + resume + logging per TUTTI gli item
## ----------------------------------------------------------------
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

      ## FIX 3: source delle defs anche qui
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
      stop("Installa 'future.apply' + 'progressr' oppure almeno 'pbapply' per la progress bar.")

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

## ---------------------------------------------------------------
## Esempio d’uso
## ---------------------------------------------------------------
itemnames_all <- c('Mobility','SelfCare','UsualActivities','PainorDiscomfort','Anxietyordepression',
                   'Gettingoutofbed','Puttingonsocks','Gettingupfromsitting',
                   'Bendingtowardthefloor_pickinguponobjectfromtheground',
                   'Twisting_pivotingontheinjuredknee','Kneeling','Squatting')

## metti in un file (es. jm_defs.R) SOLO:
##   cprobit_weibull_indexed, cumulativeprobit, .coerce_y3, weibullPH_lognormal
defs_file <- normalizePath("jm_defs.R", mustWork = TRUE)
res_item_surv <- fit_all_items_vs_surv_parallel(
  itemnames = itemnames_all,
  out_dir   = "pairwise_failure_hospost",
  resume    = TRUE,
  n_cores   = 12,
  dat_path  = "dat2_with_hosp30g.rds", #"dataknee30g_ER_first_hosp.rds",
  ctrl_joint = list(iter_EM = 0, nAGQ = 25),
  defs_file = defs_file  # <-- assoluto, visto da tutti i worker
)

# res_item_surv$summary





# for(i in 1:12){
#   ii = i
#   if(ii<10) ii=paste0("0",ii)
#   print(itemnames_all[i])
#   print(summary(readRDS(paste0("pairwise_failure_months/0",ii,"_",itemnames_all[i],"__SURV.rds"))$fit))
# }


# # piano nazionale esiti pne.genas.it





items_test <- c("SelfCare")
defs_file <- normalizePath("jm_defs.R", mustWork = TRUE)
res = fit_one_item_surv(item_name = items_test,
                          ctrl_joint = list(iter_EM = 0, nAGQ = 24, verbose = TRUE),
                          dat_path   = "dat2_with_hosp30g.rds",
                          defs_file  = defs_file)
# saveRDS(res, "/home/niccolo.cao/ord+surv/pairwise_failure_months/005_Anxietyordepression__SURV.rds")


# dataknee_ER_first_hosp=readRDS("dat2.rds")


# library(dplyr)

# # estrai stop e status a livello ID
# surv <- long %>%
#   filter(type == 1L) %>%
#   distinct(ID, stop, status)

# # conta visite longitudinali post (visit_idx 1 e 2) per ciascun ID
# post_vis <- long %>%
#   filter(type == 0L, visit_idx %in% c(1L,2L)) %>%
#   distinct(ID, visit_idx) %>%
#   count(ID, name = "n_post_vis")

# chk <- surv %>%
#   left_join(post_vis, by="ID") %>%
#   mutate(n_post_vis = coalesce(n_post_vis, 0L),
#          early_fail = (status == 1L & stop < 182.6))  # ~6 mesi in giorni

# table(chk$early_fail, chk$n_post_vis > 0)


########## Per la stabilità di stime con solo 24 eventi

# library(MASS)

# ## =========================
# ## PARAMETRI DAL TUO FIT
# ## =========================

# ## matrice D esattamente come nel tuo fit$D
# D_true <- matrix(c(
#   1.9401002,  1.13423697, -0.33304765,
#   1.1342370,  2.98343250, -0.09328921,
#  -0.3330477, -0.09328921,  0.47330520
# ), nrow = 3, byrow = TRUE)
# colnames(D_true) <- rownames(D_true) <- c("z1I0","z1Ipost","z2")

# ## phis dal fit
# phis_true <- c(
#   -4.887955,  -3.256147,  -1.649286, -0.1961103,
#   -9.440986,  -0.2163383
# )
# K             <- 5
# phis_ord_true <- phis_true[1:4]

# log_lambda_modello <- phis_true[5]   # log(lambda_model)
# log_rho_true       <- phis_true[6]

# rho_true <- exp(log_rho_true)

# ## varianza frailty (z2)
# sigma2_frailty <- D_true["z2","z2"]

# ## lambda_model = exp(phi5), poi RISCALE come hai chiesto: * exp(Var/2)
# lambda_modello    <- exp(log_lambda_modello)
# lambda_vero_true  <- lambda_modello * exp(sigma2_frailty / 2)

# ## beta ordinali SENZA prefissi is_ord:
# beta_ord_true <- c(
#   I6      = 1.87769303,
#   I12     = 1.95183369,
#   wI6     = 0.00139789,
#   wI12    = 0.07087392,
#   pI6     = -0.38095491,
#   pI12    = -0.01654053,
#   diag_oa = 0.20440917,
#   sex01   = 0.87100837,
#   loghosp = -1.30420068,
#   age     = 0.10716621
# )

# ## coefficiente di age nella survival (age:is_surv)
# beta_surv_true <- -0.32898321


# ## =========================
# ## FUNZIONE DI SIMULAZIONE
# ## =========================
# simulate_item_surv_once <- function(
#   n_id         = 679,
#   t_cens       = 3650,
#   lambda_scale = 1
# ) {
#   ## =========================
#   ## PARAMETRI FISSATI (TUO FIT)
#   ## =========================
  
#   # matrice D del fit
#   D_true <- matrix(c(
#     1.9401002,  1.13423697, -0.33304765,
#     1.1342370,  2.98343250, -0.09328921,
#    -0.3330477, -0.09328921,  0.47330520
#   ), nrow = 3, byrow = TRUE)
#   colnames(D_true) <- rownames(D_true) <- c("z1I0","z1Ipost","z2")
  
#   # phis del fit
#   phis_ord_true <- c(-4.887955,  -3.256147,  -1.649286, -0.1961103) # soglie ord
#   log_lambda_modello <- -9.440986
#   log_rho_true       <- -0.2163383
#   rho_true           <- exp(log_rho_true)
  
#   # varianza frailty
#   sigma2_frailty <- D_true["z2","z2"]
  
#   # lambda_model = exp(phi5) e poi RISCALE come vuoi tu: * exp(Var/2)
#   lambda_modello   <- exp(log_lambda_modello)
#   lambda_vero      <- lambda_modello * exp(sigma2_frailty / 2)
#   lambda           <- lambda_vero * lambda_scale
  
#   # beta ordinali (senza prefisso is_ord:)
#   beta_ord <- c(
#     I6      = 1.87769303,
#     I12     = 1.95183369,
#     wI6     = 0.00139789,
#     wI12    = 0.07087392,
#     pI6     = -0.38095491,
#     pI12    = -0.01654053,
#     diag_oa = 0.20440917,
#     sex01   = 0.87100837,
#     loghosp = -1.30420068,
#     age     = 0.10716621
#   )
  
#   # coefficiente di age nella survival (age:is_surv)
#   beta_surv <- -0.32898321
  
#   ## =========================
#   ## COVARIATE SOGGETTO
#   ## =========================
  
#   ID        <- factor(seq_len(n_id))
#   age_std   <- rnorm(n_id, 0, 1)   # analogo di scale(age)
#   loghosp   <- rnorm(n_id, 0, 1)
#   sex01     <- rbinom(n_id, 1, 0.5)
#   diag_oa   <- rbinom(n_id, 1, 0.2)
#   prost_uni <- rbinom(n_id, 1, 0.2)
#   waiting   <- runif(n_id, 0, 1)
  
#   ## =========================
#   ## RANDOM EFFECTS: b ~ N(0, D_true)
#   ## =========================
  
#   b_mat <- MASS::mvrnorm(n_id, mu = c(0, 0, 0), Sigma = D_true)
#   colnames(b_mat) <- c("b_z1I0", "b_z1Ipost", "b_z2")
  
#   ## =========================
#   ## PARTE ORDINALE (3 visite)
#   ## =========================
  
#   ord <- expand.grid(
#     ID        = ID,
#     visit_idx = 0:2
#   )
#   ord <- ord[order(ord$ID, ord$visit_idx), ]
  
#   ord$I0    <- as.integer(ord$visit_idx == 0L)
#   ord$I6    <- as.integer(ord$visit_idx == 1L)
#   ord$I12   <- as.integer(ord$visit_idx == 2L)
#   ord$Ipost <- ord$I6 + ord$I12
  
#   ord$age       <- rep(age_std,   each = 3)
#   ord$loghosp   <- rep(loghosp,   each = 3)
#   ord$sex01     <- rep(sex01,     each = 3)
#   ord$diag_oa   <- rep(diag_oa,   each = 3)
#   ord$prost_uni <- rep(prost_uni, each = 3)
#   ord$waiting   <- rep(waiting,   each = 3)
  
#   ord$wI6  <- ord$waiting * ord$I6
#   ord$wI12 <- ord$waiting * ord$I12
#   ord$pI6  <- ord$prost_uni * ord$I6
#   ord$pI12 <- ord$prost_uni * ord$I12
  
#   ord$z1I0    <- as.integer(ord$I0    == 1L)
#   ord$z1Ipost <- as.integer(ord$Ipost == 1L)
  
#   ord$b_z1I0    <- rep(b_mat[, "b_z1I0"],    each = 3)
#   ord$b_z1Ipost <- rep(b_mat[, "b_z1Ipost"], each = 3)
  
#   X_ord <- model.matrix(
#     ~ 0 + I6 + I12 + wI6 + wI12 + pI6 + pI12 + diag_oa + sex01 + loghosp + age,
#     data = ord
#   )
#   beta_ord_vec  <- beta_ord[colnames(X_ord)]
#   eta_ord_fixed <- as.vector(X_ord %*% beta_ord_vec)
#   eta_ord       <- eta_ord_fixed +
#     ord$z1I0    * ord$b_z1I0 +
#     ord$z1Ipost * ord$b_z1Ipost
  
#   cutpoints <- c(-Inf, phis_ord_true, Inf)
#   latent    <- rnorm(nrow(ord), mean = eta_ord, sd = 1)
#   y_cat     <- cut(latent, breaks = cutpoints, labels = FALSE)
  
#   ord$y_raw <- as.integer(y_cat)
#   ord$y     <- ord$y_raw
#   ord$type  <- 0L
#   ord$time  <- ord$visit_idx
#   ord$item  <- factor("Gettingupfromsitting")
  
#   ## =========================
#   ## PARTE SURVIVAL (WeibullPH + frailty z2)
#   ## =========================
  
#   surv <- data.frame(
#     ID        = ID,
#     start     = 0,
#     sex01     = sex01,
#     diag_oa   = diag_oa,
#     prost_uni = prost_uni,
#     loghosp   = loghosp,
#     age       = age_std,
#     waiting   = waiting,
#     I0        = 1L,
#     I6        = 0L,
#     I12       = 0L,
#     Ipost     = 0L
#   )
  
#   surv$wI6  <- 0
#   surv$wI12 <- 0
#   surv$pI6  <- 0
#   surv$pI12 <- 0
#   surv$z1I0    <- 0L
#   surv$z1Ipost <- 0L
#   surv$z2      <- 1L
  
#   # LP survival: age:is_surv + frailty z2
#   eta_surv <- beta_surv * surv$age + b_mat[, "b_z2"]
  
#   u       <- runif(n_id)
#   T_event <- (-log(u) / (lambda * exp(eta_surv)))^(1 / rho_true)
  
#   stop_time <- pmin(T_event, t_cens)
#   status    <- as.integer(T_event <= t_cens)
  
#   surv$stop      <- stop_time
#   surv$status    <- status
#   surv$visit_idx <- NA_integer_
#   surv$time      <- NA_integer_
#   surv$y         <- 0L
#   surv$type      <- 1L
#   surv$item      <- factor("Gettingupfromsitting")
  
#   ## =========================
#   ## COSTRUZIONE LONG
#   ## =========================
  
#   ord_long <- ord
#   ord_long$start  <- NA_real_
#   ord_long$stop   <- NA_real_
#   ord_long$status <- NA_integer_
#   ord_long$z2     <- 0L
  
#   cols <- c("ID","visit_idx","I0","I6","I12","Ipost",
#             "waiting","loghosp","age","sex01","diag_oa","prost_uni",
#             "wI6","wI12","pI6","pI12",
#             "z1I0","z1Ipost","z2",
#             "item","type","time","y","start","stop","status")
  
#   long <- rbind(
#     ord_long[, cols],
#     surv[,      cols]
#   )
  
#   long <- long[order(long$ID, long$type, long$time), ]
#   rownames(long) <- NULL
  
#   long$row_id <- seq_len(nrow(long))
#   long$Y <- cbind(
#     y   = ifelse(long$type == 0L, long$y, 0L),
#     idx = long$row_id
#   )
#   long$is_ord  <- ifelse(long$type == 0L, 1L, 0L)
#   long$is_surv <- ifelse(long$type == 1L, 1L, 0L)
  
#   list(
#     long      = long,
#     ord       = ord,
#     surv      = surv,
#     n_events  = sum(surv$status),
#     lambda    = lambda,
#     rho       = rho_true
#   )
# }



# set.seed(123)
# tuned <- tune_lambda_scale(target_events = 24)

# lambda_scale_star <- tuned$lambda_scale
# tuned$sim$n_events     # numero eventi nell'ultima simulazione



# long_sim <- tuned$sim$long

# ## numero categorie ordinale (K) dai dati simulati
# K_sim <- max(long_sim$y[long_sim$type == 0L], na.rm = TRUE)

# source("/home/niccolo.cao/ord+surv/jm_defs.R") 
# ## family indicizzata come nel tuo codice
# fam_sim <- cprobit_weibull_indexed(
#   K          = K_sim,
#   type_vec   = long_sim$type,
#   start_vec  = long_sim$start,
#   stop_vec   = long_sim$stop,
#   status_vec = long_sim$status
# )

# fixed_fml_sim <- Y ~ 0 +
#   is_ord:(I6 + I12 +
#             wI6 + wI12 +
#             pI6 + pI12 + diag_oa +
#             sex01 + loghosp +
#             age) +
#   is_surv:(age)

# random_fml_sim <- ~ 0 + z1I0 + z1Ipost + z2 | ID

# fit_sim <- GLMMadaptive::mixed_model(
#   fixed  = fixed_fml_sim,
#   random = random_fml_sim,
#   data   = long_sim,
#   family = fam_sim,
#   n_phis = (K_sim - 1L) + 2L,
#   initial_values = list(
#     phis  = phis_init,
#     D     = D_init
#   ),
#   control = list(
#     iter_EM = 0,
#     nAGQ    = 15,
#     optimizer    = "optim",
#     optim_method = "BFGS",
#     verbose      = TRUE
#   )
# )

# summary(fit_sim)
# logLik(fit_sim)
# fit_sim$D
# fit_sim$coefficients
# fit_sim$phis







# tune_lambda_scale <- function(target_events = 24,
#                               n_id         = 679,
#                               t_cens       = 3650,
#                               max_iter     = 8) {
#   lambda_scale <- 1
#   sim_last     <- NULL
  
#   for (it in seq_len(max_iter)) {
#     sim <- simulate_item_surv_once(
#       n_id         = n_id,
#       t_cens       = t_cens,
#       lambda_scale = lambda_scale
#     )
#     n_ev <- sim$n_events
#     sim_last <- sim
    
#     cat("iter", it,
#         "  lambda_scale =", signif(lambda_scale, 4),
#         "  events =", n_ev, "\n")
    
#     if (n_ev == 0) break
#     if (abs(n_ev - target_events) <= 1) break
    
#     lambda_scale <- lambda_scale * (target_events / n_ev)
#   }
  
#   list(
#     lambda_scale = lambda_scale,
#     sim          = sim_last
#   )
# }

















# # rm(list = ls())
# library(MASS)
# library(survival)
# library(GLMMadaptive)

# ## =========================================================
# ## 1. PARAMETRI "VERI" DAL FIT REALE
# ## =========================================================

# D_true <- matrix(c(
#   1.9401002,  1.13423697, -0.33304765,
#   1.1342370,  2.98343250, -0.09328921,
#  -0.3330477, -0.09328921,  0.47330520
# ), nrow = 3, byrow = TRUE)
# colnames(D_true) <- rownames(D_true) <- c("z1I0","z1Ipost","z2")

# phis_true <- c(
#   -4.887955,  -3.256147,  -1.649286, -0.1961103,
#   -9.440986,  -0.2163383
# )
# K_true           <- 5
# phis_ord_true    <- phis_true[1:4]
# log_lambda_mod   <- phis_true[5]
# log_rho_true     <- phis_true[6]
# rho_true         <- exp(log_rho_true)

# sigma2_frailty   <- D_true["z2","z2"]
# lambda_modello   <- exp(log_lambda_mod)
# lambda_vero_true <- lambda_modello * exp(sigma2_frailty / 2)

# beta_ord_true <- c(
#   I6      = 1.87769303,
#   I12     = 1.95183369,
#   wI6     = 0.00139789,
#   wI12    = 0.07087392,
#   pI6     = -0.38095491,
#   pI12    = -0.01654053,
#   diag_oa = 0.20440917,
#   sex01   = 0.87100837,
#   loghosp = -1.30420068,
#   age     = 0.10716621
# )
# beta_surv_true <- -0.32898321


# ## =========================================================
# ## 2. FUNZIONE DI SIMULAZIONE (una coppia ord + survival)
# ## =========================================================
# simulate_item_surv_once <- function(
#   n_id         = 679,
#   t_cens       = 3650,
#   lambda_scale = 1
# ) {
#   ## =========================================================
#   ## 1. PARAMETRI "VERI" INTERNI = QUELLI DEL FIT
#   ##    (copiati dalla tua versione, invariati)
#   ## =========================================================
  
#   # matrice D del fit (random effects: z1I0, z1Ipost, z2)
#   D_true <- matrix(c(
#     1.9401002,  1.13423697, -0.33304765,
#     1.1342370,  2.98343250, -0.09328921,
#    -0.3330477, -0.09328921,  0.47330520
#   ), nrow = 3, byrow = TRUE)
#   colnames(D_true) <- rownames(D_true) <- c("z1I0","z1Ipost","z2")
  
#   # soglie ordinali (phis 1:4) - INVARIATE
#   phis_ord_true <- c(-4.887955,  -3.256147,  -1.649286, -0.1961103)
  
#   # parametri Weibull dal fit - INVARIATI
#   log_lambda_mod <- -9.440986   # log(lambda_mod) usato nel modello
#   log_rho_true   <- -0.2163383
#   rho_true       <- exp(log_rho_true)
  
#   ## MODIFICA 1 (IMPORTANTE):
#   ## prima facevi:
#   ##   sigma2_frailty <- D_true["z2","z2"]
#   ##   lambda_modello <- exp(log_lambda_mod)
#   ##   lambda_vero    <- lambda_modello * exp(sigma2_frailty / 2)
#   ##   lambda         <- lambda_vero * lambda_scale
#   ##
#   ## ADESSO: uso direttamente λ_mod del modello GLMMadaptive,
#   ## senza riscalatura per exp(Var/2).
#   ## Questo rende la simulazione ESATTAMENTE coerente col modello stimato.
  
#   lambda <- exp(log_lambda_mod) * lambda_scale
  
#   # beta ordinali (senza prefisso is_ord:) - INVARIATI
#   beta_ord <- c(
#     I6      = 1.87769303,
#     I12     = 1.95183369,
#     wI6     = 0.00139789,
#     wI12    = 0.07087392,
#     pI6     = -0.38095491,
#     pI12    = -0.01654053,
#     diag_oa = 0.20440917,
#     sex01   = 0.87100837,
#     loghosp = -1.30420068,
#     age     = 0.10716621
#   )
  
#   # coefficiente di age nella survival (age:is_surv) - INVARIATO
#   beta_surv <- -0.32898321
  
#   ## =========================================================
#   ## 2. COVARIATE SOGGETTO (tue scelte, lasciate invariate)
#   ## =========================================================
  
#   ID        <- factor(seq_len(n_id))
#   age_std   <- rnorm(n_id, 0, 1)   # analogo a scale(age)
#   loghosp   <- rnorm(n_id, 0, 1)
#   sex01     <- rbinom(n_id, 1, 0.5)
#   diag_oa   <- rbinom(n_id, 1, 0.2)
#   prost_uni <- rbinom(n_id, 1, 0.2)
#   waiting   <- runif(n_id, 0, 1)
  
#   ## =========================================================
#   ## 3. RANDOM EFFECTS: b ~ N(0, D_true)
#   ## =========================================================
  
#   b_mat <- MASS::mvrnorm(n_id, mu = c(0, 0, 0), Sigma = D_true)
#   colnames(b_mat) <- c("b_z1I0", "b_z1Ipost", "b_z2")
  
#   ## =========================================================
#   ## 4. PARTE ORDINALE (3 visite: 0, 6, 12 mesi)
#   ## =========================================================
  
#   ord <- expand.grid(
#     ID        = ID,
#     visit_idx = 0:2
#   )
#   ord <- ord[order(ord$ID, ord$visit_idx), ]
  
#   ord$I0    <- as.integer(ord$visit_idx == 0L)
#   ord$I6    <- as.integer(ord$visit_idx == 1L)
#   ord$I12   <- as.integer(ord$visit_idx == 2L)
#   ord$Ipost <- ord$I6 + ord$I12
  
#   ord$age       <- rep(age_std,   each = 3)
#   ord$loghosp   <- rep(loghosp,   each = 3)
#   ord$sex01     <- rep(sex01,     each = 3)
#   ord$diag_oa   <- rep(diag_oa,   each = 3)
#   ord$prost_uni <- rep(prost_uni, each = 3)
#   ord$waiting   <- rep(waiting,   each = 3)
  
#   ord$wI6  <- ord$waiting * ord$I6
#   ord$wI12 <- ord$waiting * ord$I12
#   ord$pI6  <- ord$prost_uni * ord$I6
#   ord$pI12 <- ord$prost_uni * ord$I12
  
#   ord$z1I0    <- as.integer(ord$I0    == 1L)
#   ord$z1Ipost <- as.integer(ord$Ipost == 1L)
  
#   ord$b_z1I0    <- rep(b_mat[, "b_z1I0"],    each = 3)
#   ord$b_z1Ipost <- rep(b_mat[, "b_z1Ipost"], each = 3)
  
#   ## predittore fisso dell'item (stessa formula del tuo modello)
#   X_ord <- model.matrix(
#     ~ 0 + I6 + I12 + wI6 + wI12 + pI6 + pI12 + diag_oa + sex01 + loghosp + age,
#     data = ord
#   )
#   beta_ord_vec  <- beta_ord[colnames(X_ord)]
#   eta_ord_fixed <- as.vector(X_ord %*% beta_ord_vec)
  
#   ## aggiungo random effects come nel modello: z1I0 * b_z1I0 + z1Ipost * b_z1Ipost
#   eta_ord <- eta_ord_fixed +
#     ord$z1I0    * ord$b_z1I0 +
#     ord$z1Ipost * ord$b_z1Ipost
  
#   ## generazione risposta ordinale via variabile latente ~ N(eta_ord, 1)
#   cutpoints <- c(-Inf, phis_ord_true, Inf)
#   latent    <- rnorm(nrow(ord), mean = eta_ord, sd = 1)
#   y_cat     <- cut(latent, breaks = cutpoints, labels = FALSE)
  
#   ord$y_raw <- as.integer(y_cat)
#   ord$y     <- ord$y_raw
#   ord$type  <- 0L
#   ord$time  <- ord$visit_idx
#   ord$item  <- factor("Gettingupfromsitting")
  
#   ## =========================================================
#   ## 5. PARTE SURVIVAL (Weibull-PH + frailty lognormale)
#   ## =========================================================
  
#   surv <- data.frame(
#     ID        = ID,
#     start     = 0,
#     sex01     = sex01,
#     diag_oa   = diag_oa,
#     prost_uni = prost_uni,
#     loghosp   = loghosp,
#     age       = age_std,
#     waiting   = waiting,
#     I0        = 1L,
#     I6        = 0L,
#     I12       = 0L,
#     Ipost     = 0L
#   )
#   surv$wI6  <- 0
#   surv$wI12 <- 0
#   surv$pI6  <- 0
#   surv$pI12 <- 0
#   surv$z1I0    <- 0L
#   surv$z1Ipost <- 0L
#   surv$z2      <- 1L
  
#   ## MODIFICA 2 (CONCETTUALE, ma il codice resta uguale):
#   ## qui manteniamo eta_surv = beta_surv * age + b_z2,
#   ## esattamente come nel modello GLMMadaptive.
#   ## NON facciamo shift -sigma^2/2, perché quello sarebbe
#   ## la parametrizzazione con frailty a media 1 (altro modello).
  
#   eta_surv <- beta_surv * surv$age + b_mat[, "b_z2"]
  
#   u       <- runif(n_id)
#   T_event <- (-log(u) / (lambda * exp(eta_surv)))^(1 / rho_true)
  
#   stop_time <- pmin(T_event, t_cens)
#   status    <- as.integer(T_event <= t_cens)
  
#   surv$stop      <- stop_time
#   surv$status    <- status
#   surv$visit_idx <- NA_integer_
#   surv$time      <- NA_integer_
#   surv$y         <- 0L
#   surv$type      <- 1L
#   surv$item      <- factor("Gettingupfromsitting")
  
#   ## =========================================================
#   ## 6. COSTRUZIONE DATASET LONG (stessa struttura del tuo modello)
#   ## =========================================================
  
#   ord_long <- ord
#   ord_long$start  <- NA_real_
#   ord_long$stop   <- NA_real_
#   ord_long$status <- NA_integer_
#   ord_long$z2     <- 0L
  
#   cols <- c("ID","visit_idx","I0","I6","I12","Ipost",
#             "waiting","loghosp","age","sex01","diag_oa","prost_uni",
#             "wI6","wI12","pI6","pI12",
#             "z1I0","z1Ipost","z2",
#             "item","type","time","y","start","stop","status")
  
#   long <- rbind(
#     ord_long[, cols],
#     surv[,      cols]
#   )
  
#   long <- long[order(long$ID, long$type, long$time), ]
#   rownames(long) <- NULL
  
#   long$row_id <- seq_len(nrow(long))
#   long$Y <- cbind(
#     y   = ifelse(long$type == 0L, long$y, 0L),
#     idx = long$row_id
#   )
#   long$is_ord  <- ifelse(long$type == 0L, 1L, 0L)
#   long$is_surv <- ifelse(long$type == 1L, 1L, 0L)
  
#   list(
#     long      = long,
#     ord       = ord,
#     surv      = surv,
#     n_events  = sum(surv$status),
#     lambda    = lambda,
#     rho       = rho_true
#   )
# }

# ## =========================================================
# ## 4. SIMULAZIONE CON ~24 EVENTI
# ## =========================================================




# set.seed(123)

# tuned <- tune_lambda_scale(
#   target_events = 24,
#   n_id          = 679,
#   t_cens        = 3650
# )

# tuned$lambda_scale  # fattore di scala “giusto”
# tuned$sim$n_events  # numero eventi nell'ultima simulazione

# long_sim <- tuned$sim$long
# ord_sim  <- tuned$sim$ord
# surv_sim <- tuned$sim$surv

# ## =========================================================
# ## 5. FAMILY + FORMULE COME NEL TUO CODICE
# ## =========================================================

# source("/home/niccolo.cao/ord+surv/jm_defs.R")

# K_sim <- max(long_sim$y[long_sim$type == 0L], na.rm = TRUE)

# fam_sim <- cprobit_weibull_indexed(
#   K          = K_sim,
#   type_vec   = long_sim$type,
#   start_vec  = long_sim$start,
#   stop_vec   = long_sim$stop,
#   status_vec = long_sim$status
# )

# fixed_fml_sim <- Y ~ 0 +
#   is_ord:(I6 + I12 +
#             wI6 + wI12 +
#             pI6 + pI12 + diag_oa +
#             sex01 + loghosp +
#             age) +
#   is_surv:(age)

# random_fml_sim <- ~ 0 + z1I0 + z1Ipost + z2 | ID


# ## =========================================================
# ## 6. INIZIALI "COME IL TUO CODICE"
# ##    - univariato ordinale
# ##    - WeibullPH frailty
# ## =========================================================

# ## soglie iniziali dal polr
# zeta_init <- tryCatch({
#   fit0 <- suppressWarnings(
#     polr(ordered(ord_sim$y) ~ 1, method = "probit", Hess = FALSE)
#   )
#   as.numeric(fit0$zeta)
# }, error = function(e) {
#   qnorm(seq(1 / K_sim, (K_sim - 1) / K_sim, length.out = K_sim - 1L))
# })

# univariateGLMM <- GLMMadaptive::mixed_model(
#   fixed  = y ~ 0 + I6 + I12 + wI6 + wI12 + pI6 + pI12 +
#     loghosp + sex01 + age + diag_oa,
#   random = ~ 0 + z1I0 + z1Ipost | ID,
#   data   = ord_sim,
#   family = cumulativeprobit(K_sim),
#   n_phis = K_sim - 1L,
#   initial_values = list(phis = zeta_init),
#   control = list(iter_EM = 10, nAGQ = 9)
# )

# surv_dat <- surv_sim

# wb0 <- readRDS("/home/niccolo.cao/ord+surv/wb0fail.rds")
# bet0 <- stats::coef(wb0)[3:length(coef(wb0))]
# phi0 <- c(wb0$res.t["scale","est"], wb0$res.t["shape","est"])

# weibullPHfrailty <- GLMMadaptive::mixed_model(
#   fixed      = cbind(start, stop, status) ~ 0 + age,
#   random     = ~ 1 | ID,
#   data       = surv_dat,
#   family     = weibullPH_lognormal(),
#   n_phis     = 2L,
#   initial_values = list(phis = phi0),
#   control    = list(iter_EM = 0, nAGQ = 20)
# )

# ## betas + phis iniziali come nel tuo codice
# coefsss <- c(univariateGLMM$coefficients,
#              weibullPHfrailty$coefficients)
# names(coefsss) <- c(
#   paste0("is_ord:", names(univariateGLMM$coefficients)),
#   paste0("is_surv:", "age")
# )

# phis_initsurv <- c(univariateGLMM$phis,-9.1484155, -0.2315199 )# weibullPHfrailty$phis)
# names(phis_initsurv) <- paste0("phi_", seq_len((K_sim - 1L) + 2L))

 

# ## D0s come nel tuo D0s
# D0s <- rbind(
#   cbind(univariateGLMM$D, -0.1),
#   matrix(c(rep(-0.1, 2), 1.30), nrow = 1)
# )
# colnames(D0s)[3] <- "z2"
# rownames(D0s)[3] <- "z2"


# ## =========================================================
# ## 7. FIT JOINT SUI DATI SIMULATI
# ## =========================================================

# fit_sim <- GLMMadaptive::mixed_model(
#   fixed  = fixed_fml_sim,
#   random = random_fml_sim,
#   data   = long_sim,
#   family = fam_sim,
#   n_phis = (K_sim - 1L) + 2L,
#   initial_values = list(
#     betas = coefsss,
#     phis  = phis_initsurv,
#     D     = D0s
#   ),
#   control = list(
#     iter_EM      = 0,
#     nAGQ         = 9,
#     verbose      = TRUE
#   )
# )

# summary(fit_sim)
# fit_sim$coefficients
# fit_sim$phis
# fit_sim$D














































# tune_lambda_scale <- function(target_events = 24,
#                               n_id         = 679,
#                               t_cens       = 3650,
#                               max_iter     = 8) {
#   lambda_scale <- 1
#   sim_last     <- NULL
  
#   for (it in seq_len(max_iter)) {
#     sim <- simulate_item_surv_once(
#       n_id         = n_id,
#       t_cens       = t_cens,
#       lambda_scale = lambda_scale
#     )
#     n_ev    <- sim$n_events
#     sim_last <- sim
    
#     cat("iter", it,
#         "  lambda_scale =", signif(lambda_scale, 4),
#         "  events =", n_ev, "\n")
    
#     if (n_ev == 0) break
#     if (abs(n_ev - target_events) <= 1) break
    
#     lambda_scale <- lambda_scale * (target_events / n_ev)
#   }
  
#   list(
#     lambda_scale = lambda_scale,
#     sim          = sim_last
#   )
# }
