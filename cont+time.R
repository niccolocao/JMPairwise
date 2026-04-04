suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(GLMMadaptive); library(MASS); library(survival)
})

gaussian_weibullPH_indexed <- function(type_vec, start_vec, stop_vec, status_vec,
                                       link="identity") {
  stats <- make.link(link); eps <- 1e-12
  force(type_vec); force(start_vec); force(stop_vec); force(status_vec)

  log_dens <- function(Y, eta, mu_fun, phis, eta_zi) {
    Y <- as.matrix(Y)
    if (ncol(Y) != 2L) stop
    y_main <- as.numeric(Y[,1]); idx <- as.integer(Y[,2])

    if (length(phis) != 3L) stop
    logsigma <- phis[1]
    log_lambda <- phis[2]
    log_rho <- phis[3]

    sigma  <-  exp(logsigma)
    lambda <- exp(log_lambda)
    rho    <- exp(log_rho)

    mu <- mu_fun(eta)
    n  <- length(y_main)
    out <- if (is.matrix(mu)) matrix(0, n, ncol(mu)) else numeric(n)

    type <- type_vec[idx]

    ic <- which(type == 0L)
    if (length(ic)) {
      if (is.matrix(mu)) {
        et <- mu[ic,, drop=FALSE]
        yy <- matrix(y_main[ic], nrow=length(ic), ncol=ncol(et))
        out[ic,] <- dnorm(yy, mean=et, sd=sigma, log=TRUE)
      } else {
        et <- mu[ic]
        out[ic] <- dnorm(y_main[ic], mean=et, sd=sigma, log=TRUE)
      }
    }

    isv <- which(type == 1L)
    if (length(isv)) {
      start <- pmax(start_vec[idx[isv]], 0)
      stop  <- pmax(stop_vec[idx[isv]],  1e-18)
      stat  <- status_vec[idx[isv]]

      H0    <- lambda * ( exp(rho*log(stop)) - exp(rho*log(pmax(start, 1e-18))) )
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
    family   = "gaussian + WeibullPH (indexed)",
    link     = stats$name,
    linkfun  = stats$linkfun,
    linkinv  = stats$linkinv,
    log_dens = log_dens
  ), class="family")
}

.extract_re_summary_cont_surv <- function(fit) {
  V <- tryCatch(VarCorr(fit)$ID, error=function(e) NULL)
  if (is.null(V)) V <- tryCatch(fit$D, error=function(e) NULL)
  if (is.null(V)) {
    return(list(
      sds=c(z1I0=NA_real_, z1Ipost=NA_real_, z2=NA_real_),
      corr=c(cont_I0_vs_Ipost=NA_real_, I0_vs_frailty=NA_real_, Ipost_vs_frailty=NA_real_)
    ))
  }
  V <- as.matrix(V); nm <- rownames(V)
  getc <- function(a,b) if (!(a %in% nm) || !(b %in% nm)) NA_real_ else V[a,b]/sqrt(V[a,a]*V[b,b])
  sds <- setNames(rep(NA_real_,3), c("z1I0","z1Ipost","z2"))
  d <- diag(V); names(d) <- nm
  present <- intersect(names(sds), names(d)); sds[present] <- sqrt(d[present])
  corr <- c(
    cont_I0_vs_Ipost = getc("z1I0","z1Ipost"),
    I0_vs_frailty    = getc("z1I0","z2"),
    Ipost_vs_frailty = getc("z1Ipost","z2")
  )
  list(sds=sds, corr=corr)
}

fit_one_cont_surv <- function(dat,
                              cont_item,
                              ctrl_joint = list(iter_EM = 0, nAGQ = 15, verbose = FALSE),
                              scale100 = TRUE,
                              eps_y = 0,
                              admin_date = as.POSIXct("2025-01-01", tz = "UTC"),
                              wb0_path = "wb0fail.rds") {

  suppressPackageStartupMessages({
    library(dplyr); library(tidyr)
    library(lme4)
    library(GLMMadaptive); library(survival)
  })

 ######### data management ################
  
  base_long <- dat %>%
  group_by(ID) %>%
  arrange(times, .by_group = TRUE) %>%
  mutate(
    t_within  = times - first(times),
    visit_idx = row_number() - 1L,
    loghosp   = log(1 + hospitalizationdays),
    waiting   = as.numeric(difftime(date_intervention, date_baseline, units = "days")) / 30.4375,
    time_discrete = factor(visit_idx, levels = c(0, 1, 2), labels = c("t0", "t6", "t12")),
    I0    = as.integer(time_discrete == "t0"),
    I6    = as.integer(time_discrete == "t6"),
    I12   = as.integer(time_discrete == "t12"),
    Ipost = I6 + I12,
    loghosp_post = loghosp * Ipost
  ) %>%
  ungroup()

  surv_core <- dat %>%
    mutate(
      t_event = as.numeric(difftime(date_revision, date_intervention, units = "days"))/30.4375,
      t_admin = as.numeric(difftime(admin_date, date_intervention, units = "days"))/30.4375
    ) %>%
    group_by(ID) %>%
    summarise(
      t_event = {v <- t_event[!is.na(t_event) & t_event >= 0]; if (length(v)) min(v) else NA_real_},
      t_admin = min(t_admin, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      status = as.integer(!is.na(t_event) & (t_event <= t_admin)),
      stop   = ifelse(status == 1L, t_event, t_admin),
      start  = 0
    ) %>%
    transmute(
      ID = factor(ID),
      start = start,
      stop = pmax(as.numeric(stop), 1e-8),
      status = status
    )

  dat2 <- base_long %>%
  left_join(surv_core %>% mutate(t_event = ifelse(status == 1L, stop, NA_real_)), by = "ID") %>%
  relocate(start, stop, status, t_event, .after = ID) %>%
  tidyr::drop_na(start, stop, status, sex01, prost_uni, loghosp_post, age, diag_oa)

  dat0 <- dat2 %>%
    mutate(
      ID = as.factor(ID),
      visit_idx = coalesce(visit_idx, 0L),
      I0 = as.integer(I0), I6 = as.integer(I6), I12 = as.integer(I12), Ipost = as.integer(Ipost),
      prost_uni = as.integer(prosthesis_type == "unicompartimental"),
      age = as.numeric(scale(age))  
    )

                
  cont <- dat0 %>%
  dplyr::select(
    ID, visit_idx, I0, I6, I12, Ipost,
    waiting, loghosp_post, age, sex01, diag_oa, prost_uni,
    start, stop, status,
    dplyr::all_of(cont_item)
  ) %>%
  dplyr::filter(!is.na(.data[[cont_item]])) %>%
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
      visit_idx == 0L ~ 0L,
      visit_idx == 1L ~ 1L,
      visit_idx == 2L ~ 2L,
      TRUE ~ NA_integer_
    ),
    y_raw = as.numeric(.data[[cont_item]]),
    y = if (scale100) y_raw / 100 else y_raw,
    y = ifelse(is.na(y), NA_real_, pmin(pmax(y, 0), 1 - eps_y)),
    wI6  = waiting * I6,
    wI12 = waiting * I12,
    pI6  = prost_uni * I6,
    pI12 = prost_uni * I12,
    z1I0    = ifelse(I0 == 1L, 1, 0),
    z1Ipost = ifelse(Ipost == 1L, 1, 0)
  ) %>%
  dplyr::select(-dplyr::all_of(cont_item), -y_raw)

                 
  surv <- dat0 %>%
  group_by(ID) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    ID,
    visit_idx = NA_integer_,
    I0, I6, I12, Ipost,
    waiting,
    wI6 = waiting * I6, wI12 = waiting * I12,
    pI6 = prost_uni * I6, pI12 = prost_uni * I12,
    loghosp_post, sex01, diag_oa,
    age, prost_uni,
    type = 1L, time = NA_integer_,
    y = 0,
    start, stop, status,
    z1I0 = 0, z1Ipost = 0
  )

  long <- cont %>%
    mutate(start = NA_real_, stop = NA_real_, status = NA_integer_) %>%
    bind_rows(surv) %>%
    arrange(ID, type, time)

  long <- long %>%
    mutate(
      is_cont = ifelse(type == 0L, 1, 0),
      is_surv = ifelse(type == 1L, 1, 0),
      z2 = ifelse(type == 1L, 1, 0)
    )

  long$row_id <- seq_len(nrow(long))
  long$Y <- cbind(y = ifelse(long$type == 0L, long$y, 0), idx = long$row_id)


  ############## Univariate models for initial values ################
                 
  fam <- gaussian_weibullPH_indexed(
    type_vec = long$type,
    start_vec = long$start,
    stop_vec  = long$stop,
    status_vec = long$status
  )

                 
  fixed_fml <- Y ~ 0 +
  is_cont:(1 + I6 + I12 +
           wI6 + wI12 +
           pI6 + pI12 +
           diag_oa + sex01 + loghosp_post + age) +
  is_surv:(age)

  random_fml <- ~ 0 + z1I0 + z1Ipost + z2 | ID

                 
  fit_lmer <- lme4::lmer(
  y ~ 1 + I6 + I12 +
    wI6 + wI12 + pI6 + pI12 +
    loghosp_post + sex01 + age + diag_oa +
    (0 + z1I0 + z1Ipost | ID),
  data = cont,
  REML = FALSE
)
  b_cont <- lme4::fixef(fit_lmer)
  D_cont <- as.matrix(lme4::VarCorr(fit_lmer)$ID)
  sigma_hat <- sigma(fit_lmer)
  logsigma_init <- log(sigma_hat)

                 
  surv_dat <- long[long$type == 1L, ]
  wb0 <- readRDS(wb0_path)
  bet0 <- stats::coef(wb0)[3:length(coef(wb0))]
  phi0 <- c(wb0$res.t["scale","est"], wb0$res.t["shape","est"])

  weibullPHfrailty <- GLMMadaptive::mixed_model(
    fixed  = cbind(start, stop, status) ~ 0 + age,
    random = ~ 1 | ID,
    data   = surv_dat,
    family = weibullPH_lognormal(),
    n_phis = 2L,
    control = list(iter_EM = 0, nAGQ = 20),
    initial_values = list(betas = bet0, phis = phi0, D = matrix(0.4^2, 1, 1))
  )

                 
  X_joint <- model.matrix(fixed_fml, data = long)
  betas_init <- setNames(rep(0, ncol(X_joint)), colnames(X_joint))

  b_cont2 <- b_cont
  names(b_cont2) <- paste0("is_cont:", names(b_cont2))
  betas_init[intersect(names(b_cont2), names(betas_init))] <- b_cont2[intersect(names(b_cont2), names(betas_init))]

  b_surv <- weibullPHfrailty$coefficients
  names(b_surv) <- paste0("is_surv:", names(b_surv))
  betas_init[intersect(names(b_surv), names(betas_init))] <- b_surv[intersect(names(b_surv), names(betas_init))]

                 
  phis_init <- c(logsigma_init, as.numeric(weibullPHfrailty$phis))

                 
  var_z2 <- as.numeric(weibullPHfrailty$D)[1]
  D0s <- rbind(
    cbind(D_cont, c(0, 0)),
    c(0, 0, var_z2)
  )
  colnames(D0s) <- rownames(D0s) <- c("z1I0","z1Ipost","z2")

  stopifnot(length(betas_init) == ncol(X_joint), length(phis_init) == 3L, all(dim(D0s) == c(3L,3L)))

                 
  fit_joint <- GLMMadaptive::mixed_model(
    fixed  = fixed_fml,
    random = random_fml,
    data   = long,
    family = fam,
    n_phis = 3L,
    initial_values = list(betas = betas_init, phis = phis_init, D = D0s),
    control = ctrl_joint
  )

  re <- .extract_re_summary_cont_surv(fit_joint)

  list(
    fit = fit_joint,
    sigma  =  exp(fit_joint$phis[1]),
    lambda = exp(fit_joint$phis[2]),
    rho    = exp(fit_joint$phis[3]),
    re_sds = re$sds,
    re_corr = re$corr
  )
}

source("univariateFuns.R")

res <- fit_one_cont_surv(
  dat = readRDS("dat2_with_hosp30g.rds"),
  cont_item = "Perceivedcurrenthealthstatus",
  ctrl_joint = list(iter_EM = 0, nAGQ = 25, verbose = TRUE),
  scale100 = TRUE,
  wb0_path = "wb0fail.rds" #already fitted Weibull PH model with age as fixed effect (no frailty)
)


saveRDS(res, "Perceivedcurrenthealthstatus__SURV_CONT.rds")
