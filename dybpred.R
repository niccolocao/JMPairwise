library(MASS)
library(Matrix)
library(mvtnorm)




prep_long_all_items_surv_mixed2 <- function(
    dat,
    ord_itemnames,
    cont_itemnames = NULL,
    cont_transforms = NULL,
    surv_outcome,
    admin_date = "2025-01-01",
    rehosp_horizon_days = 30
) {
  tmp <- prep_long_all_items_surv(
    dat = dat,
    itemnames = ord_itemnames,
    surv_outcome = surv_outcome,
    admin_date = admin_date,
    rehosp_horizon_days = rehosp_horizon_days
  )
  
  if (is.null(cont_itemnames) || length(cont_itemnames) == 0) {
    if (!"is_cont" %in% names(tmp$long)) tmp$long$is_cont <- 0L
    return(tmp)
  }
  
  cont_itemnames <- as.character(cont_itemnames)
  miss_cont <- setdiff(cont_itemnames, names(dat))
  if (length(miss_cont)) {
    stop("mancano queste colonne continue in dat: ", paste(miss_cont, collapse = ", "))
  }
  
  if (is.null(cont_transforms)) {
    cont_transforms <- stats::setNames(vector("list", length(cont_itemnames)), cont_itemnames)
  }
  if (is.null(names(cont_transforms))) {
    names(cont_transforms) <- cont_itemnames[seq_along(cont_transforms)]
  }
  
  if (!"is_cont" %in% names(tmp$long)) {
    tmp$long$is_cont <- 0L
  } else {
    tmp$long$is_cont <- ifelse(is.na(tmp$long$is_cont), 0L, as.integer(tmp$long$is_cont))
  }
  
  base_long <- dat %>%
    dplyr::mutate(
      ID = factor(ID),
      sex01 = as.integer(sex == "M"),
      prost_uni = as.integer(prosthesis_type == "unicompartimental"),
      diag_oa = as.integer(diagnosis != "osteoarthritis_knee")
    ) %>%
    dplyr::group_by(ID) %>%
    dplyr::arrange(times, .by_group = TRUE) %>%
    dplyr::mutate(
      t_within = times - dplyr::first(times),
      visit_idx = dplyr::row_number() - 1L,
      post = as.integer(visit_idx >= 1L),
      time_post = ifelse(post == 1L, t_within, 0),
      loghosp = log1p(hospitalizationdays),
      waiting = as.numeric(difftime(date_intervention, date_baseline, units = "days")) / 30.4375,
      time_discrete = factor(visit_idx, levels = c(0, 1, 2), labels = c("t0", "t6", "t12")),
      I0 = as.integer(time_discrete == "t0"),
      I6 = as.integer(time_discrete == "t6"),
      I12 = as.integer(time_discrete == "t12"),
      Ipost = as.integer(I6 + I12 > 0)
    ) %>%
    dplyr::ungroup()
  
  if ("age" %in% names(base_long) &&
      !is.null(tmp$age_center) &&
      !is.null(tmp$age_scale)) {
    base_long$age <- (as.numeric(base_long$age) - tmp$age_center) / tmp$age_scale
  }
  
  cont_blocks <- lapply(cont_itemnames, function(it) {
    tfun <- cont_transforms[[it]]
    
    y_raw <- base_long[[it]]
    y_use <- if (is.null(tfun)) y_raw else tfun(y_raw)
    
    dplyr::tibble(
      ID = base_long$ID,
      visit_idx = base_long$visit_idx,
      I0 = base_long$I0,
      I6 = base_long$I6,
      I12 = base_long$I12,
      Ipost = base_long$Ipost,
      waiting = base_long$waiting,
      loghosp = base_long$loghosp,
      age = base_long$age,
      sex01 = base_long$sex01,
      diag_oa = base_long$diag_oa,
      prost_uni = base_long$prost_uni,
      item = it,
      y_raw = as.numeric(y_raw),
      type = 0L,
      time = dplyr::case_when(
        base_long$visit_idx == 0L ~ 0L,
        base_long$visit_idx == 1L ~ 1L,
        base_long$visit_idx == 2L ~ 2L,
        TRUE ~ NA_integer_
      ),
      wI6 = base_long$waiting * base_long$I6,
      wI12 = base_long$waiting * base_long$I12,
      hIpost = ifelse(base_long$Ipost == 1L, base_long$loghosp, 0),
      pI6 = base_long$prost_uni * base_long$I6,
      pI12 = base_long$prost_uni * base_long$I12,
      z1I0 = ifelse(base_long$I0 == 1L, 1L, 0L),
      z1Ipost = ifelse(base_long$Ipost == 1L, 1L, 0L),
      y = as.numeric(y_use),
      start = NA_real_,
      stop = NA_real_,
      status = NA_integer_,
      loghosp_post = base_long$loghosp * base_long$Ipost,
      is_ord = 0L,
      is_surv = 0L,
      is_cont = 1L,
      z2 = 0L,
      K_item = NA_integer_,
      is_K3 = 0L,
      is_K5 = 0L
    ) %>%
      dplyr::filter(!is.na(y), is.finite(y))
  })
  
  cont_long <- dplyr::bind_rows(cont_blocks)
  
  long2 <- dplyr::bind_rows(tmp$long, cont_long)
  
  long2$item <- factor(
    as.character(long2$item),
    levels = c(ord_itemnames, cont_itemnames, "SURV")
  )
  
  long2$is_cont <- ifelse(is.na(long2$is_cont), 0L, as.integer(long2$is_cont))
  
  long2 <- long2 %>%
    dplyr::arrange(ID, is_surv, visit_idx, item)
  
  long2$row_id <- seq_len(nrow(long2))
  long2$Y <- cbind(
    y = ifelse(long2$type == 0L, long2$y, 0),
    idx = long2$row_id
  )
  
  tmp$long <- long2
  tmp
}

prep_long_all_items_surv <- function(dat,
                                     itemnames,
                                     surv_outcome,              # OBBLIGATORIO: "failure" o "rehosp30"
                                     admin_date = "2025-01-01",  # usato SOLO se surv_outcome == "failure"
                                     rehosp_horizon_days = 30) { # finestra fissa 30gg (censura per tutti)
  
  if (missing(surv_outcome)) {
    stop("Devi passare surv_outcome = 'failure' oppure surv_outcome = 'rehosp30'.")
  }
  surv_outcome <- match.arg(surv_outcome, choices = c("failure","rehosp30"))
  
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr)
  })
  
  ## 1) ricodifiche
  dat <- dat %>%
    mutate(
      ID        = factor(ID),
      sex01     = as.integer(sex == "M"),
      prost_uni = as.integer(prosthesis_type == "unicompartimental"),
      diag_oa   = as.integer(diagnosis != "osteoarthritis_knee")
    ) %>%
    filter(!if_all(everything(), is.na))
  
  ## 2) longitudinal core
  base_long <- dat %>%
    group_by(ID) %>%
    arrange(times, .by_group = TRUE) %>%
    mutate(
      t_within   = times - first(times),
      visit_idx  = row_number() - 1L,
      post       = as.integer(visit_idx >= 1L),
      time_post  = ifelse(post == 1L, t_within, 0),
      loghosp    = log(1 + hospitalizationdays),
      waiting    = as.numeric(difftime(date_intervention, date_baseline, units = "days")) / 30.4375,
      time_discrete = factor(visit_idx, levels = c(0,1,2), labels = c("t0","t6","t12")),
      I0   = as.integer(time_discrete == "t0"),
      I6   = as.integer(time_discrete == "t6"),
      I12  = as.integer(time_discrete == "t12"),
      Ipost = I6 + I12
    ) %>%
    ungroup()
  
  ## 3) survival core (branch)
  if (surv_outcome == "failure") {
    
    # --- failure protesi: usa admin_date per censura amministrativa ---
    admin_date <- as.POSIXct(admin_date, tz = attr(dat$date_intervention, "tzone"))
    
    surv_core <- dat %>%
      mutate(
        t_event_days = as.numeric(difftime(date_revision, date_intervention, units = "days")),
        t_admin_days = as.numeric(difftime(admin_date,     date_intervention, units = "days"))
      ) %>%
      group_by(ID) %>%
      summarise(
        t_event = {
          v <- t_event_days[!is.na(t_event_days) & is.finite(t_event_days) & t_event_days >= 0]
          if (length(v)) min(v) else NA_real_
        },
        t_admin = {
          v <- t_admin_days[!is.na(t_admin_days) & is.finite(t_admin_days)]
          if (length(v)) min(v) else NA_real_
        },
        .groups = "drop"
      ) %>%
      mutate(
        status = as.integer(!is.na(t_event) & is.finite(t_admin) & (t_event <= t_admin)),
        stop   = ifelse(status == 1L, t_event, t_admin),
        start  = 0
      ) %>%
      transmute(
        ID     = factor(ID),
        start  = start,
        stop   = pmax(as.numeric(stop), 1e-8),
        status = as.integer(status)
      )
    
  } else {
    
    # --- rehosp30: censura PER TUTTI a 30 giorni (NO admin_date) ---
    if (!"date_hosp30g" %in% names(dat)) {
      stop("Manca la colonna 'date_hosp30g' nel dataset.")
    }
    
    H <- as.numeric(rehosp_horizon_days)
    
    surv_core <- dat %>%
      mutate(
        t_event_days = as.numeric(difftime(date_hosp30g, date_intervention, units = "days"))
      ) %>%
      group_by(ID) %>%
      summarise(
        t_event = {
          v <- t_event_days[!is.na(t_event_days) & is.finite(t_event_days) & t_event_days >= 0]
          if (length(v)) min(v) else NA_real_
        },
        .groups = "drop"
      ) %>%
      mutate(
        status = as.integer(!is.na(t_event) & (t_event <= H)),
        stop   = pmax(ifelse(status == 1L, t_event, H), 1e-8),
        start  = 0
      ) %>%
      transmute(
        ID     = factor(ID),
        start  = start,
        stop   = stop,
        status = status
      )
  }
  
  dat2 <- base_long %>%
    left_join(surv_core, by = "ID") %>%
    relocate(start, stop, status, .after = ID) %>%
    tidyr::drop_na(start, stop, status, sex01, prost_uni, loghosp, age, diag_oa)
  
  ## 4) dat0
  dat0 <- dat2 %>%
    mutate(
      ID        = as.factor(ID),
      visit_idx = coalesce(visit_idx, 0L),
      I0        = as.integer(I0),
      I6        = as.integer(I6),
      I12       = as.integer(I12),
      Ipost     = as.integer(Ipost),
      prost_uni = as.integer(prosthesis_type == "unicompartimental")
    )
  
  ## 5) ORD: pivot_longer su tutti gli item
  ord <- dat0 %>%
    dplyr::select(
      ID, visit_idx, I0, I6, I12, Ipost,
      waiting, loghosp, age, sex01, diag_oa, prost_uni,
      dplyr::all_of(itemnames)
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(itemnames),
      names_to = "item",
      values_to = "y_raw",
      values_transform = list(y_raw = ~{
        if (is.factor(.x)) as.integer(as.character(.x)) else as.integer(.x)
      })
    ) %>%
    dplyr::filter(!is.na(y_raw)) %>%
    dplyr::mutate(
      type = 0L,
      time = dplyr::case_when(
        visit_idx == 0 ~ 0L,
        visit_idx == 1 ~ 1L,
        visit_idx == 2 ~ 2L,
        TRUE ~ NA_integer_
      ),
      item  = factor(item, levels = itemnames),
      wI6   = waiting * I6,
      wI12  = waiting * I12, 
      hIpost = dplyr::if_else(time > 0L, loghosp, 0),
      pI6   = prost_uni * I6,
      pI12  = prost_uni * I12,
      z1I0    = ifelse(I0 == 1L, 1L, 0L),
      z1Ipost = ifelse(Ipost == 1L, 1L, 0L)
    ) %>%
    dplyr::group_by(item) %>%
    dplyr::mutate(
      y = as.integer(y_raw - min(y_raw, na.rm = TRUE) + 1L)
    ) %>%
    dplyr::ungroup()
  
  ## 6) SURV: una riga per ID
  surv <- dat0 %>%
    group_by(ID) %>% slice_head(n = 1) %>% ungroup() %>%
    transmute(
      ID,
      visit_idx = NA_integer_,
      I0, I6, I12, Ipost,
      waiting,
      wI6 = waiting*I6, wI12 = waiting*I12,
      pI6 = prost_uni*I6, pI12 = prost_uni*I12,
      sex01, diag_oa,
      age, prost_uni,
      loghosp_post = dplyr::if_else(Ipost == 1L, loghosp, 0),
      item = factor("SURV", levels = c(itemnames, "SURV")),
      type = 1L, time = NA_integer_,
      y = 0L,
      start, stop, status
    )
  
  ord2 <- ord %>%
    mutate(
      start = NA_real_, stop = NA_real_, status = NA_integer_,
      item  = factor(as.character(item), levels = c(itemnames, "SURV"))
    )
  
  long <- bind_rows(ord2, surv) %>%
    arrange(ID, type, time, item)
  
  ## 7) scaling age (su 1 riga per ID)
  age_ref <- dat0 %>% distinct(ID, age)
  m_age <- mean(age_ref$age, na.rm = TRUE)
  s_age <- sd(age_ref$age, na.rm = TRUE)
  
  long <- long %>%
    mutate(
      is_ord  = ifelse(type == 0L, 1L, 0L),
      is_surv = ifelse(type == 1L, 1L, 0L),
      z1I0    = ifelse(type == 0L & I0    == 1L, 1L, 0L),
      z1Ipost = ifelse(type == 0L & Ipost == 1L, 1L, 0L),
      z2      = ifelse(type == 1L, 1L, 0L),
      age     = (age - m_age) / s_age
    )
  
  long$row_id <- seq_len(nrow(long))
  long$Y <- cbind(
    y   = ifelse(long$type == 0L, long$y, 0L),
    idx = long$row_id
  )
  
  ## 8) K per item
  K_by_item <- ord %>%
    group_by(item) %>%
    summarise(K = max(y, na.rm = TRUE), .groups = "drop")
  
  ## join info K_item
  item_info <- tibble::tibble(
    item   = itemnames,
    K_item = vapply(dat[itemnames], function(x) {
      if (is.factor(x)) nlevels(x) else length(sort(unique(stats::na.omit(x))))
    }, integer(1))
  )
  
  long <- long %>%
    dplyr::mutate(item_chr = as.character(item)) %>%
    dplyr::left_join(item_info, by = c("item_chr" = "item")) %>%
    dplyr::mutate(
      K_item = dplyr::if_else(type == 1L, NA_integer_, K_item),
      is_K3  = dplyr::if_else(K_item == 3L, 1L, 0L, missing = 0L),
      is_K5  = dplyr::if_else(K_item == 5L, 1L, 0L, missing = 0L)
    ) %>%
    dplyr::select(-item_chr)
  
  list(
    long       = long,
    ord        = ord,
    surv       = surv,
    K_by_item  = K_by_item,
    age_center = m_age,
    age_scale  = s_age,
    surv_outcome = surv_outcome
  )
}







sanitize_name <- function(x) {
  x <- gsub("\\s+", "_", x)
  gsub("[^[:alnum:]_]+", "", x)
}

make_fullplug_object <- function(est,
                                 D,
                                 K_map,
                                 ord_items,
                                 cont_items = NULL,
                                 resid_corr = NULL,
                                 cont_scale = NULL,
                                 cont_transforms = NULL,
                                 time_unit = "months",
                                 theta_sampler = NULL,
                                 covariate_transforms = NULL) {
  
  if (is.null(names(est))) stop("est deve avere names().")
  if (is.null(rownames(D)) || is.null(colnames(D))) stop("D deve avere rownames/colnames.")
  
  sanitize_name <- function(x) {
    x <- gsub("\\s+", "_", x)
    gsub("[^[:alnum:]_]+", "", x)
  }
  
  ord_items <- as.character(ord_items)
  cont_items <- if (is.null(cont_items)) character(0) else as.character(cont_items)
  all_items <- c(ord_items, cont_items)
  
  if (is.null(cont_transforms)) {
    cont_transforms <- setNames(vector("list", length(cont_items)), cont_items)
  }
  
  get_surv_par <- function(name_plain, name_log, transform = exp) {
    if (name_plain %in% names(est)) return(unname(est[[name_plain]]))
    if (name_log %in% names(est)) return(transform(unname(est[[name_log]])))
    stop("parametro mancante: ", name_plain, " / ", name_log)
  }
  
  build_beta_item <- function(item_s) {
    pat <- paste0("^fixed_", item_s, "_")
    ii <- grep(pat, names(est))
    if (!length(ii)) return(setNames(numeric(0), character(0)))
    vv <- est[ii]
    nm <- sub(pat, "", names(vv))
    nm[nm == "Intercept"] <- "(Intercept)"
    names(vv) <- nm
    vv
  }
  
  build_thresholds <- function(item_s, K) {
    pat <- paste0("^threshold_", item_s, "_cut:")
    ii <- grep(pat, names(est))
    if (!length(ii)) stop("soglie mancanti per ", item_s)
    vv <- est[ii]
    cut_id <- as.integer(sub(pat, "", names(vv)))
    vv <- vv[order(cut_id)]
    unname(vv)
  }
  
  build_sigma_item <- function(item_s) {
    nm1 <- paste0("sigma_", item_s)
    nm2 <- paste0("log_sigma_", item_s)
    if (nm1 %in% names(est)) return(unname(est[[nm1]]))
    if (nm2 %in% names(est)) return(exp(unname(est[[nm2]])))
    stop("sigma mancante per ", item_s)
  }
  
  infer_resid_corr <- function(est, items_s) {
    R <- diag(1, length(items_s))
    dimnames(R) <- list(items_s, items_s)
    
    nm_corr <- grep("^corr_resid__", names(est), value = TRUE)
    for (nm in nm_corr) {
      ss <- sub("^corr_resid__+", "", nm)
      sp <- strsplit(ss, "__", fixed = TRUE)[[1]]
      if (length(sp) == 2L && all(sp %in% items_s)) {
        R[sp[1], sp[2]] <- est[[nm]]
        R[sp[2], sp[1]] <- est[[nm]]
      }
    }
    
    nm_atanh <- grep("^atanh_rho_resid__", names(est), value = TRUE)
    for (nm in nm_atanh) {
      ss <- sub("^atanh_rho_resid__+", "", nm)
      sp <- strsplit(ss, "__", fixed = TRUE)[[1]]
      if (length(sp) == 2L && all(sp %in% items_s)) {
        rr <- tanh(est[[nm]])
        R[sp[1], sp[2]] <- rr
        R[sp[2], sp[1]] <- rr
      }
    }
    
    R
  }
  
  surv_beta <- {
    ii <- grep("^fixed_surv_", names(est))
    vv <- est[ii]
    names(vv) <- sub("^fixed_surv_", "", names(vv))
    vv
  }
  
  outcomes <- vector("list", length(all_items))
  names(outcomes) <- sanitize_name(all_items)
  
  for (item in ord_items) {
    item_s <- sanitize_name(item)
    outcomes[[item_s]] <- list(
      name = item_s,
      original_name = item,
      type = "ordinal",
      beta = build_beta_item(item_s),
      thresholds = build_thresholds(item_s, K_map[[item]]),
      sigma = 1,
      scale_div = 1,
      transform = NULL,
      re_names = c(paste0(item_s, "_I0"), paste0(item_s, "_Ipost"))
    )
  }
  
  for (item in cont_items) {
    item_s <- sanitize_name(item)
    sc <- if (is.null(cont_scale)) 1 else if (item %in% names(cont_scale)) cont_scale[[item]] else 1
    tf <- if (item %in% names(cont_transforms)) cont_transforms[[item]] else NULL
    outcomes[[item_s]] <- list(
      name = item_s,
      original_name = item,
      type = "gaussian",
      beta = build_beta_item(item_s),
      thresholds = NULL,
      sigma = build_sigma_item(item_s),
      scale_div = sc,
      transform = tf,
      re_names = c(paste0(item_s, "_I0"), paste0(item_s, "_Ipost"))
    )
  }
  
  if (is.null(resid_corr)) {
    resid_corr <- infer_resid_corr(est, names(outcomes))
  } else {
    resid_corr <- as.matrix(resid_corr)
  }
  
  resid_corr <- resid_corr[names(outcomes), names(outcomes), drop = FALSE]
  D <- 0.5 * (as.matrix(D) + t(as.matrix(D)))
  
  list(
    est = est,
    outcomes = outcomes,
    D = D,
    resid_corr = resid_corr,
    survival = list(
      beta = surv_beta,
      lambda = get_surv_par("lambda", "log_lambda", exp),
      rho = get_surv_par("rho", "log_rho", exp),
      frailty_name = "frailty"
    ),
    theta_sampler = theta_sampler,
    time_unit = time_unit,
    covariate_transforms = covariate_transforms
  )
}

prep_newdata_fullplug <- function(newdata,
                                  object = NULL,
                                  idVar = "ID",
                                  timeVar = "visit_m",
                                  visit_idx_var = "visit_idx") {
  out <- newdata
  
  if (!timeVar %in% names(out)) {
    if (!visit_idx_var %in% names(out)) stop("manca sia ", timeVar, " che ", visit_idx_var)
    out[[timeVar]] <- c(0, 6, 12)[out[[visit_idx_var]] + 1L]
  }
  
  if (!"I0" %in% names(out)) out$I0 <- as.integer(abs(out[[timeVar]] - 0) < 1e-8)
  if (!"I6" %in% names(out)) out$I6 <- as.integer(abs(out[[timeVar]] - 6) < 1e-8)
  if (!"I12" %in% names(out)) out$I12 <- as.integer(abs(out[[timeVar]] - 12) < 1e-8)
  if (!"Ipost" %in% names(out)) out$Ipost <- as.integer(out$I6 + out$I12 > 0)
  
  if (!"prost_uni" %in% names(out) && "prosthesis_type" %in% names(out)) {
    out$prost_uni <- as.integer(out$prosthesis_type == "unicompartimental")
  }
  
  if (!"sex01" %in% names(out) && "sex" %in% names(out)) {
    out$sex01 <- as.integer(out$sex == "M")
  }
  
  if (!"diag_oa" %in% names(out) && "diagnosis" %in% names(out)) {
    out$diag_oa <- as.integer(out$diagnosis != "osteoarthritis_knee")
  }
  
  if (!"waiting" %in% names(out) &&
      all(c("date_intervention", "date_baseline") %in% names(out))) {
    out$waiting <- as.numeric(difftime(out$date_intervention, out$date_baseline, units = "days")) / 30.4375
  }
  
  if (!"loghosp" %in% names(out) && "hospitalizationdays" %in% names(out)) {
    out$loghosp <- log1p(out$hospitalizationdays)
  }
  
  if (!"loghosp_post" %in% names(out)) {
    if ("loghosp" %in% names(out)) {
      out$loghosp_post <- out$loghosp * out$Ipost
    } else if ("hospitalizationdays" %in% names(out)) {
      out$loghosp_post <- log1p(out$hospitalizationdays) * out$Ipost
    }
  }
  
  if (!"wI6" %in% names(out)) {
    if (!"waiting" %in% names(out)) stop("manca 'waiting', quindi non posso costruire wI6.")
    out$wI6 <- out$waiting * out$I6
  }
  
  if (!"wI12" %in% names(out)) {
    if (!"waiting" %in% names(out)) stop("manca 'waiting', quindi non posso costruire wI12.")
    out$wI12 <- out$waiting * out$I12
  }
  
  if (!"pI6" %in% names(out)) {
    if (!"prost_uni" %in% names(out)) stop("manca 'prost_uni', quindi non posso costruire pI6.")
    out$pI6 <- out$prost_uni * out$I6
  }
  
  if (!"pI12" %in% names(out)) {
    if (!"prost_uni" %in% names(out)) stop("manca 'prost_uni', quindi non posso costruire pI12.")
    out$pI12 <- out$prost_uni * out$I12
  }
  
  if (!is.null(object) && !is.null(object$covariate_transforms)) {
    for (nm in names(object$covariate_transforms)) {
      tfun <- object$covariate_transforms[[nm]]
      if (!is.null(tfun) && nm %in% names(out)) {
        out[[nm]] <- tfun(out[[nm]])
      }
    }
  }
  
  out <- out[order(out[[idVar]], out[[timeVar]]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.fullplug_xbeta <- function(row, beta) {
  if (!length(beta)) return(0)
  out <- 0
  for (nm in names(beta)) {
    if (nm %in% c("(Intercept)", "Intercept")) {
      out <- out + beta[[nm]]
    } else {
      if (!nm %in% names(row)) stop("covariata mancante in newdata: ", nm)
      out <- out + beta[[nm]] * as.numeric(row[[nm]][1])
    }
  }
  out
}

.fullplug_mu <- function(object, row, b, outcome_name) {
  outk <- object$outcomes[[outcome_name]]
  mu <- .fullplug_xbeta(row, outk$beta)
  rn <- outk$re_names
  if (!all(rn %in% names(b))) stop("random effects mancanti per ", outcome_name)
  mu <- mu +
    as.numeric(row$I0[1]) * b[[rn[1]]] +
    as.numeric(row$Ipost[1]) * b[[rn[2]]]
  mu
}

.fullplug_surv_lp <- function(object, row, b) {
  eta <- .fullplug_xbeta(row, object$survival$beta)
  eta + b[[object$survival$frailty_name]]
}

.fullplug_logS <- function(object, row, b, t) {
  lam <- object$survival$lambda
  rho <- object$survival$rho
  eta <- .fullplug_surv_lp(object, row, b)
  - lam * (t^rho) * exp(eta)
}

.fullplug_Sratio <- function(object, row, b, t, u) {
  lam <- object$survival$lambda
  rho <- object$survival$rho
  eta <- .fullplug_surv_lp(object, row, b)
  exp(- lam * ((u^rho) - (t^rho)) * exp(eta))
}

.fullplug_ord_bounds <- function(y, thr) {
  K <- length(thr) + 1L
  if (is.na(y)) return(c(NA_real_, NA_real_))
  y <- as.integer(y)
  if (y < 1L || y > K) stop("categoria ordinale fuori range.")
  lo <- if (y == 1L) -Inf else thr[y - 1L]
  up <- if (y == K) Inf else thr[y]
  c(lo, up)
}

.fullplug_make_obs <- function(object, row) {
  out_names <- names(object$outcomes)
  keep <- logical(length(out_names))
  vals <- vector("list", length(out_names))
  names(vals) <- out_names
  
  for (k in seq_along(out_names)) {
    nm <- out_names[k]
    ok <- object$outcomes[[nm]]
    raw_nm <- ok$original_name
    if (!raw_nm %in% names(row)) next
    
    yy <- row[[raw_nm]][1]
    if (is.na(yy)) next
    
    if (ok$type == "gaussian") {
      if (!is.null(ok$transform)) {
        yy <- ok$transform(as.numeric(yy))
      } else {
        yy <- as.numeric(yy) / ok$scale_div
      }
    } else {
      yy <- as.integer(yy)
    }
    
    vals[[nm]] <- yy
    keep[k] <- TRUE
  }
  
  vals[keep]
}

.fullplug_resid_cov <- function(object, obs_names) {
  p <- length(obs_names)
  S <- matrix(0, p, p, dimnames = list(obs_names, obs_names))
  
  for (j in seq_len(p)) {
    nmj <- obs_names[j]
    oj <- object$outcomes[[nmj]]
    S[j, j] <- if (oj$type == "ordinal") 1 else oj$sigma^2
  }
  
  if (p >= 2) {
    for (j in 2:p) {
      for (i in 1:(j - 1)) {
        ni <- obs_names[i]
        nj <- obs_names[j]
        oi <- object$outcomes[[ni]]
        oj <- object$outcomes[[nj]]
        
        rr <- 0
        if (ni %in% rownames(object$resid_corr) && nj %in% colnames(object$resid_corr)) {
          rr <- object$resid_corr[ni, nj]
        }
        
        if (oi$type == "ordinal" && oj$type == "ordinal") {
          cv <- rr
        } else if (oi$type == "ordinal" && oj$type == "gaussian") {
          cv <- rr * oj$sigma
        } else if (oi$type == "gaussian" && oj$type == "ordinal") {
          cv <- rr * oi$sigma
        } else {
          cv <- rr * oi$sigma * oj$sigma
        }
        
        S[i, j] <- cv
        S[j, i] <- cv
      }
    }
  }
  
  ee <- eigen(0.5 * (S + t(S)), symmetric = TRUE)
  ee$values[ee$values < 1e-8] <- 1e-8
  S <- ee$vectors %*% diag(ee$values, nrow = length(ee$values)) %*% t(ee$vectors)
  dimnames(S) <- list(obs_names, obs_names)
  S
}

.fullplug_visit_logdens <- function(object,
                                    row,
                                    b,
                                    pmvnorm_control = mvtnorm::GenzBretz(maxpts = 25000,
                                                                         abseps = 1e-04,
                                                                         releps = 0)) {
  obs <- .fullplug_make_obs(object, row)
  if (!length(obs)) return(0)
  
  obs_names <- names(obs)
  mu <- setNames(numeric(length(obs_names)), obs_names)
  for (nm in obs_names) mu[nm] <- .fullplug_mu(object, row, b, nm)
  
  S <- .fullplug_resid_cov(object, obs_names)
  
  ord_names <- obs_names[vapply(obs_names, function(nm) object$outcomes[[nm]]$type == "ordinal", logical(1))]
  gau_names <- obs_names[vapply(obs_names, function(nm) object$outcomes[[nm]]$type == "gaussian", logical(1))]
  
  if (!length(ord_names) && length(gau_names)) {
    yv <- unlist(obs[gau_names], use.names = FALSE)
    return(mvtnorm::dmvnorm(yv,
                            mean = mu[gau_names],
                            sigma = S[gau_names, gau_names, drop = FALSE],
                            log = TRUE))
  }
  
  if (length(ord_names) && !length(gau_names)) {
    lower <- upper <- numeric(length(ord_names))
    for (k in seq_along(ord_names)) {
      bb <- .fullplug_ord_bounds(obs[[ord_names[k]]], object$outcomes[[ord_names[k]]]$thresholds)
      lower[k] <- bb[1]
      upper[k] <- bb[2]
    }
    
    if (length(ord_names) == 1L) {
      p <- pnorm(upper - mu[ord_names]) - pnorm(lower - mu[ord_names])
      return(log(pmax(p, 1e-12)))
    }
    
    pr <- mvtnorm::pmvnorm(lower = lower,
                           upper = upper,
                           mean = mu[ord_names],
                           sigma = S[ord_names, ord_names, drop = FALSE],
                           algorithm = pmvnorm_control)[1]
    return(log(pmax(pr, 1e-12)))
  }
  
  y_c <- unlist(obs[gau_names], use.names = FALSE)
  mu_c <- mu[gau_names]
  Sig_cc <- S[gau_names, gau_names, drop = FALSE]
  
  Sig_oo <- S[ord_names, ord_names, drop = FALSE]
  Sig_oc <- S[ord_names, gau_names, drop = FALSE]
  Sig_co <- S[gau_names, ord_names, drop = FALSE]
  
  Sig_cc_inv <- solve(Sig_cc)
  mu_o_cond <- as.numeric(mu[ord_names] + Sig_oc %*% Sig_cc_inv %*% (y_c - mu_c))
  Sig_o_cond <- Sig_oo - Sig_oc %*% Sig_cc_inv %*% Sig_co
  Sig_o_cond <- 0.5 * (Sig_o_cond + t(Sig_o_cond))
  
  ee <- eigen(Sig_o_cond, symmetric = TRUE)
  ee$values[ee$values < 1e-8] <- 1e-8
  Sig_o_cond <- ee$vectors %*% diag(ee$values, nrow = length(ee$values)) %*% t(ee$vectors)
  
  lower <- upper <- numeric(length(ord_names))
  for (k in seq_along(ord_names)) {
    bb <- .fullplug_ord_bounds(obs[[ord_names[k]]], object$outcomes[[ord_names[k]]]$thresholds)
    lower[k] <- bb[1]
    upper[k] <- bb[2]
  }
  
  log_fc <- mvtnorm::dmvnorm(y_c, mean = mu_c, sigma = Sig_cc, log = TRUE)
  
  if (length(ord_names) == 1L) {
    z1 <- sqrt(Sig_o_cond[1, 1])
    p <- pnorm((upper - mu_o_cond) / z1) - pnorm((lower - mu_o_cond) / z1)
    return(log_fc + log(pmax(p, 1e-12)))
  }
  
  pr <- mvtnorm::pmvnorm(lower = lower,
                         upper = upper,
                         mean = mu_o_cond,
                         sigma = Sig_o_cond,
                         algorithm = pmvnorm_control)[1]
  
  log_fc + log(pmax(pr, 1e-12))
}

.fullplug_subject_long_logdens <- function(object,
                                           df_i,
                                           t_cut,
                                           timeVar = "visit_m",
                                           pmvnorm_control = mvtnorm::GenzBretz(maxpts = 25000,
                                                                                abseps = 1e-04,
                                                                                releps = 0),
                                           b) {
  dd <- df_i[df_i[[timeVar]] <= t_cut, , drop = FALSE]
  if (!nrow(dd)) return(0)
  sum(vapply(seq_len(nrow(dd)), function(r) {
    .fullplug_visit_logdens(object, dd[r, , drop = FALSE], b, pmvnorm_control = pmvnorm_control)
  }, numeric(1)))
}

.fullplug_logpost_b <- function(b_vec,
                                object,
                                df_i,
                                t_cut,
                                timeVar = "visit_m",
                                pmvnorm_control = mvtnorm::GenzBretz(maxpts = 25000,
                                                                     abseps = 1e-04,
                                                                     releps = 0)) {
  b <- setNames(as.numeric(b_vec), rownames(object$D))
  row0 <- df_i[1, , drop = FALSE]
  
  ll_long <- .fullplug_subject_long_logdens(object,
                                            df_i = df_i,
                                            t_cut = t_cut,
                                            timeVar = timeVar,
                                            pmvnorm_control = pmvnorm_control,
                                            b = b)
  
  ll_surv <- .fullplug_logS(object, row0, b, t_cut)
  ll_prior <- mvtnorm::dmvnorm(b_vec, mean = rep(0, length(b_vec)), sigma = object$D, log = TRUE)
  ll_long + ll_surv + ll_prior
}

.safe_hinv <- function(H) {
  H <- 0.5 * (H + t(H))
  ee <- eigen(H, symmetric = TRUE)
  ee$values[ee$values < 1e-6] <- 1e-6
  H2 <- ee$vectors %*% diag(ee$values, nrow = length(ee$values)) %*% t(ee$vectors)
  solve(H2)
}

repair_thresholds <- function(x, eps = 1e-4) {
  x <- sort(as.numeric(x))
  if (length(x) <= 1L) return(x)
  for (j in 2:length(x)) {
    if (x[j] <= x[j - 1] + eps) x[j] <- x[j - 1] + eps
  }
  x
}

project_spd <- function(S, eps = 1e-8) {
  S <- 0.5 * (S + t(S))
  ee <- eigen(S, symmetric = TRUE)
  ee$values[ee$values < eps] <- eps
  out <- ee$vectors %*% diag(ee$values, nrow = length(ee$values)) %*% t(ee$vectors)
  dimnames(out) <- dimnames(S)
  out
}

project_corr <- function(R, eps = 1e-8) {
  R <- as.matrix(R)
  R <- 0.5 * (R + t(R))
  diag(R) <- 1
  off <- row(R) != col(R)
  R[off] <- pmin(pmax(R[off], -0.999), 0.999)
  out <- as.matrix(Matrix::nearPD(R, corr = TRUE, keepDiag = TRUE)$mat)
  out <- 0.5 * (out + t(out))
  diag(out) <- 1
  out
}

theta_to_D <- function(theta, D_names) {
  D <- matrix(0, length(D_names), length(D_names), dimnames = list(D_names, D_names))
  
  for (nm in D_names) {
    vn <- paste0("var_", nm)
    if (!vn %in% names(theta)) stop("manca ", vn, " nel draw di theta.")
    D[nm, nm] <- theta[[vn]]
  }
  
  cov_names <- grep("^cov_", names(theta), value = TRUE)
  for (cn in cov_names) {
    ss <- sub("^cov_", "", cn)
    sp <- strsplit(ss, "__", fixed = TRUE)[[1]]
    if (length(sp) != 2L) next
    a <- sp[1]
    b <- sp[2]
    if (a %in% D_names && b %in% D_names) {
      D[a, b] <- theta[[cn]]
      D[b, a] <- theta[[cn]]
    }
  }
  
  D
}

theta_to_R <- function(theta, outcome_names) {
  R <- diag(1, length(outcome_names))
  dimnames(R) <- list(outcome_names, outcome_names)
  
  nm_atanh <- grep("^atanh_rho_resid__", names(theta), value = TRUE)
  for (nm in nm_atanh) {
    ss <- sub("^atanh_rho_resid__+", "", nm)
    sp <- strsplit(ss, "__", fixed = TRUE)[[1]]
    if (length(sp) == 2L && all(sp %in% outcome_names)) {
      rr <- tanh(theta[[nm]])
      R[sp[1], sp[2]] <- rr
      R[sp[2], sp[1]] <- rr
    }
  }
  
  nm_corr <- grep("^corr_resid__", names(theta), value = TRUE)
  for (nm in nm_corr) {
    ss <- sub("^corr_resid__+", "", nm)
    sp <- strsplit(ss, "__", fixed = TRUE)[[1]]
    if (length(sp) == 2L && all(sp %in% outcome_names)) {
      R[sp[1], sp[2]] <- theta[[nm]]
      R[sp[2], sp[1]] <- theta[[nm]]
    }
  }
  
  R
}

pair_raw_to_psi_exact <- function(theta_pair_pref) {
  stopifnot(!is.null(names(theta_pair_pref)))
  
  chol_names <- grep("::(log_chol_var_|chol_cov_)", names(theta_pair_pref), value = TRUE)
  nonchol_names <- setdiff(names(theta_pair_pref), chol_names)
  
  out <- theta_pair_pref[nonchol_names]
  
  if (length(chol_names)) {
    tag <- sub("::.*$", "", chol_names[1])
    
    vc <- chol_named_to_varcov(
      chol_vec = theta_pair_pref[chol_names],
      chol_names_prefixed = chol_names
    )
    
    vc <- setNames(as.numeric(vc), paste0(tag, "::", names(vc)))
    out <- c(out, vc)
  }
  
  out
}

.pair_rel_diff <- function(A, B, eps = 1e-12) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  num <- sqrt(sum((A - B)^2, na.rm = TRUE))
  den <- max(sqrt(sum(B^2, na.rm = TRUE)), eps)
  num / den
}

.pair_fill_D_in_psi <- function(psi, D, D_names) {
  D <- D[D_names, D_names, drop = FALSE]
  
  for (nm in D_names) {
    vn <- paste0("var_", nm)
    if (vn %in% names(psi)) psi[[vn]] <- D[nm, nm]
  }
  
  cov_names <- grep("^cov_", names(psi), value = TRUE)
  for (cn in cov_names) {
    ss <- sub("^cov_", "", cn)
    sp <- strsplit(ss, "__", fixed = TRUE)[[1]]
    if (length(sp) != 2L) next
    a <- sp[1]
    b <- sp[2]
    if (a %in% D_names && b %in% D_names) {
      psi[[cn]] <- D[a, b]
    }
  }
  
  psi
}

.pair_fill_R_in_psi <- function(psi, R, outcome_names) {
  R <- R[outcome_names, outcome_names, drop = FALSE]
  
  nm_atanh <- grep("^atanh_rho_resid__", names(psi), value = TRUE)
  for (nm in nm_atanh) {
    ss <- sub("^atanh_rho_resid__+", "", nm)
    sp <- strsplit(ss, "__", fixed = TRUE)[[1]]
    if (length(sp) == 2L && all(sp %in% outcome_names)) {
      rr <- R[sp[1], sp[2]]
      psi[[nm]] <- atanh(pmin(pmax(rr, -0.999999), 0.999999))
    }
  }
  
  nm_corr <- grep("^corr_resid__", names(psi), value = TRUE)
  for (nm in nm_corr) {
    ss <- sub("^corr_resid__+", "", nm)
    sp <- strsplit(ss, "__", fixed = TRUE)[[1]]
    if (length(sp) == 2L && all(sp %in% outcome_names)) {
      psi[[nm]] <- R[sp[1], sp[2]]
    }
  }
  
  psi
}

.pair_correct_draw <- function(object, psi_draw) {
  out <- list(
    psi_use = psi_draw,
    corrected_D = FALSE,
    corrected_R = FALSE,
    D_max_abs = NA_real_,
    D_rel = NA_real_,
    R_max_abs = NA_real_,
    R_rel = NA_real_
  )
  
  if (any(!is.finite(psi_draw))) {
    stop("psi_draw contiene valori non finiti.")
  }
  
  D_raw <- theta_to_D(psi_draw, rownames(object$D))
  if (any(!is.finite(D_raw))) stop("D_raw contiene valori non finiti.")
  
  D_fix <- project_spd(D_raw)
  out$D_max_abs <- max(abs(D_fix - D_raw), na.rm = TRUE)
  out$D_rel <- .pair_rel_diff(D_fix, D_raw)
  
  if (out$D_max_abs > 0) {
    out$corrected_D <- TRUE
    out$psi_use <- .pair_fill_D_in_psi(out$psi_use, D_fix, rownames(object$D))
  }
  
  R_raw <- theta_to_R(out$psi_use, names(object$outcomes))
  if (any(!is.finite(R_raw))) stop("R_raw contiene valori non finiti.")
  
  R_fix <- project_corr(R_raw)
  out$R_max_abs <- max(abs(R_fix - R_raw), na.rm = TRUE)
  out$R_rel <- .pair_rel_diff(R_fix, R_raw)
  
  if (out$R_max_abs > 0) {
    out$corrected_R <- TRUE
    out$psi_use <- .pair_fill_R_in_psi(out$psi_use, R_fix, names(object$outcomes))
  }
  
  out
}

psi_to_D <- function(psi_draw, D_names) { 
  D <- matrix(0, length(D_names), length(D_names), dimnames = list(D_names, D_names))
  for (nm in D_names) { 
    vn <- paste0("var_", nm) 
    if (!vn %in% names(psi_draw)) stop("manca ", vn)
    D[nm, nm] <- psi_draw[[vn]] 
  } 
  cov_names <- grep("^cov_", names(psi_draw), value = TRUE) 
  for (cn in cov_names) { 
    ss <- sub("^cov_", "", cn) 
    sp <- strsplit(ss, "__", fixed = TRUE)[[1]] 
    if (length(sp) != 2L) next 
    a <- sp[1]
    b <- sp[2] 
    if (a %in% D_names && b %in% D_names) { 
      D[a, b] <- psi_draw[[cn]] 
      D[b, a] <- psi_draw[[cn]] } 
  } 
  local({Df <- project_spd(D); 
    if (max(abs(D - Df), na.rm = TRUE) > 1e-10) message("iter ", getOption("fullplug_iter"), ": corrected D"); Df}) 
} 

psi_to_R <- function(psi_draw, outcome_names) { 
  R <- diag(1, length(outcome_names)) 
  dimnames(R) <- list(outcome_names, outcome_names) 
  nm_atanh <- grep("^atanh_rho_resid__", names(psi_draw), value = TRUE) 
  for (nm in nm_atanh) { 
    ss <- sub("^atanh_rho_resid__+", "", nm) 
    sp <- strsplit(ss, "__", fixed = TRUE)[[1]] 
    if (length(sp) == 2L && all(sp %in% outcome_names)) { 
      rr <- tanh(psi_draw[[nm]]) 
      R[sp[1], sp[2]] <- rr
      R[sp[2], sp[1]] <- rr
      } 
  } 
  nm_corr <- grep("^corr_resid__", names(psi_draw), value = TRUE) 
  for (nm in nm_corr) { 
    ss <- sub("^corr_resid__+", "", nm) 
    sp <- strsplit(ss, "__", fixed = TRUE)[[1]] 
    if (length(sp) == 2L && all(sp %in% outcome_names)) {
      R[sp[1], sp[2]] <- psi_draw[[nm]] 
      R[sp[2], sp[1]] <- psi_draw[[nm]] 
    } 
  } 
  local({Rf <- project_corr(R); if (max(abs(R - Rf), na.rm = TRUE) > 1e-10) message("iter ", getOption("fullplug_iter"), ": corrected R"); Rf}) 
  
  }

make_theta_sampler_pairwise_raw <- function(ex_pref_for_jk,
                                            V_theta_raw,
                                            canonicalizer = canonicalizer_all,
                                            verbose = TRUE) {
  theta_blocks <- unname(lapply(ex_pref_for_jk, function(ex) {
    th <- ex$pars[colnames(ex$H)]
    stopifnot(identical(names(th), colnames(ex$H)))
    th
  }))
  
  theta_raw_mean <- unlist(theta_blocks, use.names = TRUE)
  
  V_theta_raw <- as.matrix(V_theta_raw)
  if (is.null(colnames(V_theta_raw)) || is.null(rownames(V_theta_raw))) {
    stop("V_theta_raw deve avere rownames e colnames.")
  }
  
  if (!setequal(rownames(V_theta_raw), colnames(V_theta_raw))) {
    stop("rownames(V_theta_raw) e colnames(V_theta_raw) non coincidono come insieme.")
  }
  
  V_theta_raw <- V_theta_raw[colnames(V_theta_raw), colnames(V_theta_raw), drop = FALSE]
  
  if (!setequal(names(theta_raw_mean), colnames(V_theta_raw))) {
    miss1 <- setdiff(names(theta_raw_mean), colnames(V_theta_raw))
    miss2 <- setdiff(colnames(V_theta_raw), names(theta_raw_mean))
    stop(
      "nomi incompatibili tra theta_raw_mean e V_theta_raw.\n",
      "in theta_raw_mean ma non in V_theta_raw: ", paste(head(miss1, 20), collapse = ", "), "\n",
      "in V_theta_raw ma non in theta_raw_mean: ", paste(head(miss2, 20), collapse = ", ")
    )
  }
  
  theta_raw_mean <- theta_raw_mean[colnames(V_theta_raw)]
  stopifnot(identical(names(theta_raw_mean), colnames(V_theta_raw)))
  
  psi_blocks_mean <- unname(lapply(theta_blocks, pair_raw_to_psi_exact))
  psi_all_mean <- unlist(psi_blocks_mean, use.names = TRUE)
  
  A_all <- build_A(
    est_vec = psi_all_mean,
    canonicalizer = canonicalizer,
    make_unique_cols = FALSE
  )
  
  if (!setequal(colnames(A_all), names(psi_all_mean))) {
    miss1 <- setdiff(colnames(A_all), names(psi_all_mean))
    miss2 <- setdiff(names(psi_all_mean), colnames(A_all))
    stop(
      "nomi incompatibili tra A_all e psi_all_mean.\n",
      "in A_all ma non in psi_all_mean: ", paste(head(miss1, 20), collapse = ", "), "\n",
      "in psi_all_mean ma non in A_all: ", paste(head(miss2, 20), collapse = ", ")
    )
  }
  
  block_lengths <- vapply(theta_blocks, length, integer(1))
  block_ends <- cumsum(block_lengths)
  block_starts <- c(1L, head(block_ends, -1L) + 1L)
  
  draw_to_global <- function(theta_raw_draw) {
    theta_raw_draw <- theta_raw_draw[colnames(V_theta_raw)]
    stopifnot(identical(names(theta_raw_draw), colnames(V_theta_raw)))
    
    raw_blocks_draw <- vector("list", length(theta_blocks))
    for (i in seq_along(theta_blocks)) {
      idx <- block_starts[i]:block_ends[i]
      raw_blocks_draw[[i]] <- theta_raw_draw[idx]
    }
    
    psi_blocks_draw <- unname(lapply(raw_blocks_draw, pair_raw_to_psi_exact))
    psi_all_draw <- unlist(psi_blocks_draw, use.names = TRUE)
    
    if (!setequal(colnames(A_all), names(psi_all_draw))) {
      miss1 <- setdiff(colnames(A_all), names(psi_all_draw))
      miss2 <- setdiff(names(psi_all_draw), colnames(A_all))
      stop(
        "nomi incompatibili tra A_all e psi_all_draw.\n",
        "in A_all ma non in psi_all_draw: ", paste(head(miss1, 20), collapse = ", "), "\n",
        "in psi_all_draw ma non in A_all: ", paste(head(miss2, 20), collapse = ", ")
      )
    }
    
    psi_all_draw <- psi_all_draw[colnames(A_all)]
    psi_global_draw <- as.numeric(A_all %*% psi_all_draw)
    names(psi_global_draw) <- rownames(A_all)
    psi_global_draw
  }
  
  psi_global_mean <- as.numeric(A_all %*% psi_all_mean)
  names(psi_global_mean) <- rownames(A_all)
  
  list(
    type = "pairwise_raw",
    theta_mean = theta_raw_mean,
    V_theta = ensure_sym_pd(V_theta_raw),
    A_all = A_all,
    psi_global_mean = psi_global_mean,
    draw_to_global = draw_to_global,
    verbose = verbose
  )
}

apply_psi_draw_to_object <- function(object, psi_draw) {
  obj <- object 
  for (nm in names(obj$outcomes)) { 
    ok <- obj$outcomes[[nm]]
    pat_beta <- paste0("^fixed_", nm, "_") 
    ii_beta <- grep(pat_beta, names(psi_draw)) 
    if (length(ii_beta)) { 
      vv <- psi_draw[ii_beta] 
      nvv <- sub(pat_beta, "", names(vv)) 
      nvv[nvv == "Intercept"] <- "(Intercept)" 
      names(vv) <- nvv 
      obj$outcomes[[nm]]$beta <- vv 
    } 
    if (ok$type == "ordinal") { 
      pat_thr <- paste0("^threshold_", nm, "_cut:") 
      ii_thr <- grep(pat_thr, names(psi_draw)) 
      if (length(ii_thr)) { 
        vv <- psi_draw[ii_thr] 
        idx <- as.integer(sub(pat_thr, "", names(vv))) 
        vv <- vv[order(idx)] 
        obj$outcomes[[nm]]$thresholds <- repair_thresholds(vv)
        } 
    } 
    if (ok$type == "gaussian") {
        nm_logsig <- paste0("log_sigma_", nm) 
        nm_sig <- paste0("sigma_", nm) 
        if (nm_logsig %in% names(psi_draw)) { 
          obj$outcomes[[nm]]$sigma <- exp(psi_draw[[nm_logsig]]) 
        } else if (nm_sig %in% names(psi_draw)) {
            obj$outcomes[[nm]]$sigma <- max(psi_draw[[nm_sig]], 1e-8)
        } 
    } 
  } 
  ii_surv <- grep("^fixed_surv_", names(psi_draw))
  if (length(ii_surv)) { 
    vv <- psi_draw[ii_surv]
    names(vv) <- sub("^fixed_surv_", "", names(vv)) 
    obj$survival$beta <- vv 
    } 
  if ("log_lambda" %in% names(psi_draw)) { 
      obj$survival$lambda <- exp(psi_draw[["log_lambda"]]) 
  } else if ("lambda" %in% names(psi_draw)) {
        obj$survival$lambda <- max(psi_draw[["lambda"]], 1e-8) 
  } 
  if ("log_rho" %in% names(psi_draw)) {
    obj$survival$rho <- exp(psi_draw[["log_rho"]]) 
  } else if ("rho" %in% names(psi_draw)) {
      obj$survival$rho <- max(psi_draw[["rho"]], 1e-8) 
  } 
  obj$D <- psi_to_D(psi_draw, rownames(obj$D))
  obj$resid_corr <- project_corr(psi_to_R(psi_draw, names(obj$outcomes))) 
  obj 
} 

draw_theta_global <- function(theta_sampler) { 
  th <- as.numeric(MASS::mvrnorm( 1, mu = theta_sampler$theta_mean, Sigma = theta_sampler$V_theta )) 
  names(th) <- names(theta_sampler$theta_mean) 
  th }


draw_theta_pairwise_raw <- function(theta_sampler) {
  th <- as.numeric(MASS::mvrnorm(
    1,
    mu = theta_sampler$theta_mean,
    Sigma = theta_sampler$V_theta
  ))
  names(th) <- names(theta_sampler$theta_mean)
  th
}

draw_theta_fullplug <- function(object, theta_source = c("global", "pairwise_raw")) {
  theta_source <- match.arg(theta_source)
  options(fullplug_iter = getOption("fullplug_iter", 0L) + 1L)
  
  sampler <- object$theta_sampler
  if (is.null(sampler)) stop("object$theta_sampler mancante.")
  
  if (theta_source == "global") {
    if (sampler$type != "global") stop("theta_sampler non è di tipo 'global'.")
    
    theta_draw <- as.numeric(MASS::mvrnorm(
      1,
      mu = sampler$theta_mean,
      Sigma = ensure_sym_pd(sampler$V_theta)
    ))
    names(theta_draw) <- names(sampler$theta_mean)
    
    return(apply_psi_draw_to_object(object, theta_draw))
  }
  
  if (theta_source == "pairwise_raw") {
    if (sampler$type != "pairwise_raw") stop("theta_sampler non è di tipo 'pairwise_raw'.")
    
    draw_id <- getOption("fullplug_iter")
    
    theta_raw_draw <- as.numeric(MASS::mvrnorm(
      1,
      mu = sampler$theta_mean,
      Sigma = ensure_sym_pd(sampler$V_theta)
    ))
    names(theta_raw_draw) <- names(sampler$theta_mean)
    
    psi_global_draw <- sampler$draw_to_global(theta_raw_draw)
    chk <- .pair_correct_draw(object, psi_global_draw)
    
    if (isTRUE(sampler$verbose)) {
      cat(sprintf(
        "draw %d ACCEPTED corrected_D=%s corrected_R=%s D[max=%.3e rel=%.3e] R[max=%.3e rel=%.3e]\n",
        draw_id,
        ifelse(chk$corrected_D, "yes", "no"),
        ifelse(chk$corrected_R, "yes", "no"),
        chk$D_max_abs, chk$D_rel,
        chk$R_max_abs, chk$R_rel
      ))
    }
    
    return(apply_psi_draw_to_object(object, chk$psi_use))
  }
  
  stop("theta_source non riconosciuto.")
}

sanitize_proposal_Sigma <- function(S, eps = 1e-6, max_var = 10) {
  S <- as.matrix(S)
  storage.mode(S) <- "double"
  S <- 0.5 * (S + t(S))
  
  ee <- eigen(S, symmetric = TRUE)
  vals <- ee$values
  vals[!is.finite(vals)] <- eps
  vals[vals < eps] <- eps
  vals[vals > max_var] <- max_var
  
  S2 <- ee$vectors %*% diag(vals, nrow = length(vals)) %*% t(ee$vectors)
  S2 <- 0.5 * (S2 + t(S2))
  dimnames(S2) <- dimnames(S)
  S2
}


make_fake_hist_wide <- function(id,
                                ord_items,
                                Y_ord,
                                covariates_row,
                                cont_name = NULL,
                                Y_cont = NULL,
                                times = c(0, 6, 12),
                                times_in = c("months", "days"),
                                age_mean = NULL,
                                age_sd = NULL,
                                age_already_scaled = TRUE) {
  times_in <- match.arg(times_in)
  
  if (!is.data.frame(covariates_row) || nrow(covariates_row) != 1L) {
    stop("covariates_row deve essere un data.frame con una sola riga.")
  }
  
  if (is.null(colnames(Y_ord))) {
    stop("Y_ord deve avere colnames = ord_items.")
  }
  
  ord_items <- as.character(ord_items)
  
  if (!all(ord_items %in% colnames(Y_ord))) {
    miss <- setdiff(ord_items, colnames(Y_ord))
    stop("mancano queste colonne ordinali in Y_ord: ", paste(miss, collapse = ", "))
  }
  
  Y_ord <- as.data.frame(Y_ord[, ord_items, drop = FALSE])
  
  n_vis <- nrow(Y_ord)
  if (length(times) != n_vis) {
    stop("length(times) deve coincidere con nrow(Y_ord).")
  }
  
  if (!is.null(cont_name)) {
    if (is.null(Y_cont)) stop("se cont_name non è NULL devi passare Y_cont.")
    if (length(Y_cont) != n_vis) stop("length(Y_cont) deve coincidere con nrow(Y_ord).")
  }
  
  visit_m <- if (times_in == "months") {
    as.numeric(times)
  } else {
    as.numeric(times) / 30.4375
  }
  
  visit_idx <- rep(NA_integer_, n_vis)
  
  map_std <- function(x) {
    if (length(age_mean) == 1L && length(age_sd) == 1L && !age_already_scaled) {
      return((x - age_mean) / age_sd)
    }
    x
  }
  
  out <- covariates_row[rep(1, n_vis), , drop = FALSE]
  rownames(out) <- NULL
  
  out$ID <- id
  out$visit_m <- visit_m
  out$visit_idx <- visit_idx
  
  out$I0 <- as.integer(abs(out$visit_m - 0) < 1e-8)
  out$I6 <- as.integer(abs(out$visit_m - 6) < 1e-8)
  out$I12 <- as.integer(abs(out$visit_m - 12) < 1e-8)
  out$Ipost <- as.integer(out$I6 + out$I12 > 0)
  
  if (!"sex01" %in% names(out) && "sex" %in% names(out)) {
    out$sex01 <- as.integer(out$sex == "M")
  }
  
  if (!"diag_oa" %in% names(out) && "diagnosis" %in% names(out)) {
    out$diag_oa <- as.integer(out$diagnosis != "osteoarthritis_knee")
  }
  
  if (!"prost_uni" %in% names(out) && "prosthesis_type" %in% names(out)) {
    out$prost_uni <- as.integer(out$prosthesis_type == "unicompartimental")
  }
  
  if (!"waiting" %in% names(out)) {
    out$waiting <- 0
  }
  
  if (!"loghosp" %in% names(out)) {
    if ("hospitalizationdays" %in% names(out)) {
      out$loghosp <- log1p(out$hospitalizationdays)
    } else {
      out$loghosp <- 0
    }
  }
  
  if (!"prost_uni" %in% names(out)) {
    out$prost_uni <- 0L
  }
  
  if ("age" %in% names(out)) {
    out$age <- map_std(as.numeric(out$age))
  }
  
  out$wI6 <- out$waiting * out$I6
  out$wI12 <- out$waiting * out$I12
  out$pI6 <- out$prost_uni * out$I6
  out$pI12 <- out$prost_uni * out$I12
  out$loghosp_post <- out$loghosp * out$Ipost
  
  for (nm in ord_items) {
    out[[nm]] <- Y_ord[[nm]]
  }
  
  if (!is.null(cont_name)) {
    out[[cont_name]] <- as.numeric(Y_cont)
  }
  
  out
}

make_fake_hist_wide_baseline <- function(id,
                                         items,
                                         y_row,
                                         covariates_row,
                                         cont_name = NULL,
                                         y_cont = NULL,
                                         age_mean = NULL,
                                         age_sd = NULL,
                                         age_already_scaled = TRUE) {
  make_fake_hist_wide(
    id = id,
    ord_items = items,
    Y_ord = y_row,
    covariates_row = covariates_row,
    cont_name = cont_name,
    Y_cont = y_cont,
    times = 0,
    times_in = "months",
    age_mean = age_mean,
    age_sd = age_sd,
    age_already_scaled = age_already_scaled
  )
}


ensure_sym_pd <- function(S, eps = 1e-8) {
  S <- as.matrix(S)
  storage.mode(S) <- "double"
  
  if (nrow(S) != ncol(S)) stop("S deve essere quadrata.")
  
  S <- 0.5 * (S + t(S))
  
  ee <- eigen(S, symmetric = TRUE)
  vals <- ee$values
  vals[!is.finite(vals)] <- eps
  vals[vals < eps] <- eps
  
  S2 <- ee$vectors %*% diag(vals, nrow = length(vals)) %*% t(ee$vectors)
  S2 <- 0.5 * (S2 + t(S2))
  
  dimnames(S2) <- dimnames(S)
  S2
}
.safe_hinv <- function(H, eps = 1e-8) {
  H <- as.matrix(H)
  storage.mode(H) <- "double"
  H <- 0.5 * (H + t(H))
  
  ee <- eigen(H, symmetric = TRUE)
  vals <- ee$values
  vals[!is.finite(vals)] <- eps
  vals[vals < eps] <- eps
  
  Hpd <- ee$vectors %*% diag(vals, nrow = length(vals)) %*% t(ee$vectors)
  Hpd <- 0.5 * (Hpd + t(Hpd))
  
  solve(Hpd)
}




cd <- function(x, f, ..., eps = 1e-04) {
  x <- as.numeric(x)
  g <- numeric(length(x))
  for (j in seq_along(x)) {
    x1 <- x
    x2 <- x
    x1[j] <- x1[j] + eps
    x2[j] <- x2[j] - eps
    g[j] <- (f(x1, ...) - f(x2, ...)) / (2 * eps)
  }
  g
}

rmvt <- function (n, mu, Sigma, df) {
  p <- length(mu)
  if (is.list(Sigma)) {
    ev <- Sigma$values
    evec <- Sigma$vectors
  } else {
    ed <- eigen(Sigma, symmetric = TRUE)
    ev <- ed$values
    evec <- ed$vectors
  }
  X <- drop(mu) + tcrossprod(
    evec * rep(sqrt(pmax(ev, 0)), each = p),
    matrix(rnorm(n * p), n)
  ) / rep(sqrt(rchisq(n, df) / df), each = p)
  if (n == 1L) drop(X) else t.default(X)
}

dmvt <- function (x, mu, Sigma, df, log = FALSE) {
  if (!is.numeric(x))
    stop("'x' must be a numeric matrix or vector")
  if (!is.matrix(x))
    x <- rbind(x)
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p)) || ncol(x) != p)
    stop("incompatible arguments")
  ed <- eigen(Sigma, symmetric = TRUE)
  ev <- ed$values
  if (!all(ev >= -1e-06 * abs(ev[1])))
    stop("'Sigma' is not positive definite")
  ss <- x - rep(mu, each = nrow(x))
  inv.Sigma <- ed$vectors %*% (t(ed$vectors) / ev)
  quad <- rowSums((ss %*% inv.Sigma) * ss) / df
  fact <- lgamma((df + p)/2) - lgamma(df/2) -
    0.5 * (p * (log(pi) + log(df)) + sum(log(ev)))
  if (log)
    fact - 0.5 * (df + p) * log(1 + quad)
  else
    exp(fact) * ((1 + quad)^(- (df + p)/2))
}

.find_mode_b_fullplug <- function(object,
                                  df_i,
                                  t_cut,
                                  timeVar = "visit_m",
                                  scale = 1.6,
                                  pmvnorm_control = mvtnorm::GenzBretz(
                                    maxpts = 25000,
                                    abseps = 1e-04,
                                    releps = 0
                                  )) {
  q <- nrow(object$D)
  
  ff <- function(bb, object, df_i, t_cut, timeVar, pmvnorm_control) {
    - .fullplug_logpost_b(
      bb,
      object = object,
      df_i = df_i,
      t_cut = t_cut,
      timeVar = timeVar,
      pmvnorm_control = pmvnorm_control
    )
  }
  
  opt <- try(
    optim(
      rep(0, q),
      ff,
      object = object,
      df_i = df_i,
      t_cut = t_cut,
      timeVar = timeVar,
      pmvnorm_control = pmvnorm_control,
      method = "BFGS",
      hessian = TRUE
    ),
    TRUE
  )
  
  if (inherits(opt, "try-error")) {
    gg <- function(bb, object, df_i, t_cut, timeVar, pmvnorm_control) {
      cd(
        bb,
        ff,
        object = object,
        df_i = df_i,
        t_cut = t_cut,
        timeVar = timeVar,
        pmvnorm_control = pmvnorm_control
      )
    }
    
    opt <- optim(
      rep(0, q),
      ff,
      gg,
      object = object,
      df_i = df_i,
      t_cut = t_cut,
      timeVar = timeVar,
      pmvnorm_control = pmvnorm_control,
      method = "BFGS",
      hessian = TRUE
    )
  }
  
  list(
    mode = opt$par,
    var = scale * solve(opt$hessian),
    hessian = opt$hessian,
    opt = opt
  )
}
survfitJM_fullplug <- function(object,
                               newdata,
                               idVar = "ID",
                               timeVar = "visit_m",
                               simulate = TRUE,
                               survTimes = NULL,
                               last.time = NULL,
                               M = 200,
                               CI.levels = c(0.025, 0.975),
                               scale = 1.6,
                               theta_source = c("global", "pairwise_raw"),
                               seed = NULL,
                               repair_fake = TRUE,
                               silent_fake_repair = TRUE,
                               pmvnorm_control = mvtnorm::GenzBretz(
                                 maxpts = 25000,
                                 abseps = 1e-04,
                                 releps = 0
                               ),
                               ...) {
  theta_source <- match.arg(theta_source)
  
  if (!is.list(object))
    stop("Use only with a valid 'fullplug' object.\n")
  
  if (!is.data.frame(newdata) || nrow(newdata) == 0)
    stop("'newdata' must be a data.frame with more than one rows.\n")
  
  if (is.null(newdata[[idVar]]))
    stop("'idVar' not in 'newdata'.\n")
  
  if (is.null(object$theta_sampler))
    stop("'object$theta_sampler' is missing.\n")
  
  if (is.null(object$D) || !is.matrix(object$D))
    stop("'object$D' must be a matrix.\n")
  
  if (is.null(survTimes) || !is.numeric(survTimes))
    stop("'survTimes' must be numeric.\n")
  
  if (!is.null(seed))
    set.seed(seed)
  
  newdata <- prep_newdata_fullplug(
    newdata,
    object = object,
    idVar = idVar,
    timeVar = timeVar
  )
  
  id_chr0 <- as.character(newdata[[idVar]])
  id <- id. <- match(id_chr0, unique(id_chr0))
  
  obs.times <- split(newdata[[timeVar]], id)
  
  last.time <- if (is.null(last.time)) {
    tapply(newdata[[timeVar]], id., tail, n = 1)
  } else if (is.character(last.time) && length(last.time) == 1L) {
    tapply(newdata[[last.time]], id., tail, n = 1)
  } else if (is.numeric(last.time) && length(last.time) == length(unique(id.))) {
    last.time
  } else {
    stop("\nnot appropriate value for 'last.time' argument.")
  }
  
  times.to.pred <- lapply(last.time, function(t) survTimes[survTimes > t])
  
  
  fake_pattern <- "^FAKE"
  repair_tol <- 1e-08
  eig_floor_abs <- 1e-06
  
  split_data <- split(newdata, id.)
  ids_chr <- unique(as.character(newdata[[idVar]]))
  n_id <- length(split_data)
  q <- nrow(object$D)
  
  modes.b <- matrix(0, n_id, q)
  colnames(modes.b) <- rownames(object$D)
  rownames(modes.b) <- ids_chr
  
  Vars.b <- vector("list", n_id)
  
  proposal.repaired <- rep(FALSE, n_id)
  proposal.min_eig_raw <- rep(NA_real_, n_id)
  proposal.min_eig_fix <- rep(NA_real_, n_id)
  
  for (i in seq_len(n_id)) {
    df_i <- split_data[[i]]
    t_cut_i <- last.time[i]
    
    fit_i <- .find_mode_b_fullplug(
      object = object,
      df_i = df_i,
      t_cut = t_cut_i,
      timeVar = timeVar,
      scale = scale
    )
    
    chk <- .check_H_for_proposal(fit_i$hessian, tol = repair_tol)
    proposal.min_eig_raw[i] <- chk$min_eig
    
    is_fake_i <- .is_fake_subject(df_i, idVar = idVar, fake_pattern = fake_pattern)
    
    if (isTRUE(repair_fake) && isTRUE(is_fake_i) && !isTRUE(chk$ok)) {
      repH <- .repair_H_to_Sigma(
        H = fit_i$hessian,
        scale = 1.0,
        eig_floor_abs = 1e-06,
        floor_quantile = 0.25,
        floor_mult = 1,
        max_var = 0.05
      )
      
      fit_i$var <- repH$Sigma
      proposal.repaired[i] <- TRUE
      proposal.min_eig_fix[i] <- repH$min_eig_fix
      
      if (!isTRUE(silent_fake_repair)) {
        warning(
          "proposal repair applied for subject ",
          ids_chr[i],
          " (raw min eigenvalue of Hessian = ",
          signif(chk$min_eig, 6),
          ")."
        )
      }
    } else {
      fit_i$var <- 0.5 * (fit_i$var + t(fit_i$var))
    }
    
    modes.b[i, ] <- fit_i$mode
    Vars.b[[i]] <- sanitize_proposal_Sigma(
      fit_i$var,
      eps = 1e-6,
      max_var = 10
    )
  }
  
  if (!simulate) {
    res <- vector("list", n_id)
    
    for (i in seq_len(n_id)) {
      row0 <- split_data[[i]][1, , drop = FALSE]
      tp <- times.to.pred[[i]]
      b_i <- setNames(modes.b[i, ], rownames(object$D))
      
      S.pred <- numeric(length(tp))
      for (l in seq_along(S.pred)) {
        S.pred[l] <- .fullplug_Sratio(
          object = object,
          row = row0,
          b = b_i,
          t = last.time[i],
          u = tp[l]
        )
      }
      
      res[[i]] <- cbind(times = tp, predSurv = S.pred)
      rownames(res[[i]]) <- seq_along(S.pred)
    }
  } else {
    out <- vector("list", M)
    success.rate <- matrix(FALSE, M, n_id)
    colnames(success.rate) <- ids_chr
    
    b.old <- b.new <- modes.b
    if (n_id == 1L)
      dim(b.old) <- dim(b.new) <- c(1, q)
    
    for (m in seq_len(M)) {
      object_m <- draw_theta_fullplug(object, theta_source = theta_source)
      SS <- vector("list", n_id)
      
      for (i in seq_len(n_id)) {
        mu_i <- as.numeric(modes.b[i, ])
        Sigma_i <- Vars.b[[i]]
        
        proposed.b <- as.numeric(rmvt(1, mu = mu_i, Sigma = Sigma_i, df = 4))
        dmvt.old <- dmvt(as.numeric(b.old[i, ]), mu = mu_i, Sigma = Sigma_i, df = 4, log = TRUE)
        dmvt.proposed <- dmvt(proposed.b, mu = mu_i, Sigma = Sigma_i, df = 4, log = TRUE)
        
        loga <- .fullplug_logpost_b(
          proposed.b,
          object = object_m,
          df_i = split_data[[i]],
          t_cut = last.time[i],
          timeVar = timeVar,
          pmvnorm_control = pmvnorm_control
        ) + dmvt.old -
          .fullplug_logpost_b(
            as.numeric(b.old[i, ]),
            object = object_m,
            df_i = split_data[[i]],
            t_cut = last.time[i],
            timeVar = timeVar,
            pmvnorm_control = pmvnorm_control
          ) - dmvt.proposed
        
        a <- if (is.finite(loga)) min(exp(loga), 1) else NA_real_
        ind <- !is.na(a) && (runif(1) <= a)
        success.rate[m, i] <- ind
        
        if (ind)
          b.new[i, ] <- proposed.b
        
        row0 <- split_data[[i]][1, , drop = FALSE]
        tp <- times.to.pred[[i]]
        b_i <- setNames(as.numeric(b.new[i, ]), rownames(object_m$D))
        
        S.pred <- numeric(length(tp))
        for (l in seq_along(S.pred)) {
          S.pred[l] <- .fullplug_Sratio(
            object = object_m,
            row = row0,
            b = b_i,
            t = last.time[i],
            u = tp[l]
          )
        }
        
        SS[[i]] <- S.pred
      }
      
      b.old <- b.new
      out[[m]] <- SS
    }
    
    res <- vector("list", n_id)
    for (i in seq_len(n_id)) {
      rr <- sapply(out, "[[", i)
      if (!is.matrix(rr))
        rr <- rbind(rr)
      
      res[[i]] <- cbind(
        times = times.to.pred[[i]],
        Mean = rowMeans(rr, na.rm = TRUE),
        Median = apply(rr, 1, median, na.rm = TRUE),
        Lower = apply(rr, 1, quantile, probs = CI.levels[1], na.rm = TRUE),
        Upper = apply(rr, 1, quantile, probs = CI.levels[2], na.rm = TRUE)
      )
      rownames(res[[i]]) <- as.character(seq_len(NROW(res[[i]])))
    }
  }
  
  y_obs <- split(newdata, id.)
  fitted.times <- split(newdata[[timeVar]], factor(newdata[[idVar]]))
  
  names(res) <- names(y_obs) <- names(last.time) <- names(obs.times) <- ids_chr
  names(proposal.repaired) <- ids_chr
  names(proposal.min_eig_raw) <- ids_chr
  names(proposal.min_eig_fix) <- ids_chr
  
  out_res <- list(
    summaries = res,
    survTimes = survTimes,
    last.time = last.time,
    obs.times = obs.times,
    y = y_obs,
    fitted.times = fitted.times,
    modes.b = modes.b,
    Vars.b = Vars.b,
    proposal.repaired = proposal.repaired,
    proposal.min_eig_raw = proposal.min_eig_raw,
    proposal.min_eig_fix = proposal.min_eig_fix
  )
  
  if (simulate) {
    out_res$full.results <- out
    out_res$success.rate <- success.rate
    
  }
  
  class(out_res) <- "survfitJM_fullplug"
  out_res
}


make_theta_sampler_global <- function(theta_mean, V_theta) { 
  if (is.null(names(theta_mean))) stop("theta_mean deve avere names().") 
  V_theta <- as.matrix(V_theta) 
  if (is.null(colnames(V_theta)) || is.null(rownames(V_theta))) {
    stop("V_theta deve avere rownames e colnames.") } 
  if (!setequal(rownames(V_theta), colnames(V_theta))) { 
    stop("rownames(V_theta) e colnames(V_theta) non coincidono come insieme.") } 
  V_theta <- V_theta[colnames(V_theta), colnames(V_theta), drop = FALSE] 
  if (!setequal(names(theta_mean), colnames(V_theta))) { 
    miss1 <- setdiff(names(theta_mean), colnames(V_theta)) 
    miss2 <- setdiff(colnames(V_theta), names(theta_mean)) 
    stop( "nomi incompatibili tra theta_mean e V_theta.\n", "in theta_mean ma non in V_theta: ", 
          paste(head(miss1, 20), collapse = ", "), "\n", "in V_theta ma non in theta_mean: ",
          paste(head(miss2, 20), collapse = ", ") ) }
  theta_mean <- theta_mean[colnames(V_theta)] 
  stopifnot(identical(names(theta_mean), colnames(V_theta))) 
  list( type = "global", theta_mean = theta_mean, V_theta = V_theta ) 
  }

.check_H_for_proposal <- function(H, tol = 1e-08) {
  Hs <- 0.5 * (H + t(H))
  ee <- eigen(Hs, symmetric = TRUE, only.values = TRUE)$values
  list(
    ok = all(is.finite(ee)) && min(ee) > tol,
    min_eig = min(ee),
    eig = ee,
    Hs = Hs
  )
}

.is_fake_subject <- function(df_i, idVar = "ID", fake_pattern = "^FAKE") {
  if (!idVar %in% names(df_i)) return(FALSE)
  ids <- unique(as.character(df_i[[idVar]]))
  length(ids) == 1L && grepl(fake_pattern, ids[1])
}

.repair_H_to_Sigma <- function(H,
                               scale = 1.6,
                               eig_floor_abs = 1e-06,
                               floor_quantile = 0.25,
                               floor_mult = 1,
                               max_var = 0.05) {
  Hs <- 0.5 * (H + t(H))
  ee <- eigen(Hs, symmetric = TRUE)
  
  lam <- ee$values
  Q <- ee$vectors
  
  lam_pos <- lam[is.finite(lam) & lam > eig_floor_abs]
  
  if (length(lam_pos) == 0L) {
    lam_floor <- 1 / max_var
  } else {
    lam_ref <- unname(stats::quantile(lam_pos, probs = floor_quantile, names = FALSE))
    lam_floor <- max(eig_floor_abs, floor_mult * lam_ref, scale / max_var)
  }
  
  lam_fix <- pmax(lam, lam_floor)
  
  H_fix <- Q %*% diag(lam_fix, nrow = length(lam_fix)) %*% t(Q)
  H_fix <- 0.5 * (H_fix + t(H_fix))
  
  Sigma <- scale * solve(H_fix)
  Sigma <- 0.5 * (Sigma + t(Sigma))
  
  eeS <- eigen(Sigma, symmetric = TRUE)
  lamS <- pmin(pmax(eeS$values, eig_floor_abs), max_var)
  Sigma <- eeS$vectors %*% diag(lamS, nrow = length(lamS)) %*% t(eeS$vectors)
  Sigma <- 0.5 * (Sigma + t(Sigma))
  
  list(
    H_raw = Hs,
    H_fix = H_fix,
    Sigma = Sigma,
    eig_raw = lam,
    eig_fix = lam_fix,
    min_eig_raw = min(lam),
    min_eig_fix = min(lam_fix),
    lam_floor = lam_floor
  )
}


as_plot_out_from_fullplug <- function(pred, id = NULL) {
  if (is.null(pred$summaries) || !is.list(pred$summaries) || length(pred$summaries) == 0) {
    stop("pred deve essere un output di survfitJM_fullplug() con pred$summaries.")
  }
  
  ids0 <- names(pred$summaries)
  if (is.null(ids0) || any(ids0 == "")) {
    if (!is.null(names(pred$last.time)) && length(pred$last.time) == length(pred$summaries)) {
      ids0 <- names(pred$last.time)
    } else if (!is.null(id) && length(pred$summaries) == 1L) {
      ids0 <- as.character(id)
    } else {
      ids0 <- paste0("ID", seq_along(pred$summaries))
    }
  }
  
  if (!is.null(id) && length(pred$summaries) == 1L) {
    ids0 <- as.character(id)
  }
  
  get_cols <- function(S) {
    cn <- colnames(S)
    if (is.null(cn)) stop("Ogni elemento di pred$summaries deve avere colonne nominate.")
    
    if ("predSurv" %in% cn) {
      list(
        times  = as.numeric(S[, "times"]),
        median = as.numeric(S[, "predSurv"]),
        lower  = as.numeric(S[, "predSurv"]),
        upper  = as.numeric(S[, "predSurv"])
      )
    } else {
      med_col <- if ("Median" %in% cn) "Median" else if ("Mean" %in% cn) "Mean" else NULL
      if (is.null(med_col)) stop("Nelle summaries non trovo né 'Median' né 'Mean' né 'predSurv'.")
      lo_col <- if ("Lower" %in% cn) "Lower" else med_col
      up_col <- if ("Upper" %in% cn) "Upper" else med_col
      
      list(
        times  = as.numeric(S[, "times"]),
        median = as.numeric(S[, med_col]),
        lower  = as.numeric(S[, lo_col]),
        upper  = as.numeric(S[, up_col])
      )
    }
  }
  
  ext1 <- get_cols(pred$summaries[[1]])
  t_grid <- ext1$times
  nt <- length(t_grid)
  nid <- length(pred$summaries)
  
  M_med <- matrix(NA_real_, nrow = nid, ncol = nt, dimnames = list(ids0, NULL))
  M_lo  <- matrix(NA_real_, nrow = nid, ncol = nt, dimnames = list(ids0, NULL))
  M_up  <- matrix(NA_real_, nrow = nid, ncol = nt, dimnames = list(ids0, NULL))
  
  for (i in seq_along(pred$summaries)) {
    Si <- pred$summaries[[i]]
    ex <- get_cols(Si)
    
    if (length(ex$times) != nt || any(abs(ex$times - t_grid) > 1e-10)) {
      stop("Le griglie dei tempi nelle summaries non coincidono tra soggetti.")
    }
    
    M_med[i, ] <- ex$median
    M_lo[i, ]  <- ex$lower
    M_up[i, ]  <- ex$upper
  }
  
  list(
    ids = ids0,
    t_grid = t_grid,
    S_pointwise = list(
      median = M_med,
      lower  = M_lo,
      upper  = M_up
    ),
    S_simultaneous = list(
      median = M_med,
      lower  = M_lo,
      upper  = M_up
    )
  )
}


.make_starts_b_fullplug <- function(q,
                                    D = NULL,
                                    n_start = 12,
                                    seed = NULL,
                                    sd_small = 0.15,
                                    sd_big = 0.50) {
  if (!is.null(seed)) set.seed(seed)
  
  starts <- vector("list", n_start)
  starts[[1]] <- rep(0, q)
  
  if (!is.null(D)) {
    D <- as.matrix(D)
    sds <- sqrt(pmax(diag(0.5 * (D + t(D))), 1e-8))
    sds_small <- pmax(sd_small * sds, 0.05)
    sds_big <- pmax(sd_big * sds, 0.10)
  } else {
    sds_small <- rep(sd_small, q)
    sds_big <- rep(sd_big, q)
  }
  
  k <- 2L
  
  while (k <= n_start) {
    starts[[k]] <- rnorm(q, 0, sds_small)
    k <- k + 1L
    if (k <= n_start) {
      starts[[k]] <- rnorm(q, 0, sds_big)
      k <- k + 1L
    }
  }
  
  starts
}

.grad_info_fullplug <- function(par, ff) {
  g <- tryCatch(
    numDeriv::grad(ff, par),
    error = function(e) rep(NA_real_, length(par))
  )
  
  list(
    grad = g,
    grad_norm = sqrt(sum(g^2, na.rm = TRUE)),
    grad_maxabs = max(abs(g), na.rm = TRUE)
  )
}

.fit_one_start_fullplug <- function(start,
                                    ff,
                                    use_nm = TRUE,
                                    use_grad_refine = TRUE,
                                    nm_maxit = 400,
                                    bfgs_maxit = 400,
                                    reltol = 1e-10) {
  par0 <- as.numeric(start)
  
  if (isTRUE(use_nm)) {
    nm <- try(
      optim(
        par0,
        ff,
        method = "Nelder-Mead",
        control = list(maxit = nm_maxit, reltol = reltol)
      ),
      silent = TRUE
    )
    
    if (!inherits(nm, "try-error") && is.list(nm) && is.finite(nm$value)) {
      par0 <- nm$par
    }
  }
  
  op1 <- try(
    optim(
      par0,
      ff,
      method = "BFGS",
      hessian = TRUE,
      control = list(maxit = bfgs_maxit, reltol = reltol)
    ),
    silent = TRUE
  )
  
  if (inherits(op1, "try-error") || !is.list(op1) || !is.finite(op1$value)) {
    return(NULL)
  }
  
  best <- op1
  
  if (isTRUE(use_grad_refine)) {
    gg <- function(bb) numDeriv::grad(ff, bb)
    
    op2 <- try(
      optim(
        op1$par,
        ff,
        gr = gg,
        method = "BFGS",
        hessian = TRUE,
        control = list(maxit = bfgs_maxit, reltol = reltol)
      ),
      silent = TRUE
    )
    
    if (!inherits(op2, "try-error") && is.list(op2) && is.finite(op2$value)) {
      if (op2$value <= best$value) best <- op2
    }
  }
  
  gi <- .grad_info_fullplug(best$par, ff)
  
  best$grad <- gi$grad
  best$grad_norm <- gi$grad_norm
  best$grad_maxabs <- gi$grad_maxabs
  best
}
.find_mode_b_fullplug <- function(object,
                                  df_i,
                                  t_cut,
                                  timeVar = "visit_m",
                                  scale = 1.6,
                                  pmvnorm_control = mvtnorm::GenzBretz(
                                    maxpts = 25000,
                                    abseps = 1e-04,
                                    releps = 0
                                  ),
                                  start = NULL,
                                  incremental = TRUE,
                                  n_start = 3,
                                  seed = NULL,
                                  stage_maxit = 60,
                                  final_maxit = 120,
                                  reltol = 1e-08,
                                  check_grad = FALSE) {
  q <- nrow(object$D)
  
  ff_full <- function(bb) {
    - .fullplug_logpost_b(
      bb,
      object = object,
      df_i = df_i,
      t_cut = t_cut,
      timeVar = timeVar,
      pmvnorm_control = pmvnorm_control
    )
  }
  
  run_bfgs <- function(par0, ff, maxit, hessian = FALSE) {
    try(
      optim(
        par = par0,
        fn = ff,
        method = "BFGS",
        hessian = hessian,
        control = list(maxit = maxit, reltol = reltol)
      ),
      silent = TRUE
    )
  }
  
  make_jitter <- function(center, mult = 0.15) {
    if (is.null(object$D)) {
      sds <- rep(mult, q)
    } else {
      sds <- sqrt(pmax(diag(0.5 * (object$D + t(object$D))), 1e-8))
      sds <- pmax(mult * sds, 0.05)
    }
    as.numeric(center + rnorm(q, 0, sds))
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  if (is.null(start)) {
    start0 <- rep(0, q)
  } else {
    start0 <- as.numeric(start)
    if (length(start0) != q) stop("'start' ha lunghezza incompatibile con q.")
  }
  
  starts <- list(start0)
  
  if (isTRUE(incremental)) {
    times_i <- sort(unique(df_i[[timeVar]][df_i[[timeVar]] <= t_cut]))
    times_i <- times_i[is.finite(times_i)]
    
    if (!length(times_i)) {
      times_i <- t_cut
    }
    
    par_path <- start0
    
    for (tt in times_i) {
      df_tt <- df_i[df_i[[timeVar]] <= tt, , drop = FALSE]
      
      ff_tt <- function(bb) {
        - .fullplug_logpost_b(
          bb,
          object = object,
          df_i = df_tt,
          t_cut = tt,
          timeVar = timeVar,
          pmvnorm_control = pmvnorm_control
        )
      }
      
      op_tt <- run_bfgs(par_path, ff_tt, maxit = stage_maxit, hessian = FALSE)
      
      if (!inherits(op_tt, "try-error") && is.list(op_tt) && is.finite(op_tt$value)) {
        par_path <- op_tt$par
      }
    }
    
    starts[[length(starts) + 1L]] <- par_path
  }
  
  while (length(starts) < n_start) {
    base_center <- starts[[length(starts)]]
    starts[[length(starts) + 1L]] <- make_jitter(base_center, mult = 0.10)
  }
  
  fits0 <- vector("list", length(starts))
  vals0 <- rep(Inf, length(starts))
  
  for (j in seq_along(starts)) {
    opj <- run_bfgs(starts[[j]], ff_full, maxit = final_maxit, hessian = FALSE)
    fits0[[j]] <- opj
    
    if (!inherits(opj, "try-error") && is.list(opj) && is.finite(opj$value)) {
      vals0[j] <- opj$value
    }
  }
  
  if (!any(is.finite(vals0))) {
    stop("nessuna ottimizzazione riuscita in .find_mode_b_fullplug().")
  }
  
  jbest <- which.min(vals0)
  par_best <- fits0[[jbest]]$par
  
  opt <- run_bfgs(par_best, ff_full, maxit = final_maxit, hessian = TRUE)
  
  if (inherits(opt, "try-error") || !is.list(opt) || !is.finite(opt$value)) {
    stop("ottimizzazione finale fallita in .find_mode_b_fullplug().")
  }
  
  H <- opt$hessian
  H <- 0.5 * (H + t(H))
  
  V <- try(scale * solve(H), silent = TRUE)
  if (inherits(V, "try-error") || any(!is.finite(V))) {
    V <- scale * .safe_hinv(H)
  }
  V <- 0.5 * (V + t(V))
  
  g <- rep(NA_real_, q)
  gnorm <- NA_real_
  gmax <- NA_real_
  
  if (isTRUE(check_grad)) {
    g <- try(cd(opt$par, ff_full), silent = TRUE)
    if (!inherits(g, "try-error")) {
      gnorm <- sqrt(sum(g^2))
      gmax <- max(abs(g))
    } else {
      g <- rep(NA_real_, q)
    }
  }
  
  list(
    mode = opt$par,
    var = V,
    hessian = H,
    opt = opt,
    grad = g,
    grad_norm = gnorm,
    grad_maxabs = gmax,
    objective = opt$value,
    starts = starts,
    start_values = vals0
  )
}


.fill_D_into_theta <- function(theta, D) {
  D_names <- rownames(D)
  
  for (nm in D_names) {
    vn <- paste0("var_", nm)
    if (vn %in% names(theta)) theta[[vn]] <- D[nm, nm]
  }
  
  for (i in seq_along(D_names)) {
    for (j in seq_len(i - 1L)) {
      a <- D_names[i]
      b <- D_names[j]
      nm1 <- paste0("cov_", a, "__", b)
      nm2 <- paste0("cov_", b, "__", a)
      if (nm1 %in% names(theta)) theta[[nm1]] <- D[a, b]
      if (nm2 %in% names(theta)) theta[[nm2]] <- D[a, b]
    }
  }
  
  theta
}

.fill_R_into_theta <- function(theta, R) {
  outcome_names <- rownames(R)
  
  for (i in 2:length(outcome_names)) {
    for (j in 1:(i - 1)) {
      a <- outcome_names[i]
      b <- outcome_names[j]
      rr <- R[a, b]
      
      nm_corr_1 <- paste0("corr_resid__", a, "__", b)
      nm_corr_2 <- paste0("corr_resid__", b, "__", a)
      nm_atanh_1 <- paste0("atanh_rho_resid__", a, "__", b)
      nm_atanh_2 <- paste0("atanh_rho_resid__", b, "__", a)
      
      if (nm_corr_1 %in% names(theta)) theta[[nm_corr_1]] <- rr
      if (nm_corr_2 %in% names(theta)) theta[[nm_corr_2]] <- rr
      
      zz <- atanh(pmin(pmax(rr, -0.999999), 0.999999))
      if (nm_atanh_1 %in% names(theta)) theta[[nm_atanh_1]] <- zz
      if (nm_atanh_2 %in% names(theta)) theta[[nm_atanh_2]] <- zz
    }
  }
  
  theta
}

.repair_pairwise_global_draw <- function(psi_draw,
                                         D_names,
                                         outcome_names,
                                         repair_D = TRUE,
                                         repair_R = TRUE,
                                         eps_D = 1e-8) {
  out <- psi_draw
  
  if (isTRUE(repair_D)) {
    D <- theta_to_D(out, D_names)
    D <- project_spd(D, eps = eps_D)
    out <- .fill_D_into_theta(out, D)
  }
  
  if (isTRUE(repair_R)) {
    R <- theta_to_R(out, outcome_names)
    R <- project_corr(R)
    out <- .fill_R_into_theta(out, R)
  }
  
  out
}




