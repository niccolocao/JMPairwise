library(Matrix)
library(numDeriv)

###########################################################
####                                                   ####
###                Estimates extractor                  ###
####                                                   ####
###########################################################
extract_params_vec <- function(res,
                               obj_name = NULL,
                               return_scores = TRUE,
                               K_mode= c("map","levels","observed"),
                               K_map = NULL,
                               report_issues = TRUE,
                               stop_on_issues = TRUE) {
  
  K_mode <- match.arg(K_mode)
  if (is.null(K_map)) stop("Devi passare K_map (es. list(Mobility=3, ..., Kneeling=5, ...)).")
  
  sanitize <- function(x) {
    x <- gsub("\\s+", "_", x)
    gsub("[^[:alnum:]_]+", "", x)
  }
  
  .issues = list()
  
  .count_bad <- function(x) {
    if (is.null(x)) return(c(na=0L,nan=0L,inf=0L))
    if (is.vector(x)) x <- as.numeric(x)
    x <- as.matrix(x)
    c(
      na  = sum(is.na(x)),
      nan = sum(is.nan(x)),
      inf = sum(is.infinite(x))
    )
  }
  
  .flag = function(tag, cnt) {
    .issues[[tag]] <<- cnt
    bad <- sum(cnt)
    if (isTRUE(report_issues) && bad > 0) warning(tag, ": na=", cnt["na"], " nan=", cnt["nan"], " inf=", cnt["inf"])
    if (isTRUE(stop_on_issues) && bad > 0) stop(tag, ": na/nan/inf presenti")
  } 
  
  .parse_items_from_obj <- function(obj_name) {
    nm <- sub("\\.rds$", "", basename(obj_name))
    nm <- sub("^[0-9]+_", "", nm)
    sp <- strsplit(nm, "__", fixed = TRUE)[[1]]
    if (length(sp) < 2) return(list(item1 = sp[1], item2 = NA_character_))
    list(item1 = sp[1], item2 = sp[2])
  }
  
  # --- K helper (MAP consigliato per il bivariato) ---
  get_K <- function(dat, item_raw) {
    item_raw <- as.character(item_raw)
    
    if (K_mode == "map") {
      if (!item_raw %in% names(K_map)) stop("Item non trovato in K_map: ", item_raw)
      return(as.integer(K_map[[item_raw]]))
    }
    
    yv <- NULL
    
    if (!is.null(dat) && item_raw %in% names(dat)) {
      yv <- dat[[item_raw]]
    } else if (!is.null(dat) && "item" %in% names(dat)) {
      yv <- if ("y_cat" %in% names(dat)) dat$y_cat[dat$item == item_raw] else dat$y[dat$item == item_raw]
    } else {
      stop("Impossible to retreive K: missing column '", item_raw, "' and data$item.")
    }
    
    if (K_mode == "levels") {
      if (!is.factor(yv) && !is.ordered(yv)) stop("K_mode='levels' must be factor/ordered. Item: ", item_raw)
      return(nlevels(as.factor(yv)))
    }
    
    length(unique(yv[!is.na(yv)]))
  }
  
  .get_items_pair <- function(res, obj_name) {
    p <- res$pair
    
    if (is.data.frame(p)) {
      if (all(c("ord_item","cont_item") %in% names(p))) {
        return(list(item1 = as.character(p$ord_item[1]), item2 = as.character(p$cont_item[1])))
      }
      if (all(c("item1","item2") %in% names(p))) {
        return(list(item1 = as.character(p$item1[1]), item2 = as.character(p$item2[1])))
      }
    } else if (is.list(p)) {
      if (!is.null(p$ord_item) && !is.null(p$cont_item)) {
        return(list(item1 = as.character(p$ord_item), item2 = as.character(p$cont_item)))
      }
      if (!is.null(p$item1) && !is.null(p$item2)) {
        return(list(item1 = as.character(p$item1), item2 = as.character(p$item2)))
      }
    }
    
    nm <- sub("\\.rds$", "", basename(obj_name))
    nm <- sub("^[0-9]+_", "", nm)
    sp <- strsplit(nm, "__", fixed = TRUE)[[1]]
    if (length(sp) < 2) stop("file name not valid: ", obj_name)
    list(item1 = sp[1], item2 = sp[2])
  }
  
  pair_get <- function(pair, nm) {
    if (is.null(pair)) return(NULL)
    
    if (is.data.frame(pair)) {
      if (nm %in% names(pair)) return(pair[[nm]][1])
      return(NULL)
    }
    
    nms <- names(pair)
    if (is.null(nms) || !(nm %in% nms)) return(NULL)
    pair[[nm]]
  }
  
  .get_pair_items_any <- function(res) {
    p <- res$pair
    
    o <- pair_get(p, "ord_item")
    c <- pair_get(p, "cont_item")
    if (!is.null(o) && !is.null(c) && !is.na(o) && !is.na(c)) {
      return(list(item1 = as.character(o), item2 = as.character(c)))
    }
    
    i1 <- pair_get(p, "item1")
    i2 <- pair_get(p, "item2")
    if (!is.null(i1) && !is.null(i2) && !is.na(i1) && !is.na(i2)) {
      return(list(item1 = as.character(i1), item2 = as.character(i2)))
    }
    
    dat <- res$fit$data
    cont_item <- if (!is.null(dat) && CONT_OUTCOME %in% names(dat)) CONT_OUTCOME else NULL
    if (is.null(cont_item)) stop("impossible to find cont_item (missing pair e missing ", CONT_OUTCOME, " in data).")
    
    cand <- intersect(names(K_map), names(dat))
    if (length(cand) == 0) stop("impossible to find ord_item from data/K_map.")
    ord_item <- cand[1]
    
    list(item1 = ord_item, item2 = cont_item)
  }
  
  .is_ord_norm_rho_fit <- function(res, obj_name) {
    pp <- .get_items_pair(res, obj_name)
    is_cont <- !is.na(pp$item2) && !(pp$item2 %in% names(K_map))
    isTRUE(is_cont && !is.null(res$fit$gammas) && length(res$fit$gammas) > 0)
  }
  
  .make_phi_names_ord_norm_rho <- function(ord_item_raw, cont_item_raw, K_ord) {
    c(
      paste0("threshold_", sanitize(ord_item_raw), "_cut:", seq_len(K_ord - 1L)),
      paste0("log_sigma_", sanitize(cont_item_raw)),
      paste0("atanh_rho_resid__", sanitize(ord_item_raw), "__", sanitize(cont_item_raw))
    )
  }
  
  get_item_surv = function(res) {
    it <- pair_get(res$pair, "item")
    if (!is.null(it) && length(it) >= 1 && !is.na(it[1])) return(as.character(it[1]))
    
    dat <- res$fit$data
    if (!is.null(dat) && "item" %in% names(dat)) {
      u <- unique(as.character(dat$item))
      if (length(u) >= 1) return(u[1])
    }
    
    stop("Impossible to find item for fit SURV (missing pair$item e data$item).")
  }
  
  has_SURV_fit <- function(res, obj_name) {
    if (!is.null(obj_name) && grepl("SURV", obj_name, ignore.case = TRUE)) return(TRUE)
    
    bnm <- names(res$fit$coefficients)
    if (!is.null(bnm) && any(grepl("(^|:)is_surv(:|$)", bnm))) return(TRUE)
    
    rnD <- rownames(res$fit$D)
    if (!is.null(rnD) && "z2" %in% rnD) return(TRUE)   # frailty tipica del SURV
    
    FALSE
  }
  
  CONT_OUTCOME <- "Perceivedcurrenthealthstatus"
  
  get_cont_outcome <- function(dat, obj_name = NULL) {
    if (!is.null(dat) && CONT_OUTCOME %in% names(dat)) return(CONT_OUTCOME)
    if (!is.null(obj_name)) {
      p <- .parse_items_from_obj(obj_name)
      if (!is.na(p$item1) && nzchar(p$item1)) return(p$item1)
    }
    CONT_OUTCOME
  }
  
  is_CONT_SURV_fit <- function(res, obj_name) {
    has_SURV_fit(res, obj_name) &&
      any(grepl("(^|:)is_cont(:|$)", names(res$fit$coefficients)))
  }
  
  .is_ord_beta_pair <- function(fit) {
    dat <- fit$data
    if (is.null(dat)) return(FALSE)
    if (!all(c("item","y_cat","y01") %in% names(dat))) return(FALSE)
    has_ord  <- any(!is.na(dat$y_cat))
    has_cont <- any(is.finite(dat$y01))
    nphi     <- length(fit$phis)
    isTRUE(has_ord && has_cont && nphi >= 2)
  }
  
  .ord_beta_items <- function(fit) {
    dat <- fit$data
    levs <- levels(dat$item)
    if (is.null(levs)) levs <- unique(as.character(dat$item))
    
    ord_counts  <- vapply(levs, function(it) sum(dat$item == it & !is.na(dat$y_cat)), integer(1))
    cont_counts <- vapply(levs, function(it) sum(dat$item == it &  is.finite(dat$y01)), integer(1))
    
    ord_item  <- levs[which.max(ord_counts)]
    cont_item <- levs[which.max(cont_counts)]
    list(ord_item = ord_item, cont_item = cont_item)
  }

                          
  .parse_Dij <- function(x) {
    tok <- sub("^D_", "", x)
    pr  <- strsplit(tok, "_", fixed = TRUE)[[1]]
    if (length(pr) == 2) {
      i <- as.integer(pr[1]); j <- as.integer(pr[2])
    } else {
      m <- regexec("^([0-9]+)([0-9]+)$", tok)
      r <- regmatches(tok, m)[[1]]
      if (length(r) != 3) stop("Formato D non riconosciuto: ", x)
      i <- as.integer(r[2]); j <- as.integer(r[3])
    }
    c(i = i, j = j)
  }

                          
  parse_z_pair_any <- function(z, items_map) {
    z0 <- sub("^zi_", "", z)  
    m  <- regexec("^z([12])(I0|Ipost)$", z0)
    r  <- regmatches(z0, m)[[1]]
    if (length(r) != 3) stop("Nome RE non riconosciuto (atteso z1I0/z1Ipost/z2I0/z2Ipost con opzionale zi_): ", z)
    idx  <- as.integer(r[2])
    when <- r[3]
    item <- items_map[as.character(idx)]
    list(item = item, when = when, idx = idx)
  }
  
  D_raw_from_Dij <- function(D_name, D_mat) {
    ij <- .parse_Dij(D_name)
    ii <- min(ij[1], ij[2])
    jj <- max(ij[1], ij[2])
    rn <- rownames(D_mat)
    q  <- length(rn)
    
    # SURV branch 
    if (q == 3 && all(c("z1I0","z1Ipost","z2") %in% rn)) {
      item_raw <- if (is_CONT_SURV_fit(res, obj_name)) {
        get_cont_outcome(res$fit$data, obj_name)
      } else {
        get_item_surv(res)
      }
      item <- sanitize(item_raw)
      
      re_label <- function(z) {
        if (z == "z2") return("frailty")
        m <- regexec("^z1(I0|Ipost)$", z)
        r <- regmatches(z, m)[[1]]
        paste0(item, "_", r[2])
      }
      labs <- vapply(rn, re_label, character(1))
      
      if (ii == jj) {
        lab <- labs[ii]
        if (lab == "frailty") "log_chol_var_frailty" else paste0("log_chol_var_", lab)
      } else {
        li <- labs[ii]; lj <- labs[jj]
        if ("frailty" %in% c(li, lj)) {
          lab_it <- if (li == "frailty") lj else li
          when <- sub(".*_", "", lab_it)
          it0  <- sub("_(I0|Ipost)$", "", lab_it)
          paste0("chol_cov_", it0, "_frailty_", when)
        } else {
          wi <- sub(".*_", "", li); wj <- sub(".*_", "", lj)
          ord <- c(wi, wj); ord <- ord[order(match(ord, c("I0","Ipost")))]
          paste0("chol_cov_", item, "_", ord[2], "_", ord[1])
        }
      }
    } else {
      # PAIR branch
      pp <- .get_items_pair(res, obj_name)
      items_raw <- c(pp$item1, pp$item2)
      items_map <- setNames(sanitize(items_raw), c("1","2"))
      
      info <- lapply(rn, parse_z_pair_any, items_map = items_map)
      labs <- vapply(info, function(u) paste0(u$item, "_", u$when), character(1))
      
      if (ii == jj) {
        paste0("log_chol_var_", labs[ii])
      } else {
        li <- info[[ii]]; lj <- info[[jj]]
        if (li$idx == lj$idx) {
          ord <- c(li$when, lj$when)
          ord <- ord[order(match(ord, c("I0","Ipost")))]
          paste0("chol_cov_", li$item, "_", ord[2], "_", ord[1])
        } else {
          paste0("chol_cov_cross_", li$item, "_", lj$item, "_", li$when, "_", lj$when)
        }
      }
    }
  }
  
  
  # CHOL vector da D
  chol_from_D_SURV <- function(D, item_sanit) {
    dn <- rownames(D)
    
    re_label <- function(z) {
      if (z == "z2") return("frailty")
      m <- regexec("^z1(I0|Ipost)$", z)
      r <- regmatches(z, m)[[1]]
      paste0(item_sanit, "_", r[2])
    }
    labs <- vapply(dn, re_label, character(1))
    
    L <- t(chol(D))
    out <- numeric(0)
    
    diag_names <- ifelse(labs == "frailty", "log_chol_var_frailty", paste0("log_chol_var_", labs))
    out <- c(out, setNames(log(diag(L)), diag_names))
    
    off_vals <- numeric(0); off_names <- character(0)
    for (j in seq_len(nrow(L))) {
      if (j == 1) next
      for (i in seq_len(j - 1)) {
        li <- labs[i]; lj <- labs[j]
        nm <- if ("frailty" %in% c(li, lj)) {
          it_side <- if (li == "frailty") lj else li
          when <- sub(".*_", "", it_side)
          it0  <- sub("_(I0|Ipost)$", "", it_side)
          paste0("chol_cov_", it0, "_frailty_", when)
        } else {
          wi <- sub(".*_", "", li); wj <- sub(".*_", "", lj)
          ord <- c(wi, wj); ord <- ord[order(match(ord, c("I0","Ipost")))]
          paste0("chol_cov_", item_sanit, "_", ord[2], "_", ord[1])
        }
        off_vals <- c(off_vals, L[j,i])
        off_names <- c(off_names, nm)
      }
    }
    c(out, setNames(off_vals, off_names))
  }

                     
  chol_from_D_PAIR <- function(D, items_raw) {
    dn <- rownames(D)
    items_map <- setNames(sanitize(items_raw), as.character(seq_along(items_raw)))
    
    info <- lapply(dn, parse_z_pair_any, items_map = items_map)
    
    L <- t(chol(D))
    out <- numeric(0)
    
    diag_names <- vapply(seq_len(nrow(L)), function(i) {
      ri <- info[[i]]
      paste0("log_chol_var_", ri$item, "_", ri$when)
    }, character(1))
    out <- c(out, setNames(log(diag(L)), diag_names))
    
    off_vals <- numeric(0); off_names <- character(0)
    for (j in seq_len(nrow(L))) {
      if (j == 1) next
      for (i in seq_len(j - 1)) {
        ri <- info[[i]]; rj <- info[[j]]
        nm <- if (ri$idx == rj$idx) {
          ord <- c(ri$when, rj$when)
          ord <- ord[order(match(ord, c("I0","Ipost")))]
          paste0("chol_cov_", ri$item, "_", ord[2], "_", ord[1])
        } else {
          paste0("chol_cov_cross_", ri$item, "_", rj$item, "_", ri$when, "_", rj$when)
        }
        off_vals <- c(off_vals, L[j,i])
        off_names <- c(off_names, nm)
      }
    }
    c(out, setNames(off_vals, off_names))
  }

                     
  # rename betas (stacked)
  rename_betas <- function(betas, has_SURV, item_raw) {
    bnames <- names(betas)
    
    if (has_SURV) {
      token_non <- if (any(grepl("(^|:)is_cont(:|$)", bnames))) "is_cont" else "is_ord"
      
      base_raw <- if (token_non == "is_cont") get_cont_outcome(res$fit$data, obj_name) else item_raw
      base <- sanitize(base_raw)
      
      is_non  <- grepl(paste0("(^|:)", token_non, "(:|$)"), bnames)
      is_surv <- grepl("(^|:)is_surv(:|$)", bnames)
      
      strip_token <- function(nm, token) {
        nm <- gsub(paste0("(^|:)", token, "(:|$)"), ":", nm)
        gsub("^:|:$", "", nm)
      }
      
      non_cov  <- sanitize(gsub(":", "_", strip_token(bnames[is_non],  token_non)))
      surv_cov <- sanitize(gsub(":", "_", strip_token(bnames[is_surv], "is_surv")))
      
      bnew <- character(length(betas))
      bnew[is_non]  <- paste0("fixed_", base, "_", non_cov)
      bnew[is_surv] <- paste0("fixed_surv_", surv_cov)
      
      other <- !(is_non | is_surv)
      bare_intercept <- other & bnames == "(Intercept)"
      
      if (any(bare_intercept) && token_non == "is_cont") {
        bnew[bare_intercept] <- paste0("fixed_", base, "_Intercept")
      }
      
      other_rest <- other & !bare_intercept
      if (any(other_rest)) {
        bnew[other_rest] <- paste0("fixed_", sanitize(bnames[other_rest]))
      }
      
      names(betas) <- bnew
      return(betas)
    }
    
    bnew <- sanitize(gsub(":", "_", sub("^item", "fixed_", bnames)))
    names(betas) <- bnew
    betas
  }

                     
  rename_betas_bivar_pair <- function(vec, item_raw) {
    nm <- names(vec)
    if (is.null(nm)) nm <- paste0("V", seq_along(vec))
    names(vec) <- paste0("fixed_", sanitize(item_raw), "_", sanitize(gsub(":", "_", nm)))
    vec
  }
  
 
  # phis names helpers (e.g. thresholds)
  make_phi_names_SURV <- function(item_raw, K) {
    item <- sanitize(item_raw)
    c(paste0("threshold_", item, "_cut:", seq_len(K - 1L)),
      "log_lambda", "log_rho")
  }
  
  make_phi_names_CONT_SURV <- function(cont_raw = CONT_OUTCOME) {
    c(paste0("log_sigma_", sanitize(cont_raw)),
      "log_lambda", "log_rho")
  }
  
  make_phi_names_ordbeta <- function(ord_item_raw, cont_item_raw, K_ord) {
    c(paste0("threshold_", sanitize(ord_item_raw), "_cut:", seq_len(K_ord - 1L)),
      paste0("beta_log_precision_", sanitize(cont_item_raw)))
  }
  
  make_phi_names_pairord <- function(item1_raw, item2_raw, K1, K2, with_resid_corr = FALSE) {
    base <- c(
      paste0("threshold_", sanitize(item1_raw), "_cut:", seq_len(K1 - 1L)),
      paste0("threshold_", sanitize(item2_raw), "_cut:", seq_len(K2 - 1L))
    )
    if (with_resid_corr) {
      base <- c(base, paste0("atanh_rho_resid__", sanitize(item1_raw), "__", sanitize(item2_raw)))
    }
    base
  }
  
  
  #  pair ordinal + rho residual (Fisher-z)
  make_phi_names_pairord_rho <- function(item1_raw, item2_raw, K1, K2) {
    make_phi_names_pairord(item1_raw, item2_raw, K1, K2, with_resid_corr = TRUE)
  }
  
  
  # CORE: pars_chol
  core_extract <- function() {
    has_SURV <- has_SURV_fit(res, obj_name)
    
    if (has_SURV) {
      if (is_CONT_SURV_fit(res, obj_name)) {
        base_raw <- get_cont_outcome(res$fit$data, obj_name)
        betas <- rename_betas(res$fit$coefficients, has_SURV = TRUE, item_raw = base_raw)
        
        phis <- res$fit$phis
        stopifnot(length(phis) == 3L)
        names(phis) <- make_phi_names_CONT_SURV(base_raw)
        
        D_chol <- chol_from_D_SURV(res$fit$D, sanitize(base_raw))
        return(c(betas, phis, D_chol))
      }
      
      item_raw <- get_item_surv(res)
      betas <- rename_betas(res$fit$coefficients, has_SURV = TRUE, item_raw = item_raw)
      
      K <- get_K(res$fit$data, item_raw)
      phis <- res$fit$phis
      stopifnot(length(phis) == (K - 1L) + 2L)
      names(phis) <- make_phi_names_SURV(item_raw, K)
      
      D_chol <- chol_from_D_SURV(res$fit$D, sanitize(item_raw))
      return(c(betas, phis, D_chol))
    }
    
    if (.is_ord_norm_rho_fit(res, obj_name)) {
      pp <- .get_items_pair(res, obj_name)
      ord_item_raw  <- pp$item1
      cont_item_raw <- pp$item2
      
      K_ord <- get_K(res$fit$data, ord_item_raw)
      
      bet1 <- rename_betas_bivar_pair(res$fit$coefficients, ord_item_raw)
      bet2 <- rename_betas_bivar_pair(res$fit$gammas,       cont_item_raw)
      
      phis <- res$fit$phis
      stopifnot(length(phis) == (K_ord - 1L) + 2L)
      names(phis) <- c(
        paste0("threshold_", sanitize(ord_item_raw), "_cut:", seq_len(K_ord - 1L)),
        paste0("log_sigma_", sanitize(cont_item_raw)),
        paste0("atanh_rho_resid__", sanitize(ord_item_raw), "__", sanitize(cont_item_raw))
      )
      
      D_chol <- chol_from_D_PAIR(res$fit$D, c(ord_item_raw, cont_item_raw))
      return(c(bet1, bet2, phis, D_chol))
    }
    
    if (.is_ord_beta_pair(res$fit)) {
      
      pp <- .get_items_pair(res, obj_name)
      items_raw <- c(pp$item1, pp$item2)
      if (length(items_raw) != 2 || any(is.na(items_raw))) {
        stop("Fit PAIR: pair$item1/item2 missing.")
      }
      if (length(items_raw) != 2) items_raw <- unique(as.character(res$fit$data$item))[1:2]
      
      betas <- rename_betas(res$fit$coefficients, has_SURV = FALSE, item_raw = NULL)
      
      info <- .ord_beta_items(res$fit)
      ord_item_raw  <- info$ord_item
      cont_item_raw <- info$cont_item
      
      K_ord <- get_K(res$fit$data, ord_item_raw)
      phis <- res$fit$phis
      stopifnot(length(phis) == K_ord)
      names(phis) <- make_phi_names_ordbeta(ord_item_raw, cont_item_raw, K_ord)
      
      D_chol <- chol_from_D_PAIR(res$fit$D, items_raw)
      return(c(betas, phis, D_chol))
    }
    
    
    pp <- .get_items_pair(res, obj_name)
    items_raw <- c(pp$item1, pp$item2)
    
    K1 <- get_K(res$fit$data, items_raw[1])
    K2 <- get_K(res$fit$data, items_raw[2])
    
    phis <- res$fit$phis
    
    # if gammas are available (i.e., residual correlation)
    if (!is.null(res$fit$gammas) && length(res$fit$gammas) > 0) {
      bet1 <- rename_betas_bivar_pair(res$fit$coefficients, items_raw[1])
      bet2 <- rename_betas_bivar_pair(res$fit$gammas,       items_raw[2])

      
      stopifnot(length(phis) == (K1 - 1L) + (K2 - 1L) + 1L)
      names(phis) <- make_phi_names_pairord_rho(items_raw[1], items_raw[2], K1, K2)
      
      D_chol <- chol_from_D_PAIR(res$fit$D, items_raw)
      return(c(bet1, bet2, phis, D_chol))
    }

    # without residual correlation
    betas <- rename_betas(res$fit$coefficients, has_SURV = FALSE, item_raw = NULL)
    stopifnot(length(phis) == (K1 - 1L) + (K2 - 1L))
    names(phis) <- make_phi_names_pairord(items_raw[1], items_raw[2], K1, K2)
    
    D_chol <- chol_from_D_PAIR(res$fit$D, items_raw)
    c(betas, phis, D_chol)
  }
  
  pars_chol <- core_extract()
  .flag("pars", .count_bad(pars_chol))
  if (!isTRUE(return_scores)) return(pars_chol)
  
  # Hessian rename and reorder
  H_pretty <- res$fit$Hessian
  .flag("H", .count_bad(H_pretty))
  
  rename_phi <- function(vec) {
    has_SURV <- has_SURV_fit(res, obj_name)
    
    if (has_SURV) {
      if (is_CONT_SURV_fit(res, obj_name)) {
        cont_raw <- get_cont_outcome(res$fit$data, obj_name)
        newn <- make_phi_names_CONT_SURV(cont_raw)
        for (j in seq_len(3)) vec[vec == paste0("phi_", j)] <- newn[j]
        return(vec)
      }
      
      dat <- res$fit$data
      item_raw <- get_item_surv(res)
      K <- get_K(dat, item_raw)
      nth <- K - 1L
      
      if (nth > 0) {
        for (j in seq_len(nth)) {
          vec[vec == paste0("phi_", j)] <- paste0("threshold_", sanitize(item_raw), "_cut:", j)
        }
      }
      vec[vec == paste0("phi_", nth + 1L)] <- "log_lambda"
      vec[vec == paste0("phi_", nth + 2L)] <- "log_rho"
      return(vec)
    }
    
    if (.is_ord_norm_rho_fit(res, obj_name)) {
      pp <- .get_items_pair(res, obj_name)
      ord_item_raw  <- pp$item1
      cont_item_raw <- pp$item2
      
      K_ord <- get_K(res$fit$data, ord_item_raw)
      
      nth <- K_ord - 1L
      if (nth > 0) {
        for (j in seq_len(nth)) {
          vec[vec == paste0("phi_", j)] <- paste0("threshold_", sanitize(ord_item_raw), "_cut:", j)
        }
      }
      vec[vec == paste0("phi_", nth + 1L)] <- paste0("log_sigma_", sanitize(cont_item_raw))
      vec[vec == paste0("phi_", nth + 2L)] <- paste0("atanh_rho_resid__", sanitize(ord_item_raw), "__", sanitize(cont_item_raw))
      return(vec)
    }
    
    if (.is_ord_beta_pair(res$fit)) {
      info <- .ord_beta_items(res$fit)
      ord_item_raw  <- info$ord_item
      cont_item_raw <- info$cont_item
      K_ord <- get_K(res$fit$data, ord_item_raw)
      kth <- K_ord - 1L
      if (kth > 0) {
        for (j in seq_len(kth)) {
          vec[vec == paste0("phi_", j)] <- paste0("threshold_", sanitize(ord_item_raw), "_cut:", j)
        }
      }
      vec[vec == paste0("phi_", K_ord)] <- paste0("beta_log_precision_", sanitize(cont_item_raw))
      return(vec)
    }
    
    # pair ordinal
    pp <- .get_items_pair(res, obj_name)
    items_raw <- c(pp$item1, pp$item2)
    if (length(items_raw) != 2 || any(is.na(items_raw))) return(vec)
    K1 <- get_K(res$fit$data, items_raw[1])
    K2 <- get_K(res$fit$data, items_raw[2])
    k1m1 <- K1 - 1L
    k2m1 <- K2 - 1L
    
    if (k1m1 > 0) {
      for (j in seq_len(k1m1)) {
        vec[vec == paste0("phi_", j)] <- paste0("threshold_", sanitize(items_raw[1]), "_cut:", j)
      }
    }
    if (k2m1 > 0) {
      for (j in seq_len(k2m1)) {
        vec[vec == paste0("phi_", k1m1 + j)] <- paste0("threshold_", sanitize(items_raw[2]), "_cut:", j)
      }
    }

    
    idx_rho <- k1m1 + k2m1 + 1L
    old <- paste0("phi_", idx_rho)
    new <- paste0("atanh_rho_resid__", sanitize(items_raw[1]), "__", sanitize(items_raw[2]))
    if (any(vec == old)) {
      vec[vec == old] <- new
    }
    
    vec
  }
  
  if (!is.null(dim(H_pretty))) {
    cn <- colnames(H_pretty); rn <- rownames(H_pretty)
    
    vc <- colnames(vcov(res$fit))
    Dij_vcov <- vc[grepl("^D_", vc)]
    target_D_raw <- vapply(Dij_vcov, function(x) D_raw_from_Dij(x, res$fit$D), character(1))
    
    D_idx_c <- grepl("^D_", cn)
    D_idx_r <- grepl("^D_", rn)
    
    if (any(D_idx_c)) cn[D_idx_c] <- vapply(cn[D_idx_c], function(x) D_raw_from_Dij(x, res$fit$D), character(1))
    if (any(D_idx_r)) rn[D_idx_r] <- vapply(rn[D_idx_r], function(x) D_raw_from_Dij(x, res$fit$D), character(1))
    
    colnames(H_pretty) <- cn
    rownames(H_pretty) <- rn

    cn <- colnames(H_pretty)
    
    phi_cols  <- which(grepl("^phi_\\d+$", cn))
    
    d_cols <- match(target_D_raw, cn)
    d_cols <- d_cols[!is.na(d_cols)]
    
    
                                            
    b_old <- names(res$fit$coefficients)
    beta_cols <- match(b_old, cn); beta_cols <- beta_cols[!is.na(beta_cols)]
    
   
                                            
    gamma_cols <- integer(0)
    g_old <- if (!is.null(res$fit$gammas)) names(res$fit$gammas) else NULL
    if (!is.null(g_old) && length(g_old)) {
      g_candidates <- ifelse(paste0("zi_", g_old) %in% cn, paste0("zi_", g_old), g_old)
      gamma_cols <- match(g_candidates, cn)
      gamma_cols <- gamma_cols[!is.na(gamma_cols)]
      gamma_cols <- setdiff(gamma_cols, beta_cols)
    }
    
    keep_core <- c(beta_cols, gamma_cols, d_cols, phi_cols)
    keep_core <- keep_core[!is.na(keep_core)]
    other_cols <- setdiff(seq_along(cn), keep_core)
    
    perm <- c(beta_cols, gamma_cols, d_cols, phi_cols, other_cols)
    perm <- perm[!is.na(perm)]
    perm <- unique(perm)
    H_pretty <- H_pretty[perm, perm, drop = FALSE]
    
    # rename betas e gammas in H 
    # betas: fixed_<item1>_<cov>
    pp <- .get_items_pair(res, obj_name)
    it1_raw <- pp$item1
    it2_raw <- pp$item2
    
    if (!is.null(res$fit$gammas) && length(res$fit$gammas) > 0) {
      
      it1 <- sanitize(it1_raw)
      it2 <- sanitize(it2_raw)
      
      cn <- colnames(H_pretty); rn <- rownames(H_pretty)
      
      # map betas
      b_map <- setNames(paste0("fixed_", it1, "_", sanitize(b_old)), b_old)
      
      # map gammas
      gH_old <- ifelse(paste0("zi_", g_old) %in% cn, paste0("zi_", g_old), g_old)
      g_map  <- setNames(paste0("fixed_", it2, "_", sanitize(g_old)), gH_old)
      
      replace_names <- function(v, mp) {
        hit <- match(v, names(mp))
        ok  <- hit > 0 & !is.na(hit)
        v[ok] <- unname(mp[v[ok]])
        v
      }
      
      cn <- replace_names(cn, b_map)
      rn <- replace_names(rn, b_map)
      cn <- replace_names(cn, g_map)
      rn <- replace_names(rn, g_map)
      
      colnames(H_pretty) <- cn
      rownames(H_pretty) <- rn
      
    } else {
      b_old <- names(res$fit$coefficients)
      if (!is.null(b_old)) {
        has_SURV <- has_SURV_fit(res, obj_name)
        
        if (has_SURV) {
          item_raw <- if (is_CONT_SURV_fit(res, obj_name)) get_cont_outcome(res$fit$data) else get_item_surv(res)
          
          
          b_new <- names(rename_betas(setNames(res$fit$coefficients, b_old),
                                      has_SURV = TRUE, item_raw = item_raw))
        } else {
          b_new <- names(rename_betas(setNames(res$fit$coefficients, b_old),
                                      has_SURV = FALSE, item_raw = NULL))
        }
        
        cn <- colnames(H_pretty); rn <- rownames(H_pretty)
        m  <- match(cn, b_old, nomatch = 0); if (any(m)) cn[m > 0] <- b_new[m[m > 0]]
        m  <- match(rn, b_old, nomatch = 0); if (any(m)) rn[m > 0] <- b_new[m[m > 0]]
        colnames(H_pretty) <- cn
        rownames(H_pretty) <- rn
      }
    }
    
    colnames(H_pretty) <- rename_phi(colnames(H_pretty))
    rownames(H_pretty) <- rename_phi(rownames(H_pretty))
  }
  
  # gradient (same as in H_pretty)
  IDfac <- res$fit$id$ID
  Sb <- res$fit$score_vect_contributions$score.betas
  Sg <- if (!is.null(res$fit$score_vect_contributions$score.gammas)) res$fit$score_vect_contributions$score.gammas else NULL
  Sp <- res$fit$score_vect_contributions$score.phis
  SD <- res$fit$score_vect_contributions$score.D
  
  ids <- {
    pm <- res$fit$post_modes
    if (!is.null(pm) && !is.null(rownames(pm))) rownames(pm)
    else unique(as.character(IDfac))
  }
  
  align_rows <- function(M, ids) {
    if (is.null(M)) return(NULL)
    M <- as.matrix(M)
    if (is.null(rownames(M))) stop("Blocco senza rownames: impossibile allineare.")
    out <- matrix(0, nrow = length(ids), ncol = ncol(M),
                  dimnames = list(ids, colnames(M)))
    rr <- match(rownames(M), ids)
    out[rr[!is.na(rr)], ] <- M[!is.na(rr), , drop = FALSE]
    if (!all(ids %in% rownames(M))) {
      miss <- setdiff(ids, rownames(M))
      .issues[[paste0("align_missing_rows_", paste(colnames(M), collapse="|"))]] <<- head(miss, 20)
      if (isTRUE(report_issues)) warning("align_rows: missing ", length(miss), " ids in this block")
      if (isTRUE(stop_on_issues)) stop("align_rows: missing ids in score block")
    }
    out
  }
  
  # ---- betas ----
  Gb <- if (!is.null(Sb)) {
    .flag("score.betas", .count_bad(Sb))
    tmp <- rowsum(Sb, group = IDfac, na.rm = FALSE)
    .flag("Gb_after_rowsum", .count_bad(tmp))
    
    # fixed_<item1>_<cov>
    pp <- .get_items_pair(res, obj_name)
    it1_raw <- pp$item1
    it2_raw <- pp$item2
    
    if (!is.null(res$fit$gammas) && length(res$fit$gammas) > 0) {
      b_old <- names(res$fit$coefficients)
      b_new <- paste0("fixed_", sanitize(it1_raw), "_", sanitize(b_old))
      if (ncol(tmp) == length(b_new)) colnames(tmp) <- b_new
    } else {
       has_SURV <- has_SURV_fit(res, obj_name)
      if (has_SURV) {
        item_raw <- if (is_CONT_SURV_fit(res, obj_name)) get_cont_outcome(res$fit$data)
        else if (!is.null(res$pair$item)) res$pair$item
        else as.character(unique(res$fit$data$item))[1]
        
        b_new <- names(rename_betas(setNames(res$fit$coefficients, names(res$fit$coefficients)),
                                    has_SURV = TRUE, item_raw = item_raw))
      } else {
        b_new <- names(rename_betas(setNames(res$fit$coefficients, names(res$fit$coefficients)),
                                    has_SURV = FALSE, item_raw = NULL))
      }
      if (ncol(tmp) == length(b_new)) colnames(tmp) <- b_new
    }
    
    tmp
  } else NULL
  
  # gammas 
  Gg = if (!is.null(Sg)) {
    .flag("score.gammas", .count_bad(Sg))
    tmp <- rowsum(Sg, group = IDfac, na.rm = FALSE)
    .flag("Gg_after_rowsum", .count_bad(tmp))
    pp <- .get_items_pair(res, obj_name)
    it1_raw <- pp$item1
    it2_raw <- pp$item2
    
    if (!is.null(res$fit$gammas) && length(res$fit$gammas) > 0) {
      g_old <- names(res$fit$gammas)
      pp <- .get_items_pair(res, obj_name)
      g_new <- paste0("fixed_", sanitize(pp$item2), "_", sanitize(g_old))
      if (ncol(tmp) == length(g_new)) colnames(tmp) <- g_new
    }
    tmp
  } else NULL
  
  #phis 
  Gp <- if (!is.null(Sp)) {
    .flag("score.phis", .count_bad(Sp))
    tmp <- rowsum(Sp, group = IDfac, na.rm = FALSE)
    .flag("Gp_after_rowsum", .count_bad(tmp))
    pnames <- paste0("phi_", seq_len(ncol(tmp)))
    pnames <- rename_phi(pnames)
    colnames(tmp) <- pnames
    tmp
  } else NULL
  
  # D 
  GD <- if (!is.null(SD)) {
    .flag("score.D", .count_bad(SD))
    tmp <- as.matrix(SD)
    if (is.null(rownames(tmp))) {
      pm <- res$fit$post_modes
      if (!is.null(pm) && !is.null(rownames(pm)) && nrow(pm) == nrow(tmp)) {
        rownames(tmp) <- rownames(pm)
      } else {
        rownames(tmp) <- unique(as.character(IDfac))
      }
    }
    
    D_names <- colnames(H_pretty)[grepl("^(log_chol_var_|chol_cov_)", colnames(H_pretty))]
    if (length(D_names) == 0) D_names <- colnames(H_pretty)[grepl("^D_", colnames(H_pretty))]
    stopifnot(length(D_names) == ncol(tmp))
    colnames(tmp) <- D_names
    tmp
  } else NULL
  
  Gb <- align_rows(Gb, ids)
  Gg <- align_rows(Gg, ids)
  GD <- align_rows(GD, ids)
  Gp <- align_rows(Gp, ids)
  
  blocks <- Filter(function(x) !is.null(x) && ncol(x) > 0, list(Gb, Gg, GD, Gp))
  
  G_all  <- if (length(blocks)) do.call(cbind, blocks) else NULL
  if (is.null(G_all)) {
    .issues[["grad_is_null"]] <- TRUE
    if (isTRUE(report_issues)) warning("grad is NULL")
    if (isTRUE(stop_on_issues)) stop("grad is NULL")
  } else {
    .flag("grad", .count_bad(G_all))
  }
  
  if (!is.null(G_all) && !is.null(H_pretty) && length(colnames(H_pretty))) {
    target <- colnames(H_pretty)
    cnG <- colnames(G_all)
    
    keep <- target[target %in% cnG]
    extra <- setdiff(cnG, keep)
    
    ord <- c(keep, extra)
    ord <- ord[!is.na(ord)]
    ord <- ord[ord %in% cnG]
    
    if (length(ord)) G_all <- G_all[, ord, drop = FALSE]
  }  
  if (!is.null(G_all) && !is.null(H_pretty)) {
    miss_in_grad <- setdiff(colnames(H_pretty), colnames(G_all))
    extra_in_grad <- setdiff(colnames(G_all), colnames(H_pretty))
    .issues[["grad_missing_cols"]] <- miss_in_grad
    .issues[["grad_extra_cols"]] <- extra_in_grad
    if (isTRUE(report_issues) && (length(miss_in_grad) || length(extra_in_grad))) {
      warning("mismatch cols: missing_in_grad=", length(miss_in_grad), " extra_in_grad=", length(extra_in_grad))
    }
    if (isTRUE(stop_on_issues) && (length(miss_in_grad) || length(extra_in_grad))) {
      stop("mismatch colonne tra H e grad")
    }
  }
  list(pars = pars_chol, grad = G_all, H = H_pretty, issues = .issues)
}


process_one <- function(fp, K_map) {
  res <- readRDS(fp)
  
  if (is.data.frame(res$pair)) {
    if (!"item1" %in% names(res$pair) && "ord_item" %in% names(res$pair)) res$pair$item1 <- res$pair$ord_item
    if (!"item2" %in% names(res$pair) && "cont_item" %in% names(res$pair)) res$pair$item2 <- res$pair$cont_item
  }
  
  ex <- tryCatch(
    extract_params_vec(
      res,
      obj_name      = basename(fp),
      return_scores = TRUE,
      K_mode        = "map",
      K_map         = K_map
    ),
    error = function(e) stop(basename(fp), " -> ", conditionMessage(e))
  )
  
  list(file = fp, ex = ex, pars = ex$pars, grad = ex$grad, H = ex$H, issues = ex$issues)
}

####################################
####################################
######                       #######
##      J and K matrices          ##
######                       #######
####################################

JK_all_pairs <- function(ex_list, hessian_of = c("neglogLik", "logLik")) {
  hessian_of <- match.arg(hessian_of)

  
  ids_all <- Reduce(union, lapply(ex_list, function(ex) rownames(ex$grad)))
  if (is.null(ids_all)) stop("Ogni ex$grad deve avere rownames=ID.")
  N <- length(ids_all)
  if (N < 2) stop("Serve N>=2.")
  
  J_blocks <- vector("list", length(ex_list))
  U_blocks <- vector("list", length(ex_list))
  
  for (p in seq_along(ex_list)) {
    ex <- ex_list[[p]]
    U <- as.matrix(ex$grad)
    H <- as.matrix(ex$H)
    
    U[is.na(U)] <- 0
    H <- 0.5 * (H + t(H))

    
    if (!is.null(colnames(H)) && !is.null(colnames(U))) {
      tgt <- colnames(H)
      miss <- setdiff(tgt, colnames(U))
      if (length(miss)) stop("Pair ", p, ": grad manca colonne in H: ", paste(miss, collapse=", "))
      U <- U[, tgt, drop = FALSE]
      H <- H[tgt, tgt, drop = FALSE]
    } else {
      if (ncol(U) != ncol(H)) stop("Pair ", p, ": mismatch grad/H.")
    }
    
    U_big <- matrix(0, nrow = N, ncol = ncol(U), dimnames = list(ids_all, colnames(U)))
    rr <- match(rownames(U), ids_all)
    U_big[rr, ] <- U
    
     Jp <- if (hessian_of == "neglogLik") (1/N) * H else -(1/N) * H
    
    J_blocks[[p]] <- Jp
    U_blocks[[p]] <- U_big
  }
  
  # concatenate score blocks: U_all = [U_1 | U_2 | ...]
  U_all <- do.call(cbind, U_blocks)
  
  theta_names <- colnames(U_all)
  
  # J
  J <- as.matrix(Matrix::bdiag(lapply(J_blocks, Matrix::Matrix)))
  dimnames(J) <- list(theta_names, theta_names)
  
  # K 
  K <- crossprod(U_all) / N
  dimnames(K) <- list(theta_names, theta_names)
  
  J_inv <- solve(J)
  Sigma <- J_inv %*% K %*% t(J_inv)
  
  V_theta <- Sigma / N
  dimnames(V_theta) <- list(theta_names, theta_names)
  
  list(N = N, J = J, K = K, Sigma = Sigma, vcov_theta = V_theta)
}




########################################################
######                                           #######
##      from  LL' and diag(L))=logdiag(L)             ##
##        to D  for each pair                         ##
######                                           #######
########################################################

                                  
prefix_ex <- function(ex, tag) {
  stopifnot(!is.null(ex$pars), !is.null(ex$H), !is.null(ex$grad))
  
  pfx <- function(nm) paste0(tag, "::", nm)
  
  # pars
  names(ex$pars) <- pfx(names(ex$pars))
  
  # Hessian
  colnames(ex$H) <- pfx(colnames(ex$H))
  rownames(ex$H) <- pfx(rownames(ex$H))
  
  # grad: SOLO le colonne (le righe sono ID)
  colnames(ex$grad) <- pfx(colnames(ex$grad))
  
  ex
}

theta_in_H_order <- function(ex) {
  theta <- ex$pars[colnames(ex$H)]
  stopifnot(identical(names(theta), colnames(ex$H)))
  theta
}
                                  
strip_tag <- function(nm) sub("^.*::", "", nm)

chol_named_to_varcov <- function(chol_vec, chol_names_prefixed) {
  nms0 <- strip_tag(chol_names_prefixed)
  vals <- setNames(as.numeric(chol_vec), nms0)
  
  diag_names <- names(vals)[grepl("^log_chol_var_", names(vals))]
  off_names  <- names(vals)[grepl("^chol_cov_",     names(vals))]
  
  q <- length(diag_names)
  stopifnot(q >= 1)
  
  re_labels <- sub("^log_chol_var_", "", diag_names)
  
  # L con dimnames = re_labels
  L <- matrix(0, q, q, dimnames = list(re_labels, re_labels))
  
  # diag(L)
  for (nm in diag_names) {
    lab <- sub("^log_chol_var_", "", nm)
    L[lab, lab] <- exp(vals[[nm]])
  }
  
  # parser nomi off-diagonal -> due label
  re_from_cov <- function(nm, re_labels) {
    s <- sub("^chol_cov_", "", nm)
    
    # frailty: chol_cov_<item>_frailty_<T>
    if (grepl("_frailty_", s, fixed = TRUE)) {
      parts <- strsplit(s, "_frailty_", fixed = TRUE)[[1]]
      stopifnot(length(parts) == 2)
      item <- parts[1]; when <- parts[2]
      return(c(paste0(item, "_", when), "frailty"))
    }
    
    # cross: chol_cov_cross_<ItemA>_<ItemB>_<TA>_<TB>
    if (startsWith(s, "cross_")) {
      core <- sub("^cross_", "", s)
      toks <- strsplit(core, "_", fixed = TRUE)[[1]]
      stopifnot(length(toks) >= 4)
      tA <- toks[length(toks) - 1]
      tB <- toks[length(toks)]
      mid <- toks[1:(length(toks) - 2)]
      
      for (k in 1:(length(mid) - 1)) {
        itemA <- paste(mid[1:k], collapse = "_")
        itemB <- paste(mid[(k + 1):length(mid)], collapse = "_")
        labA <- paste0(itemA, "_", tA)
        labB <- paste0(itemB, "_", tB)
        if (labA %in% re_labels && labB %in% re_labels) return(c(labA, labB))
      }
      stop("Impossibile parsare cross: ", nm)
    }
    
    # same item different times: chol_cov_<Item>_<T1>_<T2>
    toks <- strsplit(s, "_", fixed = TRUE)[[1]]
    stopifnot(length(toks) >= 3)
    t1 <- toks[length(toks) - 1]
    t2 <- toks[length(toks)]
    item <- paste(toks[1:(length(toks) - 2)], collapse = "_")
    c(paste0(item, "_", t1), paste0(item, "_", t2))
  }
  
  # off-diagonal: L[j,i] con j>i
  for (nm in off_names) {
    labs2 <- re_from_cov(nm, re_labels)
    i1 <- match(labs2[1], re_labels)
    i2 <- match(labs2[2], re_labels)
    stopifnot(!anyNA(c(i1, i2)))
    i <- min(i1, i2); j <- max(i1, i2)
    L[j, i] <- vals[[nm]]
  }
  
  D <- L %*% t(L)
  
  var <- diag(D); names(var) <- paste0("var_", re_labels)
  cov_vals <- c(); cov_nms <- c()
  for (j in 2:q) for (i in 1:(j - 1)) {
    cov_vals <- c(cov_vals, D[j, i])
    cov_nms  <- c(cov_nms, paste0("cov_", re_labels[j], "__", re_labels[i]))
  }
  cov <- setNames(cov_vals, cov_nms)
  
  c(var, cov)
}


#### First delta method from chol+log to var-cov
make_G_pair <- function(ex_pref) {
  theta <- theta_in_H_order(ex_pref)
  
  chol_names <- grep("::(log_chol_var_|chol_cov_)", names(theta), value = TRUE)
  nonchol_names <- setdiff(names(theta), chol_names)
  
  stopifnot(length(chol_names) > 0)
  
  varcov0 <- chol_named_to_varcov(theta[chol_names], chol_names)
  varcov_names <- paste0(tag, "::", names(varcov0))
  
  psi_names <- c(nonchol_names, varcov_names)
  
  gD <- function(x) {
    out <- chol_named_to_varcov(x, chol_names)
    as.numeric(out)
  }
  JD <- jacobian(gD, x = as.numeric(theta[chol_names]))
  rownames(JD) <- names(varcov0)
  colnames(JD) <- strip_tag(chol_names)  
  
  G <- Matrix(0, nrow = length(psi_names), ncol = length(theta),
              dimnames = list(psi_names, names(theta)), sparse = TRUE)
  
  
  
  idx_r <- match(nonchol_names, psi_names)
  idx_c <- match(nonchol_names, names(theta))
  G[cbind(idx_r, idx_c)] <- 1
  
  r_vc <- match(varcov_names, psi_names)
  c_ch <- match(chol_names, names(theta))

  
  G[r_vc, c_ch] <- JD
  
  list(G = G, psi_names = psi_names)
}

predelta_prepare <- function(ex_list, pair_tags, jk_vcov_theta) {
  stopifnot(length(ex_list) == length(pair_tags))

  
  ex_pref <- Map(prefix_ex, ex_list, pair_tags)

  
  ex_pref <- lapply(ex_pref, function(ex) {
    ex$pars <- ex$pars[colnames(ex$H)]
    stopifnot(identical(names(ex$pars), colnames(ex$H)))
    ex
  })
  
  G_list <- vector("list", length(ex_pref))
  psi_list <- vector("list", length(ex_pref))
  
  for (i in seq_along(ex_pref)) {
    tmp <- make_G_pair(ex_pref[[i]])
    G_list[[i]] <- tmp$G
    
    theta <- theta_in_H_order(ex_pref[[i]])
    chol_names <- grep("::(log_chol_var_|chol_cov_)", names(theta), value = TRUE)
    
    nonchol_names <- setdiff(names(theta), chol_names)
    tag <- pair_tags[i]
    
    varcov0 <- chol_named_to_varcov(theta[chol_names], chol_names)
    varcov  <- setNames(as.numeric(varcov0), paste0(tag, "::", names(varcov0)))
    
    psi <- c(theta[nonchol_names], varcov)
    psi_list[[i]] <- psi[tmp$psi_names]
  }


  
  G_all <- as.matrix(Matrix::bdiag(G_list))

  
  psi_all <- unlist(psi_list, use.names = TRUE)

  
  stopifnot(ncol(jk_vcov_theta) == ncol(G_all))
  V_psi <- G_all %*% jk_vcov_theta %*% t(G_all)
  V_psi <- 0.5 * (V_psi + t(V_psi))

  
  colnames(V_psi) <- rownames(V_psi) <- names(psi_all)
  
  list(ex_pref = ex_pref, G_all = G_all, psi_all = psi_all, V_psi = V_psi)
}





###########################################################
####                                                   ####
###                A average matrix                     ###
####                                                   ####
###########################################################



build_A <- function(est_vec,
                    canonicalizer   = identity,
                    element_weights = NULL,
                    make_unique_cols = TRUE,
                    tol_row_sum = 1e-10) {
  if (is.null(names(est_vec))) stop("est_vec deve avere names().")
  raw_names <- names(est_vec)
  
  canon_names <- canonicalizer(raw_names)
  if (anyNA(canon_names)) stop("Ci sono nomi NA dopo canonicalizer().")

  
  levels <- unique(canon_names)

  
  row_index <- match(canon_names, levels)
  col_index <- seq_along(est_vec)

  
  if (is.null(element_weights)) {
    element_weights <- rep(1, length(est_vec))
  } else {
    stopifnot(length(element_weights) == length(est_vec))
    if (any(!is.finite(element_weights)) || any(element_weights < 0))
      stop("element_weights must be finite and >= 0.")
  }

  
  denom <- numeric(length(levels))
  denom_by_row <- tapply(element_weights, row_index, sum)
  denom[as.integer(names(denom_by_row))] <- denom_by_row

  
  
  weights <- element_weights / denom[row_index]
  
  # nomi colonna univoci (consigliato)
  col_names <- if (make_unique_cols) make.unique(raw_names, sep = "__") else raw_names
  
  A <- Matrix::sparseMatrix(
    i = row_index,
    j = col_index,
    x = weights,
    dims = c(length(levels), length(est_vec)),
    dimnames = list(levels, col_names)
  )
  
  rs <- as.numeric(A %*% rep(1, length(est_vec)))
  bad <- which(abs(rs - 1) > tol_row_sum)
  if (length(bad)) {
    warning(sprintf("Unnormalized rows (first 10): %s",
                    paste(head(rownames(A)[bad], 10), collapse = ", ")))
  }
  
  A
}


theta_pair_allpars <- function(ex, do_check = TRUE) {
  p <- ex$pars
  stopifnot(!is.null(names(p)))
  if (!is.numeric(p)) stop("ex$pars must be numeric.")

  
  p_nochol <- p[!grepl("^(log_chol_var_|chol_cov_)", names(p))]

  
  vc <- cholpars_to_D_per_pair(ex, do_check = do_check)$varcov
  if (!is.numeric(vc)) stop("varcov not numeric.")

  
  theta <- c(p_nochol, vc)

  
  if (anyDuplicated(names(theta))) {
    dup <- unique(names(theta)[duplicated(names(theta))])
    stop("Duplicated names for pair: ", paste(dup, collapse=", "))
  }
  
  theta
}

                                  
canonicalizer_all <- function(nms) {
  nms <- sub("^.*::", "", nms)

  is_cov <- grepl("^cov_", nms)
  if (any(is_cov)) {
    for (k in which(is_cov)) {
      xy <- strsplit(sub("^cov_", "", nms[k]), "__", fixed = TRUE)[[1]]
      if (length(xy) == 2) {
        ord <- sort(xy)
        nms[k] <- paste0("cov_", ord[1], "__", ord[2])
      }
    }
  }
  nms
}

##############################################################
# Transformation in the scale of interest for the parameters #
##############################################################

derive_interest <- function(psi) {
  out <- c()
  
  # SD da var_*
  var_names <- grep("^var_", names(psi), value = TRUE)
  for (vn in var_names) {
    v <- psi[[vn]]
    v <- max(v, 0) 
    out[paste0("sd_", sub("^var_", "", vn))] <- sqrt(v)
  }
  
  #corr from cov_*
  cov_names <- grep("^cov_", names(psi), value = TRUE)
  for (cn in cov_names) {
    ab <- strsplit(sub("^cov_", "", cn), "__", fixed = TRUE)[[1]]
    A <- ab[1]; B <- ab[2]
    vA <- psi[[paste0("var_", A)]]
    vB <- psi[[paste0("var_", B)]]
    denom <- sqrt(max(vA,0) * max(vB,0))
    out[paste0("corr_", A, "__", B)] <- if (denom > 0) psi[[cn]] / denom else NA_real_
  }
  
  # Fisher z -> rho for residual correlation
  atanh_names <- grep("^atanh_rho_resid__", names(psi), value = TRUE)
  for (nm in atanh_names) {
    pair <- sub("^atanh_rho_resid__", "", nm)
    out[paste0("corr_resid__", pair)] <- tanh(psi[[nm]])
  }
  
  
  #Weibull
  if ("log_lambda" %in% names(psi)) out["lambda"] <- exp(psi[["log_lambda"]])
  if ("log_rho"    %in% names(psi)) out["rho"]    <- exp(psi[["log_rho"]])
  
  # HR- betas survival 
  bsurv <- grep("^fixed_surv_", names(psi), value = TRUE)
  for (bn in bsurv) out[paste0("HR_", sub("^fixed_surv_", "", bn))] <- exp(psi[[bn]])
  
  # precision for Beta 
  blogphi <- grep("^beta_log_precision_", names(psi), value = TRUE)
  for (pn in blogphi) out[paste0("phi_", sub("^beta_log_precision_", "", pn))] <- exp(psi[[pn]])
  
  logsig_names <- grep("^log_sigma_", names(psi), value = TRUE)
  for (sn in logsig_names) {
    out[paste0("sigma_", sub("^log_sigma_", "", sn))] <- exp(psi[[sn]])
  }
  
  out
}

# Jacobian numeric
g <- function(x) {
  tmp <- psi_global
  tmp[names(x)] <- x
  as.numeric(derive_interest(tmp))
}

#### Reconstruct vector with trasfomed and untrasfomed estimates with corresponding VCOV matrix
build_interest_full <- function(psi_global, V_global,
                                est_trans, J_trans,
                                keep_untransformed = NULL,
                                drop_untransformed_regex = NULL,
                                out_order = NULL,
                                symmetrize = TRUE) {
  stopifnot(is.numeric(psi_global), !is.null(names(psi_global)))
  
  V_global <- as.matrix(V_global)
  
  stopifnot(identical(rownames(V_global), names(psi_global)))
  stopifnot(identical(colnames(V_global), names(psi_global)))
  
  stopifnot(is.numeric(est_trans), !is.null(names(est_trans)))
  
  J_trans <- as.matrix(J_trans)
  stopifnot(is.matrix(J_trans))
  stopifnot(identical(colnames(J_trans), names(psi_global)))
  
  trans_names <- intersect(names(est_trans), rownames(J_trans))
 
  missing_in_J <- setdiff(names(est_trans), rownames(J_trans))

  est_trans <- est_trans[trans_names]
  J_trans   <- J_trans[trans_names, , drop = FALSE]

  
  if (is.null(keep_untransformed)) keep_untransformed <- names(psi_global)
  keep_untransformed <- intersect(keep_untransformed, names(psi_global))

  if (!is.null(drop_untransformed_regex)) {
    keep_untransformed <- keep_untransformed[!grepl(drop_untransformed_regex, keep_untransformed)]
  }

  keep_untransformed <- setdiff(keep_untransformed, trans_names)
  

  out_names <- c(keep_untransformed, trans_names)
  
  if (!is.null(out_order)) {
    miss <- setdiff(out_order, out_names)
    if (length(miss)) stop("out_order no names available: ", paste(miss, collapse=", "))
    out_names <- out_order
  }

  
  est_full <- setNames(rep(NA_real_, length(out_names)), out_names)
  est_full[keep_untransformed] <- psi_global[keep_untransformed]
  est_full[trans_names]        <- est_trans[trans_names]

  
  J_full <- matrix(0, nrow = length(out_names), ncol = length(psi_global),
                   dimnames = list(out_names, names(psi_global)))

  
  idx_r <- match(keep_untransformed, out_names)
  idx_c <- match(keep_untransformed, names(psi_global))
  J_full[cbind(idx_r, idx_c)] <- 1

  
  J_full[trans_names, ] <- J_trans[trans_names, , drop = FALSE]
  

  
  V_full <- J_full %*% V_global %*% t(J_full)
  if (isTRUE(symmetrize)) V_full <- 0.5 * (V_full + t(V_full))
  
  se_full <- sqrt(pmax(0, diag(V_full)))
  names(se_full) <- out_names
  
  list(est = est_full, V = V_full, se = se_full, J = J_full)
}



#####
# Var covar radnom effects
#################################
# Build D (var-cov matrix) from psi (var_*, cov_*)
#########################################
re_varcov_matrix2 <- function(psi,
                              strip_tag = TRUE,
                              fill_missing = "NA",
                              order = c("appearance", "sorted"),
                              label_filter = NULL,   # es: "(_I0|_Ipost|frailty)$"
                              make_psd = FALSE) {
  fill_missing <- match.arg(fill_missing)
  order <- match.arg(order)
  
  if (is.null(names(psi))) stop("psi deve avere names().")
  
  nms <- names(psi)
  if (strip_tag) nms <- sub("^.*::", "", nms)
  psi2 <- setNames(as.numeric(psi), nms)
  
  var_nms <- grep("^var_", names(psi2), value = TRUE)
  cov_nms <- grep("^cov_", names(psi2), value = TRUE)
  
  if (length(var_nms) == 0 && length(cov_nms) == 0)
    stop("Nessun var_* o cov_* trovato in psi.")

  
  labs_var <- sub("^var_", "", var_nms)

  
  cov_pairs <- lapply(cov_nms, function(cn) {
    s <- sub("^cov_", "", cn)
    parts <- strsplit(s, "__", fixed = TRUE)[[1]]
    if (length(parts) < 2) return(NULL)
    A <- parts[1]
    B <- paste(parts[-1], collapse = "__")
    c(A = A, B = B)
  })
  cov_pairs <- Filter(Negate(is.null), cov_pairs)
  labs_cov <- unique(unlist(lapply(cov_pairs, unname)))

  
  labs <- c(labs_var, labs_cov)
  if (order == "appearance") {
    labs <- unique(labs)
  } else {
    labs <- sort(unique(labs))
  }
  

  if (!is.null(label_filter)) {
    keep <- grepl(label_filter, labs)
    labs <- labs[keep]
  }
  
  q <- length(labs)
  if (q == 0) stop("Dopo label_filter non resta nessuna dimensione.")
  
  
  D <- matrix(NA, q, q, dimnames = list(labs, labs))
  
  #set diagonale (var)
  missing_var <- character(0)
  for (lab in labs) {
    vn <- paste0("var_", lab)
    if (vn %in% names(psi2)) {
      D[lab, lab] <- psi2[[vn]]
    } else {
      missing_var <- c(missing_var, vn)
    }
  }
  
  #set off-diagonal (cov)
  unused_cov <- NA
  for (k in seq_along(cov_pairs)) {
    A <- cov_pairs[[k]]["A"]
    B <- cov_pairs[[k]]["B"]
    cn <- cov_nms[k]
    
    if (!(A %in% labs && B %in% labs)) {
      unused_cov <- c(unused_cov, cn)
      next
    }
    
    val <- psi2[[cn]]
    D[A, B] <- val
    D[B, A] <- val
  }
  
  # simmetrizza
  # D <- 0.5 * (D + t(D))
  
  # # opzional PSD
  # if (make_psd) {
  #   if (!requireNamespace("Matrix", quietly = TRUE))
  #     stop("Per make_psd=TRUE serve il pacchetto Matrix.")
  #   D <- as.matrix(Matrix::nearPD(D, corr = FALSE)$mat)
  #   dimnames(D) <- list(labs, labs)
  # }
  
  # correlazioni
  sdv <- sqrt(pmax(diag(D), 0))
  Cor <- D / (sdv %o% sdv)
  diag(Cor) <- 1
  
  list(
    D = D,
    Cor = Cor,
    missing_var = missing_var,   
    unused_cov = unused_cov    
  )
}


### residual error correlation matrix
corr_resid_matrix_simple <- function(psi, base = "corr_resid", sep = "__") {
  if (is.null(names(psi))) stop("psi deve avere names().")
  
  nms <- names(psi)
  keep <- startsWith(nms, base)
  if (!any(keep)) stop("Nessun parametro che inizi con '", base, "'.")
  
  v <- as.numeric(psi[keep])
  nm <- names(psi[keep])

  
  s <- sub(paste0("^", base, "_+"), "", nm)
  
  parts <- strsplit(s, sep, fixed = TRUE)
  if (any(lengths(parts) != 2)) {
    bad <- nm[lengths(parts) != 2]
    stop("Nomi non parsabili (atteso: ", base, "__A__B): ", paste(bad, collapse = ", "))
  }
  
  A <- vapply(parts, `[`, character(1), 1)
  B <- vapply(parts, `[`, character(1), 2)
  
  items <- unique(c(A, B))
  R <- matrix(NA_real_, length(items), length(items), dimnames = list(items, items))
  diag(R) <- 1
  
  for (k in seq_along(v)) {
    R[A[k], B[k]] <- v[k]
    R[B[k], A[k]] <- v[k]
  }
  
  R
}

### update correlations and sd in D after adjustment for SPD 
update_est_interest_from_D <- function(est_interest, Dadj, Corradj = NULL,
                                       representation = c("corr","cov"),
                                       prefix_sd = "sd_", prefix_corr = "corr_",
                                       prefix_var = "var_", prefix_cov = "cov_",
                                       var_names = NULL, cov_names = NULL,
                                       add_missing = FALSE,
                                       keep_if_missing = TRUE) {
  representation <- match.arg(representation)
  
  if (is.null(names(est_interest))) stop("est_interest deve essere un vettore nominato.")
  if (is.null(rownames(Dadj)) || is.null(colnames(Dadj))) {
    stop("Dadj deve avere dimnames coerenti (rownames/colnames).")
  }
  
  out <- est_interest
  
  ensure_names <- function(x, nm) {
    if (!isTRUE(add_missing) || is.null(nm)) return(x)
    miss <- setdiff(nm, names(x))
    if (length(miss)) x <- c(x, setNames(rep(NA_real_, length(miss)), miss))
    x
  }

  
  if (representation == "corr") {
    
    if (is.null(Corradj)) Corradj <- cov2cor(Dadj)
    
    nms <- names(out)
    
    # SD
    is_sd <- startsWith(nms, prefix_sd)
    if (any(is_sd)) {
      idx  <- which(is_sd)
      labs <- sub(paste0("^", prefix_sd), "", nms[idx])
      ok   <- labs %in% rownames(Dadj)
      
      vals <- rep(NA_real_, length(labs))
      vals[ok] <- sqrt(pmax(diag(Dadj)[labs[ok]], 0))
      
      if (!keep_if_missing) out[idx[!ok]] <- NA_real_
      out[idx[ok]] <- vals[ok]
    }
    
    # CORR
    nms <- names(out)
    is_corr <- startsWith(nms, prefix_corr)
    if (any(is_corr)) {
      idx <- which(is_corr)
      labs12 <- sub(paste0("^", prefix_corr), "", nms[idx])
      spl <- strsplit(labs12, "__", fixed = TRUE)
      
      vals <- rep(NA_real_, length(spl))
      ok <- logical(length(spl))
      
      for (i in seq_along(spl)) {
        if (length(spl[[i]]) != 2) next
        a <- spl[[i]][1]; b <- spl[[i]][2]
        if (a %in% rownames(Corradj) && b %in% colnames(Corradj)) {
          vals[i] <- Corradj[a, b]; ok[i] <- TRUE
        } else if (b %in% rownames(Corradj) && a %in% colnames(Corradj)) {
          vals[i] <- Corradj[b, a]; ok[i] <- TRUE
        }
      }
      
      if (!keep_if_missing) out[idx[!ok]] <- NA_real_
      out[idx[ok]] <- vals[ok]
    }
    
    return(out)
  }

  
  if (representation == "cov") {

    
    out <- ensure_names(out, var_names)
    out <- ensure_names(out, cov_names)
    
    nms <- names(out)
    
    # VAR
    is_var <- startsWith(nms, prefix_var)
    if (any(is_var)) {
      idx <- which(is_var)
      labs <- sub(paste0("^", prefix_var), "", nms[idx])
      ok <- labs %in% rownames(Dadj)
      
      vals <- rep(NA_real_, length(labs))
      vals[ok] <- diag(Dadj)[labs[ok]]
      
      if (!keep_if_missing) out[idx[!ok]] <- NA_real_
      out[idx[ok]] <- vals[ok]
    }
    
    # COV
    nms <- names(out)
    is_cov <- startsWith(nms, prefix_cov)
    if (any(is_cov)) {
      idx <- which(is_cov)
      labs12 <- sub(paste0("^", prefix_cov), "", nms[idx])
      spl <- strsplit(labs12, "__", fixed = TRUE)
      
      vals <- rep(NA_real_, length(spl))
      ok <- logical(length(spl))
      
      for (i in seq_along(spl)) {
        if (length(spl[[i]]) != 2) next
        a <- spl[[i]][1]; b <- spl[[i]][2]
        if (a %in% rownames(Dadj) && b %in% colnames(Dadj)) {
          vals[i] <- Dadj[a, b]; ok[i] <- TRUE
        } else if (b %in% rownames(Dadj) && a %in% colnames(Dadj)) {
          vals[i] <- Dadj[b, a]; ok[i] <- TRUE
        }
      }
      
      if (!keep_if_missing) out[idx[!ok]] <- NA_real_
      out[idx[ok]] <- vals[ok]
    }
    
    return(out)
  }
}








                                  
