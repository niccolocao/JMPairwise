

# latex util
latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%#$&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x
}

############# wald for all + FDR for correlations only
wald_table_raw <- function(est, se,
                           alpha = 0.05,
                           null = NULL,       
                           bounds_corr = c(-1, 1),
                           set_p_na_if_out_of_bounds = FALSE,
                           add_fdr_corr = TRUE,
                           fdr_method = "fdr",
                           sort_by = c("p", "p_fdr_corr", "none")) {
  
  sort_by <- match.arg(sort_by)
  
  stopifnot(is.numeric(est), is.numeric(se))
  if (is.null(names(est)) || is.null(names(se)))
    stop("est e se devono essere *named vectors*.")
  
  common <- intersect(names(est), names(se))
  if (!length(common)) stop("Nessun nome in comune tra est e se.")
  
  est <- est[common]
  se  <- se[common]
  
  # type solo da prefisso (non parsing)
  type <- ifelse(grepl("^corr_", common), "corr",
                 ifelse(grepl("^sd_",   common), "sd",
                        ifelse(grepl("^HR_",   common), "HR", "other")))
  
  # H0
  default_h0 <- ifelse(type == "HR", 1, 0)
  
  if (is.null(null)) {
    h0 <- default_h0
  } else if (is.function(null)) {
    h0 <- mapply(null, param = common, type = type)
  } else if (length(null) == 1 && is.numeric(null)) {
    h0 <- rep(null, length(common))
  } else {
    if (is.null(names(null))) stop("Se 'null' non è scalare deve essere un *named vector*.")
    h0 <- null[common]
    miss <- is.na(h0)
    if (any(miss)) h0[miss] <- default_h0[miss]
  }
  
  z <- (as.numeric(est) - as.numeric(h0)) / as.numeric(se)
  p <- 2 * pnorm(-abs(z))
  
  tab <- data.frame(
    param = common,
    type  = type,
    est   = as.numeric(est),
    se    = as.numeric(se),
    null  = as.numeric(h0),
    z     = as.numeric(z),
    p     = as.numeric(p),
    sig_0.05 = as.numeric(p) < alpha,
    stringsAsFactors = FALSE
  )
  
  tab$out_of_bounds <- with(tab, type == "corr" & (est < bounds_corr[1] | est > bounds_corr[2]))
  
  if (set_p_na_if_out_of_bounds) {
    idx <- which(tab$out_of_bounds)
    tab$z[idx] <- NA_real_
    tab$p[idx] <- NA_real_
    tab$sig_0.05[idx] <- NA
  }
  
  # FDR solo corr
  tab$p_fdr_corr <- NA_real_
  tab$sig_fdr_0.05 <- NA
  if (add_fdr_corr) {
    idxc <- which(tab$type == "corr" & is.finite(tab$p))
    tab$p_fdr_corr[idxc] <- p.adjust(tab$p[idxc], method = fdr_method)
    tab$sig_fdr_0.05[idxc] <- tab$p_fdr_corr[idxc] < alpha
  }
  
  if (sort_by == "p") {
    tab <- tab[order(tab$p), ]
  } else if (sort_by == "p_fdr_corr") {
    # mette prima corr con FDR definito
    key <- ifelse(is.na(tab$p_fdr_corr), Inf, tab$p_fdr_corr)
    tab <- tab[order(key), ]
  }
  
  rownames(tab) <- NULL
  tab
}

wald_fisher_from_r <- function(est_r, se_r,
                               alpha = 0.05,
                               eps = 1e-6,
                               add_fdr = TRUE,
                               fdr_method = "fdr") {
  stopifnot(is.numeric(est_r), is.numeric(se_r),
            !is.null(names(est_r)), !is.null(names(se_r)))
  
  nm <- intersect(names(est_r), names(se_r))
  if (!length(nm)) stop("nessun nome in comune tra est_r e se_r")
  
  r  <- est_r[nm]
  sr <- se_r[nm]
  
  r_clamp <- pmin(pmax(r, -1 + eps), 1 - eps)
  
  zhat <- atanh(r_clamp)
  sez  <- sr / (1 - r_clamp^2)
  
  w <- zhat / sez
  p <- 2 * pnorm(-abs(w))
  
  out <- data.frame(
    param = nm,
    r = as.numeric(r),
    se_r = as.numeric(sr),
    fisher_z = as.numeric(zhat),
    se_fisher_z = as.numeric(sez),
    wald_z = as.numeric(w),
    p = as.numeric(p),
    sig_0.05 = as.numeric(p) < alpha,
    stringsAsFactors = FALSE
  )
  
  if (isTRUE(add_fdr)) {
    out$p_fdr <- p.adjust(out$p, method = fdr_method)
    out$sig_fdr_0.05 <- out$p_fdr < alpha
  }
  
  out[order(out$p), ]
}

get_name_anyorder <- function(vnames, prefix, A, B) {
  n1 <- paste0(prefix, A, "__", B)
  n2 <- paste0(prefix, B, "__", A)
  if (n1 %in% vnames) return(n1)
  if (n2 %in% vnames) return(n2)
  NA_character_
}

corrdiff_frailty_table_from_interest <- function(est_interest_adj, V_interest,
                                                 items = NULL,
                                                 t_pre = "I0",
                                                 t_post = "Ipost",
                                                 frailty = "frailty",
                                                 alpha = 0.05,
                                                 add_fdr = TRUE,
                                                 fdr_method = "fdr",
                                                 sort_by = c("p","p_fdr","none")) {
  sort_by <- match.arg(sort_by)
  
  stopifnot(is.numeric(est_interest_adj), !is.null(names(est_interest_adj)))
  V_interest <- as.matrix(V_interest)
  stopifnot(identical(rownames(V_interest), colnames(V_interest)))
  
  vnm <- names(est_interest_adj)
  
  if (is.null(items)) {
    sdp  <- grep(paste0("^sd_(.+)_", t_pre,  "$"),  vnm, value = TRUE)
    sdpo <- grep(paste0("^sd_(.+)_", t_post, "$"), vnm, value = TRUE)
    ipre  <- sub(paste0("^sd_(.+)_", t_pre,  "$"), "\\1", sdp)
    ipost <- sub(paste0("^sd_(.+)_", t_post, "$"), "\\1", sdpo)
    items <- sort(intersect(ipre, ipost))
  }
  
  one_item <- function(it) {
    lab_pre  <- paste0(it, "_", t_pre)
    lab_post <- paste0(it, "_", t_post)
    
    sdP_nm <- paste0("sd_", lab_post)
    sd0_nm <- paste0("sd_", lab_pre)
    
    rPf_nm <- get_name_anyorder(vnm, "corr_", frailty, lab_post)
    r0f_nm <- get_name_anyorder(vnm, "corr_", frailty, lab_pre)
    rP0_nm <- get_name_anyorder(vnm, "corr_", lab_post, lab_pre)
    
    need <- c(sdP_nm, sd0_nm, rPf_nm, r0f_nm, rP0_nm)
    
    if (anyNA(need) || any(!need %in% vnm) || any(!need %in% rownames(V_interest))) {
      return(data.frame(
        param = paste0("corrdiff_frailty__", it, "__", t_post, "_minus_", t_pre),
        type="corrdiff", est=NA_real_, se=NA_real_, null=0,
        z=NA_real_, p=NA_real_, sig_0.05=NA,
        stringsAsFactors=FALSE
      ))
    }
    
    sP  <- as.numeric(est_interest_adj[[sdP_nm]])
    s0  <- as.numeric(est_interest_adj[[sd0_nm]])
    rPf <- as.numeric(est_interest_adj[[rPf_nm]])
    r0f <- as.numeric(est_interest_adj[[r0f_nm]])
    rP0 <- as.numeric(est_interest_adj[[rP0_nm]])
    
    b <- sP^2 + s0^2 - 2*rP0*sP*s0
    d <- rPf*sP - r0f*s0
    
    est <- if (is.finite(b) && b > 0) as.numeric(d / sqrt(b)) else NA_real_
    
    se <- NA_real_
    if (is.finite(est) && is.finite(b) && b > 0) {
      g_sdP <- rPf/sqrt(b) - d*(sP - rP0*s0)/(b^(3/2))
      g_sd0 <- -r0f/sqrt(b) - d*(s0 - rP0*sP)/(b^(3/2))
      g_rPf <- sP/sqrt(b)
      g_r0f <- -s0/sqrt(b)
      g_rP0 <- d*(sP*s0)/(b^(3/2))
      
      gvec <- c(g_sdP, g_sd0, g_rPf, g_r0f, g_rP0)
      names(gvec) <- need
      
      Vsub <- V_interest[need, need, drop = FALSE]
      v <- as.numeric(t(gvec) %*% Vsub %*% gvec)
      se <- sqrt(max(v, 0))
    }
    
    z <- if (is.finite(est) && is.finite(se) && se > 0) est/se else NA_real_
    p <- if (is.finite(z)) 2*pnorm(-abs(z)) else NA_real_
    
    data.frame(
      param = paste0("corrdiff_frailty__", it, "__", t_post, "_minus_", t_pre),
      type  = "corrdiff",
      est   = est,
      se    = se,
      null  = 0,
      z     = z,
      p     = p,
      sig_0.05 = if (is.finite(p)) p < alpha else NA,
      stringsAsFactors = FALSE
    )
  }
  
  tab <- do.call(rbind, lapply(items, one_item))
  
  tab$p_fdr <- NA_real_
  tab$sig_fdr_0.05 <- NA
  if (isTRUE(add_fdr)) {
    ok <- which(is.finite(tab$p))
    tab$p_fdr[ok] <- p.adjust(tab$p[ok], method = fdr_method)
    tab$sig_fdr_0.05[ok] <- tab$p_fdr[ok] < alpha
  }
  
  if (sort_by == "p") tab <- tab[order(tab$p), ]
  if (sort_by == "p_fdr") {
    key <- ifelse(is.na(tab$p_fdr), Inf, tab$p_fdr)
    tab <- tab[order(key), ]
  }
  
  rownames(tab) <- NULL
  tab
}


########################################################
##########       Plots for correlations      ###########
########################################################
library(RColorBrewer)

scale_corr_colors <- function() {
  scale_color_gradientn(
    colours = rev(brewer.pal(11, "RdBu")),
    limits = c(-1, 1),
    name = "$r$"
  )
}

make_pretty_label_vec <- function(x, map_base = short_map_base) {
  x <- as.character(x)
  x[x == "frailty"] <- "Frailty"
  base <- sub("_(I0|Ipost)$", "", x)
  
  out <- unname(map_base[base])
  out[is.na(out)] <- base[is.na(out)]
  
  out
}

update_corr_est_from_matrix <- function(est, R, prefix = "corr_") {
  stopifnot(is.numeric(est), !is.null(names(est)))
  stopifnot(is.matrix(R), !is.null(rownames(R)), !is.null(colnames(R)))
  
  out <- est
  nm <- names(out)
  
  if (prefix == "corr_") {
    target <- nm[startsWith(nm, "corr_") & !startsWith(nm, "corr_resid__")]
  } else {
    target <- nm[startsWith(nm, prefix)]
  }
  
  for (p in target) {
    s <- sub(paste0("^", prefix), "", p)
    ab <- strsplit(s, "__", fixed = TRUE)[[1]]
    if (length(ab) != 2L) next
    
    a <- ab[1]
    b <- ab[2]
    
    if (a %in% rownames(R) && b %in% colnames(R)) {
      out[p] <- R[a, b]
    } else if (b %in% rownames(R) && a %in% colnames(R)) {
      out[p] <- R[b, a]
    }
  }
  
  out
}

make_corr_table_fisher_fdr <- function(est, se,
                                       family = c("all","std","resid"),
                                       alpha = 0.05,
                                       fdr_method = "fdr",
                                       eps = 1e-6) {
  family <- match.arg(family)
  
  stopifnot(is.numeric(est), is.numeric(se),
            !is.null(names(est)), !is.null(names(se)))
  
  nm <- intersect(names(est), names(se))
  est <- est[nm]
  se  <- se[nm]
  
  cand_std <- nm[startsWith(nm, "corr_") & !startsWith(nm, "corr_resid__")]
  cand_res <- nm[startsWith(nm, "corr_resid__")]
  
  pars <- switch(
    family,
    all = c(cand_std, cand_res),
    std = cand_std,
    resid = cand_res
  )
  
  if (!length(pars)) return(data.frame())
  
  fam <- ifelse(startsWith(pars, "corr_resid__"), "resid", "std")
  raw <- sub("^corr_resid__|^corr_", "", pars)
  spl <- strsplit(raw, "__", fixed = TRUE)
  
  ok <- lengths(spl) == 2L
  pars <- pars[ok]
  fam  <- fam[ok]
  spl  <- spl[ok]
  
  A <- vapply(spl, `[`, "", 1)
  B <- vapply(spl, `[`, "", 2)
  
  r  <- as.numeric(est[pars])
  sr <- as.numeric(se[pars])
  
  valid <- is.finite(r) & is.finite(sr) & sr > 0 & abs(r) < 1
  r_clamp <- pmin(pmax(r, -1 + eps), 1 - eps)
  
  fisher_z <- rep(NA_real_, length(pars))
  se_fisher_z <- rep(NA_real_, length(pars))
  z <- rep(NA_real_, length(pars))
  p <- rep(NA_real_, length(pars))
  
  fisher_z[valid] <- atanh(r_clamp[valid])
  se_fisher_z[valid] <- sr[valid] / (1 - r_clamp[valid]^2)
  z[valid] <- fisher_z[valid] / se_fisher_z[valid]
  p[valid] <- 2 * pnorm(-abs(z[valid]))
  
  p_fdr <- rep(NA_real_, length(pars))
  if (any(valid)) p_fdr[valid] <- p.adjust(p[valid], method = fdr_method)
  
  data.frame(
    param = pars,
    type = "corr",
    corr_family = fam,
    A = A,
    B = B,
    est = r,
    se = sr,
    fisher_z = fisher_z,
    se_fisher_z = se_fisher_z,
    z = z,
    p = p,
    p_fdr = p_fdr,
    sig_0.05 = !is.na(p) & p < alpha,
    sig_fdr_0.05 = !is.na(p_fdr) & p_fdr < alpha,
    out_of_bounds = !valid,
    stringsAsFactors = FALSE
  )
}

parse_corr_tab <- function(tab, which_corr = c("std", "resid", "all")) {
  which_corr <- match.arg(which_corr)
  
  tt <- tab[tab$type == "corr", , drop = FALSE]
  if (!nrow(tt)) return(tt)
  
  if (!all(c("A", "B") %in% names(tt))) {
    raw <- as.character(tt$param)
    
    fam <- ifelse(
      grepl("^corr_resid__", raw),
      "resid",
      ifelse(grepl("^corr_", raw), "std", NA_character_)
    )
    
    ab <- sub("^corr_resid__|^corr_", "", raw)
    ss <- strsplit(ab, "__", fixed = TRUE)
    ok <- lengths(ss) == 2L

     
    tt  <- tt[ok, , drop = FALSE]
    fam <- fam[ok]
    ss  <- ss[ok]
    
    tt$corr_family <- fam
    tt$A <- vapply(ss, `[`, character(1), 1)
    tt$B <- vapply(ss, `[`, character(1), 2)
  } else {
    raw <- as.character(tt$param)
    tt$corr_family <- ifelse(grepl("^corr_resid__", raw), "resid", "std")
  }
  
  if (which_corr == "std") {
    tt <- tt[tt$corr_family == "std", , drop = FALSE]
  } else if (which_corr == "resid") {
    tt <- tt[tt$corr_family == "resid", , drop = FALSE]
  }
  
  tt
}


make_RP_from_tab <- function(tab, labs, p_col = c("p_fisher","p"),
                             set_oob_to_NA = TRUE,
                             which_corr = c("std", "resid", "all")) {
  
  which_corr <- match.arg(which_corr)
  
  p_col <- p_col[p_col %in% names(tab)][1]
  if (length(p_col) == 0 || is.na(p_col)) stop("tab non ha né p_fisher né p.")
  
  n <- length(labs)
  R <- matrix(NA_real_, n, n, dimnames = list(labs, labs))
  diag(R) <- 1
  
  P <- matrix(NA_real_, n, n, dimnames = list(labs, labs))
  diag(P) <- NA_real_
  
  tt <- parse_corr_tab(tab, which_corr = which_corr)
  if (!nrow(tt)) return(list(R = R, P = P, tab_used = tt))
  
  keep <- tt$A %in% labs & tt$B %in% labs
  tt <- tt[keep, , drop = FALSE]
  
  if (set_oob_to_NA && "out_of_bounds" %in% names(tt)) {
    tt$est[tt$out_of_bounds] <- NA_real_
    tt[[p_col]][tt$out_of_bounds] <- NA_real_
  } else {
    bad <- !is.na(tt$est) & abs(tt$est) > 1
    tt$est[bad] <- NA_real_
    tt[[p_col]][bad] <- NA_real_
  }
  
  for (k in seq_len(nrow(tt))) {
    a <- tt$A[k]
    b <- tt$B[k]
    R[a, b] <- R[b, a] <- tt$est[k]
    P[a, b] <- P[b, a] <- tt[[p_col]][k]
  }
  
  list(R = R, P = P, tab_used = tt)
}


make_cross_RP_from_tab <- function(tab, rows, cols, p_col = c("p_fisher","p"),
                                   set_oob_to_NA = TRUE,
                                   which_corr = c("std", "resid", "all")) {
  
  which_corr <- match.arg(which_corr)
  
  p_col <- p_col[p_col %in% names(tab)][1]
  if (length(p_col) == 0 || is.na(p_col)) stop("tab non ha né p_fisher né p.")
  
  R <- matrix(NA_real_, length(rows), length(cols),
              dimnames = list(rows, cols))
  
  P <- matrix(NA_real_, length(rows), length(cols),
              dimnames = list(rows, cols))
  
  tt <- parse_corr_tab(tab, which_corr = which_corr)
  if (!nrow(tt)) return(list(R = R, P = P, tab_used = tt))
  
  keep_dir <- tt$A %in% rows & tt$B %in% cols
  keep_rev <- tt$A %in% cols & tt$B %in% rows
  tt <- tt[keep_dir | keep_rev, , drop = FALSE]
  
  if (set_oob_to_NA && "out_of_bounds" %in% names(tt)) {
    tt$est[tt$out_of_bounds] <- NA_real_
    tt[[p_col]][tt$out_of_bounds] <- NA_real_
  } else {
    bad <- !is.na(tt$est) & abs(tt$est) > 1
    tt$est[bad] <- NA_real_
    tt[[p_col]][bad] <- NA_real_
  }
  
  for (k in seq_len(nrow(tt))) {
    a <- tt$A[k]
    b <- tt$B[k]
    r <- tt$est[k]
    p <- tt[[p_col]][k]
    
    if (a %in% rows && b %in% cols) {
      R[a, b] <- r
      P[a, b] <- p
    } else if (a %in% cols && b %in% rows) {
      R[b, a] <- r
      P[b, a] <- p
    }
  }
  
  list(R = R, P = P, tab_used = tt)
}

library(ggplot2)
library(RColorBrewer)

scale_corr_colors <- function() {
  scale_color_gradientn(
    colours = rev(brewer.pal(11, "RdBu")),
    limits = c(-1, 1),
    name = "$r$"
  )
}

theme_corr_text <- function(
    axis_x_size = 9,
    axis_y_size = 9,
    axis_title_x_size = 11,
    axis_title_y_size = 11,
    plot_title_size = 12,
    legend_title_size = 10,
    legend_text_size = 9
) {
  theme(
    axis.text.x = element_text(size = axis_x_size, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = axis_y_size),
    axis.text.y.right = element_text(size = axis_y_size),
    axis.title.x = element_text(size = axis_title_x_size),
    axis.title.y = element_text(size = axis_title_y_size),
    axis.title.y.right = element_text(size = axis_title_y_size),
    plot.title = element_text(size = plot_title_size),
    legend.title = element_text(size = legend_title_size),
    legend.text = element_text(size = legend_text_size)
  )
}

plot_corr_circles_lower <- function(R, P = NULL, alpha = 0.05,
                                    sd_vec = NULL, title = NULL,
                                    circle_max = 12, size_gamma = 1.7,
                                    digits_sd = 2, digits_r = 2,
                                    label_fun = make_pretty_label_vec,
                                    grid_color = "grey85",
                                    r_text_size = 2.8, sd_text_size = 3.0,
                                    block_split = 5,
                                    eq_label = "EQ-5D",
                                    koos_label = "KOOS",
                                    block_text_size = 5,
                                    mask_block_upper = TRUE,
                                    block_line_color = "grey40",
                                    block_line_width = 0.6,
                                    size_breaks = c(0.2, 0.5, 0.8)) {
  stopifnot(is.matrix(R))
  rn <- rownames(R)
  cn <- colnames(R)
  stopifnot(identical(rn, cn))
  n <- length(rn)
  
  if (!is.null(P)) {
    stopifnot(
      is.matrix(P),
      identical(rownames(P), rn),
      identical(colnames(P), cn)
    )
  }
  
  df <- expand.grid(row = rn, col = cn, stringsAsFactors = FALSE)
  df$i <- match(df$row, rn)
  df$j <- match(df$col, cn)
  df$r <- as.numeric(R[cbind(df$i, df$j)])
  df$p <- if (!is.null(P)) as.numeric(P[cbind(df$i, df$j)]) else NA_real_
  
  df$is_lower <- df$i >= df$j
  df$is_diag <- df$i == df$j
  df <- df[df$is_lower, , drop = FALSE]
  
  df$x <- df$j
  df$y <- n - df$i + 1
  
  diag_df <- df[df$is_diag, , drop = FALSE]
  if (!is.null(sd_vec)) {
    diag_df$sd <- sd_vec[diag_df$row]
    diag_df$sd_lab <- ifelse(
      is.na(diag_df$sd), "",
      formatC(diag_df$sd, digits = digits_sd, format = "f")
    )
  } else {
    diag_df$sd_lab <- ""
  }
  
  off_df <- df[!df$is_diag & is.finite(df$r), , drop = FALSE]
  off_df$absr <- abs(off_df$r)
  off_df$size_val <- off_df$absr^size_gamma
  off_df$star <- ifelse(!is.na(off_df$p) & off_df$p < alpha, "*", "")
  off_df$r_lab <- paste0(formatC(off_df$r, digits = digits_r, format = "f"), off_df$star)
  
  grid_df <- unique(df[, c("x", "y")])
  
  size_breaks <- sort(size_breaks)
  size_breaks <- size_breaks[size_breaks > 0 & size_breaks <= 1]
  size_breaks_trans <- size_breaks^size_gamma
  
  x_labs <- label_fun(cn)
  y_labs <- label_fun(rev(rn))
  
  has_blocks <- !is.null(block_split) && is.numeric(block_split) &&
    block_split >= 2 && block_split < n
  m <- if (has_blocks) as.integer(block_split) else NA_integer_
  
  if (has_blocks) {
    eq_tx <- (0.5 + (m + 0.5) + (m + 0.5)) / 3
    eq_ty <- ((n + 0.5) + (n + 0.5) + (n - m + 0.5)) / 3
    
    koos_tx <- ((m + 0.5) + (n + 0.5) + (n + 0.5)) / 3
    koos_ty <- ((n - m + 0.5) + (n - m + 0.5) + 0.5) / 3
  }
  
  txt_sz_off <- r_text_size * min(1, 12 / n)
  txt_sz_off <- max(1.2, txt_sz_off)
  
  txt_sz_diag <- sd_text_size * min(1, 12 / n)
  txt_sz_diag <- max(1.2, txt_sz_diag)
  
  p <- ggplot() +
    {if (has_blocks && isTRUE(mask_block_upper))
      annotate(
        "rect",
        xmin = 0.5, xmax = m + 0.5,
        ymin = n - m + 0.5, ymax = n + 0.5,
        fill = "white", color = NA
      )
    } +
    {if (has_blocks && isTRUE(mask_block_upper))
      annotate(
        "rect",
        xmin = m + 0.5, xmax = n + 0.5,
        ymin = 0.5, ymax = n - m + 0.5,
        fill = "white", color = NA
      )
    } +
    {if (has_blocks && !is.null(eq_label) && nzchar(eq_label))
      annotate(
        "text",
        x = eq_tx, y = eq_ty,
        label = eq_label, angle = -45,
        fontface = "bold", size = block_text_size, colour = "grey30"
      )
    } +
    {if (has_blocks && !is.null(koos_label) && nzchar(koos_label))
      annotate(
        "text",
        x = koos_tx, y = koos_ty,
        label = koos_label, angle = -45,
        fontface = "bold", size = block_text_size, colour = "grey30"
      )
    } +
    geom_tile(
      data = grid_df,
      aes(x = x, y = y),
      fill = "white",
      color = grid_color,
      linewidth = 0.3
    ) +
    geom_point(
      data = off_df,
      aes(x = x, y = y, size = size_val, color = r, alpha = absr)
    ) +
    geom_text(
      data = off_df,
      aes(x = x, y = y, label = r_lab),
      size = txt_sz_off,
      color = "black"
    ) +
    geom_text(
      data = diag_df,
      aes(x = x, y = y, label = sd_lab),
      size = txt_sz_diag,
      color = "black"
    ) +
    scale_size_continuous(
      name = "$| r |$",
      breaks = size_breaks_trans,
      labels = formatC(size_breaks, digits = 1, format = "f"),
      range = c(1, circle_max),
      guide = guide_legend(
        override.aes = list(alpha = 1, color = "grey30"),
        order = 2
      )
    ) +
    scale_alpha(range = c(0.15, 1), guide = "none") +
    scale_corr_colors() +
    coord_fixed(
      xlim = c(0.5, n + 0.5),
      ylim = c(0.5, n + 0.5),
      expand = FALSE
    ) +
    labs(x = NULL, y = NULL, title = title) +
    scale_x_continuous(breaks = 1:n, labels = x_labs) +
    scale_y_continuous(breaks = 1:n, labels = y_labs) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = element_blank(),
      legend.box = "vertical",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  if (has_blocks) {
    cut <- m + 0.5
    ycut <- (n - m) + 0.5
    
    p <- p +
      annotate(
        "segment",
        x = cut, xend = cut, y = 0.5, yend = n + 0.5,
        color = block_line_color, linewidth = block_line_width
      ) +
      annotate(
        "segment",
        x = 0.5, xend = n + 0.5, y = ycut, yend = ycut,
        color = block_line_color, linewidth = block_line_width
      )
  }
  
  p
}

plot_corr_circles_rect <- function(R, P = NULL, alpha = 0.05,
                                   title = NULL,
                                   x_title = "Post-intervention",
                                   y_title = "Baseline",
                                   y_axis_position = c("left", "right"),
                                   circle_max = 12, size_gamma = 1.7,
                                   digits_r = 2,
                                   label_fun = make_pretty_label_vec,
                                   grid_color = "grey85",
                                   r_text_size = 2.8,
                                   row_split = 5, col_split = 5,
                                   block_line_color = "grey20",
                                   block_line_width = 1.4,
                                   size_breaks = c(0.2, 0.5, 0.8),
                                   oob_cells = NULL,
                                   oob_mark_size = 3.2,
                                   show_y_labels = FALSE,
                                   right_title_size = 5.2,
                                   right_title_x = 0.95) {
  stopifnot(is.matrix(R))
  rn <- rownames(R)
  cn <- colnames(R)
  nr <- length(rn)
  nc <- length(cn)
  
  y_axis_position <- match.arg(y_axis_position)
  
  if (!is.null(P)) {
    stopifnot(
      is.matrix(P),
      identical(rownames(P), rn),
      identical(colnames(P), cn)
    )
  }
  
  grid_df <- expand.grid(row = rn, col = cn, stringsAsFactors = FALSE)
  grid_df$i <- match(grid_df$row, rn)
  grid_df$j <- match(grid_df$col, cn)
  grid_df$x <- grid_df$j
  grid_df$y <- nr - grid_df$i + 1
  grid_df <- unique(grid_df[, c("x", "y")])
  
  df <- expand.grid(row = rn, col = cn, stringsAsFactors = FALSE)
  df$i <- match(df$row, rn)
  df$j <- match(df$col, cn)
  df$r <- as.numeric(R[cbind(df$i, df$j)])
  df$p <- if (!is.null(P)) as.numeric(P[cbind(df$i, df$j)]) else NA_real_
  df <- df[is.finite(df$r), , drop = FALSE]
  
  df$x <- df$j
  df$y <- nr - df$i + 1
  df$absr <- abs(df$r)
  df$size_val <- df$absr^size_gamma
  df$star <- ifelse(!is.na(df$p) & df$p < alpha, "*", "")
  df$r_lab <- paste0(formatC(df$r, digits = digits_r, format = "f"), df$star)
  
  size_breaks <- sort(size_breaks)
  size_breaks <- size_breaks[size_breaks > 0 & size_breaks <= 1]
  size_breaks_trans <- size_breaks^size_gamma
  
  x_labs <- label_fun(cn)
  y_labs <- label_fun(rev(rn))
  
  txt_sz <- r_text_size * min(1, 12 / max(nr, nc))
  txt_sz <- max(1.2, txt_sz)
  
  extra_right_space <- if (!is.null(y_title) && nzchar(y_title)) 1.2 else 0
  
  p <- ggplot() +
    geom_tile(
      data = grid_df,
      aes(x = x, y = y),
      fill = "white",
      color = grid_color,
      linewidth = 0.3
    ) +
    geom_point(
      data = df,
      aes(x = x, y = y, size = size_val, color = r, alpha = absr)
    ) +
    geom_text(
      data = df,
      aes(x = x, y = y, label = r_lab),
      size = txt_sz,
      color = "black"
    ) +
    {if (!is.null(oob_cells) && nrow(oob_cells) > 0) {
      tmp <- oob_cells
      tmp$i <- match(tmp$A, rn)
      tmp$j <- match(tmp$B, cn)
      tmp <- tmp[is.finite(tmp$i) & is.finite(tmp$j), , drop = FALSE]
      tmp$x <- tmp$j
      tmp$y <- nr - tmp$i + 1
      
      geom_point(
        data = tmp,
        aes(x = x, y = y),
        shape = 4,
        size = oob_mark_size,
        stroke = 1.1,
        colour = "black"
      )
    }} +
    scale_size_continuous(
      name = "$| r |$",
      breaks = size_breaks_trans,
      labels = formatC(size_breaks, digits = 1, format = "f"),
      range = c(1, circle_max),
      guide = guide_legend(
        override.aes = list(alpha = 1, color = "grey30")
      )
    ) +
    scale_alpha(range = c(0.15, 1), guide = "none") +
    scale_corr_colors() +
    coord_fixed(
      xlim = c(0.5, nc + 0.5 + extra_right_space),
      ylim = c(0.5, nr + 0.5),
      expand = FALSE
    ) +
    labs(x = x_title, y = NULL, title = title) +
    scale_x_continuous(breaks = 1:nc, labels = x_labs) +
    scale_y_continuous(
      breaks = 1:nr,
      labels = if (show_y_labels) y_labs else rep("", nr)
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = if (show_y_labels) element_text() else element_blank(),
      axis.ticks.y = if (show_y_labels) element_line() else element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(5.5, 24, 5.5, 5.5)
    )
  
  if (!is.null(y_title) && nzchar(y_title)) {
    p <- p + annotate(
      "text",
      x = nc + right_title_x,
      y = (nr + 1) / 2,
      label = y_title,
      angle = -90,
      size = right_title_size
    )
  }
  
  if (!is.null(col_split) && col_split >= 1 && col_split < nc) {
    xcut <- col_split + 0.5
    p <- p + annotate(
      "segment",
      x = xcut, xend = xcut,
      y = 0.5, yend = nr + 0.5,
      linewidth = block_line_width,
      colour = block_line_color
    )
  }
  
  if (!is.null(row_split) && row_split >= 1 && row_split < nr) {
    ycut <- (nr - row_split) + 0.5
    p <- p + annotate(
      "segment",
      x = 0.5, xend = nc + 0.5,
      y = ycut, yend = ycut,
      linewidth = block_line_width,
      colour = block_line_color
    )
  }
  
  p
}


make_P_resid_from_tab <- function(tab, labs,
                                  base = "corr_resid",
                                  sep = "__",
                                  p_col = c("p_fdr", "p_fisher", "p")) {
  p_name <- p_col[p_col %in% names(tab)][1]
  if (length(p_name) == 0 || is.na(p_name)) {
    stop("nessuna colonna p trovata in tab")
  }
  
  par <- if ("label" %in% names(tab)) {
    tab$label
  } else if ("param" %in% names(tab)) {
    tab$param
  } else if ("parameter" %in% names(tab)) {
    tab$parameter
  } else if ("term" %in% names(tab)) {
    tab$term
  } else {
    rownames(tab)
  }
  
  P <- matrix(NA_real_, length(labs), length(labs), dimnames = list(labs, labs))
  
  hit <- startsWith(par, paste0(base, sep))
  
  for (k in which(hit)) {
    rest <- sub(paste0("^", base, sep), "", par[k])
    ij <- strsplit(rest, sep, fixed = TRUE)[[1]]
    if (length(ij) != 2) next
    
    i <- ij[1]
    j <- ij[2]
    
    if (i %in% labs && j %in% labs) {
      P[i, j] <- tab[[p_name]][k]
      P[j, i] <- tab[[p_name]][k]
    }
  }
  
  diag(P) <- NA_real_
  P
}


      
make_labs_time <- function(est, time = c("I0","Ipost"),
                           include_frailty = FALSE,
                           base_order = .default_base_order) {
  time <- match.arg(time)
  sd_names <- names(est)[grepl("^sd_", names(est))]
  vars <- sub("^sd_", "", sd_names)
  
  want <- paste0(setdiff(base_order, "frailty"), "_", time)
  if (isTRUE(include_frailty)) want <- c(want, "frailty")
  want[want %in% vars]
}


panel_title <- function(tag, title, add = add_abc) {
  if (!add) return(title)
  sprintf("(%s) %s", tag, title)
}



      
#### wald su cov
make_covtab_frailty <- function(psi, se = NULL, times = c("I0","Ipost")) {
  nm <- names(psi)
  if (is.null(nm)) stop("psi deve essere un named numeric vector.")
  
  # due formati ammessi:
  # 1) cov_frailty__Mobility_I0
  # 2) cov_Mobility_I0__frailty
  pat1 <- "^cov_frailty__(.+)$"
  pat2 <- "^cov_(.+)__frailty$"
  
  cov_names <- unique(c(grep(pat1, nm, value = TRUE),
                        grep(pat2, nm, value = TRUE)))
  
  if (length(cov_names) == 0) {
    warning("Nessuna cov frailty–item trovata (pattern cov_frailty__* o cov_*__frailty).")
    return(data.frame())
  }
  
  get_item_time <- function(pname) {
    it <- NA_character_
    if (grepl(pat1, pname)) it <- sub(pat1, "\\1", pname)
    if (grepl(pat2, pname)) it <- sub(pat2, "\\1", pname)
    
    # time = suffisso _I0 / _Ipost
    time <- sub("^.*_(I0|Ipost)$", "\\1", it)
    if (!time %in% times) time <- NA_character_
    
    base_item <- sub("_(I0|Ipost)$", "", it)
    
    list(item_time = it, item = base_item, time = time)
  }
  
  info <- lapply(cov_names, get_item_time)
  item_time <- vapply(info, `[[`, character(1), "item_time")
  item      <- vapply(info, `[[`, character(1), "item")
  time      <- vapply(info, `[[`, character(1), "time")
  
  out <- data.frame(
    param     = cov_names,                    # <-- nome originale completo
    item_time = item_time,
    item      = item,
    time      = time,
    est       = unname(psi[cov_names]),
    se        = if (!is.null(se)) unname(se[cov_names]) else NA_real_,
    stringsAsFactors = FALSE
  )
  
  # tieni solo tempi richiesti se presenti (I0/Ipost); se time NA, li lasciamo fuori
  out <- out[!is.na(out$time), , drop = FALSE]
  
  out <- out[order(out$time, out$item_time), ]
  rownames(out) <- NULL
  out
}


wald_joint <- function(psi, params, V = NULL, se = NULL, rank_tol = 1e-10) {
  nm <- names(psi)
  if (is.null(nm)) stop("psi deve essere un named numeric vector.")
  
  miss <- setdiff(params, nm)
  if (length(miss)) stop("Parametri non trovati in psi: ", paste(miss, collapse = ", "))
  
  c_hat <- as.numeric(psi[params])
  
  if (!is.null(V)) {
    # usa la VCOV vera
    V <- V[nm, nm, drop = FALSE]
    V_sub <- V[params, params, drop = FALSE]
  } else {
    # fallback: VCOV diagonale da se (approssimazione)
    if (is.null(se)) stop("Serve V (vcov completa) oppure se (per VCOV diagonale).")
    miss_se <- setdiff(params, names(se))
    if (length(miss_se)) stop("SE mancanti per: ", paste(miss_se, collapse = ", "))
    V_sub <- diag(as.numeric(se[params])^2, nrow = length(params))
    rownames(V_sub) <- colnames(V_sub) <- params
  }
  
  # simmetrizza per sicurezza numerica
  V_sub <- 0.5 * (V_sub + t(V_sub))
  
  # df effettivi (rank)
  eig <- eigen(V_sub, symmetric = TRUE, only.values = TRUE)$values
  df_eff <- sum(eig > max(eig) * rank_tol)
  
  if (!requireNamespace("MASS", quietly = TRUE)) stop("Installa MASS per ginv().")
  Vinv <- tryCatch(solve(V_sub), error = function(e) MASS::ginv(V_sub))
  
  W <- as.numeric(t(c_hat) %*% Vinv %*% c_hat)
  p <- pchisq(W, df = df_eff, lower.tail = FALSE)
  
  list(W = W, df = df_eff, p = p, params = params, V_sub = V_sub, c_hat = c_hat)
}


subset_covtab <- function(cov_tab, items = NULL, time = NULL) {
  out <- cov_tab
  if (!is.null(items)) out <- out[out$item %in% items, , drop = FALSE]
  if (!is.null(time))  out <- out[out$time %in% time,  , drop = FALSE]
  out
}

make_test_row_cov <- function(label, cov_tab_sub, psi, V = NULL, se = NULL) {
  if (is.null(cov_tab_sub) || nrow(cov_tab_sub) == 0) {
    return(data.frame(test = label, n = 0L, W = NA_real_, df = NA_real_, p = NA_real_))
  }
  r <- wald_joint(psi = psi, params = cov_tab_sub$param, V = V, se = se)
  data.frame(test = label, n = nrow(cov_tab_sub), W = r$W, df = r$df, p = r$p)
}


wald_frailty_cov_suite <- function(psi_global, se_global = NULL, V_global = NULL,
                                   EQ_items, KOOS_items,
                                   times = c("I0","Ipost")) {
  cov_tab <- make_covtab_frailty(psi = psi_global, se = se_global, times = times)
  
  cov_I0    <- subset_covtab(cov_tab, time = "I0")
  cov_Ipost <- subset_covtab(cov_tab, time = "Ipost")
  
  cov_EQ_all<- subset_covtab(cov_tab, items = EQ_items)
  cov_EQ_I0 <- subset_covtab(cov_tab, items = EQ_items, time = "I0")
  cov_EQ_Ipost<- subset_covtab(cov_tab, items = EQ_items, time = "Ipost")
  
  cov_KOOS_all<- subset_covtab(cov_tab, items = KOOS_items)
  cov_KOOS_I0 <- subset_covtab(cov_tab, items = KOOS_items, time = "I0")
  cov_KOOS_Ipost<- subset_covtab(cov_tab, items = KOOS_items, time = "Ipost")
  
  tests_df <- do.call(rbind, list(
    make_test_row_cov("All_frailty_cov",       cov_tab,       psi_global, V_global, se_global),
    make_test_row_cov("I0_frailty_cov",        cov_I0,        psi_global, V_global, se_global),
    make_test_row_cov("Ipost_frailty_cov",     cov_Ipost,     psi_global, V_global, se_global),
    
    make_test_row_cov("EQ_all_frailty_cov",    cov_EQ_all,    psi_global, V_global, se_global),
    make_test_row_cov("EQ_I0_frailty_cov",     cov_EQ_I0,     psi_global, V_global, se_global),
    make_test_row_cov("EQ_Ipost_frailty_cov",  cov_EQ_Ipost,  psi_global, V_global, se_global),
    
    make_test_row_cov("KOOS_all_frailty_cov",  cov_KOOS_all,  psi_global, V_global, se_global),
    make_test_row_cov("KOOS_I0_frailty_cov",   cov_KOOS_I0,   psi_global, V_global, se_global),
    make_test_row_cov("KOOS_Ipost_frailty_cov",cov_KOOS_Ipost,psi_global, V_global, se_global)
  ))
  
  rownames(tests_df) <- NULL
  
  list(
    cov_tab = cov_tab,   # <-- qui hai param = nome originale intero
    tests   = tests_df
  )
}

                   
#####################################################
##############     PCA and plots     ################
#####################################################               
get_pca_xy <- function(M) {
  out <- pca_loadings_df(M)$df
  list(x = out$PC1, y = out$PC2)
}

strip_D_summaries_from_est <- function(est) {
  keep <- grepl(
    "^(fixed_|threshold_|log_sigma_|sigma_|fixed_surv_|log_lambda$|log_rho$|lambda$|rho$|HR_|atanh_rho_resid__|corr_resid__)",
    names(est)
  )
  est[keep]
}
                   
pca_loadings_df <- function(M) {
  M <- as.matrix(M)
  eg <- eigen(M, symmetric = TRUE)
  ord <- order(eg$values, decreasing = TRUE)
  vals <- eg$values[ord]
  vecs <- eg$vectors[, ord, drop = FALSE]
  load <- vecs %*% diag(sqrt(pmax(vals, 0)))
  x <- load[, 1]
  y <- load[, 2]
  if (sum(x) < 0) x <- -x
  if (sum(y) < 0) y <- -y
  list(
    df = data.frame(name = rownames(M), PC1 = x, PC2 = y, stringsAsFactors = FALSE),
    prop = vals / sum(vals)
  )
}


pca_plot_grouped <- function(M,
                             main = "",
                             short_map,
                             eq_base,
                             koos_base,
                             mode = c("timed", "resid"),
                             show_suffix = FALSE,
                             xlim_fixed = NULL,
                             ylim_fixed = NULL,
                             label_size = 4.0,
                             axis_title_size = 10,
                             axis_text_size = 11,
                             plot_title_size = 12,
                             legend_text_size = 9.5,
                             force = 0.10,
                             force_pull = 5,
                             box_padding = 0.12,
                             point_padding = 0.02,
                             escape_tex = FALSE,
                             show_legend = TRUE,
                             legend_inside = TRUE,
                             legend_pos = c(0.02, 0.98),
                             legend_reserve_left = 0.24,
                             legend_reserve_top = 0.02,
                             legend_bg_alpha = 0.0,
                             cols_main = c(
                               "EQ-5D_I0"    = "#0F766E",
                               "EQ-5D_Ipost" = "#d1a00d",
                               "KOOS_I0"     = "#5B21B6",
                               "KOOS_Ipost"  = "#52340b",
                               "Frailty"     = "#C2410C"
                             ),
                             cols_resid = c(
                               "EQ"   = "#7a7a7a",
                               "KOOS" = "#111111"
                             )) {
  
  mode <- match.arg(mode)
  
  esc_tex <- function(x) {
    x <- gsub("\\\\", "\\\\textbackslash{}", x)
    x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
    x
  }
  
  pca_loadings_df <- function(M) {
    M <- as.matrix(M)
    eg <- eigen(M, symmetric = TRUE)
    ord <- order(eg$values, decreasing = TRUE)
    vals <- eg$values[ord]
    vecs <- eg$vectors[, ord, drop = FALSE]
    load <- vecs %*% diag(sqrt(pmax(vals, 0)))
    x <- load[, 1]
    y <- load[, 2]
    if (sum(x) < 0) x <- -x
    if (sum(y) < 0) y <- -y
    list(
      df = data.frame(
        name = rownames(M),
        PC1 = x,
        PC2 = y,
        stringsAsFactors = FALSE
      ),
      prop = vals / sum(vals)
    )
  }
  
  make_k <- function(n) {
    out <- numeric(0)
    j <- 0
    while (length(out) < n) {
      if (j == 0) {
        out <- c(out, 0)
      } else {
        out <- c(out, j, -j)
      }
      j <- j + 1
    }
    out[seq_len(n)]
  }
  
  short_label <- function(x, short_map, show_suffix = FALSE) {
    if (identical(x, "frailty")) return("Frailty")
    base <- sub("(_I0|_Ipost)$", "", x, perl = TRUE)
    suff <- sub(paste0("^", base), "", x)
    lab <- if (base %in% names(short_map)) unname(short_map[base]) else base
    if (show_suffix) paste0(lab, suff) else lab
  }
  
  group_timed <- function(x, eq_base, koos_base) {
    if (identical(x, "frailty")) return("Frailty")
    base <- sub("(_I0|_Ipost)$", "", x, perl = TRUE)
    if (base %in% eq_base   && grepl("_I0$", x))    return("EQ-5D_I0")
    if (base %in% eq_base   && grepl("_Ipost$", x)) return("EQ-5D_Ipost")
    if (base %in% koos_base && grepl("_I0$", x))    return("KOOS_I0")
    if (base %in% koos_base && grepl("_Ipost$", x)) return("KOOS_Ipost")
    NA_character_
  }
  
  group_resid <- function(x, eq_base, koos_base) {
    base <- sub("(_I0|_Ipost)$", "", x, perl = TRUE)
    if (base %in% eq_base)   return("EQ")
    if (base %in% koos_base) return("KOOS")
    NA_character_
  }
  
  add_repulsion_columns <- function(df) {
    df$r <- sqrt(df$PC1^2 + df$PC2^2)
    df$r[df$r < 1e-8] <- 1e-8
    df$ux <- df$PC1 / df$r
    df$uy <- df$PC2 / df$r
    df$nx <- -df$uy
    df$ny <-  df$ux
    df$theta <- atan2(df$PC2, df$PC1)
    
    ref_lim <- max(abs(c(df$PC1, df$PC2, 0)))
    grp_quad <- paste0(ifelse(df$PC1 >= 0, "R", "L"),
                       ifelse(df$PC2 >= 0, "U", "D"))
    
    df$k <- 0
    for (g in unique(grp_quad)) {
      id <- which(grp_quad == g)
      oo <- order(df$theta[id], df$r[id])
      df$k[id[oo]] <- make_k(length(id))
    }
    
    push0 <- 0.012 * ref_lim
    perp0 <- 0.008 * ref_lim
    
    df$nudge_x <- df$ux * push0 + df$nx * (df$k * perp0)
    df$nudge_y <- df$uy * push0 + df$ny * (df$k * perp0)
    df$hjust <- ifelse(df$PC1 >= 0, 0, 1)
    df
  }
  
  pct_lab <- function(txt, p) sprintf("%s (%.1f%%)", txt, 100 * p)
  pct_lab_tex <- function(txt, p) sprintf("%s (%.1f\\%%)", txt, 100 * p)
  
  M <- as.matrix(M)
  if (nrow(M) != ncol(M)) stop("R deve essere quadrata.")
  if (anyNA(M)) stop("R contiene NA.")
  if (max(abs(M - t(M))) > 1e-10) stop("R deve essere simmetrica.")
  if (is.null(rownames(M))) stop("servono i rownames di R")
  
  out <- pca_loadings_df(M)
  df <- out$df
  prop <- out$prop
  
  df$label <- vapply(
    df$name, short_label, character(1),
    short_map = short_map,
    show_suffix = show_suffix
  )
  
  if (mode == "timed") {
    df$group <- vapply(
      df$name, group_timed, character(1),
      eq_base = eq_base,
      koos_base = koos_base
    )
    if (anyNA(df$group)) {
      stop(paste("variabili non classificate:",
                 paste(df$name[is.na(df$group)], collapse = ", ")))
    }
    pal <- cols_main
    legend_labels <- c(
      "EQ-5D_I0"    = "EQ-5D_I0",
      "EQ-5D_Ipost" = "EQ-5D_Ipost",
      "KOOS_I0"     = "KOOS_I0",
      "KOOS_Ipost"  = "KOOS_Ipost",
      "Frailty"     = "Frailty"
    )
  } else {
    df$group <- vapply(
      df$name, group_resid, character(1),
      eq_base = eq_base,
      koos_base = koos_base
    )
    if (anyNA(df$group)) {
      stop(paste("variabili non classificate:",
                 paste(df$name[is.na(df$group)], collapse = ", ")))
    }
    pal <- cols_resid
    legend_labels <- c(
      "EQ"   = "EQ",
      "KOOS" = "KOOS"
    )
  }
  
  used_levels <- names(pal)[names(pal) %in% unique(df$group)]
  df$group <- factor(df$group, levels = used_levels)
  
  if (escape_tex) {
    df$label <- esc_tex(df$label)
    main <- esc_tex(main)
    legend_labels <- setNames(esc_tex(unname(legend_labels)), names(legend_labels))
    x_lab <- pct_lab_tex("First Principal Component", prop[1])
    y_lab <- pct_lab_tex("Second Principal Component", prop[2])
  } else {
    x_lab <- pct_lab("First Principal Component", prop[1])
    y_lab <- pct_lab("Second Principal Component", prop[2])
  }
  
  df <- add_repulsion_columns(df)
  
  ref_lim <- max(abs(c(df$PC1, df$PC2, 0)))
  if (is.null(xlim_fixed) || is.null(ylim_fixed)) {
    xr <- range(c(0, df$PC1), na.rm = TRUE)
    yr <- range(c(0, df$PC2), na.rm = TRUE)
    xpad <- 0.10 * diff(xr)
    ypad <- 0.10 * diff(yr)
    if (!is.finite(xpad) || xpad == 0) xpad <- 0.10 * ref_lim
    if (!is.finite(ypad) || ypad == 0) ypad <- 0.10 * ref_lim
    xlim <- xr + c(-xpad, xpad)
    ylim <- yr + c(-ypad, ypad)
  } else {
    xlim <- xlim_fixed
    ylim <- ylim_fixed
  }
  
  if (show_legend && legend_inside) {
    xspan <- diff(xlim)
    yspan <- diff(ylim)
    xlim <- c(xlim[1] - legend_reserve_left * xspan, xlim[2])
    ylim <- c(ylim[1], ylim[2] + legend_reserve_top * yspan)
  }
  
  p <- ggplot(df, aes(PC1, PC2, color = group)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey88") +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey88") +
    geom_segment(
      aes(x = 0, y = 0, xend = PC1, yend = PC2),
      linewidth = 0.85,
      arrow = arrow(type = "closed", length = grid::unit(0.18, "cm")),
      show.legend = show_legend
    ) +
    ggrepel::geom_text_repel(
      aes(label = label),
      nudge_x = df$nudge_x,
      nudge_y = df$nudge_y,
      hjust = df$hjust,
      size = label_size,
      fontface = "bold",
      direction = "both",
      force = force,
      force_pull = force_pull,
      box.padding = box_padding,
      point.padding = point_padding,
      min.segment.length = 0,
      segment.color = "grey45",
      segment.size = 0.25,
      max.overlaps = Inf,
      seed = 123,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = pal[used_levels],
      breaks = used_levels,
      labels = legend_labels[used_levels],
      drop = TRUE,
      guide = if (show_legend) {
        guide_legend(
          ncol = 1,
          byrow = TRUE,
          override.aes = list(linewidth = 1.4)
        )
      } else {
        "none"
      }
    ) +
    coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE, clip = "on") +
    labs(
      title = main,
      x = x_lab,
      y = y_lab,
      color = NULL
    ) +
    theme_minimal(base_size = 12) +
   theme(
  plot.title = element_text(hjust = 0.5, face = "bold", size = plot_title_size),
  axis.title = element_text(face = "bold", size = axis_title_size),
  axis.text = element_text(size = axis_text_size),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(linewidth = 0.3, color = "grey92"),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
  plot.margin = margin(8, 10, 8, 10),
  legend.position = if (legend_inside) legend_pos else "bottom",
  legend.justification = c(0, 1),
  legend.direction = "vertical",
  legend.background = element_blank(),
  legend.box.background = element_blank(),
  legend.key = element_rect(fill = "transparent", color = NA),
  legend.key.height = grid::unit(0.32, "cm"),
  legend.key.width = grid::unit(0.55, "cm"),
  legend.spacing.y = grid::unit(0.03, "cm"),
  legend.margin = margin(1, 2, 1, 2),
  legend.text = element_text(size = legend_text_size)
)
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  p
}

  
extract_frailty_corr_wide_fisher_fdr <- function(est, se,
                                                 alpha = 0.05,
                                                 item_levels = NULL,
                                                 corr_levels = c("I0","Ipost"),
                                                 corr_prefix = "corrFrailty",
                                                 fdr_method = "fdr",
                                                 eps = 1e-6) {
  stopifnot(is.numeric(est), is.numeric(se),
            !is.null(names(est)), !is.null(names(se)))
  
  nm <- intersect(names(est), names(se))
  est <- est[nm]
  se  <- se[nm]
  
  cand <- nm[startsWith(nm, "corr_") & !startsWith(nm, "corr_resid__")]
  cand <- cand[grepl("frailty", cand, fixed = TRUE)]
  
  if (!length(cand)) {
    return(data.frame(item = character(0), stringsAsFactors = FALSE))
  }
  
  parse_one <- function(nm1) {
    s <- sub("^corr_", "", nm1)
    ab <- strsplit(s, "__", fixed = TRUE)[[1]]
    if (length(ab) != 2L) return(NULL)
    
    a <- ab[1]
    b <- ab[2]
    
    if (a == "frailty" && b != "frailty") {
      other <- b
    } else if (b == "frailty" && a != "frailty") {
      other <- a
    } else {
      return(NULL)
    }
    
    if (!grepl("_(I0|Ipost)$", other)) return(NULL)
    
    lev  <- sub("^.*_(I0|Ipost)$", "\\1", other)
    item <- sub("_(I0|Ipost)$", "", other)
    
    list(item = item, re_level = lev)
  }
  
  parsed <- lapply(cand, parse_one)
  ok <- !vapply(parsed, is.null, logical(1))
  
  cand <- cand[ok]
  parsed <- parsed[ok]
  
  if (!length(cand)) {
    return(data.frame(item = character(0), stringsAsFactors = FALSE))
  }
  
  df <- data.frame(
    item = vapply(parsed, `[[`, "", "item"),
    re_level = vapply(parsed, `[[`, "", "re_level"),
    estimate = as.numeric(est[cand]),
    se = as.numeric(se[cand]),
    stringsAsFactors = FALSE
  )
  
  df <- df[df$re_level %in% corr_levels, , drop = FALSE]
  if (!nrow(df)) return(data.frame(item = character(0), stringsAsFactors = FALSE))
  
  r_clamp <- pmin(pmax(df$estimate, -1 + eps), 1 - eps)
  valid <- is.finite(df$estimate) & is.finite(df$se) & df$se > 0 & abs(df$estimate) < 1
  
  df$fisher_z <- NA_real_
  df$se_fisher_z <- NA_real_
  df$z <- NA_real_
  df$p_value <- NA_real_
  
  df$fisher_z[valid] <- atanh(r_clamp[valid])
  df$se_fisher_z[valid] <- df$se[valid] / (1 - r_clamp[valid]^2)
  df$z[valid] <- df$fisher_z[valid] / df$se_fisher_z[valid]
  df$p_value[valid] <- 2 * pnorm(-abs(df$z[valid]))
  
  df$p_fdr <- NA_real_
  if (any(valid)) {
    df$p_fdr[valid] <- p.adjust(df$p_value[valid], method = fdr_method)
  }
  
  df$sig_0.05 <- !is.na(df$p_value) & df$p_value < alpha
  df$sig_fdr_0.05 <- !is.na(df$p_fdr) & df$p_fdr < alpha
  
  if (!is.null(item_levels)) df$item <- factor(df$item, levels = item_levels)
  
  out <- tidyr::pivot_wider(
    df[, c("item","re_level","estimate","se","p_value","p_fdr","sig_0.05","sig_fdr_0.05")],
    names_from = re_level,
    values_from = c(estimate, se, p_value, p_fdr, sig_0.05, sig_fdr_0.05),
    names_glue = paste0(corr_prefix, "_{re_level}_{.value}")
  )
  
  if (!is.null(item_levels)) {
    out$item <- factor(out$item, levels = item_levels)
    out <- out[order(out$item), , drop = FALSE]
  } else {
    out <- out[order(out$item), , drop = FALSE]
  }
  
  rownames(out) <- NULL
  out
}



  extract_fixed_wald_wide_fdr <- function(psi, se,
                                        alpha = 0.05,
                                        predictors = c("I6","I12","wI6","wI12","pI6","pI12",
                                                       "diag_oa","sex01","loghosp_post","age"),
                                        item_levels = NULL,
                                        exclude_regex = "^fixed_surv_",
                                        fdr_method = "fdr") {
  stopifnot(is.numeric(psi), is.numeric(se),
            !is.null(names(psi)), !is.null(names(se)))
  
  nm_all <- intersect(names(psi), names(se))
  psi <- psi[nm_all]
  se  <- se[nm_all]
  
  keep <- grepl("^fixed_", nm_all) & !grepl(exclude_regex, nm_all)
  psi_fix <- psi[keep]
  se_fix  <- se[names(psi_fix)]
  
  nm <- names(psi_fix)
  if (!length(nm)) return(data.frame(item = character(0), stringsAsFactors = FALSE))
  
  pred_pat <- paste(predictors, collapse = "|")
  re <- paste0("^fixed_(.*)_(", pred_pat, ")$")
  
  m <- regexec(re, nm)
  mm <- regmatches(nm, m)
  ok <- lengths(mm) == 3L
  
  nm_ok <- nm[ok]
  mm_ok <- mm[ok]
  
  if (!length(nm_ok)) {
    return(data.frame(item = character(0), stringsAsFactors = FALSE))
  }
  
  item <- vapply(mm_ok, `[`, "", 2)
  pred <- vapply(mm_ok, `[`, "", 3)
  
  df <- data.frame(
    item = item,
    predictor = pred,
    estimate = as.numeric(psi_fix[nm_ok]),
    se = as.numeric(se_fix[nm_ok]),
    stringsAsFactors = FALSE
  )
  
  df$z <- df$estimate / df$se
  df$p_value <- 2 * pnorm(-abs(df$z))
  df$p_fdr <- p.adjust(df$p_value, method = fdr_method)
  df$sig_0.05 <- !is.na(df$p_value) & df$p_value < alpha
  df$sig_fdr_0.05 <- !is.na(df$p_fdr) & df$p_fdr < alpha
  
  if (!is.null(item_levels)) df$item <- factor(df$item, levels = item_levels)
  
  out <- tidyr::pivot_wider(
    df[, c("item","predictor","estimate","se","p_value","p_fdr","sig_0.05","sig_fdr_0.05")],
    names_from = predictor,
    values_from = c(estimate, se, p_value, p_fdr, sig_0.05, sig_fdr_0.05),
    names_glue = "{predictor}_{.value}"
  )
  
  if (!is.null(item_levels)) {
    out$item <- factor(out$item, levels = item_levels)
    out <- out[order(out$item), , drop = FALSE]
  } else {
    out <- out[order(out$item), , drop = FALSE]
  }
  
  rownames(out) <- NULL
  out
}


  make_wald_table_3cols_fdr <- function(psi, se, est_corr,
                                      alpha = 0.05,
                                      predictors = c("I6","I12","wI6","wI12","pI6","pI12",
                                                     "diag_oa","sex01","loghosp_post","age"),
                                      item_levels = NULL,
                                      exclude_regex = "^fixed_surv_",
                                      include_corr_frailty = TRUE,
                                      corr_levels = c("I0","Ipost"),
                                      corr_prefix = "corrFrailty",
                                      fdr_method_fixed = "fdr",
                                      fdr_method_corr = "fdr") {
  
  wide_fix <- extract_fixed_wald_wide_fdr(
    psi = psi,
    se = se,
    alpha = alpha,
    predictors = predictors,
    item_levels = item_levels,
    exclude_regex = exclude_regex,
    fdr_method = fdr_method_fixed
  )
  
  if (!isTRUE(include_corr_frailty)) return(wide_fix)
  
  wide_corr <- extract_frailty_corr_wide_fisher_fdr(
    est = est_corr,
    se = se,
    alpha = alpha,
    item_levels = item_levels,
    corr_levels = corr_levels,
    corr_prefix = corr_prefix,
    fdr_method = fdr_method_corr
  )
  
  if (!nrow(wide_corr)) return(wide_fix)
  
  out <- merge(wide_fix, wide_corr, by = "item", all.x = TRUE, sort = FALSE)
  
  if (!is.null(item_levels)) {
    out$item <- factor(out$item, levels = item_levels)
    out <- out[order(out$item), , drop = FALSE]
  }
  
  rownames(out) <- NULL
  out
}
make_latex_beta_table_fdr <- function(tab_items3,
                                      predictors = c("I6","I12","wI6","wI12","pI6","pI12",
                                                     "diag_oa","sex01","loghosp_post","age"),
                                      extra_terms = NULL,
                                      item_name_map = NULL,
                                      alpha = 0.05,
                                      digits_est = 2,
                                      digits_se  = 2,
                                      p_col_suffix = "p_fdr",
                                      caption = NULL,
                                      label = NULL,
                                      use_resizebox = TRUE,
                                      font_size_cmd = "\\small") {
  
  stopifnot(is.data.frame(tab_items3), "item" %in% names(tab_items3))
  
  term_prefixes <- predictors
  
  col_headers_pred <- vapply(predictors, function(pr) {
    pr2 <- escape_for_texttt(pr)
    paste0("$\\\\widehat{\\\\beta_{\\\\texttt{", pr2, "}}}$")
  }, character(1))
  
  term_headers <- col_headers_pred
  
  if (!is.null(extra_terms)) {
    extra_prefix <- vapply(extra_terms, `[[`, "", "prefix")
    extra_head   <- vapply(extra_terms, `[[`, "", "header")
    term_prefixes <- c(term_prefixes, extra_prefix)
    term_headers  <- c(term_headers,  extra_head)
  }
  
  need <- unlist(lapply(term_prefixes, function(pr) {
    c(
      paste0(pr, "_estimate"),
      paste0(pr, "_se"),
      paste0(pr, "_", p_col_suffix)
    )
  }))
  
  miss <- setdiff(need, names(tab_items3))
  if (length(miss) > 0) stop("colonne mancanti in tab_items3: ", paste(miss, collapse = ", "))
  
  item_raw <- as.character(tab_items3$item)
  item_disp <- item_raw
  
  if (!is.null(item_name_map)) {
    if (is.null(names(item_name_map))) stop("item_name_map deve essere named")
    hit <- item_raw %in% names(item_name_map)
    item_disp[hit] <- unname(item_name_map[item_raw[hit]])
  }
  
  item_disp <- latex_escape(item_disp)
  
  fmt <- function(x, d) ifelse(is.na(x), "", formatC(x, digits = d, format = "f"))
  
  n <- nrow(tab_items3)
  m <- length(term_prefixes)
  cells <- matrix("", nrow = n, ncol = m)
  
  for (j in seq_along(term_prefixes)) {
    pr <- term_prefixes[j]
    est <- tab_items3[[paste0(pr, "_estimate")]]
    se  <- tab_items3[[paste0(pr, "_se")]]
    pv  <- tab_items3[[paste0(pr, "_", p_col_suffix)]]
    
    star <- ifelse(!is.na(pv) & pv < alpha, "$^{*}$", "")
    
    cells[, j] <- ifelse(
      is.na(est) | is.na(se),
      "",
      paste0(fmt(est, digits_est), star, " (", fmt(se, digits_se), ")")
    )
  }
  
  header_line <- paste0("Item & ", paste(term_headers, collapse = " & "), " \\\\")
  body_lines <- vapply(seq_len(n), function(i) {
    paste0(item_disp[i], " & ", paste(cells[i, ], collapse = " & "), " \\\\")
  }, character(1))
  
  align <- paste0("l", paste(rep("r", m), collapse = ""))
  
  tabular <- paste0(
    font_size_cmd, "\n",
    "\\begin{tabular}{", align, "}\n",
    "\\toprule\n",
    header_line, "\n",
    "\\midrule\n",
    paste(body_lines, collapse = "\n"), "\n",
    "\\bottomrule\n",
    "\\end{tabular}\n"
  )
  
  out <- "\\begin{table}[htbp]\n\\centering\n"
  if (use_resizebox) out <- paste0(out, "\\resizebox{\\textwidth}{!}{%\n")
  out <- paste0(out, tabular)
  if (use_resizebox) out <- paste0(out, "}\n")
  if (!is.null(caption)) out <- paste0(out, "\\caption{", caption, "}\n")
  if (!is.null(label))   out <- paste0(out, "\\label{", label, "}\n")
  out <- paste0(out, "\\end{table}\n")
  
  out
}


escape_for_texttt <- function(x) {
  # dentro \texttt{...} gli underscore vanno escapati
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%#$&{}_])", "\\\\\\1", x, perl = TRUE)
  x
}



#### Wald congiunto per H0 : Cov(b_{ik2}-b_{ik1},s_i)

wald_linear <- function(psi, L, V = NULL, se = NULL, c0 = NULL, rank_tol = 1e-10) {
  stopifnot(is.numeric(psi), !is.null(names(psi)))
  
  L <- as.matrix(L)
  if (is.null(colnames(L))) stop("L deve avere colnames uguali ai nomi di psi.")
  
  miss <- setdiff(colnames(L), names(psi))
  if (length(miss)) stop("colonne di L non trovate in psi: ", paste(miss, collapse = ", "))
  
  psi <- psi[colnames(L)]
  
  if (is.null(c0)) c0 <- rep(0, nrow(L))
  est <- as.numeric(L %*% psi - c0)
  names(est) <- rownames(L)
  
  if (!is.null(V)) {
    V <- as.matrix(V)
    V <- V[colnames(L), colnames(L), drop = FALSE]
    Vc <- L %*% V %*% t(L)
  } else {
    if (is.null(se)) stop("serve V oppure se")
    se <- se[colnames(L)]
    Vc <- L %*% diag(se^2, nrow = length(se)) %*% t(L)
  }
  
  Vc <- 0.5 * (Vc + t(Vc))
  
  eig <- eigen(Vc, symmetric = TRUE, only.values = TRUE)$values
  df_eff <- sum(eig > max(eig) * rank_tol)
  
  if (!requireNamespace("MASS", quietly = TRUE)) stop("serve MASS")
  Vinv <- tryCatch(solve(Vc), error = function(e) MASS::ginv(Vc))
  
  W <- as.numeric(t(est) %*% Vinv %*% est)
  p <- pchisq(W, df = df_eff, lower.tail = FALSE)
  
  list(W = W, df = df_eff, p = p, est = est, V = Vc, L = L)
}

make_L_shiftcov_frailty <- function(psi_names, items,
                                    t_pre = "I0",
                                    t_post = "Ipost",
                                    frailty = "frailty") {
  L <- matrix(0, nrow = length(items), ncol = length(psi_names),
              dimnames = list(paste0("covdiff_frailty__", items), psi_names))
  
  for (k in seq_along(items)) {
    it <- items[k]
    
    nm_pre <- get_name_anyorder(psi_names, "cov_", paste0(it, "_", t_pre), frailty)
    nm_post <- get_name_anyorder(psi_names, "cov_", paste0(it, "_", t_post), frailty)
    
    if (is.na(nm_pre) || is.na(nm_post)) {
      stop("covarianze con frailty mancanti per item: ", it)
    }
    
    L[k, nm_post] <- 1
    L[k, nm_pre] <- -1
  }
  
  L
}

wald_shiftcov_frailty <- function(psi_global, V_global = NULL, se_global = NULL,
                                  items,
                                  t_pre = "I0",
                                  t_post = "Ipost",
                                  frailty = "frailty") {
  L <- make_L_shiftcov_frailty(
    psi_names = names(psi_global),
    items = items,
    t_pre = t_pre,
    t_post = t_post,
    frailty = frailty
  )
  
  wald_linear(
    psi = psi_global,
    L = L,
    V = V_global,
    se = se_global
  )
}

         

