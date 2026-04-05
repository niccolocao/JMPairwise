

# tabella stime complete 
make_table_est_se_items_rows <- function(est, se,
                                         item_order = c(names(K_map), "Perceivedcurrenthealthstatus", "Global"),
                                         digits = 3) {
  stopifnot(!is.null(names(est)), !is.null(names(se)))
  
  se <- se[names(est)]
  
  keep <- grepl("^threshold_", names(est)) |
    grepl("^sigma_", names(est)) |
    grepl("^phi_", names(est)) |
    names(est) %in% c("lambda", "rho")
  
  df <- data.frame(
    name = names(est)[keep],
    est = as.numeric(est[keep]),
    se  = as.numeric(se[keep]),
    stringsAsFactors = FALSE
  )
  
  parse_name <- function(x) {
    if (grepl("^threshold_", x)) {
      item <- sub("^threshold_(.*)_cut:[0-9]+$", "\\1", x)
      cut  <- as.integer(sub("^threshold_.*_cut:([0-9]+)$", "\\1", x))
      return(c(item = item, parameter = paste0("threshold_cut", cut)))
    }
    
    if (grepl("^sigma_", x)) {
      item <- sub("^sigma_", "", x)
      return(c(item = item, parameter = "sigma"))
    }
    
    if (grepl("^phi_", x)) {
      item <- sub("^phi_", "", x)
      return(c(item = item, parameter = "phi"))
    }
    
    if (x %in% c("lambda", "rho")) {
      return(c(item = "Global", parameter = x))
    }
    
    c(item = NA_character_, parameter = NA_character_)
  }
  
  parsed <- t(vapply(df$name, parse_name, character(2)))
  df$item <- parsed[, "item"]
  df$parameter <- parsed[, "parameter"]
  df$cell <- sprintf(paste0("%.", digits, "f (%.", digits, "f)"), df$est, df$se)
  
  df <- df[!is.na(df$item) & !is.na(df$parameter), , drop = FALSE]
  
  thr_present <- df$parameter[grepl("^threshold_cut[0-9]+$", df$parameter)]
  thr_idx <- sort(unique(as.integer(sub("^threshold_cut", "", thr_present))))
  
  param_order <- c(
    paste0("threshold_cut", thr_idx),
    "sigma",
    "phi",
    "lambda",
    "rho"
  )
  
  rest <- setdiff(unique(df$parameter), param_order)
  param_order <- c(param_order, rest)
  
  row_order <- c(intersect(item_order, unique(df$item)),
                 setdiff(unique(df$item), item_order))
  
  df$parameter <- factor(df$parameter, levels = param_order)
  df$item <- factor(df$item, levels = row_order)
  
  df <- df |>
    dplyr::arrange(item, parameter)
  
  tab <- df |>
    dplyr::select(item, parameter, cell) |>
    tidyr::pivot_wider(
      names_from = parameter,
      values_from = cell,
      values_fill = ""
    ) |>
    dplyr::arrange(item) |>
    dplyr::mutate(item = as.character(item))
  
  names(tab)[1] <- "item"
  
  keep_cols <- c("item", intersect(param_order, names(tab)))
  tab[, keep_cols, drop = FALSE]
}
latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([_%&#$])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x
}

pretty_col_latex <- function(x) {
  out <- x
  
  thr <- grepl("^threshold_cut[0-9]+$", out)
  out[thr] <- paste0("Threshold ", sub("^threshold_cut", "", out[thr]))
  
  out[out == "sigma"] <- "$\\sigma$"
  out[out == "phi"] <- "$\\phi$"
  out[out == "lambda"] <- "$\\lambda$"
  out[out == "rho"] <- "$\\rho$"
  
  out
}

export_tab_items_rows_latex <- function(tab,
                                        file = "tab_stime_items_rows.tex",
                                        caption = "Parameter estimates with standard errors",
                                        label = "tab:stime_items_rows",
                                        resize = TRUE,
                                        font_size = "\\scriptsize") {
  
  tab <- as.data.frame(tab, stringsAsFactors = FALSE)
  
  head_tex <- names(tab)
  head_tex[1] <- "Item"
  head_tex[-1] <- pretty_col_latex(head_tex[-1])
  head_tex[1] <- latex_escape(head_tex[1])
  
  body <- tab
  body[[1]] <- latex_escape(body[[1]])
  for (j in 2:ncol(body)) {
    body[[j]] <- latex_escape(body[[j]])
  }
  
  rows_tex <- apply(body, 1, function(z) paste(z, collapse = " & "))
  rows_tex <- paste0(rows_tex, " \\\\")
  
  align <- paste0("l", paste(rep("c", ncol(tab) - 1), collapse = ""))
  
  tex <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    font_size
  )
  
  if (resize) tex <- c(tex, "\\resizebox{\\textwidth}{!}{%")
  
  tex <- c(
    tex,
    paste0("\\begin{tabular}{", align, "}"),
    "\\hline",
    paste(head_tex, collapse = " & "),
    " \\\\",
    "\\hline",
    rows_tex,
    "\\hline",
    "\\end{tabular}"
  )
  
  if (resize) tex <- c(tex, "}")
  
  tex <- c(
    tex,
    paste0("\\caption{", latex_escape(caption), "}"),
    paste0("\\label{", label, "}"),
    "\\end{table}"
  )
  
  writeLines(tex, con = file)
  invisible(tex)
}

############# wald for all + FDR for correlations only
wald_table_raw <- function(est, se,
                           alpha = 0.05,
                           null = NULL,          # NULL = HR->1, altrimenti 0; oppure scalare; named vector; o function(param,type)->H0
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


############
## Plots correlations
library(RColorBrewer)

split_corr_param <- function(param) {
  # param tipo: corr_<A>__<B>
  s <- sub("^corr_", "", param)
  parts <- strsplit(s, "__", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(list(A = NA_character_, B = NA_character_))
  list(A = parts[1], B = parts[2])
}

scale_corr_colors <- function() {
  scale_color_gradientn(
    colours = rev(brewer.pal(11, "RdBu")),
    limits = c(-1, 1),
    name = "$r$"
  )
}


make_pretty_label_vec <- function(x, map_base = short_map_base) {
  x <- as.character(x)
  
  # caso frailty
  x[x == "frailty"] <- "Frailty"
  
  # togli suffisso temporale prima di mappare
  base <- sub("_(I0|Ipost)$", "", x)
  
  # mappa -> short
  out <- unname(map_base[base])
  
  # fallback (se non mappato)
  out[is.na(out)] <- base[is.na(out)]
  
  out
}


make_labs_time <- function(est, time = c("I0","Ipost"), include_frailty = TRUE,
                           base_order = .default_base_order) {
  time <- match.arg(time)
  # variabili disponibili dagli sd_
  sd_names <- names(est)[grepl("^sd_", names(est))]
  vars <- sub("^sd_", "", sd_names)
  
  # costruisci ordine atteso
  want <- c(paste0(setdiff(base_order, "frailty"), "_", time))
  if (include_frailty) want <- c(want, "frailty")
  
  # tieni solo quelli presenti, nell'ordine desiderato
  want[want %in% vars]
}

make_sd_vec_for_labs <- function(est, labs) {
  sd_names <- paste0("sd_", labs)
  out <- unname(est[sd_names])
  names(out) <- labs
  out
}

#parser comune
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
    
    if (any(!ok)) {
      warning("param corr malformati ignorati: ",
              paste(raw[!ok], collapse = ", "))
    }
    
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


# Ordine "standard" (EQ-5D poi KOOS). Modifica se vuoi.
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

make_cross_RP_from_tab <- function(tab, rows, cols,
                                   p_col = c("p_fisher","p"),
                                   set_oob_to_NA = TRUE) {
  p_col <- p_col[p_col %in% names(tab)][1]
  if (is.na(p_col)) stop("tab non ha né p_fisher né p.")
  
  # prendi solo correlazioni
  tt <- tab[tab$type == "corr", , drop = FALSE]
  
  # se A/B non ci sono, li ricaviamo da param (solo per riempire matrici)
  if (!all(c("A","B") %in% names(tt))) {
    AB <- do.call(rbind, strsplit(sub("^corr_", "", tt$param), "__", fixed = TRUE))
    tt$A <- AB[,1]; tt$B <- AB[,2]
  }
  
  nr <- length(rows); nc <- length(cols)
  R <- matrix(NA_real_, nr, nc, dimnames = list(rows, cols))
  P <- matrix(NA_real_, nr, nc, dimnames = list(rows, cols))
  
  # tieni solo coppie che stanno nel rettangolo (anche invertite)
  keep <- (tt$A %in% rows & tt$B %in% cols) | (tt$A %in% cols & tt$B %in% rows)
  tt <- tt[keep, , drop = FALSE]
  
  # identifica out-of-bounds
  oob <- rep(FALSE, nrow(tt))
  if ("out_of_bounds" %in% names(tt)) oob <- oob | isTRUE(tt$out_of_bounds)
  oob <- oob | (!is.na(tt$est) & abs(tt$est) > 1)
  
  # riempi le matrici
  for (k in seq_len(nrow(tt))) {
    a <- tt$A[k]; b <- tt$B[k]
    r <- tt$est[k]
    p <- tt[[p_col]][k]
    
    # se stava invertita, ribalta
    if (a %in% cols && b %in% rows) {
      tmp <- a; a <- b; b <- tmp
    }
    
    if (set_oob_to_NA && oob[k]) {
      r <- NA_real_
      p <- NA_real_
    }
    
    R[a, b] <- r
    P[a, b] <- p
  }
  
  # celle OOB (per disegnare una X anche se mettiamo NA)
  oob_cells <- tt[oob, c("A","B","est",p_col), drop = FALSE]
  # normalizza orientamento (riga=rows, col=cols)
  flip <- oob_cells$A %in% cols & oob_cells$B %in% rows
  if (any(flip)) {
    tmp <- oob_cells$A[flip]
    oob_cells$A[flip] <- oob_cells$B[flip]
    oob_cells$B[flip] <- tmp
  }
  oob_cells <- oob_cells[oob_cells$A %in% rows & oob_cells$B %in% cols, , drop = FALSE]
  names(oob_cells)[names(oob_cells)==p_col] <- "p_used"
  
  list(R = R, P = P, oob_cells = oob_cells)
}


panel_title <- function(tag, title, add = add_abc) {
  if (!add) return(title)
  sprintf("(%s) %s", tag, title)
}


########## wald su corr
make_frailty_corr_tab_from_est <- function(est_g,
                                           times = c("I0","Ipost"),
                                           frailty = "frailty",
                                           drop_oob = TRUE) {
  nm <- names(est_g)
  if (is.null(nm)) stop("est_g deve essere un named numeric.")
  
  # prendi solo parametri corr_*
  corr_names <- grep("^corr_", nm, value = TRUE)
  if (length(corr_names) == 0) {
    warning("Nessun parametro che inizia con 'corr_' in est_g.")
    return(data.frame())
  }
  
  # split corr_<A>__<B>
  parse_one <- function(pn) {
    s <- sub("^corr_", "", pn)
    ab <- strsplit(s, "__", fixed = TRUE)[[1]]
    if (length(ab) != 2) return(NULL)
    
    A <- ab[1]; B <- ab[2]
    if (A == frailty && B != frailty) {
      other <- B
    } else if (B == frailty && A != frailty) {
      other <- A
    } else {
      return(NULL) # non è una corr con frailty
    }
    
    # time = suffisso dell'altro (se c'è)
    time <- if (grepl("_(I0|Ipost)$", other)) sub("^.*_(I0|Ipost)$", "\\1", other) else NA_character_
    if (!is.na(time) && !(time %in% times)) return(NULL)
    
    item_base <- if (!is.na(time)) sub("_(I0|Ipost)$", "", other) else other
    
    est <- unname(est_g[pn])
    
    oob <- is.na(est) || !is.finite(est) || abs(est) > 1
    if (drop_oob && oob) {
      est <- NA_real_
    }
    
    data.frame(
      param     = pn,          # NOME ORIGINALE COMPLETO
      type      = "corr_frailty",
      item      = item_base,
      time      = time,
      other     = other,       # es. Mobility_I0
      est       = est,
      out_of_bounds = oob,
      stringsAsFactors = FALSE
    )
  }
  
  out <- do.call(rbind, lapply(corr_names, parse_one))
  if (is.null(out) || nrow(out) == 0) {
    warning("Nessuna correlazione con frailty trovata (pattern corr_*__frailty o corr_frailty__*).")
    return(data.frame())
  }
  
  # ordina bene: time poi item
  out <- out[order(out$time, out$item, out$param), , drop = FALSE]
  rownames(out) <- NULL
  out
}


wald_global_corr0 <- function(corr_tab, est_g, V_g,
                              transform = c("fisher_z","r"),
                              rank_tol = 1e-10) {
  transform <- match.arg(transform)
  
  if (is.null(corr_tab) || nrow(corr_tab) == 0) {
    return(list(W = NA_real_, df = 0L, p = NA_real_,
                params = character(0), est = numeric(0), V = NULL))
  }
  
  stopifnot(is.matrix(V_g), !is.null(rownames(V_g)), !is.null(colnames(V_g)))
  nm <- names(est_g)
  if (is.null(nm)) stop("est_g deve essere named numeric.")
  
  # allinea V_g a est_g (nomi)
  V_g <- V_g[nm, nm, drop = FALSE]
  
  pars <- corr_tab$param
  pars <- unique(pars)
  miss <- setdiff(pars, nm)
  if (length(miss)) stop("Questi parametri non sono in est_g: ", paste(miss, collapse = ", "))
  
  r_hat <- as.numeric(est_g[pars])
  ok <- is.finite(r_hat)
  
  # se vuoi test su Fisher z, devi stare in (-1,1)
  if (transform == "fisher_z") {
    ok <- ok & (abs(r_hat) < 1)
  }
  
  pars <- pars[ok]
  r_hat <- r_hat[ok]
  
  if (length(pars) == 0) {
    return(list(W = NA_real_, df = 0L, p = NA_real_,
                params = character(0), est = numeric(0), V = NULL))
  }
  
  V_r <- V_g[pars, pars, drop = FALSE]
  V_r <- 0.5 * (V_r + t(V_r))
  
  if (transform == "fisher_z") {
    # z = atanh(r), J = diag(1/(1-r^2))
    z_hat <- atanh(r_hat)
    J <- diag(1 / (1 - r_hat^2), nrow = length(r_hat))
    V_t <- J %*% V_r %*% t(J)
    t_hat <- z_hat
  } else {
    V_t <- V_r
    t_hat <- r_hat
  }
  
  V_t <- 0.5 * (V_t + t(V_t))
  
  eig <- eigen(V_t, symmetric = TRUE, only.values = TRUE)$values
  df_eff <- sum(eig > max(eig) * rank_tol)
  
  if (!requireNamespace("MASS", quietly = TRUE)) stop("Installa MASS.")
  Vinv <- tryCatch(solve(V_t), error = function(e) MASS::ginv(V_t))
  
  W <- as.numeric(t(t_hat) %*% Vinv %*% t_hat)
  p <- pchisq(W, df = df_eff, lower.tail = FALSE)
  
  list(W = W, df = df_eff, p = p, params = pars, est = t_hat, V = V_t)
}

subset_frailty_corrtab <- function(corr_tab, items = NULL, time = NULL) {
  out <- corr_tab
  if (!is.null(items)) out <- out[out$item %in% items, , drop = FALSE]
  if (!is.null(time))  out <- out[out$time %in% time, , drop = FALSE]
  out
}

make_test_row <- function(label, corr_tab_sub, est_g, V_g,
                          transform = c("fisher_z","r")) {
  transform <- match.arg(transform)
  r <- tryCatch(
    wald_global_corr0(corr_tab_sub, est_g = est_g, V_g = V_g, transform = transform),
    error = function(e) NULL
  )
  if (is.null(r)) {
    return(data.frame(test = label, n = NA_integer_, W = NA_real_, df = NA_real_, p = NA_real_))
  }
  data.frame(
    test = label,
    n    = length(r$params),
    W    = r$W,
    df   = r$df,
    p    = r$p,
    stringsAsFactors = FALSE
  )
}

wald_tests_frailty_corr <- function(est_g, V_g,
                                    EQ_items, KOOS_items,
                                    times = c("I0","Ipost"),
                                    transform = c("fisher_z","r"),
                                    drop_oob = TRUE) {
  transform <- match.arg(transform)
  
  corr_tab <- make_frailty_corr_tab_from_est(
    est_g = est_g, times = times, frailty = "frailty", drop_oob = drop_oob
  )
  
  # split per tempo
  corr_I0    <- subset_frailty_corrtab(corr_tab, time = "I0")
  corr_Ipost <- subset_frailty_corrtab(corr_tab, time = "Ipost")
  
  # split per questionario
  corr_EQ_all    <- subset_frailty_corrtab(corr_tab, items = EQ_items)
  corr_KOOS_all  <- subset_frailty_corrtab(corr_tab, items = KOOS_items)
  
  corr_EQ_I0     <- subset_frailty_corrtab(corr_tab, items = EQ_items,   time = "I0")
  corr_EQ_Ipost  <- subset_frailty_corrtab(corr_tab, items = EQ_items,   time = "Ipost")
  corr_KOOS_I0   <- subset_frailty_corrtab(corr_tab, items = KOOS_items, time = "I0")
  corr_KOOS_Ipost<- subset_frailty_corrtab(corr_tab, items = KOOS_items, time = "Ipost")
  
  tests_df <- do.call(rbind, list(
    make_test_row("All (I0+Ipost)", corr_tab,     est_g, V_g, transform),
    
    make_test_row("I0 (all items)",    corr_I0,    est_g, V_g, transform),
    make_test_row("Ipost (all items)", corr_Ipost, est_g, V_g, transform),
    
    make_test_row("EQ (I0+Ipost)",   corr_EQ_all,   est_g, V_g, transform),
    make_test_row("EQ (I0)",         corr_EQ_I0,    est_g, V_g, transform),
    make_test_row("EQ (Ipost)",      corr_EQ_Ipost, est_g, V_g, transform),
    
    make_test_row("KOOS (I0+Ipost)", corr_KOOS_all,   est_g, V_g, transform),
    make_test_row("KOOS (I0)",       corr_KOOS_I0,    est_g, V_g, transform),
    make_test_row("KOOS (Ipost)",    corr_KOOS_Ipost, est_g, V_g, transform)
  ))
  rownames(tests_df) <- NULL
  
  list(
    corr_tab = corr_tab,
    tests_df = tests_df
  )
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

make_test_row <- function(label, cov_tab_sub, psi, V = NULL, se = NULL) {
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
    make_test_row("All_frailty_cov",       cov_tab,       psi_global, V_global, se_global),
    make_test_row("I0_frailty_cov",        cov_I0,        psi_global, V_global, se_global),
    make_test_row("Ipost_frailty_cov",     cov_Ipost,     psi_global, V_global, se_global),
    
    make_test_row("EQ_all_frailty_cov",    cov_EQ_all,    psi_global, V_global, se_global),
    make_test_row("EQ_I0_frailty_cov",     cov_EQ_I0,     psi_global, V_global, se_global),
    make_test_row("EQ_Ipost_frailty_cov",  cov_EQ_Ipost,  psi_global, V_global, se_global),
    
    make_test_row("KOOS_all_frailty_cov",  cov_KOOS_all,  psi_global, V_global, se_global),
    make_test_row("KOOS_I0_frailty_cov",   cov_KOOS_I0,   psi_global, V_global, se_global),
    make_test_row("KOOS_Ipost_frailty_cov",cov_KOOS_Ipost,psi_global, V_global, se_global)
  ))
  
  rownames(tests_df) <- NULL
  
  list(
    cov_tab = cov_tab,   # <-- qui hai param = nome originale intero
    tests   = tests_df
  )
}

### 
make_abbr2 <- function(x) {
  x0 <- gsub("[^A-Za-z]", "", x)
  caps <- regmatches(x0, gregexpr("[A-Z]", x0, perl = TRUE))[[1]]
  caps <- paste(caps, collapse = "")
  
  if (nchar(caps) >= 2) return(substr(caps, 1, 2))
  if (nchar(x0) >= 2) return(toupper(substr(x0, 1, 2)))
  toupper(x0)
}

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

pca_plot_cor_nice <- function(R, main = "", labels = NULL,
                              short_map_base = NULL,
                              drop_suffix = FALSE,
                              suffix_regex = "(_I0|_Ipost)$",
                              common_from = NULL,
                              xlim_fixed = NULL,
                              ylim_fixed = NULL,
                              manual_nudge = NULL,
                              label_size = 3.9,
                              label_style = c("abbr2", "full"),
                              eq_items = NULL,
                              koos_items = NULL,
                              eq_col = "#6D28D9",
                              koos_col = "#475569",
                              other_col = "#6B7280",
                              leader_col = "grey45",
                              force = 0.08,
                              force_pull = 6,
                              box_padding = 0.10,
                              point_padding = 0.02,
                              axis_title_size = 11.5,
                              axis_text_size = 10.2,
                              plot_title_size = 12.5) {
  
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("serve il pacchetto 'ggrepel'")
  }
  
  label_style <- match.arg(label_style)
  
  R <- as.matrix(R)
  if (nrow(R) != ncol(R)) stop("R deve essere quadrata.")
  if (anyNA(R)) stop("R contiene NA.")
  if (max(abs(R - t(R))) > 1e-10) stop("R deve essere simmetrica.")
  if (is.null(rownames(R))) rownames(R) <- seq_len(nrow(R))
  if (is.null(colnames(R))) colnames(R) <- rownames(R)
  
  raw_labels <- rownames(R)
  raw_base <- sub(suffix_regex, "", raw_labels, perl = TRUE)
  
  eg <- eigen(R, symmetric = TRUE)
  ord <- order(eg$values, decreasing = TRUE)
  vals <- eg$values[ord]
  vecs <- eg$vectors[, ord, drop = FALSE]
  
  load <- vecs %*% diag(sqrt(pmax(vals, 0)))
  x <- load[, 1]
  y <- load[, 2]
  
  if (sum(x) < 0) x <- -x
  if (sum(y) < 0) y <- -y
  
  prop <- vals / sum(vals)
  
  if (is.null(labels)) labels <- raw_labels
  labels <- as.character(labels)
  
  if (!is.null(short_map_base)) {
    map_one <- function(lbl) {
      base <- sub(suffix_regex, "", lbl, perl = TRUE)
      if (base %in% names(short_map_base)) unname(short_map_base[base]) else base
    }
    labels <- vapply(labels, map_one, character(1))
  } else {
    labels <- sub(suffix_regex, "", labels, perl = TRUE)
  }
  
  if (isTRUE(drop_suffix)) {
    labels <- sub(suffix_regex, "", labels, perl = TRUE)
  }
  
  display_labels <- if (label_style == "abbr2") {
    vapply(labels, make_abbr2, character(1))
  } else {
    labels
  }
  
  group <- rep("Other", length(raw_base))
  if (!is.null(eq_items)) group[raw_base %in% eq_items] <- "EQ"
  if (!is.null(koos_items)) group[raw_base %in% koos_items] <- "KOOS"
  group <- factor(group, levels = c("EQ", "KOOS", "Other"))
  
  df <- data.frame(
    raw = raw_base,
    label_full = labels,
    label = display_labels,
    PC1 = x,
    PC2 = y,
    group = group,
    stringsAsFactors = FALSE
  )
  
  df$r <- sqrt(df$PC1^2 + df$PC2^2)
  df$r[df$r < 1e-8] <- 1e-8
  df$ux <- df$PC1 / df$r
  df$uy <- df$PC2 / df$r
  df$nx <- -df$uy
  df$ny <-  df$ux
  df$theta <- atan2(df$PC2, df$PC1)
  
  ref_lim <- max(abs(c(df$PC1, df$PC2, 0)))
  
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
  
  grp <- paste0(ifelse(df$PC1 >= 0, "R", "L"), ifelse(df$PC2 >= 0, "U", "D"))
  df$k <- 0
  for (g in unique(grp)) {
    id <- which(grp == g)
    oo <- order(df$theta[id], df$r[id])
    df$k[id[oo]] <- make_k(length(id))
  }
  
  push0 <- 0.012 * ref_lim
  perp0 <- 0.008 * ref_lim
  
  df$nudge_x <- df$ux * push0 + df$nx * (df$k * perp0)
  df$nudge_y <- df$uy * push0 + df$ny * (df$k * perp0)
  df$hjust <- ifelse(df$PC1 >= 0, 0, 1)
  
  if (!is.null(manual_nudge)) {
    for (nm in names(manual_nudge)) {
      id <- which(df$label == nm | df$label_full == nm | df$raw == nm)
      if (length(id) == 1) {
        nud <- manual_nudge[[nm]]
        if (length(nud) == 2) {
          df$nudge_x[id] <- df$nudge_x[id] + nud[1]
          df$nudge_y[id] <- df$nudge_y[id] + nud[2]
        }
      }
    }
  }
  
  if (!is.null(xlim_fixed) && !is.null(ylim_fixed)) {
    xlim <- xlim_fixed
    ylim <- ylim_fixed
  } else if (!is.null(common_from)) {
    all_x <- 0
    all_y <- 0
    
    for (M in common_from) {
      M <- as.matrix(M)
      ee <- eigen(M, symmetric = TRUE)
      oo <- order(ee$values, decreasing = TRUE)
      vv <- ee$values[oo]
      UU <- ee$vectors[, oo, drop = FALSE]
      LL <- UU %*% diag(sqrt(pmax(vv, 0)))
      xx <- LL[, 1]
      yy <- LL[, 2]
      if (sum(xx) < 0) xx <- -xx
      if (sum(yy) < 0) yy <- -yy
      all_x <- c(all_x, xx)
      all_y <- c(all_y, yy)
    }
    
    ref_all <- max(abs(c(all_x, all_y)))
    xlim <- range(all_x, 0) + c(-0.06, 0.16) * ref_all
    ylim <- range(all_y, 0) + c(-0.10, 0.10) * ref_all
  } else {
    xr <- range(c(0, df$PC1), na.rm = TRUE)
    yr <- range(c(0, df$PC2), na.rm = TRUE)
    xpad <- 0.08 * diff(xr)
    ypad <- 0.08 * diff(yr)
    if (!is.finite(xpad) || xpad == 0) xpad <- 0.08 * ref_lim
    if (!is.finite(ypad) || ypad == 0) ypad <- 0.08 * ref_lim
    xlim <- xr + c(-xpad, xpad)
    ylim <- yr + c(-ypad, ypad)
  }
  
  pct_lab <- function(txt, p) sprintf("%s (%.1f\\%%)", txt, 100 * p)
  
  cols <- c(EQ = eq_col, KOOS = koos_col, Other = other_col)
  
  library(ggplot2)
  
  ggplot(df, aes(PC1, PC2, color = group)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey88") +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey88") +
    geom_segment(
      aes(x = 0, y = 0, xend = PC1, yend = PC2),
      linewidth = 0.8,
      arrow = arrow(type = "closed", length = grid::unit(0.18, "cm")),
      show.legend = FALSE
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
      segment.color = leader_col,
      segment.size = 0.28,
      max.overlaps = Inf,
      seed = 123,
      show.legend = FALSE
    ) +
    scale_color_manual(values = cols) +
    coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE, clip = "on") +
    labs(
      title = main,
      x = pct_lab("First Principal Component", prop[1]),
      y = pct_lab("Second Principal Component", prop[2])
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = plot_title_size),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, color = "grey92"),
      axis.title = element_text(face = "bold", size = axis_title_size),
      axis.text = element_text(size = axis_text_size),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      legend.position = "none",
      plot.margin = margin(8, 10, 8, 10)
    )
}


pca_plot_cor_full <- function(R,
                              main = "PCA random effects",
                              short_map,
                              eq_base,
                              koos_base,
                              label_size = 4.0,
                              axis_title_size = 12,
                              axis_text_size = 11,
                              plot_title_size = 13,
                              legend_text_size = 10.5,
                              force = 0.10,
                              force_pull = 5,
                              box_padding = 0.12,
                              point_padding = 0.02,
                              escape_tex = FALSE) {
  
  esc_tex <- function(x) {
    gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  }
  
  make_short_label_with_suffix <- function(x, short_map, suffix_regex = "(_I0|_Ipost)$") {
    if (identical(x, "frailty")) return("Frailty")
    base <- sub(suffix_regex, "", x, perl = TRUE)
    suff <- sub(paste0("^", base), "", x)
    short <- if (base %in% names(short_map)) unname(short_map[base]) else base
    paste0(short, suff)
  }
  
  get_group_full <- function(x, eq_base, koos_base) {
    if (identical(x, "frailty")) return("Frailty")
    base <- sub("(_I0|_Ipost)$", "", x, perl = TRUE)
    
    if (base %in% eq_base && grepl("_I0$", x)) return("EQ-5D_I0")
    if (base %in% eq_base && grepl("_Ipost$", x)) return("EQ-5D_Ipost")
    if (base %in% koos_base && grepl("_I0$", x)) return("KOOS_I0")
    if (base %in% koos_base && grepl("_Ipost$", x)) return("KOOS_Ipost")
    
    NA_character_
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
  
  pct_lab <- function(txt, p) sprintf("%s (%.1f%%)", txt, 100 * p)
  pct_lab_tex <- function(txt, p) sprintf("%s (%.1f\\%%)", txt, 100 * p)
  
  R <- as.matrix(R)
  
  if (nrow(R) != ncol(R)) stop("R deve essere quadrata.")
  if (anyNA(R)) stop("R contiene NA.")
  if (max(abs(R - t(R))) > 1e-10) stop("R deve essere simmetrica.")
  if (is.null(rownames(R))) stop("servono i rownames di R")
  
  eg <- eigen(R, symmetric = TRUE)
  ord <- order(eg$values, decreasing = TRUE)
  vals <- eg$values[ord]
  vecs <- eg$vectors[, ord, drop = FALSE]
  
  load <- vecs %*% diag(sqrt(pmax(vals, 0)))
  x <- load[, 1]
  y <- load[, 2]
  
  if (sum(x) < 0) x <- -x
  if (sum(y) < 0) y <- -y
  
  prop <- vals / sum(vals)
  
  labs_raw <- rownames(R)
  
  labs_short <- vapply(
    labs_raw,
    make_short_label_with_suffix,
    character(1),
    short_map = short_map
  )
  
  groups <- vapply(
    labs_raw,
    get_group_full,
    character(1),
    eq_base = eq_base,
    koos_base = koos_base
  )
  
  if (anyNA(groups)) {
    stop(
      paste(
        "variabili non classificate:",
        paste(labs_raw[is.na(groups)], collapse = ", ")
      )
    )
  }
  
  if (escape_tex) {
    labs_short <- esc_tex(labs_short)
    main <- esc_tex(main)
    legend_labels <- c(
      "EQ-5D_I0"    = "EQ-5D\\_I0",
      "EQ-5D_Ipost" = "EQ-5D\\_Ipost",
      "KOOS_I0"     = "KOOS\\_I0",
      "KOOS_Ipost"  = "KOOS\\_Ipost",
      "Frailty"     = "Frailty"
    )
    x_lab <- pct_lab_tex("First Principal Component", prop[1])
    y_lab <- pct_lab_tex("Second Principal Component", prop[2])
  } else {
    legend_labels <- c(
      "EQ-5D_I0"    = "EQ-5D_I0",
      "EQ-5D_Ipost" = "EQ-5D_Ipost",
      "KOOS_I0"     = "KOOS_I0",
      "KOOS_Ipost"  = "KOOS_Ipost",
      "Frailty"     = "Frailty"
    )
    x_lab <- pct_lab("First Principal Component", prop[1])
    y_lab <- pct_lab("Second Principal Component", prop[2])
  }
  
  used_levels <- c("EQ-5D_I0", "EQ-5D_Ipost", "KOOS_I0", "KOOS_Ipost", "Frailty")
  used_levels <- used_levels[used_levels %in% unique(groups)]
  
  df <- data.frame(
    raw = labs_raw,
    label = labs_short,
    group = factor(groups, levels = used_levels),
    PC1 = x,
    PC2 = y,
    stringsAsFactors = FALSE
  )
  
  df$r <- sqrt(df$PC1^2 + df$PC2^2)
  df$r[df$r < 1e-8] <- 1e-8
  df$ux <- df$PC1 / df$r
  df$uy <- df$PC2 / df$r
  df$nx <- -df$uy
  df$ny <-  df$ux
  df$theta <- atan2(df$PC2, df$PC1)
  
  ref_lim <- max(abs(c(df$PC1, df$PC2, 0)))
  
  grp_quad <- paste0(ifelse(df$PC1 >= 0, "R", "L"), ifelse(df$PC2 >= 0, "U", "D"))
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
  
  xr <- range(c(0, df$PC1), na.rm = TRUE)
  yr <- range(c(0, df$PC2), na.rm = TRUE)
  xpad <- 0.10 * diff(xr)
  ypad <- 0.10 * diff(yr)
  
  if (!is.finite(xpad) || xpad == 0) xpad <- 0.10 * ref_lim
  if (!is.finite(ypad) || ypad == 0) ypad <- 0.10 * ref_lim
  
  xlim <- xr + c(-xpad, xpad)
  ylim <- yr + c(-ypad, ypad)
  
  cols <- c(
    "EQ-5D_I0"    = "#0F766E",
    "EQ-5D_Ipost" = "#d1a00d",
    "KOOS_I0"     = "#5B21B6",
    "KOOS_Ipost"  = "#52340b",
    "Frailty"     = "#C2410C"
  )
  
  ggplot(df, aes(PC1, PC2, color = group)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey88") +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey88") +
    geom_segment(
      aes(x = 0, y = 0, xend = PC1, yend = PC2),
      linewidth = 0.85,
      arrow = arrow(type = "closed", length = unit(0.18, "cm"))
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
      values = cols[used_levels],
      breaks = used_levels,
      labels = legend_labels[used_levels],
      drop = TRUE
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
      legend.position = "bottom",
      legend.text = element_text(size = legend_text_size),
      plot.margin = margin(8, 10, 8, 10)
    ) +
    guides(
      color = guide_legend(
        nrow = 2,
        byrow = TRUE,
        override.aes = list(linewidth = 1.3)
      )
    )
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
  if (base %in% eq_base && grepl("_I0$", x)) return("EQ-5D_I0")
  if (base %in% eq_base && grepl("_Ipost$", x)) return("EQ-5D_Ipost")
  if (base %in% koos_base && grepl("_I0$", x)) return("KOOS_I0")
  if (base %in% koos_base && grepl("_Ipost$", x)) return("KOOS_Ipost")
  NA_character_
}

group_resid <- function(x, eq_base, koos_base) {
  base <- sub("(_I0|_Ipost)$", "", x, perl = TRUE)
  if (base %in% eq_base) return("EQ")
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
  grp_quad <- paste0(ifelse(df$PC1 >= 0, "R", "L"), ifelse(df$PC2 >= 0, "U", "D"))
  
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

##### Table fixed effects longitudinal (no Surv)

override_named <- function(base, add) {
  base[names(add)] <- add
  base
}

shiftcov_to_wide <- function(tab_shiftcov_all,
                             item_levels = NULL,
                             prefix = "covDiffFrailty") {
  stopifnot(is.data.frame(tab_shiftcov_all))
  stopifnot(all(c("param","est","se","p") %in% names(tab_shiftcov_all)))
  
  out <- data.frame(
    item = sub("^covdiff_frailty__", "", tab_shiftcov_all$param),
    estimate = tab_shiftcov_all$est,
    se = tab_shiftcov_all$se,
    p_value = tab_shiftcov_all$p,
    sig_0.05 = tab_shiftcov_all$sig_0.05,
    stringsAsFactors = FALSE
  )
  
  names(out)[names(out) == "estimate"] <- paste0(prefix, "_estimate")
  names(out)[names(out) == "se"]       <- paste0(prefix, "_se")
  names(out)[names(out) == "p_value"]  <- paste0(prefix, "_p_value")
  names(out)[names(out) == "sig_0.05"] <- paste0(prefix, "_sig_0.05")
  
  if (!is.null(item_levels)) {
    out$item <- factor(out$item, levels = item_levels)
    out <- out[order(out$item), , drop = FALSE]
  }
  
  rownames(out) <- NULL
  out
}

extract_frailty_corr_wide <- function(psi, se,
                                      alpha = 0.05,
                                      item_levels = NULL,
                                      corr_levels = c("I0","Ipost"),
                                      corr_prefix = "corrFrailty") {
  stopifnot(is.numeric(psi), is.numeric(se),
            !is.null(names(psi)), !is.null(names(se)))
  
  nm <- intersect(names(psi), names(se))
  psi <- psi[nm]; se <- se[nm]
  
  cand <- grep("^corr_", nm, value = TRUE)
  cand <- cand[grepl("frailty", cand, fixed = TRUE)]
  if (length(cand) == 0) {
    # ritorna solo item vuoto coerente
    out <- data.frame(item = character(0), stringsAsFactors = FALSE)
    return(out)
  }
  
  parse_one <- function(nm1) {
    s <- sub("^corr_", "", nm1)
    parts <- strsplit(s, "__", fixed = TRUE)[[1]]
    if (length(parts) != 2) return(NULL)
    
    a <- parts[1]; b <- parts[2]
    if (a == "frailty") other <- b
    else if (b == "frailty") other <- a
    else return(NULL)
    
    # other deve finire con _I0 o _Ipost
    if (!grepl("_(I0|Ipost)$", other)) return(NULL)
    lev  <- sub("^.*_(I0|Ipost)$", "\\1", other)
    item <- sub("_(I0|Ipost)$", "", other)
    
    list(item = item, re_level = lev)
  }
  
  parsed <- lapply(cand, parse_one)
  ok <- !vapply(parsed, is.null, logical(1))
  cand <- cand[ok]
  parsed <- parsed[ok]
  if (length(cand) == 0) {
    return(data.frame(item = character(0), stringsAsFactors = FALSE))
  }
  
  df <- data.frame(
    item = vapply(parsed, `[[`, "", "item"),
    re_level = vapply(parsed, `[[`, "", "re_level"),
    estimate = as.numeric(psi[cand]),
    se = as.numeric(se[cand]),
    stringsAsFactors = FALSE
  )
  
  # filtra livelli richiesti
  df <- df[df$re_level %in% corr_levels, , drop = FALSE]
  
  df$z <- df$estimate / df$se
  df$p_value <- 2 * pnorm(-abs(df$z))
  df$sig_0.05 <- !is.na(df$p_value) & df$p_value < alpha
  
  if (!is.null(item_levels)) df$item <- factor(df$item, levels = item_levels)
  
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Serve il pacchetto 'tidyr'.")
  
  tidyr::pivot_wider(
    df[, c("item","re_level","estimate","se","p_value","sig_0.05")],
    names_from = re_level,
    values_from = c(estimate, se, p_value, sig_0.05),
    names_glue = paste0(corr_prefix, "_{re_level}_{.value}")
  )
}

make_wald_table_3cols <- function(psi, se,
                                  alpha = 0.05,
                                  predictors = c("I6","I12","wI6","wI12","pI6","pI12",
                                                 "diag_oa","sex01","loghosp","age"),
                                  item_levels = NULL,
                                  exclude_regex = "^fixed_surv_",
                                  include_corr_frailty = TRUE,
                                  corr_levels = c("I0","Ipost"),
                                  corr_prefix = "corrFrailty") {
  stopifnot(is.numeric(psi), is.numeric(se), !is.null(names(psi)), !is.null(names(se)))
  
  nm_all <- intersect(names(psi), names(se))
  psi <- psi[nm_all]
  se  <- se[nm_all]
  
  # 1) solo fixed_* e (di default) ESCLUDI fixed_surv_*
  keep <- grepl("^fixed_", nm_all) & !grepl(exclude_regex, nm_all)
  psi_fix <- psi[keep]
  nm <- names(psi_fix)
  
  se_fix <- se[nm]
  
  # 2) parse robusto: fixed_<ITEM>_<PRED>
  pred_pat <- paste(predictors, collapse = "|")
  re <- paste0("^fixed_(.*)_(", pred_pat, ")$")
  
  m <- regexec(re, nm)
  mm <- regmatches(nm, m)
  ok <- lengths(mm) == 3
  
  nm_ok <- nm[ok]
  mm_ok <- mm[ok]
  
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
  df$sig_0.05 <- !is.na(df$p_value) & df$p_value < alpha
  
  if (!is.null(item_levels)) df$item <- factor(df$item, levels = item_levels)
  
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Serve il pacchetto 'tidyr'.")
  
  wide_fix <- tidyr::pivot_wider(
    df[, c("item","predictor","estimate","se","p_value","sig_0.05")],
    names_from = predictor,
    values_from = c(estimate, se, p_value, sig_0.05),
    names_glue = "{predictor}_{.value}"
  )
  
  # 3) aggiungi corr frailty–random effects (I0/Ipost) per item
  if (isTRUE(include_corr_frailty)) {
    wide_corr <- extract_frailty_corr_wide(
      psi = psi, se = se,
      alpha = alpha,
      item_levels = item_levels,
      corr_levels = corr_levels,
      corr_prefix = corr_prefix
    )
    
    # merge by item
    if (nrow(wide_corr) > 0) {
      wide <- merge(wide_fix, wide_corr, by = "item", all.x = TRUE, sort = FALSE)
    } else {
      wide <- wide_fix
    }
  } else {
    wide <- wide_fix
  }
  
  # 4) ordina righe se factor
  if (!is.null(item_levels)) {
    wide$item <- factor(wide$item, levels = item_levels)
    wide <- wide[order(wide$item), , drop = FALSE]
  } else {
    wide <- wide[order(wide$item), , drop = FALSE]
  }
  rownames(wide) <- NULL
  
  wide
}

make_latex_beta_table <- function(tab_items3,
                                  predictors = c("I6","I12","wI6","wI12","pI6","pI12",
                                                 "diag_oa","sex01","loghosp","age"),
                                  extra_terms = NULL,      # lista: list(list(prefix="...", header="..."), ...)
                                  item_name_map = NULL,
                                  alpha = 0.05,
                                  digits_est = 2,
                                  digits_se  = 2,
                                  caption = NULL,
                                  label = NULL,
                                  use_resizebox = TRUE,
                                  font_size_cmd = "\\small") {
  
  stopifnot(is.data.frame(tab_items3), "item" %in% names(tab_items3))
  
  # termini = predictors + extra_terms
  term_prefixes <- predictors
  term_headers  <- NULL
  
  # header predictors
  col_headers_pred <- vapply(predictors, function(pr) {
    pr2 <- escape_for_texttt(pr)
    paste0("$\\\\widehat{\\\\beta_{\\\\texttt{", pr2, "}}}$")
  }, character(1))
  
  term_headers <- col_headers_pred
  
  if (!is.null(extra_terms)) {
    stopifnot(is.list(extra_terms))
    extra_prefix <- vapply(extra_terms, `[[`, "", "prefix")
    extra_head   <- vapply(extra_terms, `[[`, "", "header")
    
    term_prefixes <- c(term_prefixes, extra_prefix)
    term_headers  <- c(term_headers,  extra_head)
  }
  
  # controlla colonne richieste
  need <- unlist(lapply(term_prefixes, function(pr)
    c(paste0(pr,"_estimate"), paste0(pr,"_se"), paste0(pr,"_p_value"))
  ))
  miss <- setdiff(need, names(tab_items3))
  if (length(miss) > 0) stop("Colonne mancanti in tab_items3: ", paste(miss, collapse = ", "))
  
  # nomi item (mappati se richiesto)
  item_raw <- as.character(tab_items3$item)
  item_disp <- item_raw
  if (!is.null(item_name_map)) {
    if (!is.null(names(item_name_map))) {
      hit <- item_raw %in% names(item_name_map)
      item_disp[hit] <- unname(item_name_map[item_raw[hit]])
    } else {
      stop("item_name_map deve essere un named vector: names=original item, values=label.")
    }
  }
  item_disp <- latex_escape(item_disp)
  
  fmt <- function(x, d) ifelse(is.na(x), "", formatC(x, digits = d, format = "f"))
  
  n <- nrow(tab_items3)
  m <- length(term_prefixes)
  cells <- matrix("", nrow = n, ncol = m)
  
  for (j in seq_along(term_prefixes)) {
    pr <- term_prefixes[j]
    est <- tab_items3[[paste0(pr,"_estimate")]]
    se  <- tab_items3[[paste0(pr,"_se")]]
    pv  <- tab_items3[[paste0(pr,"_p_value")]]
    
    star <- ifelse(!is.na(pv) & pv < alpha, "$^{*}$", "")
    cells[, j] <- ifelse(is.na(est) | is.na(se), "",
                         paste0(fmt(est, digits_est), star, " (", fmt(se, digits_se), ")"))
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


latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%#$&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x
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

# test singolo 
shiftcov_table <- function(res_wald, alpha = 0.05) {
  est <- res_wald$est
  Vc <- res_wald$V
  se <- sqrt(pmax(diag(Vc), 0))
  z <- est / se
  p <- 2 * pnorm(-abs(z))
  
  out <- data.frame(
    param = names(est),
    est = as.numeric(est),
    se = as.numeric(se),
    z = as.numeric(z),
    p = as.numeric(p),
    sig_0.05 = as.numeric(p) < alpha,
    stringsAsFactors = FALSE
  )
  
  rownames(out) <- NULL
  out
}

##### joint per prost_uni
get_beta_params_prostuni <- function(psi,
                                     predictors = c("pI6", "pI12"),
                                     items = NULL) {
  stopifnot(is.numeric(psi), !is.null(names(psi)))
  
  nm <- names(psi)
  pred_pat <- paste(predictors, collapse = "|")
  pars <- grep(paste0("^fixed_.*_(", pred_pat, ")$"), nm, value = TRUE)
  
  if (!is.null(items)) {
    item_pat <- paste(items, collapse = "|")
    pars <- pars[grepl(paste0("^fixed_(", item_pat, ")_(", pred_pat, ")$"), pars)]
  }
  
  if (!length(pars)) {
    stop("nessun parametro trovato per i predittori: ", paste(predictors, collapse = ", "))
  }
  
  pars
}

wald_prostuni_betas <- function(psi, V = NULL, se = NULL,
                                predictors = c("pI6", "pI12"),
                                items = NULL,
                                alpha = 0.05) {
  pars <- get_beta_params_prostuni(
    psi = psi,
    predictors = predictors,
    items = items
  )
  
  res_global <- wald_joint(
    psi = psi,
    params = pars,
    V = V,
    se = se
  )
  
  if (!is.null(V)) {
    Vsub <- V[pars, pars, drop = FALSE]
    se_sub <- sqrt(pmax(diag(Vsub), 0))
  } else {
    se_sub <- se[pars]
  }
  
  est <- psi[pars]
  z <- est / se_sub
  p <- 2 * pnorm(-abs(z))
  
  tab <- data.frame(
    param = pars,
    est = as.numeric(est),
    se = as.numeric(se_sub),
    z = as.numeric(z),
    p = as.numeric(p),
    sig_0.05 = p < alpha,
    stringsAsFactors = FALSE
  )
  
  list(
    global = data.frame(
      test = "H0: tutti i beta associati a prost_uni (pI6,pI12) uguali a zero",
      n = length(pars),
      W = res_global$W,
      df = res_global$df,
      p = res_global$p,
      stringsAsFactors = FALSE
    ),
    table = tab[order(tab$p), , drop = FALSE],
    params = pars
  )
}

##### test per beta_{k,pI12}+beta_{k,pI6}=0 for all k
make_L_postmean_prostuni <- function(psi_names, items,
                                     pred1 = "pI6",
                                     pred2 = "pI12",
                                     average = TRUE) {
  w <- if (average) 0.5 else 1
  
  L <- matrix(0, nrow = length(items), ncol = length(psi_names),
              dimnames = list(paste0("postmean_prostuni__", items), psi_names))
  
  for (k in seq_along(items)) {
    it <- items[k]
    n1 <- paste0("fixed_", it, "_", pred1)
    n2 <- paste0("fixed_", it, "_", pred2)
    
    if (!(n1 %in% psi_names)) stop("parametro mancante: ", n1)
    if (!(n2 %in% psi_names)) stop("parametro mancante: ", n2)
    
    L[k, n1] <- w
    L[k, n2] <- w
  }
  
  L
}

wald_postmean_prostuni <- function(psi, V = NULL, se = NULL, items,
                                   pred1 = "pI6",
                                   pred2 = "pI12",
                                   average = TRUE) {
  L <- make_L_postmean_prostuni(
    psi_names = names(psi),
    items = items,
    pred1 = pred1,
    pred2 = pred2,
    average = average
  )
  
  wald_linear(
    psi = psi,
    L = L,
    V = V,
    se = se
  )
}

wald_contrast_table <- function(res_wald, alpha = 0.05) {
  est <- res_wald$est
  Vc <- res_wald$V
  se <- sqrt(pmax(diag(Vc), 0))
  z <- est / se
  p <- 2 * pnorm(-abs(z))
  
  out <- data.frame(
    param = names(est),
    est = as.numeric(est),
    se = as.numeric(se),
    z = as.numeric(z),
    p = as.numeric(p),
    sig_0.05 = p < alpha,
    stringsAsFactors = FALSE
  )
  
  rownames(out) <- NULL
  out
}


#### post-mean
make_L_posteq <- function(psi_names, items,
                          term6 = "I6", term12 = "I12") {
  L <- matrix(0, nrow = length(items), ncol = length(psi_names),
              dimnames = list(paste0("posteq__", items), psi_names))
  
  for (k in seq_along(items)) {
    it <- items[k]
    n6  <- paste0("fixed_", it, "_", term6)
    n12 <- paste0("fixed_", it, "_", term12)
    
    if (!(n6 %in% psi_names)) stop("parametro mancante: ", n6)
    if (!(n12 %in% psi_names)) stop("parametro mancante: ", n12)
    
    L[k, n12] <- 1
    L[k, n6]  <- -1
  }
  
  L
}

wald_posteq <- function(psi, V = NULL, se = NULL, items,
                        term6 = "I6", term12 = "I12") {
  L <- make_L_posteq(
    psi_names = names(psi),
    items = items,
    term6 = term6,
    term12 = term12
  )
  
  wald_linear(
    psi = psi,
    L = L,
    V = V,
    se = se
  )
}


############# joint log-likelihood for BIC #########################################################################################
#install.packages(c("mvtnorm","Matrix","dplyr"))
## =========================
## joint full likelihood + mc
## =========================

#dipendenze
library(dplyr)
library(Matrix)
library(MASS)

logmeanexp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

near_pd_chol <- function(S) {
  S2 <- as.matrix(nearPD(S, corr = FALSE)$mat)
  chol(S2)
}

#1) prepara dati survival con 
prepare_surv_from_dates <- function(dat, id_var="ID",
                                    admin_censor="2025-01-01",
                                    month_unit=30.4375,
                                    use_min_lastvisit=TRUE) {
  
  stopifnot(all(c("date_baseline","date_intervention","date_revision") %in% names(dat)))
  
  d <- dat %>%
    dplyr::mutate(
      ID_chr = trimws(as.character(.data[[id_var]])),
      base   = as.POSIXct(date_baseline),
      interv = as.POSIXct(date_intervention),
      rev    = as.POSIXct(date_revision)
    )
  
  admin_cens <- as.POSIXct(admin_censor, tz="UTC")
  
  surv_tbl <- d %>%
    dplyr::group_by(ID_chr) %>%
    dplyr::summarise(
      base0   = base[which(!is.na(base))[1]],
      interv0 = interv[which(!is.na(interv))[1]],
      rev0    = rev[which(!is.na(rev))[1]],
      
      wait_m = dplyr::if_else(!is.na(base0) & !is.na(interv0),
                              as.numeric(difftime(interv0, base0, units="days"))/month_unit,
                              NA_real_),
      
      last_date = suppressWarnings(max(c(base, interv), na.rm=TRUE)),
      .groups="drop"
    ) %>%
    dplyr::mutate(
      delta = as.integer(!is.na(rev0)),
      censor_date = dplyr::if_else(delta==1L, rev0, admin_cens),
      
      #se non sei sicuro che tutti siano seguiti fino a 01/2025, questo evita di “inventare” follow-up
      censor_date = dplyr::if_else(use_min_lastvisit & delta==0L & is.finite(as.numeric(last_date)),
                                   pmin(censor_date, last_date),
                                   censor_date),
      
      T = dplyr::if_else(!is.na(interv0) & !is.na(censor_date),
                         as.numeric(difftime(censor_date, interv0, units="days"))/month_unit,
                         NA_real_)
    )
  
  checks <- list(
    n_id = nrow(surv_tbl),
    n_event = sum(surv_tbl$delta==1L, na.rm=TRUE),
    n_T_bad = sum(!is.finite(surv_tbl$T) | surv_tbl$T <= 0, na.rm=TRUE),
    n_rev_before_interv = sum(!is.na(surv_tbl$rev0) & !is.na(surv_tbl$interv0) & surv_tbl$rev0 < surv_tbl$interv0)
  )
  
  list(surv_tbl=surv_tbl, checks=checks)
}

#2) parsing soglie (robusto)
build_thresholds <- function(psi, outcomes) {
  thr <- vector("list", length(outcomes)); names(thr) <- outcomes
  for (y in outcomes) {
    pat <- paste0("^threshold_", gsub("([\\^\\$\\.|\\+\\(\\)\\[\\]\\{\\}\\\\])","\\\\\\1", y), "_cut:(\\d+)$")
    idx <- grep(pat, names(psi))
    if (length(idx) == 0) stop(paste0("missing thresholds for ", y))
    cuts <- psi[idx]
    k <- as.integer(sub(paste0("threshold_", y, "_cut:"), "", names(cuts)))
    cuts <- cuts[order(k)]
    thr[[y]] <- as.numeric(cuts)
  }
  thr
}

#3) parsing fixed effects e costruzione covariate mancanti (wI6,wI12,pI6,pI12)
build_fixed <- function(psi, outcomes, dat) {
  fixed <- vector("list", length(outcomes)); names(fixed) <- outcomes
  
  for (y in outcomes) {
    pref <- paste0("^fixed_", gsub("([\\^\\$\\.|\\+\\(\\)\\[\\]\\{\\}\\\\])","\\\\\\1", y), "_")
    idx <- grep(pref, names(psi))
    if (length(idx) == 0) fixed[[y]] <- numeric(0) else {
      coefs <- psi[idx]
      covn <- sub(paste0("fixed_", y, "_"), "", names(coefs))
      names(coefs) <- covn
      fixed[[y]] <- as.numeric(coefs); names(fixed[[y]]) <- covn
    }
  }
  
  #covariate richieste globali
  need_cov <- sort(unique(unlist(lapply(fixed, names))))
  need_cov <- need_cov[need_cov != "(Intercept)" & !is.na(need_cov)]
  
  d <- dat
  #creazioni standard attese nel tuo modello
  if (!("wI6" %in% names(d)) && ("waiting" %in% names(d)) && ("I6" %in% names(d))) d$wI6 <- d$waiting * d$I6
  if (!("wI12" %in% names(d)) && ("waiting" %in% names(d)) && ("I12" %in% names(d))) d$wI12 <- d$waiting * d$I12
  
  if (!("pI6" %in% names(d)) && ("prost_uni" %in% names(d)) && ("I6" %in% names(d))) d$pI6 <- d$prost_uni * d$I6
  if (!("pI12" %in% names(d)) && ("prost_uni" %in% names(d)) && ("I12" %in% names(d))) d$pI12 <- d$prost_uni * d$I12
  
  missing_cov <- setdiff(need_cov, names(d))
  if (length(missing_cov) > 0) {
    stop(paste0("missing covariates in dat: ", paste(missing_cov, collapse=", ")))
  }
  
  list(fixed=fixed, dat=d, need_cov=need_cov)
}

#4) matrice di correlazione residua R da atanh_rho_resid__A__B
build_resid_R <- function(psi, outcomes) {
  K <- length(outcomes)
  R <- diag(K); rownames(R) <- outcomes; colnames(R) <- outcomes
  
  idx <- grep("^atanh_rho_resid__", names(psi))
  if (length(idx) > 0) {
    for (nm in names(psi)[idx]) {
      s <- sub("^atanh_rho_resid__", "", nm)
      parts <- strsplit(s, "__", fixed=TRUE)[[1]]
      if (length(parts) != 2) next
      a <- parts[1]; b <- parts[2]
      if (!(a %in% outcomes) || !(b %in% outcomes)) next
      rho <- tanh(as.numeric(psi[[nm]]))
      R[a,b] <- rho
      R[b,a] <- rho
    }
  }
  #assicura pd
  as.matrix(nearPD(R, corr=TRUE)$mat)
}

#5) matrice cov random effects: frailty + (y_I0,y_Ipost) per ogni outcome
build_Sigma_b <- function(psi, outcomes) {
  re_names <- c("frailty", as.vector(rbind(paste0(outcomes,"_I0"), paste0(outcomes,"_Ipost"))))
  q <- length(re_names)
  S <- matrix(0, q, q, dimnames=list(re_names, re_names))
  
  #var frailty
  if (!("var_frailty" %in% names(psi))) stop("missing var_frailty")
  S["frailty","frailty"] <- as.numeric(psi[["var_frailty"]])
  
  #var longitudinali
  for (y in outcomes) {
    v0 <- paste0("var_", y, "_I0")
    vp <- paste0("var_", y, "_Ipost")
    if (!(v0 %in% names(psi))) stop(paste0("missing ", v0))
    if (!(vp %in% names(psi))) stop(paste0("missing ", vp))
    S[paste0(y,"_I0"), paste0(y,"_I0")] <- as.numeric(psi[[v0]])
    S[paste0(y,"_Ipost"), paste0(y,"_Ipost")] <- as.numeric(psi[[vp]])
  }
  
  #cov: supporta sia cov_A__B che cov_frailty__X
  cov_idx <- grep("^cov_", names(psi))
  if (length(cov_idx) > 0) {
    for (nm in names(psi)[cov_idx]) {
      rest <- sub("^cov_", "", nm)
      
      if (startsWith(rest, "frailty__")) {
        x <- sub("^frailty__", "", rest)
        if (x %in% re_names) {
          S["frailty", x] <- as.numeric(psi[[nm]])
          S[x, "frailty"] <- as.numeric(psi[[nm]])
        }
        next
      }
      
      parts <- strsplit(rest, "__", fixed=TRUE)[[1]]
      if (length(parts) != 2) next
      a <- parts[1]; b <- parts[2]
      if (!(a %in% re_names) || !(b %in% re_names)) next
      S[a,b] <- as.numeric(psi[[nm]])
      S[b,a] <- as.numeric(psi[[nm]])
    }
  }
  
  #assicura simmetria/pd
  S <- (S + t(S))/2
  as.matrix(nearPD(S, corr=FALSE)$mat)
}

#6) loglik survival (weibull) condizionata su frailty u
loglik_surv_i <- function(T, delta, age0, u, psi) {
  beta_age <- as.numeric(psi[["fixed_surv_age"]])
  lambda <- exp(as.numeric(psi[["log_lambda"]]))
  rho    <- exp(as.numeric(psi[["log_rho"]]))
  
  eta <- beta_age * age0 + u
  H <- lambda * (T^rho) * exp(eta)
  if (delta==1) {
    logh <- log(lambda) + log(rho) + (rho-1)*log(T) + eta
    logh - H
  } else {
    -H
  }
}

#8) loglik long condizionata su b (frailty non entra qui direttamente)
loglik_long_i <- function(data_i, outcomes, thr, fixed, R, b_named, Rz=80, antithetic=TRUE) {
  
  if (nrow(data_i)==0) return(0)
  
  #design random: I0 e Ipost devono esistere
  if (!("I0" %in% names(data_i)) || !("Ipost" %in% names(data_i))) {
    stop("missing I0/Ipost in data")
  }
  
  ll <- 0
  for (r in seq_len(nrow(data_i))) {
    
    #eta per tutti gli outcome
    eta_row <- numeric(length(outcomes)); names(eta_row) <- outcomes
    
    for (y in outcomes) {
      beta <- fixed[[y]]
      if (length(beta) > 0) {
        x <- rep(0, length(beta)); names(x) <- names(beta)
        for (nm in names(beta)) x[nm] <- data_i[[nm]][r]
        eta_row[y] <- sum(beta * x)
      }
      
      #random part
      b0 <- b_named[paste0(y,"_I0")]
      bp <- b_named[paste0(y,"_Ipost")]
      eta_row[y] <- eta_row[y] + b0 * data_i$I0[r] + bp * data_i$Ipost[r]
    }
    
    #y osservati
    y_row <- sapply(outcomes, function(y) data_i[[y]][r])
    names(y_row) <- outcomes
    
    ll <- ll + row_logprob_loglik__genz(y_row, eta_row, thr, R)
    if (!is.finite(ll)) return(-Inf)
  }
  
  ll
}


logmeanexp_safe <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  #se tutto è NA o -Inf -> integral ~ 0 -> log = -Inf
  if (all((is.na(x)) | (is.infinite(x) & x < 0))) return(-Inf)
  
  x2 <- x[is.finite(x)]
  if (length(x2) == 0) return(-Inf)
  
  m <- max(x2)
  m + log(mean(exp(x2 - m)))
}

logmeanexp_safe <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  if (all((is.na(x)) | (is.infinite(x) & x < 0))) return(-Inf)
  x2 <- x[is.finite(x)]
  if (length(x2) == 0) return(-Inf)
  m <- max(x2)
  m + log(mean(exp(x2 - m)))
}

rows_with_obs <- function(df, cols) {
  if (nrow(df) == 0 || length(cols) == 0) return(rep(FALSE, nrow(df)))
  tmp <- df[, cols, drop = FALSE]
  Reduce(`|`, lapply(tmp, function(x) !is.na(x)))
}

extract_outcomes_from_psi <- function(psi) {
  thr_names <- grep("^threshold_", names(psi), value = TRUE)
  sort(unique(sub("^threshold_(.+)_cut:.*$", "\\1", thr_names)))
}

extract_need_cov_from_psi <- function(psi, outcomes) {
  fx_names <- grep("^fixed_", names(psi), value = TRUE)
  covs <- character(0)
  for (y in outcomes) {
    pat <- paste0("^fixed_", y, "_")
    covs <- c(covs, sub(pat, "", fx_names[grepl(pat, fx_names)]))
  }
  sort(unique(covs))
}

ensure_derived_covariates <- function(dat, need_cov) {
  # crea solo se servono davvero
  if ("wI6" %in% need_cov && !"wI6" %in% names(dat)) {
    if (!all(c("waiting", "I6") %in% names(dat))) stop("manca waiting o I6 per costruire wI6")
    dat$wI6 <- dat$waiting * dat$I6
  }
  if ("wI12" %in% need_cov && !"wI12" %in% names(dat)) {
    if (!all(c("waiting", "I12") %in% names(dat))) stop("manca waiting o I12 per costruire wI12")
    dat$wI12 <- dat$waiting * dat$I12
  }
  
  if ("pI6" %in% need_cov && !"pI6" %in% names(dat)) {
    if ("time_post" %in% names(dat) && "I6" %in% names(dat)) {
      dat$pI6 <- dat$time_post * dat$I6
    } else if (all(c("post", "I6") %in% names(dat))) {
      dat$pI6 <- dat$post * dat$I6
    } else stop("manca time_post (o post) o I6 per costruire pI6")
  }
  if ("pI12" %in% need_cov && !"pI12" %in% names(dat)) {
    if ("time_post" %in% names(dat) && "I12" %in% names(dat)) {
      dat$pI12 <- dat$time_post * dat$I12
    } else if (all(c("post", "I12") %in% names(dat))) {
      dat$pI12 <- dat$post * dat$I12
    } else stop("manca time_post (o post) o I12 per costruire pI12")
  }
  
  dat
}

.lemon_squeeze <- function(y, n) {
  (y * (n - 1) + 0.5) / n
}

library(mvtnorm)

mc_loglik_total_by_id <- function(
    dat,
    psi_adj_cov,
    id_var = "ID",
    M = 100,
    Rz = 80,
    antithetic = TRUE,
    include_R = TRUE, #opzione: usa o no la correlazione residua tra ordinali
    scale100_cont = TRUE,
    squeeze_cont  = TRUE,
    seed = 123,
    month_unit = 30.4375,
    admin_censor = "2025-01-01",
    cont_outcomes = NULL,
    cont_phi_prefix = "logphi_",
    cont_phi_global = "logphi",
    y_eps = 1e-12,
    mu_eps = 1e-12,
    phi_min = 1e-8,
    cont_transform = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  psi <- psi_adj_cov
  if (is.null(names(psi))) stop("psi_adj_cov deve avere names()")
  
  d <- dat
  d$ID_chr <- trimws(as.character(d[[id_var]]))
  
  #ord outcomes da threshold_
  thr_names <- grep("^threshold_", names(psi), value = TRUE)
  outcomes_ord <- if (length(thr_names) == 0) character(0) else
    sort(unique(sub("^threshold_(.+)_cut:.*$", "\\1", thr_names)))
  outcomes_ord <- outcomes_ord[outcomes_ord %in% names(d)]
  
  #cont outcomes: o passati, o da logphi_<name>
  if (is.null(cont_outcomes)) {
    nmphi <- grep(paste0("^", cont_phi_prefix), names(psi), value = TRUE)
    outcomes_cont <- sub(paste0("^", cont_phi_prefix), "", nmphi)
  } else {
    outcomes_cont <- cont_outcomes
  }
  outcomes_cont <- sort(unique(outcomes_cont))
  outcomes_cont <- outcomes_cont[outcomes_cont %in% names(d)]
  
  if (length(outcomes_ord) == 0 && length(outcomes_cont) == 0) {
    stop("nessun outcome: né threshold_ (ordinali) né cont_outcomes/logphi_ (beta) presenti anche in dat")
  }
  
  outcomes_all <- c(outcomes_ord, outcomes_cont)
  
  #assicura factor solo per ordinali
  for (o in outcomes_ord) {
    if (!is.factor(d[[o]])) {
      lv <- sort(unique(d[[o]][!is.na(d[[o]])]))
      d[[o]] <- factor(d[[o]], levels = lv)
    }
  }
  
  #soglie solo per ordinali
  build_thresholds <- function(psi, outcomes) {
    thr <- vector("list", length(outcomes)); names(thr) <- outcomes
    for (o in outcomes) {
      pat <- paste0("^threshold_", o, "_cut:")
      nm <- names(psi)[grepl(pat, names(psi))]
      if (length(nm) == 0) stop(paste0("missing thresholds for ", o))
      k <- as.integer(sub(".*_cut:(\\d+)$", "\\1", nm))
      ord <- order(k)
      thr[[o]] <- as.numeric(psi[nm[ord]])
    }
    thr
  }
  
  #fixed per tutti gli outcome (ord+cont), se presenti
  build_fixed <- function(psi, outcomes) {
    fixed <- vector("list", length(outcomes)); names(fixed) <- outcomes
    for (o in outcomes) {
      pat <- paste0("^fixed_", o, "_")
      nm <- names(psi)[grepl(pat, names(psi))]
      if (length(nm) == 0) fixed[[o]] <- numeric(0) else {
        covs <- sub(paste0("^fixed_", o, "_"), "", nm)
        b <- as.numeric(psi[nm]); names(b) <- covs
        fixed[[o]] <- b
      }
    }
    fixed
  }
  
  #corr residua solo per ordinali
  build_R <- function(psi, outcomes) {
    J <- length(outcomes)
    R <- diag(1, J); dimnames(R) <- list(outcomes, outcomes)
    nm <- names(psi)[grepl("^atanh_rho_resid__", names(psi))]
    if (length(nm) > 0) {
      for (s in nm) {
        rest <- sub("^atanh_rho_resid__", "", s)
        ab <- strsplit(rest, "__", fixed = TRUE)[[1]]
        if (length(ab) != 2) next
        a <- ab[1]; b <- ab[2]
        if (a %in% outcomes && b %in% outcomes) {
          rho <- tanh(as.numeric(psi[s]))
          R[a,b] <- rho; R[b,a] <- rho
        }
      }
    }
    R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
    R
  }
  
  #D random: (b0,bpost) per ogni outcome_all + frailty
  build_D <- function(psi, outcomes_all) {
    re <- c(as.vector(rbind(paste0(outcomes_all, "_I0"),
                            paste0(outcomes_all, "_Ipost"))),
            "frailty")
    q <- length(re)
    D <- matrix(0, q, q, dimnames = list(re, re))
    
    for (o in outcomes_all) {
      v0 <- paste0("var_", o, "_I0")
      vp <- paste0("var_", o, "_Ipost")
      if (!(v0 %in% names(psi))) stop(paste0("missing ", v0))
      if (!(vp %in% names(psi))) stop(paste0("missing ", vp))
      D[paste0(o,"_I0"),   paste0(o,"_I0")]   <- as.numeric(psi[[v0]])
      D[paste0(o,"_Ipost"),paste0(o,"_Ipost")]<- as.numeric(psi[[vp]])
    }
    
    if ("var_frailty" %in% names(psi)) {
      D["frailty","frailty"] <- as.numeric(psi[["var_frailty"]])
    } else stop("missing var_frailty")
    
    nmc <- names(psi)[grepl("^cov_", names(psi))]
    for (s in nmc) {
      rest <- sub("^cov_", "", s)
      
      if (startsWith(rest, "frailty__")) {
        x <- sub("^frailty__", "", rest)
        if (x %in% re) {
          cv <- as.numeric(psi[s])
          D["frailty", x] <- cv; D[x, "frailty"] <- cv
        }
        next
      }
      
      ab <- strsplit(rest, "__", fixed = TRUE)[[1]]
      if (length(ab) != 2) next
      a <- ab[1]; b <- ab[2]
      if (a %in% re && b %in% re) {
        cv <- as.numeric(psi[s])
        D[a,b] <- cv; D[b,a] <- cv
      }
    }
    
    D <- (D + t(D)) / 2
    D <- as.matrix(Matrix::nearPD(D, corr = FALSE)$mat)
    D
  }
  
  build_surv <- function(psi) {
    if (!("log_lambda" %in% names(psi))) stop("missing log_lambda")
    if (!("log_rho" %in% names(psi))) stop("missing log_rho")
    lam <- exp(as.numeric(psi[["log_lambda"]]))
    rho <- exp(as.numeric(psi[["log_rho"]]))
    nm <- names(psi)[grepl("^fixed_surv_", names(psi))]
    b <- numeric(0)
    if (length(nm) > 0) {
      covs <- sub("^fixed_surv_", "", nm)
      b <- as.numeric(psi[nm]); names(b) <- covs
    }
    list(lambda = lam, rho = rho, beta = b)
  }
  
  thr <- if (length(outcomes_ord) > 0) build_thresholds(psi, outcomes_ord) else list()
  fixed <- build_fixed(psi, outcomes_all)
  
  Rfull <- NULL
  if (isTRUE(include_R) && length(outcomes_ord) > 1) {
    Rfull <- build_R(psi, outcomes_ord)
  }
  
  D <- build_D(psi, outcomes_all)
  surv_par <- build_surv(psi)
  
  #cov richieste
  need_cov_long <- sort(unique(unlist(lapply(fixed, names))))
  need_cov_long <- need_cov_long[!is.na(need_cov_long) & need_cov_long != "(Intercept)"]
  need_cov_surv <- names(surv_par$beta)
  need_cov_surv <- need_cov_surv[!is.na(need_cov_surv) & need_cov_surv != "(Intercept)"]
  need_cov <- sort(unique(c(need_cov_long, need_cov_surv)))
  
  need_basic <- c(id_var, "date_baseline", "date_intervention", "date_revision", "I0", "Ipost")
  miss_basic <- setdiff(need_basic, names(d))
  if (length(miss_basic) > 0) stop(paste("missing columns in dat:", paste(miss_basic, collapse = ", ")))
  
  #derivate standard se servono
  if ("wI6" %in% need_cov && !"wI6" %in% names(d)) {
    if (!all(c("waiting","I6") %in% names(d))) stop("manca waiting o I6 per wI6")
    d$wI6 <- as.numeric(d$waiting) * as.numeric(d$I6)
  }
  if ("wI12" %in% need_cov && !"wI12" %in% names(d)) {
    if (!all(c("waiting","I12") %in% names(d))) stop("manca waiting o I12 per wI12")
    d$wI12 <- as.numeric(d$waiting) * as.numeric(d$I12)
  }
  
  tp <- if ("time_post" %in% names(d)) as.numeric(d$time_post) else NA_real_
  if (all(!is.finite(tp)) && ("t_within" %in% names(d))) tp <- as.numeric(d$t_within)
  if (all(!is.finite(tp)) && ("post" %in% names(d))) tp <- as.numeric(d$post)
  
  if ("pI6" %in% need_cov && !"pI6" %in% names(d)) {
    if (!("I6" %in% names(d))) stop("manca I6 per pI6")
    d$pI6 <- as.numeric(tp) * as.numeric(d$I6)
  }
  if ("pI12" %in% need_cov && !"pI12" %in% names(d)) {
    if (!("I12" %in% names(d))) stop("manca I12 per pI12")
    d$pI12 <- as.numeric(tp) * as.numeric(d$I12)
  }
  
  miss_cov <- setdiff(need_cov, names(d))
  if (length(miss_cov) > 0) stop(paste("missing covariates in dat:", paste(miss_cov, collapse = ", ")))
  
  #cont transform
  if (!is.null(cont_transform)) {
    for (o in intersect(names(cont_transform), outcomes_cont)) {
      d[[o]] <- cont_transform[[o]](d[[o]])
    }
  }
  
  
  # --- PREP outcome beta (bounded) come nel tuo script ---
  if (length(outcomes_cont) > 0) {
    for (o in outcomes_cont) {
      
      yraw <- as.numeric(d[[o]])
      
      # 1) scala in [0,1] se era 0-100
      y01 <- if (isTRUE(scale100_cont)) yraw / 100 else yraw
      
      # 2) clamp in [0,1]
      y01 <- pmin(pmax(y01, 0), 1)
      
      # 3) squeeze in (0,1)
      n_eff <- sum(is.finite(y01))
      if (n_eff <= 1) stop(paste0("troppi pochi valori finiti per outcome cont: ", o))
      
      if (isTRUE(squeeze_cont)) {
        y01 <- .lemon_squeeze(y01, n_eff)
      } else {
        y01 <- pmin(pmax(y01, y_eps), 1 - y_eps)
      }
      
      d[[o]] <- y01
    }
  }
  
  #date
  base <- as.POSIXct(d$date_baseline, tz = "UTC")
  interv <- as.POSIXct(d$date_intervention, tz = "UTC")
  rev <- as.POSIXct(d$date_revision, tz = "UTC")
  admin_cens <- as.POSIXct(admin_censor, tz = "UTC")
  
  #visit_date
  postv <- if ("post" %in% names(d)) as.numeric(d$post) else NA_real_
  visit_date <- rep(as.POSIXct(NA, tz="UTC"), nrow(d))
  for (i in seq_len(nrow(d))) {
    if (!is.na(base[i]) && (is.na(postv[i]) || postv[i] == 0)) {
      visit_date[i] <- base[i]
    } else if (!is.na(interv[i]) && is.finite(tp[i])) {
      visit_date[i] <- as.POSIXct(as.numeric(interv[i]) + tp[i] * month_unit * 86400,
                                  origin = "1970-01-01", tz = "UTC")
    } else if (!is.na(interv[i])) {
      visit_date[i] <- interv[i]
    } else {
      visit_date[i] <- base[i]
    }
  }
  d$visit_date <- visit_date
  
  ids_all <- unique(d$ID_chr)
  ids_all <- ids_all[!is.na(ids_all) & nzchar(ids_all)]
  n_id <- length(ids_all)
  
  #chol D
  D <- (D + t(D)) / 2
  D <- as.matrix(Matrix::nearPD(D, corr = FALSE)$mat)
  U_D <- chol(D)
  
  draw_u <- function(M, antithetic, U_D) {
    q <- ncol(U_D)
    M2 <- if (antithetic) ceiling(M / 2) else M
    Z <- matrix(rnorm(M2 * q), nrow = M2, ncol = q)
    U <- Z %*% U_D
    if (antithetic) {
      U <- rbind(U, -U)
      U <- U[seq_len(M), , drop = FALSE]
    }
    U
  }
  
  y_to_int0 <- function(x) {
    if (is.factor(x)) return(as.integer(x) - 1L)
    xi <- as.integer(x)
    if (!is.finite(xi)) return(xi)
    if (xi >= 1L) return(xi - 1L)
    xi
  }
  
  any_obs_row <- function(row, outs) {
    for (o in outs) if (!is.na(row[[o]])) return(TRUE)
    FALSE
  }
  
  cache_C <- new.env(parent = emptyenv())
  get_C_sub <- function(obs_names) {
    obs_names <- unname(obs_names)
    key <- paste(obs_names, collapse = ",")
    if (exists(key, envir = cache_C, inherits = FALSE)) {
      return(get(key, envir = cache_C, inherits = FALSE))
    }
    Rsub <- Rfull[obs_names, obs_names, drop = FALSE]
    Rsub <- (Rsub + t(Rsub)) / 2
    Rsub <- as.matrix(Matrix::nearPD(Rsub, corr = TRUE)$mat)
    L <- t(chol(Rsub))
    L[upper.tri(L)] <- 0
    Lvec <- as.numeric(L[lower.tri(L, diag = TRUE)])
    C <- mvtnorm::ltMatrices(Lvec, diag = TRUE)
    assign(key, C, envir = cache_C)
    C
  }
  
  #eta fixed per outcomes_all
  make_eta_fixed <- function(data_pre) {
    n <- nrow(data_pre)
    J <- length(outcomes_all)
    X_eta <- matrix(0, nrow = n, ncol = J, dimnames = list(NULL, outcomes_all))
    if (n == 0) return(X_eta)
    for (j in seq_along(outcomes_all)) {
      o <- outcomes_all[j]
      b <- fixed[[o]]
      if (length(b) == 0) next
      eta <- rep(0, n)
      if ("(Intercept)" %in% names(b)) eta <- eta + as.numeric(b["(Intercept)"])
      for (cv in setdiff(names(b), "(Intercept)")) {
        eta <- eta + as.numeric(b[cv]) * as.numeric(data_pre[[cv]])
      }
      X_eta[, j] <- eta
    }
    X_eta
  }
  
  #phi per cont outcome
  get_phi <- function(o) {
    nm <- paste0(cont_phi_prefix, o)
    if (nm %in% names(psi)) {
      phi <- exp(as.numeric(psi[[nm]]))
    } else if (cont_phi_global %in% names(psi)) {
      phi <- exp(as.numeric(psi[[cont_phi_global]]))
    } else {
      stop(paste0("missing ", nm, " e missing ", cont_phi_global))
    }
    if (!is.finite(phi) || phi < phi_min) phi <- phi_min
    phi
  }
  
  log_phi_diff <- function(b, a) {
    lb <- pnorm(b, log.p = TRUE)
    la <- pnorm(a, log.p = TRUE)
    if (is.infinite(la) && la < 0) return(lb)
    if (lb <= la) return(-Inf)
    lb + log1p(-exp(la - lb))
  }
  
  #loglik long: ordinal mvn + beta indipendente condizionata a u
  loglik_long_given_u <- function(data_pre, u, eta_fixed_mat, I0v, Ipostv) {
    if (nrow(data_pre) == 0) return(0)
    Jall <- length(outcomes_all)
    b0_all <- u[seq(1, 2 * Jall, by = 2)]
    bp_all <- u[seq(2, 2 * Jall, by = 2)]
    idx_all <- setNames(seq_len(Jall), outcomes_all)
    
    ll <- 0
    
    for (r in seq_len(nrow(data_pre))) {
      
      #ordinal part
      if (length(outcomes_ord) > 0) {
        obs <- which(!is.na(data_pre[r, outcomes_ord, drop = TRUE]))
        if (length(obs) > 0) {
          obs_names <- outcomes_ord[obs]
          jj <- idx_all[obs_names]
          mu <- eta_fixed_mat[r, obs_names] + I0v[r] * b0_all[jj] + Ipostv[r] * bp_all[jj]
          
          a <- numeric(length(obs_names))
          b <- numeric(length(obs_names))
          for (k in seq_along(obs_names)) {
            o <- obs_names[k]
            Klev <- nlevels(data_pre[[o]])
            cuts <- thr[[o]]
            if (length(cuts) != (Klev - 1L)) stop(paste0("soglie incoerenti per ", o))
            thrfull <- c(-Inf, cuts, Inf)
            yc <- y_to_int0(data_pre[[o]][r])
            lo <- thrfull[yc + 1L]
            hi <- thrfull[yc + 2L]
            a[k] <- lo - mu[k]
            b[k] <- hi - mu[k]
          }
          
          if (isTRUE(include_R) && length(obs_names) > 1) {
            C <- get_C_sub(obs_names)
            
            lp <- mvtnorm::lpmvnorm(
              lower  = matrix(a, ncol = 1),
              upper  = matrix(b, ncol = 1),
              chol   = C,
              M      = Rz,
              logLik = TRUE
            )
            lp <- as.numeric(lp)
            if (!is.finite(lp)) return(-Inf)
            ll <- ll + lp
          } else {
            lp <- 0
            for (k in seq_along(a)) {
              lp <- lp + log_phi_diff(b[k], a[k])
              if (!is.finite(lp)) return(-Inf)
            }
            ll <- ll + lp
          }
        }
      }
      
      #beta part
      if (length(outcomes_cont) > 0) {
        for (o in outcomes_cont) {
          y <- data_pre[[o]][r]
          if (is.na(y)) next
          
          j <- idx_all[[o]]
          et <- eta_fixed_mat[r, o] + I0v[r] * b0_all[j] + Ipostv[r] * bp_all[j]
          mu <- plogis(et)
          mu <- pmin(pmax(mu, mu_eps), 1 - mu_eps)
          
          y2 <- pmin(pmax(as.numeric(y), y_eps), 1 - y_eps)
          phi <- get_phi(o)
          a <- mu * phi
          bb <- (1 - mu) * phi
          lp <- dbeta(y2, a, bb, log = TRUE)
          if (!is.finite(lp)) return(-Inf)
          ll <- ll + lp
        }
      }
      
      if (!is.finite(ll)) return(-Inf)
    }
    
    ll
  }
  
  loglik_surv_given_u <- function(data_i, u, T, delta) {
    if (!is.finite(T) || T <= 0) return(-Inf)
    Jall <- length(outcomes_all)
    frail <- u[2 * Jall + 1L]
    
    lp <- frail
    b <- surv_par$beta
    if (length(b) > 0) {
      if ("(Intercept)" %in% names(b)) lp <- lp + as.numeric(b["(Intercept)"])
      for (cv in setdiff(names(b), "(Intercept)")) {
        x0 <- data_i[[cv]][1]
        if (is.na(x0)) x0 <- 0
        lp <- lp + b[cv] * as.numeric(x0)
      }
    }
    
    lam <- surv_par$lambda
    rho <- surv_par$rho
    lt <- log(T)
    logH <- log(lam) + rho * lt + lp
    H <- exp(logH)
    
    if (delta == 1L) {
      logh <- log(lam) + log(rho) + (rho - 1) * lt + lp
      logh - H
    } else {
      -H
    }
  }
  
  get_subj_surv <- function(id) {
    ii <- which(d$ID_chr == id)
    i0 <- interv[ii][which(!is.na(interv[ii]))[1]]
    r0 <- rev[ii][which(!is.na(rev[ii]))[1]]
    delta <- as.integer(!is.na(r0))
    censor_date <- if (delta == 1L) r0 else admin_cens
    T <- if (is.na(i0) || is.na(censor_date)) NA_real_
    else as.numeric(difftime(censor_date, i0, units = "days")) / month_unit
    list(delta = delta, T = T, censor_date = censor_date, rev0 = r0)
  }
  surv_map <- lapply(ids_all, get_subj_surv)
  names(surv_map) <- ids_all
  
  by_id <- vector("list", n_id)
  
  for (ii in seq_len(n_id)) {
    id <- ids_all[ii]
    idx <- which(d$ID_chr == id)
    data_i <- d[idx, , drop = FALSE]
    
    sm <- surv_map[[id]]
    T <- sm$T
    delta <- sm$delta
    censor_date <- sm$censor_date
    rev0 <- sm$rev0
    
    pre_mask <- rep(TRUE, nrow(data_i))
    if (!is.na(censor_date)) {
      pre_mask <- !is.na(data_i$visit_date) & (data_i$visit_date <= censor_date)
    }
    data_pre <- data_i[pre_mask, , drop = FALSE]
    
    n_post <- 0
    max_excess <- 0
    if (!is.na(rev0)) {
      post_mask <- !is.na(data_i$visit_date) & (data_i$visit_date > rev0)
      n_post <- sum(post_mask)
      if (n_post > 0) {
        max_excess <- max(as.numeric(difftime(data_i$visit_date[post_mask], rev0, units = "days")) / month_unit)
      }
    }
    
    n_long_obs_pre <- 0
    if (nrow(data_pre) > 0) {
      for (r in seq_len(nrow(data_pre))) {
        if (any_obs_row(data_pre[r, , drop = FALSE], outcomes_all)) n_long_obs_pre <- n_long_obs_pre + 1L
      }
    }
    
    X_eta <- make_eta_fixed(data_pre)
    I0v <- as.numeric(data_pre$I0)
    Ipostv <- as.numeric(data_pre$Ipost)
    
    U <- draw_u(M, antithetic, U_D)
    
    logw <- rep(NA_real_, M)
    for (m in seq_len(M)) {
      u <- U[m, ]
      ll_long <- if (n_long_obs_pre == 0) 0 else loglik_long_given_u(data_pre, u, X_eta, I0v, Ipostv)
      if (!is.finite(ll_long)) {
        logw[m] <- -Inf
      } else {
        ll_surv <- loglik_surv_given_u(data_i, u, T, delta)
        logw[m] <- ll_long + ll_surv
      }
    }
    
    n_finite <- sum(is.finite(logw))
    n_na <- sum(is.na(logw))
    n_nonfinite <- sum(!is.finite(logw) & !is.na(logw))
    
    if (n_finite == 0) {
      loglik_i <- NaN
      status <- "all_nonfinite"
    } else {
      mx <- max(logw[is.finite(logw)])
      loglik_i <- mx + log(mean(exp(logw[is.finite(logw)] - mx)))
      status <- if (is.finite(loglik_i)) "ok" else "nonfinite_logmean"
    }
    cat(sprintf("[%d/%d] id=%s  loglik_i=%.6f  status=%s\n",
                ii, n_id, id, loglik_i, status))
    
    by_id[[ii]] <- data.frame(
      ID_chr = id,
      delta = delta,
      T = T,
      n_rows_total = nrow(data_i),
      n_rows_pre = nrow(data_pre),
      n_long_obs_pre = n_long_obs_pre,
      n_post_rows = n_post,
      max_excess_m = max_excess,
      loglik_i = loglik_i,
      mc_n_finite = n_finite,
      mc_n_na = n_na,
      mc_n_nonfinite = n_nonfinite,
      status = status,
      stringsAsFactors = FALSE,
      row.names = NULL,
      check.names = FALSE
    )
  }
  
  by_id <- do.call(rbind, by_id)
  loglik_total <- if (all(is.finite(by_id$loglik_i))) sum(by_id$loglik_i) else NaN
  
  list(
    loglik = loglik_total,
    outcomes_ord = outcomes_ord,
    outcomes_cont = outcomes_cont,
    outcomes_all = outcomes_all,
    checks = list(
      n_id = n_id,
      n_event = sum(by_id$delta == 1L, na.rm = TRUE),
      n_post_ids = sum(by_id$n_post_rows > 0, na.rm = TRUE),
      n_T_bad = sum(!is.finite(by_id$T) | by_id$T <= 0, na.rm = TRUE)
    ),
    by_id = by_id
  )
}


