library(Matrix)
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(tikzDevice)
library(knitr)
library(kableExtra)
library(ggrepel)
library(grid)
source("pairaveFuns.R")



K_map <- c(
  setNames(rep(3L, 5L), c('Mobility','SelfCare','UsualActivities','PainorDiscomfort','Anxietyordepression')),
  setNames(rep(5L, 7L), c('Gettingoutofbed','Puttingonsocks','Gettingupfromsitting',
                          'Bendingtowardthefloor_pickinguponobjectfromtheground',
                          'Twisting_pivotingontheinjuredknee','Kneeling','Squatting'))
)

get_K_from_itemname <- function(item) {
  item <- as.character(item)
  if (!item %in% names(K_map)) {
    stop("Item non trovato in K_map: ", item)
  }
  as.integer(K_map[[item]])
}

files <- list.files("pairwise_9_bivar/", pattern=".rds", full.names=TRUE)


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


# Extract bivariate fit estimates
out_list <- lapply(files, process_one, K_map = K_map)
names(out_list) <- basename(files)
ex_list <- lapply(out_list, `[[`, "ex")         
pair_tags <- names(out_list)                 

# Estimated J and K matrices, with Sigma(\hat theta) for estimates on original scale 
# Cholesky + log on diagonal for D and log for weibull and beta
ex_pref_for_jk <- Map(prefix_ex, ex_list, pair_tags)
jk_all <- JK_all_pairs(ex_pref_for_jk, hessian_of = "neglogLik")
V_theta <- jk_all$vcov_theta   

# transformation of D from cholesky + log on diagonal to var covar first with delta method on the elements os Sigma(\hat theta)
pre <- predelta_prepare(ex_list, pair_tags, jk_vcov_theta = V_theta)
psi_all <- pre$psi_all
V_psi_all <- pre$V_psi

# Averaging over var covar for D and the other parameters on original scale
A_all <- build_A(psi_all, canonicalizer = canonicalizer_all)
psi_global <- as.numeric(A_all %*% psi_all)
names(psi_global) <- rownames(A_all)

V_global <- A_all %*% V_psi_all %*% t(A_all)
V_global <- 0.5 * (V_global + t(V_global)) #symmetry for numeric stability

se_global <- sqrt(diag(V_global))


# Trasfomation of the other parameters on scale of interest (e.g., sd and corr) + delta method
est_interest <- derive_interest(psi_global)

Jg <- jacobian(func = g, x = psi_global)
colnames(Jg) <- names(psi_global)
rownames(Jg) <- names(derive_interest(psi_global))

V_interest <- Jg %*% V_global %*% t(Jg)
V_interest <- 0.5 * (V_interest + t(V_interest)) #symmetry for numeric stability

se_interest <- sqrt(diag(V_interest))


# D var covar and sd corr
outD = re_varcov_matrix2(
  psi_global,
  label_filter = "(_I0|_Ipost|frailty)$"       # tieni solo questi
)



### correction  for non SPD
eigenvalues <- eigen(outD$D)$values
eigenvectors = eigen(outD$D)$vectors
eigenvalues[eigenvalues<0]<-1e-12

Dadj = eigenvectors%*%diag(eigenvalues)%*%t(eigenvectors)
rownames(Dadj) = rownames(outD$D)
colnames(Dadj) = colnames(outD$D)


itemnames <- c('Mobility','SelfCare','UsualActivities','PainorDiscomfort','Anxietyordepression',
               'Gettingoutofbed','Puttingonsocks','Gettingupfromsitting',
               'Bendingtowardthefloor_pickinguponobjectfromtheground',
               'Twisting_pivotingontheinjuredknee','Kneeling','Squatting',"Perceivedcurrenthealthstatus"
               )
ordD = c(paste0(itemnames,"_I0"),paste0(itemnames,"_Ipost"),'frailty')
Dadj <- Dadj[ordD,ordD]

Corradj <- cov2cor(Dadj)


# update est_interest for sd_ e corr_
est_interest_adj <- update_est_interest_from_D(est_interest, Dadj, Corradj,representation = "corr")




##### only correlations
res_full <- build_interest_full(
  psi_global = psi_global,
  V_global   = V_global,
  est_trans  = est_interest_adj,
  J_trans    = Jg,
  keep_untransformed = names(psi_global),
  drop_untransformed_regex = "^(var_|cov_)",   # <<<<<< toglie var/cov dal blocco "as-is"
  out_order = NULL
)

est_full_nocov <- res_full$est
V_full_nocov   <- res_full$V
se_full_nocov  <- res_full$se



# Vector of estimates for dynamic prediction with  var cov and untrasformed estiamtes, except for SPD correction to D
psi_adj_cov <- update_est_interest_from_D(
  psi_global, Dadj, Corradj,
  representation = "cov"
)





######### Table for all fixed effects (longitudinal submodel ############################################

item_order <- setdiff(.default_base_order, "frailty")

tab_items3 <- make_wald_table_3cols(
  psi = psi_global,
  se  = se_global,
  alpha = 0.05,
  predictors = c("I6","I12","wI6","wI12","pI6","pI12","diag_oa","sex01","loghosp_post","age"),
  item_levels = item_order
)


# esempio mappa (opzionale) per rinominare le righe:
item_map <- c(
  Mobility = "Mobility",
  SelfCare = "SelfCare",
  UsualActivities = "UsualAct",
  PainorDiscomfort = "PainDisc",
  Anxietyordepression = "AnxDepr",   
  Gettingoutofbed = "OutBed",
  Puttingonsocks = "PutSocks",
  Gettingupfromsitting = "UpSit", 
  Bendingtowardthefloor_pickinguponobjectfromtheground="Bend",
  Twisting_pivotingontheinjuredknee = "Twist",
  Kneeling = "Kneel",
  Squatting = "Squat"
)


preds <- c("I6","I12","wI6","wI12","pI6","pI12","diag_oa","sex01","loghosp_post","age")

latex_tab <- make_latex_beta_table(
  tab_items3,
  predictors = preds,
  item_name_map = item_map,   # oppure NULL
  alpha = 0.05,
  digits_est = 3,
  digits_se  = 3,
  caption = "Fixed effects by item (Wald test; $^{*}$ p<0.05).",
  label   = "tab:fixed_items",
  use_resizebox = TRUE
)

cat(latex_tab)





#################################### Wald tests ########################################################
tab <- wald_table_raw(est_interest_adj, se_interest, add_fdr_corr = TRUE)
tab_use = tab

# Wald on fisher z for correlations Cor(b_ik1, s_i) and Cor(b_ik2, s_i)
corr_frail <- est_interest_adj[grepl("frailty", names(est_interest_adj))]
se_frail   <- se_interest[names(corr_frail)]

tab_fisher <- wald_fisher_from_r(corr_frail, se_frail)



################# Simultaneous Wald test on covariances  Cov(b_ik1, s_i) = Cov(b_ik2, s_i) = 0 for different subgroups of items
EQ_items <- c("Mobility","SelfCare","UsualActivities","PainorDiscomfort","Anxietyordepression", "Perceivedcurrenthealthstatus")
KOOS_items <- c("Gettingoutofbed","Gettingupfromsitting",
                "Bendingtowardthefloor_pickinguponobjectfromtheground",
                "Twisting_pivotingontheinjuredknee","Kneeling","Squatting","Puttingonsocks")

res_wald <- wald_frailty_cov_suite(psi_global, se_global = se_global, V_global = V_global,
                              EQ_items = EQ_items, KOOS_items = KOOS_items)


res_wald$tests = res_wald$tests[,-2]
kable(res_wald$tests,
      format = "latex", booktabs = TRUE,
      caption = "Wald tests congiunti: H0 corr(item, frailty)=0",
      digits = c(0,3,0,3)) |>
  kable_styling(latex_options = c("hold_position"))




# Correlations Cor(b_ik2 - b_ik1, s_i)
tab_cd2 <- corrdiff_frailty_table_from_interest(
  est_interest_adj = est_interest_adj,
  V_interest = V_interest,
  sort_by = "p_fdr"
)
tab_cd2


###### Simultaneous Wald test on covariances Cov(b_ik2 - b_ik1, s_i) = 0 
EQ_items <- c(
  "Mobility","SelfCare","UsualActivities",
  "PainorDiscomfort","Anxietyordepression",
  "Perceivedcurrenthealthstatus"
)

KOOS_items <- c(
  "Gettingoutofbed","Gettingupfromsitting",
  "Bendingtowardthefloor_pickinguponobjectfromtheground",
  "Twisting_pivotingontheinjuredknee","Kneeling",
  "Squatting","Puttingonsocks"
)

all_items <- c(EQ_items, KOOS_items)

res_shiftcov_all <- wald_shiftcov_frailty(
  psi_global = psi_global,
  V_global = V_global,
  se_global = se_global,
  items = all_items
)

res_shiftcov_EQ <- wald_shiftcov_frailty(
  psi_global = psi_global,
  V_global = V_global,
  se_global = se_global,
  items = EQ_items
)

res_shiftcov_KOOS <- wald_shiftcov_frailty(
  psi_global = psi_global,
  V_global = V_global,
  se_global = se_global,
  items = KOOS_items
)

res_shiftcov_all$W
res_shiftcov_all$df
res_shiftcov_all$p

res_shiftcov_EQ$W
res_shiftcov_EQ$df
res_shiftcov_EQ$p

res_shiftcov_KOOS$W
res_shiftcov_KOOS$df
res_shiftcov_KOOS$p





######################################################## Correlation plots ######################################################################

itemnamesEQ   <- c('Mobility','SelfCare','UsualAct','PainDisc','AnxDepr',"Healt")
itemnamesKOOS <- c('OutBed','PutSocks','UpSit','Bend','Twist','Kneel','Squat')

origEQ   <- c('Mobility','SelfCare','UsualActivities','PainorDiscomfort','Anxietyordepression',"Perceivedcurrenthealthstatus")
origKOOS <- c('Gettingoutofbed','Puttingonsocks','Gettingupfromsitting',
              'Bendingtowardthefloor_pickinguponobjectfromtheground',
              'Twisting_pivotingontheinjuredknee','Kneeling','Squatting')

short_map_base <- setNames(c(itemnamesEQ, itemnamesKOOS), c(origEQ, origKOOS))

add_abc <- TRUE 
labs_I0 <- make_labs_time(est_interest, time = "I0", include_frailty = TRUE)
sd_I0   <- make_sd_vec_for_labs(est_interest, labs_I0)
RP_I0   <- make_RP_from_tab(tab, labs_I0, p_col = c("p_fisher","p"))


I0frail <- plot_corr_circles_lower(
  R = RP_I0$R,
  P = RP_I0$P,
  sd_vec = sd_I0,
  title = panel_title("a","RE at baseline + Frailty"),
  label_fun = make_pretty_label_vec,
  block_split = 6,
  circle_max = 6,
  r_text_size = 3.5,
  sd_text_size = 3.5,
  block_text_size = 4.8
) + theme_corr_text(
  axis_x_size = 12,
  axis_y_size = 12,
  plot_title_size = 15,
  legend_title_size = 10,
  legend_text_size = 9
)

# --- post + frailty
labs_Ipost <- make_labs_time(est_interest, time = "Ipost", include_frailty = TRUE)
sd_Ipost   <- make_sd_vec_for_labs(est_interest, labs_Ipost)
RP_Ipost   <- make_RP_from_tab(tab, labs_Ipost, p_col = c("p_fisher","p"))


Ipostfrail <- plot_corr_circles_lower(
  R = RP_Ipost$R,
  P = RP_Ipost$P,
  sd_vec = sd_Ipost,
  title = panel_title("b","RE at post-intervention + Frailty"),
  label_fun = make_pretty_label_vec,
  block_split = 6,
  circle_max = 6,
  r_text_size = 3.5,
  sd_text_size = 3.5,
  block_text_size = 4.8
) + theme_corr_text(
  axis_x_size = 12,
  axis_y_size = 12,
  plot_title_size = 15,
  legend_title_size = 10,
  legend_text_size = 9
)

I0frail
Ipostfrail



rows_I0   <- make_labs_time(est_interest, time = "I0",   include_frailty = FALSE)
cols_Ipost<- make_labs_time(est_interest, time = "Ipost",include_frailty = FALSE)

cross <- make_cross_RP_from_tab(
  tab = tab,
  rows = rows_I0,
  cols = cols_Ipost,
  p_col = c("p_fisher","p"),
  set_oob_to_NA = TRUE
)

p_cross <- plot_corr_circles_rect(
  R = cross$R,
  P = cross$P,
  title = panel_title("c","Baseline × Post-intervention"),
  x_title = "Post-intervention",
  y_title = "Baseline",
  y_axis_position = "right",
  row_split = 6,
  col_split = 6,
  block_line_width = 1.8,
  circle_max = 6,
  oob_cells = cross$oob_cells,
  r_text_size = 3.5,
  show_y_labels = FALSE,
  right_title_size = 5.6,
  right_title_x = 0.98
) + theme_corr_text(
  axis_x_size = 12,
  axis_y_size = 12,
  axis_title_y_size = 12,
  axis_title_x_size = 14,
  plot_title_size = 15,
  legend_title_size = 10,
  legend_text_size = 9
)



p_all <- (I0frail + guides(size = "none", colour = "none") |
            Ipostfrail + guides(size = "none", colour = "none") |
            p_cross) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


 

no_y_left <- theme(
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)

no_y_right <- theme(
  axis.text.y.right  = element_blank(),
  axis.ticks.y.right = element_blank()
)

p_all <- (
  I0frail + guides(size = "none", colour = "none") |
    (Ipostfrail + guides(size = "none", colour = "none") + no_y_left) |
    p_cross
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

tikz("/home/nc/MEGA/BioPsychometrics/biomlatex/corrplots.tex", width = 14, height = 6.5, standAlone = FALSE,engine = "xetex")
print(p_all)
dev.off()


labs_resid <- c(origEQ, origKOOS)

R_resid <- corr_resid_matrix_simple(
  est_interest_adj,
  base = "corr_resid",
  sep = "__"
)

labs_resid <- labs_resid[labs_resid %in% rownames(R_resid)]
R_resid <- R_resid[labs_resid, labs_resid, drop = FALSE]
R_resid <- 0.5 * (R_resid + t(R_resid))

P_resid <- make_P_resid_from_tab(
  tab = tab_use,
  labs = labs_resid,
  base = "corr_resid",
  sep = "__",
  p_col = c("p_fdr", "p_fisher", "p")
)

p_resid <- plot_corr_circles_lower(
  R = R_resid,
  P = P_resid,
  sd_vec = setNames(rep(1, length(labs_resid)), labs_resid),
  title = "Residual correlations",
  label_fun = make_pretty_label_vec,
  block_split = 6,
  circle_max = 6,
  r_text_size = 3.5,
  sd_text_size = 3.5,
  
)+ theme_corr_text(
  axis_x_size = 12,
  axis_y_size = 12,
  plot_title_size = 15,
  legend_title_size = 10,
  legend_text_size = 9
)

p_resid


tikz("/home/nc/MEGA/BioPsychometrics/biomlatex/residcor.tex", width = 6, height = 6, standAlone = FALSE,engine = "xetex")
print(p_resid)
dev.off()


######################## Principal component analysis on correlation matrices ##############################################
itemnamesEQ   <- c("Mobility","SelfCare","UsualAct","PainDisc","AnxDepr","Health")
itemnamesKOOS <- c('OutBed','PutSocks','UpSit','Bend','Twist','Kneel','Squat')

origEQ_full <- c(
  "Mobility",
  "SelfCare",
  "UsualActivities",
  "PainorDiscomfort",
  "Anxietyordepression",
  "Perceivedcurrenthealthstatus"
)

origKOOS_full <- c(
  "Gettingoutofbed",
  "Puttingonsocks",
  "Gettingupfromsitting",
  "Bendingtowardthefloor_pickinguponobjectfromtheground",
  "Twisting_pivotingontheinjuredknee",
  "Kneeling",
  "Squatting"
)

short_map_full <- c(
  setNames(itemnamesEQ, origEQ_full),
  setNames(itemnamesKOOS, origKOOS_full)
)

cols_main <- c(
  "EQ-5D_I0"    = "#0F766E",
  "EQ-5D_Ipost" = "#d1a00d",
  "KOOS_I0"     = "#5B21B6",
  "KOOS_Ipost"  = "#52340b",
  "Frailty"     = "#C2410C"
)

cols_resid <- c(
  "EQ"   = "grey40",
  "KOOS" = "#111111"
)

R_I0 <- Corradj[grepl("I0", rownames(Corradj)), grepl("I0", colnames(Corradj)), drop = FALSE]
R_Ipost <- Corradj[grepl("Ipost", rownames(Corradj)), grepl("Ipost", colnames(Corradj)), drop = FALSE]

R_resid <- corr_resid_matrix_simple(
  est_interest_adj,
  base = "corr_resid",
  sep = "__"
)

ord_resid <- c(origEQ_full, origKOOS_full)
ord_resid <- ord_resid[ord_resid %in% rownames(R_resid)]
R_resid <- R_resid[ord_resid, ord_resid, drop = FALSE]
R_resid <- 0.5 * (R_resid + t(R_resid))
R_resid <- as.matrix(nearPD(R_resid, corr = TRUE)$mat)
rownames(R_resid) <- colnames(R_resid) <- ord_resid


xy_I0 <- get_pca_xy(R_I0)
xy_Ipost <- get_pca_xy(R_Ipost)
xy_resid <- get_pca_xy(R_resid)

all_x <- c(0, xy_I0$x, xy_Ipost$x, xy_resid$x)
all_y <- c(0, xy_I0$y, xy_Ipost$y, xy_resid$y)

ref_all <- max(abs(c(all_x, all_y)))
xlim_common <- range(all_x) + c(-0.06, 0.30) * ref_all
ylim_common <- range(all_y) + c(-0.10, 0.10) * ref_all



R_full <- Corradj[!grepl("frail", rownames(Corradj)), !grepl("frail", colnames(Corradj)), drop = FALSE]

xy_full  <- get_pca_xy(R_full)
xy_I0    <- get_pca_xy(R_I0)
xy_Ipost <- get_pca_xy(R_Ipost)
xy_resid <- get_pca_xy(R_resid)

all_x <- c(0, xy_full$x, xy_I0$x, xy_Ipost$x, xy_resid$x)
all_y <- c(0, xy_full$y, xy_I0$y, xy_Ipost$y, xy_resid$y)

ref_all <- max(abs(c(all_x, all_y)))

xlim_common_all <- range(all_x) + c(-0.06, 0.30) * ref_all
ylim_common_all <- range(all_y) + c(-0.10, 0.10) * ref_all


p_full <- pca_plot_grouped(
  R_full,
  main = "(a)",
  short_map = short_map_full,
  eq_base = origEQ_full,
  koos_base = origKOOS_full,
  mode = "timed",
  show_suffix = FALSE,
  xlim_fixed = xlim_common_all,
  ylim_fixed = ylim_common_all,
  label_size = 6,
  axis_title_size = 10,
  axis_text_size = 11,
  plot_title_size = 19,
  escape_tex = TRUE,
  legend_text_size = 12
)

p_pca_I0 <- pca_plot_grouped(
  R_I0,
  main = "(b)",
  short_map = short_map_full,
  eq_base = origEQ_full,
  koos_base = origKOOS_full,
  mode = "timed",
  show_suffix = FALSE,
  xlim_fixed = xlim_common_all,
  ylim_fixed = ylim_common_all,
  label_size = 6,
  axis_title_size = 9,
  axis_text_size = 12,
  plot_title_size = 19,
  escape_tex = TRUE,
  legend_text_size = 14
)

p_pca_Ipost <- pca_plot_grouped(
  R_Ipost,
  main = "(c)",
  short_map = short_map_full,
  eq_base = origEQ_full,
  koos_base = origKOOS_full,
  mode = "timed",
  show_suffix = FALSE,
  xlim_fixed = xlim_common_all,
  ylim_fixed = ylim_common_all,
  label_size = 6,
  axis_title_size = 9,
  axis_text_size = 12,
  plot_title_size = 19,
  escape_tex = TRUE,
  legend_text_size = 14
)

p_pca_resid <- pca_plot_grouped(
  R_resid,
  main = "(d)",
  short_map = short_map_full,
  eq_base = origEQ_full,
  koos_base = origKOOS_full,
  mode = "resid",
  show_suffix = FALSE,
  xlim_fixed = xlim_common_all,
  ylim_fixed = ylim_common_all,
  label_size = 6,
  axis_title_size = 9,
  axis_text_size = 12,
  plot_title_size = 19,
  escape_tex = TRUE,
  legend_text_size = 14
)


no_y_numbers <- theme(
  axis.text.y = element_blank(),
  axis.ticks.y = element_line(color = "black"),
  axis.ticks.length.y = grid::unit(2.2, "mm")
)

p_all <- wrap_plots(
  p_full,
  p_pca_I0 + no_y_numbers,
  p_pca_Ipost + no_y_numbers,
  p_pca_resid + no_y_numbers,
  ncol = 2,
  byrow = TRUE
) &
  theme(
    axis.title = element_text(size = 19),
    legend.text = element_text(size = 18)
  )

out_tex <- "/home/nc/MEGA/BioPsychometrics/biomlatex/pca.tex"

tikz(out_tex, width = 20, height = 20, standAlone = FALSE,engine = "xetex")
print(p_all)

dev.off()






############################ Dynamic Predictions ##########################################
source("dynpred.R")
source("dynplots.R")

dat2 = readRDS("dat2_with_hosp30g.rds")
dat_pred <- dat2
dat_pred$visit_m <- c(0, 6, 12)[dat_pred$visit_idx + 1L]



itemnames <- c('Mobility','SelfCare','UsualActivities','PainorDiscomfort','Anxietyordepression',
               'Gettingoutofbed','Puttingonsocks','Gettingupfromsitting','Bendingtowardthefloor_pickinguponobjectfromtheground','Twisting_pivotingontheinjuredknee','Kneeling','Squatting')

q_map <- list(
  QoL = c("Mobility","SelfCare","UsualActivities","PainorDiscomfort","Anxietyordepression"),
  Mob = c("Gettingoutofbed","Puttingonsocks","Gettingupfromsitting","Bendingtowardthefloor_pickinguponobjectfromtheground","Twisting_pivotingontheinjuredknee","Kneeling","Squatting"))

item_labels <- c(Mobility = "Mobility", SelfCare = "SelfCare", UsualActivities = "UsualAct", PainorDiscomfort = "PainDisc", Anxietyordepression = "AnxDepr",
  Gettingoutofbed = "OutofBed", Puttingonsocks = "PutSocks", Gettingupfromsitting = "UpSit", Bendingtowardthefloor_pickinguponobjectfromtheground = "Bend",
  Twisting_pivotingontheinjuredknee = "Twist", Kneeling = "Kneel", Squatting = "Squat")



prep <- prep_long_all_items_surv_mixed2( dat = dat, ord_itemnames = itemnames, cont_itemnames = "Perceivedcurrenthealthstatus", cont_transforms = cont_tf, surv_outcome = "failure",  admin_date = "2025-01-01")

long_dat <- prep$long
K_by_item  <- prep$K_by_item

outcome_info <- data.frame(item = c(itemnames,"Perceivedcurrenthealthstatus"), type = c(rep("ordinal", length(itemnames)),"continuous"), stringsAsFactors = FALSE)


itemnames_ord <- c( "Mobility","SelfCare","UsualActivities","PainorDiscomfort","Anxietyordepression",
  "Gettingoutofbed","Puttingonsocks","Gettingupfromsitting", "Bendingtowardthefloor_pickinguponobjectfromtheground",
  "Twisting_pivotingontheinjuredknee","Kneeling","Squatting")


theta_sampler_global <- make_theta_sampler_global(theta_mean = psi_adj_cov, V_theta = V_global)



est_model <- strip_D_summaries_from_est(est_full_nocov)

fit_fullplug_global <- make_fullplug_object(est = est_model,D = Dadj, K_map = K_map, ord_items = itemnames_ord, cont_items = "Perceivedcurrenthealthstatus",cont_scale = c(Perceivedcurrenthealthstatus = 100),
  cont_transforms = cont_tf,
  theta_sampler = theta_sampler_global,
  covariate_transforms = list(
    age = function(x) (as.numeric(x) - prep$age_center) / prep$age_scale
  )
)



### new (fake) subjects

fake1_id <- "FAKE1"
fake2_id <- "FAKE2"
fake3_id <- "FAKE3"

times_m <- c(0, 6, 12)

# longitudinal profiles
Ymat1 <- rbind(
  c(1,1,1,1,1, 1,1,1,1,1,1,1),
  c(2,1,2,1,2, 2,3,3,2,2,1,1),
  c(2,2,2,2,1, 1,2,4,2,3,1,1)
)
colnames(Ymat1) <- itemnames

Ymat2 <- rbind(
  c(2,1,2,1,1, 2,2,2,1,1,2,1),
  c(2,3,2,2,2, 3,4,4,3,2,3,2),
  c(2,3,3,2,2, 4,4,5,4,3,4,2)
)
colnames(Ymat2) <- itemnames

Ymat3 <- rbind(
  c(1,2,1,2,2, 2,3,2,3,1,2,1),
  c(3,2,3,1,3, 3,3,4,3,2,2,3),
  c(3,3,3,2,2, 5,4,4,5,3,5,5)
)
colnames(Ymat3) <- itemnames


Ycont1 <- c(30, 40, 45)
Ycont2 <- c(40, 60, 70)
Ycont3 <- c(70, 80, 85)



# covariates
dat0_fake1 <- data.frame(
  ID = fake1_id,
  age = mean(long_dat$age, na.rm = TRUE),
  sex01 = 0L,
  diag_oa = 0L,
  loghosp = mean(long_dat$loghosp, na.rm = TRUE),
  prost_uni = 0L,
  waiting = mean(long_dat$waiting, na.rm = TRUE)
)

dat0_fake2 <- data.frame(
  ID = fake2_id,
  age = mean(long_dat$age, na.rm = TRUE),
  sex01 = 0L,
  diag_oa = 0L,
  loghosp = mean(long_dat$loghosp, na.rm = TRUE),
  prost_uni = 0L,
  waiting = mean(long_dat$waiting, na.rm = TRUE)
)

dat0_fake3 <- data.frame(
  ID = fake3_id,
  age = mean(long_dat$age, na.rm = TRUE),
  sex01 = 0L,
  diag_oa = 0L,
  loghosp = mean(long_dat$loghosp, na.rm = TRUE),
  prost_uni = 0L,
  waiting = mean(long_dat$waiting, na.rm = TRUE)
)

# history construction
hist_fake1_wide <- make_fake_hist_wide(
  id = fake1_id,
  ord_items = itemnames,
  Y_ord = Ymat1,
  covariates_row = dat0_fake1,
  cont_name = "Perceivedcurrenthealthstatus",
  Y_cont = Ycont1,
  times = times_m,
  times_in = "months"
)

hist_fake2_wide <- make_fake_hist_wide(
  id = fake2_id,
  ord_items = itemnames,
  Y_ord = Ymat2,
  covariates_row = dat0_fake2,
  cont_name = "Perceivedcurrenthealthstatus",
  Y_cont = Ycont2,
  times = times_m,
  times_in = "months"
)

hist_fake3_wide <- make_fake_hist_wide(
  id = fake3_id,
  ord_items = itemnames,
  Y_ord = Ymat3,
  covariates_row = dat0_fake3,
  cont_name = "Perceivedcurrenthealthstatus",
  Y_cont = Ycont3,
  times = times_m,
  times_in = "months"
)


newdata_fake <- dplyr::bind_rows(
  hist_fake1_wide,
  hist_fake2_wide,
  hist_fake3_wide
)

######################## parallel computing for predictions ######################################
library(parallel)
library(doParallel)
library(foreach)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

n_cores <- 9

tasks <- data.frame(
  label = c("dm1_0","dm1_6","dm1_12",
            "dm2_0","dm2_6","dm2_12",
            "dm3_0","dm3_6","dm3_12"),
  obj_name = c("hist_fake1_wide","hist_fake1_wide","hist_fake1_wide",
               "hist_fake2_wide","hist_fake2_wide","hist_fake2_wide",
               "hist_fake3_wide","hist_fake3_wide","hist_fake3_wide"),
  last_time = c(0,6,12, 0,6,12, 0,6,12),
  seed = c(28963,28963,28963,
           28543,28543,28543,
           98043,98043,98043),
  outfile = c("pred_fake1_0.rds","pred_fake1_6.rds","pred_fake1_12.rds",
              "pred_fake2_0.rds","pred_fake2_6.rds","pred_fake2_12.rds",
              "pred_fake3_0_hist.rds","pred_fake3_6_hist.rds","pred_fake3_12_hist.rds"),
  stringsAsFactors = FALSE
)

cl <- makeCluster(n_cores, type = "PSOCK", outfile = "")

clusterEvalQ(cl, {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    BLIS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
  setwd("~/dynpred")
  library(Matrix)
  library(MASS)
  library(mvtnorm)
  library(doParallel)
  library(foreach)
  source("dynpred.R")
  NULL
})

clusterExport(cl, c(
  "tasks",
  "fit_fullplug_global",
  "hist_fake1_wide",
  "hist_fake2_wide",
  "hist_fake3_wide"
))

registerDoParallel(cl)

pred_summary <- foreach(i = 1:nrow(tasks),
                        .combine = rbind,
                        .errorhandling = "pass") %dopar% {
  tt <- tasks[i, , drop = FALSE]
  
  hist_obj <- get(tt$obj_name)
  
  cat(tt$label, "\n")
  
  res <- tryCatch({
    pred <- survfitJM_fullplug(
      object = fit_fullplug_global,
      newdata = hist_obj,
      idVar = "ID",
      timeVar = "visit_m",
      survTimes = seq(0, 12 * 15, by = 1),
      last.time = tt$last_time,
      simulate = TRUE,
      silent_fake_repair = TRUE,
      repair_fake = FALSE,
      theta_source = "global",
      M = 1000,
      scale = 1.6,
      seed = tt$seed
    )
    
    saveRDS(pred, tt$outfile)
    
    data.frame(
      label = tt$label,
      outfile = tt$outfile,
      success_rate = mean(pred$success.rate, na.rm = TRUE),
      ok = TRUE,
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      label = tt$label,
      outfile = tt$outfile,
      success_rate = NA_real_,
      ok = FALSE,
      error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
  
  res
}

stopCluster(cl)
##################################################################################


pred_fake1_0 = readRDS("pred_fake1_0.rds")
pred_fake2_0 = readRDS("pred_fake2_0.rds")
pred_fake3_0 = readRDS("pred_fake3_0_hist.rds")

pred_fake1_6 = readRDS("pred_fake1_6.rds")
pred_fake2_6 = readRDS("pred_fake2_6.rds")
pred_fake3_6 = readRDS("pred_fake3_6_hist.rds")

pred_fake1_12 = readRDS("pred_fake1_12.rds")
pred_fake2_12 = readRDS("pred_fake2_12.rds")
pred_fake3_12 = readRDS("pred_fake3_12_hist.rds")



histories <- list(
  FAKE1 = hist_fake1_wide,
  FAKE2 = hist_fake2_wide,
  FAKE3 = hist_fake3_wide
)

outs <- list(
  FAKE1 = as_plot_out_from_fullplug(pred_fake1_12, id = "FAKE1"),
  FAKE2 = as_plot_out_from_fullplug(pred_fake2_12, id = "FAKE2"),
  FAKE3 = as_plot_out_from_fullplug(pred_fake3_12, id = "FAKE3")
)

titles <- c(
  FAKE1 = "(a) poor profile",
  FAKE2 = "(b) average profile",
  FAKE3 = "(c) good profile"
)


t0_days <- 0
t0_x <- 0

q_map$Cont <- "Perceivedcurrenthealthstatus"

item_labels <- c(
  item_labels,
  Perceivedcurrenthealthstatus = "Health"
)



outs_33 <- list(
  FAKE1 = make_out_triplet_survfit(
    pred0 = pred_fake1_0,
    pred6 = pred_fake1_6,
    pred12 = pred_fake1_12,
    id = "FAKE1",
    t_followup = c(0, 6, 12),
    relative = "auto"
  ),
  FAKE2 = make_out_triplet_survfit(
    pred0 = pred_fake2_0,
    pred6 = pred_fake2_6,
    pred12 = pred_fake2_12,
    id = "FAKE2",
    t_followup = c(0, 6, 12),
    relative = "auto"
  ),
  FAKE3 = make_out_triplet_survfit(
    pred0 = pred_fake3_0,
    pred6 = pred_fake3_6,
    pred12 = pred_fake3_12,
    id = "FAKE3",
    t_followup = c(0, 6, 12),
    relative = "auto"
  )
)

histories_33 <- histories[c("FAKE1","FAKE2", "FAKE3")]


p_33 <- plot_3_subjects_3x3_compact(
  outs = outs_33,
  histories = histories_33,
  ids = c("FAKE1","FAKE2", "FAKE3"),
  q_map = q_map,
  labels_map = item_labels,
  row_titles = c("poor profile","average profile", "good profile"),
  cond_labels = c("baseline", "6 months", "12 months"),
  time_in = "months",
  out_time_in = "months",
  ci = "pointwise",
  unit_out = "months",
  horizon = 12 * 5,
  pt_dodge = 0.6,
  y_jitter = 0.05,
  pt_size = 4.5,
  pt_stroke = 1.1,
  line_size = 1.2,
  item_linetype = 1,
  surv_size = 1.2,
  surv_linetype = 1,
  surv_ci_size = 0.9,
  surv_ci_linetype = 2,
  sep_size = 0.8,
  sep_colour = "grey75",
  sep_linetype = 1,
  sep_extend_to = "cond",
  sep_extend = 0,
  cont_sep_pad = 0.15,
  cond_vline_size = 0.8,
  cond_vline_linetype = 3,
  cond_vline_colour = "grey35",
  cex = 2.3,
  left_axis_title = "Item support",
  cont_limits = c(0, 100),
  cont_band = 0.85,
  title_size = 14,
  axis_title_size = 13,
  axis_text_size = 12,
  axis_line_size = 0.9,
  axis_tick_size = 0.9,
  panel_border_size = 0.9,
  axis_ticks_length_pt = 7,
  x_right_pad = 0,
  x_left_pad = 0.05,
  legend_height = 0.10,
  legend_title_size = 28,
  legend_text_size = 25,
  legend_key_pt = 25,
  
  top_row_height = 1.20,
  middle_row_height = 1.20,
  bottom_row_height = 1.38,
  legend_key_width_pt = 4,
  legend_spacing_x_pt = 5,
  legend_text_margin_left_pt = 5,
  legend_text_margin_right_pt = 15,
  
  
  question_col_width = 0.1,
  panel_widths = c(1, 1, 1),
  plot_margin_left_pt = 0.8,
  plot_margin_right_title_pt = 1.0,
  plot_margin_right_numbers_pt = 0.35,
  plot_margin_right_none_pt = 0.1,
  legend_nrow = 2,
  legend_byrow = TRUE
)

tikzDevice::tikz("/home/nc/MEGA/BioPsychometrics/biomlatex/pred3x3.tex", width=25, height=18)
print(p_33) 
dev.off()



