library(Matrix)
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(tikzDevice)
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




# Complete table of results
tab_stime <- make_table_est_se_items_rows(est_full_nocov, se_full_nocov, digits = 3)
ù

# Wald tests
tab <- wald_table_raw(est_interest_adj, se_interest, add_fdr_corr = TRUE)
tab_use = tab

# Wald on fisher z for correlations Cor(b_ik1, s_i) and Cor(b_ik2, s_i)
corr_frail <- est_interest_adj[grepl("frailty", names(est_interest_adj))]
se_frail   <- se_interest[names(corr_frail)]

tab_fisher <- wald_fisher_from_r(corr_frail, se_frail)



# Correlations Cor(b_ik2 - b_ik1, s_i)
tab_cd2 <- corrdiff_frailty_table_from_interest(
  est_interest_adj = est_interest_adj,
  V_interest = V_interest,
  sort_by = "p_fdr"
)
tab_cd2


###### Test on covariances Cor(b_ik2 - b_ik1, s_i) = 0 simultaneously
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

