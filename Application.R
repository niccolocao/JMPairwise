library(Matrix)
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)
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


# Wald on fisher z for correlations Cor(b_ik1, s_i) and Cor(b_ik2, s_i)
corr_frail <- est_interest_adj[grepl("frailty", names(est_interest_adj))]
se_frail   <- se_interest[names(corr_frail)]

tab_fisher <- wald_fisher_from_r(corr_frail, se_frail)
tab_fisher


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

