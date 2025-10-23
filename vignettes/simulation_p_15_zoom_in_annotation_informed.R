# ROAR
# .libPaths("/storage/home/mgc5166/R/x86_64-pc-linux-gnu-library/4.3")

############# Helper functions #######################################
add_susine_diagnostics = function(results_df, X, i, susine_output, label, beta){
  #Store data
  coef = colSums(susine_output$b_hat) / compute_colSds(X) + colMeans(X)
  fitted_y = X %*% coef

  results_df[[paste0("L2_", label)]][i] = get_coef_L2_error(coef, beta)
  results_df[[paste0("R2_", label)]][i] = get_r2_adj(coef, fitted_y, y)
  results_df[[paste0("auc_", label)]][i] = get_auc(susine_output$alpha, beta)
  results_df[[paste0("pr_auc_", label)]][i] = get_pr_auc(susine_output$alpha, beta)
  results_df[[paste0("sigma_2_", label)]][i] = get_sigma_2_from_output(susine_output)
  results_df[[paste0("mu_0_", label)]][i] = paste0(round(susine_output$mu_0,3), collapse="_")
  results_df[[paste0("sigma_0_2_", label)]][i] = paste0(signif(susine_output$sigma_0_2,3), collapse="_")

  return(results_df)
}

get_sigma_0_2_inf = function(L, sigma_star_1, sigma_star_2, star_1_idx, y){
  sigma_0_2_inf = c(rep(sigma_star_1^2, star_1_idx),rep(sigma_star_2^2, L - star_1_idx))
  sigma_0_2_inf = sigma_0_2_inf / rep(var(y), L)
  sigma_0_2_inf = matrix(sigma_0_2_inf,nrow = L, ncol = 1)
}

get_sigma_2_from_output = function(susine_output){
  output_sigma_2 = susine_output$sigma_2[!is.na(susine_output$sigma_2)]
  output_sigma_2 = output_sigma_2[length(output_sigma_2)]
  output_sigma_2 = round(output_sigma_2,3)
  return(output_sigma_2)
}

add_row_to_vary_p_df = function(results_df, i, L, X, y, oracle_inclusion_weights, beta){
  p = dim(X)[2]
  mu_0_inf = rep(0,p)

  susine_fit_var = susine(L, X, y,
                          prior_update_method = "var",
                          mu_0 = mu_0_inf,
                          max_iter = 1000)
  susine_fit_mean = susine(L, X, y,
                           prior_update_method = "mean",
                           mu_0 = mu_0_inf,
                           max_iter = 1000)
  susine_fit_both = susine(L, X, y,
                           prior_update_method = "both",
                           mu_0 = mu_0_inf,
                           max_iter = 1000)
  susine_fit_none = susine(L, X, y,
                           prior_update_method = "none",
                           mu_0 = mu_0_inf,
                           max_iter = 1000)
  susine_fit_oracle = susine(L, X, y,
                             prior_update_method = "none",
                             prior_inclusion_weights = oracle_inclusion_weights,
                             mu_0 = mu_0_inf,
                             max_iter = 1000)

  #Evaluation
  results_df = add_susine_diagnostics(results_df, X, i, susine_fit_var, "var", beta)
  results_df = add_susine_diagnostics(results_df, X, i, susine_fit_mean, "mean", beta)
  results_df = add_susine_diagnostics(results_df, X, i, susine_fit_both, "both", beta)
  results_df = add_susine_diagnostics(results_df, X, i, susine_fit_none, "naive", beta)
  results_df = add_susine_diagnostics(results_df, X, i, susine_fit_oracle, "full-oracle", beta)

  #Run susine random
  susine_output <- susine(
    L=L,
    X=X[sample(1:nrow(X)), sample(1:ncol(X))],
    y=y,
    prior_update_params = "none"
  )
  results_df = add_susine_diagnostics(results_df, X, i, susine_output, "random", beta)

  return(results_df)
}

#Shuffle the columns of the data
get_X_sample = function(seed){
  X <- X[sample(n), sample(p_data)]
  return(X)
}

get_v = function(p_annots,p){
  v = rep(0, p)
  v_snp_select = sample(1:p, p_annots, replace=FALSE)
  v[v_snp_select] = sample(annots, p_annots, replace = FALSE)
  return(v)
}

get_beta = function(v,p,p_star,beta_noise){
  beta = rep(0, p)
  star_select = sample(1:p, p_star, replace=FALSE)
  sigma_star = sqrt(var(v)/(1/beta_noise -1))

  effects = rnorm(p, v, sigma_star)
  beta[star_select] = effects[star_select]
  return(beta)
}

############# Dependencies #######################################

library(devtools)
# library(susieR)
library(matrixStats)
library(ROCR)
library(PRROC)
library(pushoverr)
library(VGAM) #laplace

# load_all()
set_pushover_user(user = "ukxcksez6rmyh4g9pnez13iq6nuvfg")
set_pushover_app(token = "ahycq5kiyytvd4fae25sz862xsbat6")


#Annotations data (annots)
#Annotations data (annots)
#folder_path = "/storage/work/mgc5166/susine/data/" #ROAR
#subfolder = "Annotations" #ROAR

folder_path = "~/School/Research/Genetics/cTWAS Generalization/Code/Other data/SydhK562Irf1Ifng6h" #Edit
subfolder = "1.sannot" #Edit

file_name = "1.sannot"
df_file_path <- file.path(folder_path, subfolder,file_name)

data <- read.table(df_file_path, header = TRUE, sep = "", stringsAsFactors = FALSE)
annots <- data[,6]
rm(data)
####### TESTING

# summary(lm(beta[beta!=0] ~ v[beta!=0] ))
# which(beta!=0)
# beta[beta!=0]
# v[beta!=0]
#
# inform_susine_output <- susine(
#   L=10,
#   X=X,
#   y=y,
#   mu_0=mu_0_inf,
#   prior_update_params = "var"
# )
# inform_susine_output$sigma_0_2[,1:2]
# pips = aggregate_pips(inform_susine_output$alpha)
# sum(pips[v!=0])
#
# sum_annotated_beta_nonzero <- sum((pips)[(v != 0) & (beta != 0)])
# sum_nonannotated_beta_nonzero <- sum((pips)[(v == 0) & (beta != 0)])
# sum_annotated_beta_zero <- sum((pips)[(v != 0) & (beta == 0)])
# sum_nonannotated_beta_zero <- sum((pips)[(v == 0) & (beta == 0)])
#
# # Create the 2x2 table
# table_2x2 <- matrix(
#   c(
#     sum_annotated_beta_nonzero,
#     sum_nonannotated_beta_nonzero,
#     sum_annotated_beta_zero,
#     sum_nonannotated_beta_zero
#   ),
#   nrow = 2,
#   byrow = TRUE,
#   dimnames = list(
#     c("Beta != 0", "Beta == 0"),
#     c("Annotated SNPs (v != 0)", "Non-Annotated SNPs (v == 0)")
#   )
# )
#
# # Print the table
# print(table_2x2)
#
# which(pips > 0.1)
# pips[pips > 0.1]
# v[pips > 0.1]
# beta[pips > 0.1]
# inform_susine_output$alpha[,797]
# which(inform_susine_output$alpha[1,]>0.1)
# inform_susine_output$alpha[1,inform_susine_output$alpha[1,]>0.1]
#
# susine_output <- susine(
#   L=10,
#   X=X,
#   y=y,
#   mu_0=0,
#   prior_update_params = "var"
# )
# susine_output$sigma_0_2[,1:2]
# pips = aggregate_pips(susine_output$alpha)
#
#
# sum((pips)[(v!=0) & (beta!=0)]) #selection err on annotated SNPs
# sum((pips)[(v==0)& (beta!=0)]) #selection err on non-annotated SNPs
#
# sum((pips)[(v!=0) & (beta==0)]) #selection err on annotated SNPs
# sum((pips)[(v==0)& (beta==0)]) #selection err on non-annotated SNPs
#
# sum(pips[v!=0])
# which(pips > 0.1)
# pips[pips > 0.1]
# v[pips > 0.1]
# beta[pips > 0.1]
# susine_output$alpha[2,797]
# susine_output$alpha[1,c(180, 669)]
# which(susine_output$alpha[1,]>0.1)
# susine_output$alpha[1,susine_output$alpha[1,]>0.1]


############# Main - setup ####################################################

#Experimental loop vars
#Number of effects
L_min = 5
L_max = 10 #20
L_seq = seq(L_min, L_max, by = 5)

#Effect groups
p_star = 10
p_annots = seq(100, 500, by = 200)
beta_noise = seq(0.1, 0.9, by=0.2)

#Initialization
seed_max = 25
seeds <- seq(1, seed_max, by = 1)

#Experimental steady vars
noise = seq(0.3, 0.7, by = 0.4)

#Generate experimental parameter DF
loop_df <- expand.grid(L_seq = L_seq, beta_noise = beta_noise)

# Create a data frame with all combinations of L, p_star -- one row = one task
experiment_df <- expand.grid(noise = noise, seed = seeds, p_annots=p_annots)
n_exp = nrow(experiment_df)

experiment_ID = 36 #36 is teh trouble maker
# args = commandArgs(trailingOnly=TRUE) #TASK ID -- ROAR
# experiment_ID = as.numeric(args[1]) #TASK ID -- ROAR

p_annots <-experiment_df$p_annots[experiment_ID]
noise <-experiment_df$noise[experiment_ID]
seed <- experiment_df$seed[experiment_ID]

set.seed(seed)

data("SuSiE_N3_X") #This comes from the SuSiE pkg under 'data/'N3finemapping.rData'

#Get n, p from X data
X = as.matrix(SuSiE_N3_X)
X = standardize_X(X)
rm(SuSiE_N3_X)
p = ncol(X)
n =  nrow(X)

############# Main - loop #######################################################

metric_list = c("L2_", "R2_", "auc_", "pr_auc_", "sigma_2_", "mu_0_", "sigma_0_2_")
model_list = c("var", "both", "naive", "inform", "inform-EB","full-oracle",
               "random")

metric_model_list <- expand.grid(metric_list, model_list)
metric_model_list <- apply(metric_model_list, 1, paste0, collapse = "")


results_df = expand.grid(noise = noise, p = p, L = L_seq, seed = seed,
                        p_star = p_star, p_annots = p_annots,
                        beta_noise = beta_noise)

results_df[metric_model_list] = NA

########
i=0
for (k in 1:nrow(loop_df)){
  i=i+1
  beta_noise <-loop_df$beta_noise[k]
  L <-loop_df$L[k]

  v = get_v(p_annots,p)
  beta = get_beta(v,p,p_star,beta_noise)
  y = simulate_y_noise(X, beta, noise, seed+1)

  #Model settings -- not dependent on L
  active_betas = 1*(beta!=0)
  oracle_inclusion_weights = matrix(active_betas / sum(active_betas),nrow = 1,ncol = p)
  mu_0_inf = v
  mu_0_oracle = beta

  results_df = add_row_to_vary_p_df(results_df, i, L, X, y, oracle_inclusion_weights, beta)
}

############## Output ######################################################

output_folder_path = "/mgc5166/work/output..." #Edit
output_folder_path = "~/School/Research/Genetics/cTWAS Generalization/Code/susine/data" #Edit

# main_output_folder_path = "/storage/work/mgc5166/simulation_R_output" #Edit -- ROAR
# subDir = paste0('Job_ID_',args[2]) #Edit -- ROAR
# dir.create(file.path(main_output_folder_path, subDir), showWarnings = FALSE) #Edit -- ROAR
# output_folder_path = paste0(main_output_folder_path, '/',subDir) #Edit -- ROAR

file_name = paste0("sim_annotations_exp_id_",experiment_ID, ".Rdata")
df_output_file_path <- file.path(output_folder_path, file_name)
save(results_df, file=df_output_file_path)

if (experiment_ID == nrow(experiment_df)){
  pushover(message = "Finished running!")
}

print("Completed.")

