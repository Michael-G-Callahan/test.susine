standardize_X = function(X){
  #Center X
  X_col_means <- colMeans(X)
  X_centered <- sweep(X, 2, X_col_means, "-")

  #Scale X
  X_scaled <- apply(X_centered, 2, scale)

  return(X_scaled)
}

simulate_b = function(p, L, mu_0, sigma_0_2, prior_inclusion_weights, seed=1){
  set.seed(seed)
  beta = matrix(0, nrow = L, ncol = p)

  for (l in 1:L){
    b_l = rnorm(p, mean = mu_0[l], sd = sqrt(sigma_0_2[l]))
    gamma_l = rmultinom(1, 1, prob = prior_inclusion_weights)
    beta[l,] = b_l*gamma_l
  }

  return(colSums(beta))
}

simulate_y_noise = function(X, beta, noise, seed){
  set.seed(seed)
  y = X %*% beta #LABEL: MATRIX COMPUTATION
  sigma_2 = (1/(1-noise) - 1)*var(y)
  if (noise==1){
    y=0
    sigma_2=1
    }
  y = y + rnorm(nrow(X), mean = 0, sd = sqrt(sigma_2))
  return(y)
}

simulate_y = function(X, beta, sigma){
  set.seed(1)
  y = X %*% beta + rnorm(nrow(X), mean = 0, sd = sigma) #LABEL: MATRIX COMPUTATION
  return(y)
}

aggregate_pips = function(alpha){
  complement = 1 - alpha
  column_product = apply(complement, 2, prod)
  PIP = 1 - column_product

  return(PIP)
}

# For finding SusieR coef
compute_colSds = function(X) {
  if (is.matrix(X))
    return(colSds(X))
  else {
    n = nrow(X)
    Y = apply_nonzeros(X,function (u) u^2)
    d = colMeans(Y) - colMeans(X)^2
    return(sqrt(d*n/(n-1)))
  }
}

extract_susie_coef = function(alpha, mu, X){
  coef = colSums(alpha * mu) / compute_colSds(X)

  return(round(coef, 3))
}

get_coef_L2_error = function(fitted_coef, true_coef){
  err= fitted_coef - true_coef
  L2 = 1 - var(err) / var(true_coef)
  return(L2)
}

get_r2_adj = function(fitted_coef, y_hat, y){
  R_squared <- 1 - sum((y - y_hat)^2) / sum((y - mean(y))^2)
  n <- length(y)
  p_coef = sum(round(fitted_coef, 3) != 0)

  adjusted_R_squared <- 1 - ((1 - R_squared) * (n - 1)) / (n - p_coef - 1)
  return(adjusted_R_squared)
}

get_auc = function(alpha, beta){
  true_labels = as.integer(beta != 0)
  PIPs = aggregate_pips(alpha)

  # Create prediction objects
  predictions <- prediction(PIPs, true_labels)

  #Get AUC
  auc <- performance(predictions, "auc")@y.values[[1]]
  return(auc)
}

get_pr_auc = function(alpha, beta){
  true_labels = as.integer(beta != 0)
  PIPs = aggregate_pips(alpha)

  true_PIP = PIPs[true_labels == 1]
  false_PIP = PIPs[true_labels != 1]

  pr <- pr.curve(scores.class0 = true_PIP, scores.class1 = false_PIP)

  auc = pr$auc.davis.goadrich
  return(auc)
}

add_row_to_results_df = function(results_df, i, L, X, y, mu_0, prior_inclusion_weights, beta){

  #Run susine
  susine_output <- susine(
    L=L,
    X=X,
    y=y,
    mu_0=mu_0,
    prior_inclusion_weights = prior_inclusion_weights,
    prior_update_params = "var"
  )
  coef = colSums(susine_output$b_hat) / compute_colSds(X) + colMeans(X)
  fitted_y = X %*% coef
  results_df$L2_susine[i] = get_coef_L2_error(coef, beta)
  results_df$r2_susine[i] = get_r2_adj(coef, fitted_y, y)
  results_df$auc_susine[i] = get_auc(susine_output$alpha, beta)
  results_df$pr_auc_susine[i] = get_pr_auc(susine_output$alpha, beta)

  #Run susie
  susie_output <- susie(
    X=X,
    y=y,
    L=L,
    estimate_residual_variance = TRUE,
    estimate_prior_variance = TRUE
  )

  coef = extract_susie_coef(susie_output$alpha, susie_output$mu, X)
  fitted_y = X %*% coef + mean(y)
  results_df$L2_susie[i] = get_coef_L2_error(coef, beta)
  results_df$r2_susie[i] = get_r2_adj(coef, fitted_y, y)
  results_df$auc_susie[i] = get_auc(susie_output$alpha, beta)
  results_df$pr_auc_susie[i] = get_pr_auc(susie_output$alpha, beta)

  return(results_df)
}
