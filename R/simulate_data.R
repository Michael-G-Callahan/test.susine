# Simulation data generation helpers ----------------------------------------

#' Column-wise centering and scaling of a numeric matrix.
#'
#' This mirrors the behaviour of `scale()` but keeps the input type and avoids
#' copying when possible. Columns with zero variance are left unchanged.
#'
#' @param X Numeric matrix.
#' @param center Logical; subtract column means when TRUE.
#' @param scale Logical; divide by column standard deviations when TRUE.
#'
#' @return Matrix with standardized columns.
#' @keywords internal
standardize_x <- function(X, center = TRUE, scale = TRUE) {
  stopifnot(is.matrix(X))
  if (!center && !scale) {
    return(X)
  }
  cm <- if (center) matrixStats::colMeans2(X) else rep(0, ncol(X))
  if (scale) {
    csd <- matrixStats::colSds(X)
    csd[csd == 0] <- 1
  } else {
    csd <- rep(1, ncol(X))
  }
  sweep(sweep(X, 2, cm, `-`), 2, csd, `/`)
}

#' Simulate sparse effect sizes for p SNPs.
#'
#' @param p Number of SNPs.
#' @param p_star Number of causal SNPs to activate.
#' @param effect_sd Standard deviation of true effects.
#' @param seed Optional seed for reproducibility.
#'
#' @return List with `beta` vector and `causal_idx` integer indices.
#' @keywords internal
simulate_effect_sizes <- function(p,
                                  p_star,
                                  effect_sd = 1,
                                  seed = NULL) {
  stopifnot(p > 0, p_star >= 0)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  causal_idx <- if (p_star > 0) {
    sort(sample.int(p, size = min(p_star, p), replace = FALSE))
  } else integer(0)
  beta <- numeric(p)
  if (length(causal_idx)) {
    beta[causal_idx] <- stats::rnorm(length(causal_idx), mean = 0, sd = effect_sd)
  }
  list(beta = beta, causal_idx = causal_idx)
}

#' Construct mu_0 and sigma_0_2 priors from effect sizes.
#'
#' @param beta Numeric vector of true effects.
#' @param noise_causal Variance fraction for mu_0 - beta on causal SNPs.
#' @param noise_nonc Variance fraction for mu_0 - beta on non-causal SNPs.
#' @param base_sigma2 Optional baseline prior variance.
#'
#' @return List with `mu_0`, `sigma_0_2`, and `prior_inclusion_weights`.
#' @keywords internal
simulate_priors <- function(beta,
                            noise_causal,
                            noise_nonc,
                            base_sigma2 = NULL) {
  p <- length(beta)
  causal_idx <- which(beta != 0)
  noncausal_idx <- setdiff(seq_len(p), causal_idx)

  beta_causal <- beta[causal_idx]
  effect_var <- stats::var(beta_causal)
  if (is.na(effect_var) || effect_var == 0) {
    effect_var <- max(stats::var(beta), 1e-3)
  }
  if (is.na(effect_var) || effect_var == 0) {
    effect_var <- 1
  }

  sd_causal <- sqrt(max(noise_causal, 0)) * sqrt(effect_var)
  sd_noncausal <- sqrt(max(noise_nonc, 0)) * sqrt(effect_var)

  mu_0 <- numeric(p)
  if (length(causal_idx)) {
    mu_0[causal_idx] <- beta_causal + stats::rnorm(length(causal_idx), sd = sd_causal)
  }
  if (length(noncausal_idx)) {
    mu_0[noncausal_idx] <- stats::rnorm(length(noncausal_idx), mean = 0, sd = sd_noncausal)
  }

  base_sigma2 <- base_sigma2 %||% effect_var
  if (is.na(base_sigma2) || base_sigma2 <= 0) {
    base_sigma2 <- effect_var
  }
  sigma_0_2 <- rep(base_sigma2, p)
  if (length(causal_idx)) {
    sigma_0_2[causal_idx] <- pmax(base_sigma2 * (1 - noise_causal), base_sigma2 * 0.25)
  }
  if (length(noncausal_idx)) {
    sigma_0_2[noncausal_idx] <- base_sigma2 * (1 + noise_nonc)
  }

  scores <- abs(mu_0)
  scores <- scores + 1e-6
  prior_weights <- scores / sum(scores)

  list(
    mu_0 = mu_0,
    sigma_0_2 = sigma_0_2,
    prior_inclusion_weights = prior_weights
  )
}

#' Simulate phenotype with a target noise fraction.
#'
#' @param X Design matrix.
#' @param beta Effect sizes.
#' @param noise_fraction Fraction of variance attributed to noise (between 0 and 1).
#' @param seed Optional seed.
#'
#' @return List with `y` and implied residual variance `sigma2`.
#' @keywords internal
simulate_phenotype <- function(X, beta, noise_fraction, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  signal <- as.vector(X %*% beta)
  signal_var <- stats::var(signal)
  if (is.na(signal_var) || signal_var == 0) {
    signal_var <- max(stats::var(beta), 1e-3)
  }
  noise_fraction <- max(min(noise_fraction, 0.999), 0)
  if (noise_fraction == 1) {
    sigma2 <- 1
    y <- stats::rnorm(nrow(X), mean = 0, sd = 1)
  } else {
    sigma2 <- signal_var * (noise_fraction / (1 - noise_fraction))
    y <- signal + stats::rnorm(nrow(X), mean = 0, sd = sqrt(sigma2))
  }
  list(y = as.numeric(y), sigma2 = sigma2)
}

#' Generate a single simulation dataset for a run specification.
#'
#' @param spec Named list or data.frame row containing simulation controls.
#'   Required fields: `seed`, `p_star`, `y_noise`, `prior_noise_causal`,
#'   `prior_noise_nonc`. Optional: `effect_sd`, `standardize_X`.
#' @param base_X Optional matrix to reuse instead of loading the package data.
#'
#' @return List with design matrix, phenotype, ground-truth beta, priors, and
#'   metadata required by the modelling pipeline.
#' @export
generate_simulation_data <- function(spec,
                                     base_X = NULL) {
  if (is.null(base_X)) {
    data_env <- new.env(parent = emptyenv())
    utils::data("SuSiE_N3_X", package = "test_susine", envir = data_env)
    base_X <- get("SuSiE_N3_X", envir = data_env)
  }
  X <- as.matrix(base_X)
  if (resolve_flag(spec$standardize_X, FALSE)) {
    X <- standardize_x(X)
  }

  seed <- spec$seed %||% 1L
  effects <- simulate_effect_sizes(
    p = ncol(X),
    p_star = spec$p_star %||% 5L,
    effect_sd = spec$effect_sd %||% 1,
    seed = seed
  )
  phenotype <- simulate_phenotype(
    X = X,
    beta = effects$beta,
    noise_fraction = spec$y_noise %||% 0.5,
    seed = seed + 1L
  )
  priors <- simulate_priors(
    beta = effects$beta,
    noise_causal = spec$prior_noise_causal %||% 0,
    noise_nonc = spec$prior_noise_nonc %||% 1,
    base_sigma2 = stats::var(phenotype$y)
  )

  list(
    X = X,
    y = phenotype$y,
    beta = effects$beta,
    sigma2 = phenotype$sigma2,
    mu_0 = priors$mu_0,
    sigma_0_2 = priors$sigma_0_2,
    prior_inclusion_weights = priors$prior_inclusion_weights,
    causal_idx = effects$causal_idx
  )
}
