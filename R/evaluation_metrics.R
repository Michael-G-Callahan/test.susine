# evaluation_metrics.R
# Comprehensive evaluation for SuSiE/SuSiNE-like model fits.
# Produces BOTH purity-FILTERED (primary, as in SuSiE) and UNFILTERED metrics.
#
# Inputs:
#   fit : model object with
#         - fit$effect_fits$alpha  (L x p matrix of per-effect PIPs α_lj)
#         - fit$settings$L         (integer)
#         - optional: fit$model_fit$fitted_y, fit$model_fit$sigma_2, fit$model_fit$coef
#   X   : n x p matrix (columns align with α / variables)
#   y   : length-n response vector (only needed for h_g^2)
#   causal_idx : integer indices of true causal variables (for simulation-based metrics)
#
# Key definitions (SuSiE-style):
#   - ρ-level CS: sort α decreasing, take smallest prefix with cumulative sum ≥ ρ.
#   - Purity: min absolute correlation among all pairs inside a CS.
#   - Coverage (per-effect): 1 if CS contains ≥ 1 true causal, else 0.
#   - Power (per-model): fraction of true causals captured by ≥ 1 CS.
#   - Combined PIP per variable: 1 - prod_l (1 - α_lj).
#
# References:
#   - Credible sets & purity filtering: SuSiE §3.2.2. :contentReference[oaicite:1]{index=1}
#   - Power/FDR vs PIP thresholds (PR-curve relation). :contentReference[oaicite:2]{index=2}

# ---------- helpers ----------

# ρ-level credible set by cumulative PIP (SuSiE definition)
get_credible_set <- function(alpha, rho = 0.95) {
  if (!is.numeric(alpha)) stop("alpha must be numeric PIPs")
  if (!is.finite(rho) || rho <= 0 || rho >= 1) stop("rho must be in (0,1)")

  a <- alpha
  a[!is.finite(a)] <- 0
  if (sum(a) <= 0) return(integer(0))  # nothing to select

  # order by decreasing PIP; break ties by original index for determinism
  o  <- order(-a, seq_along(a))
  cs <- cumsum(a[o])
  k  <- which(cs >= rho)[1]
  if (is.na(k)) integer(0) else o[seq_len(k)]
}

# Purity = min absolute pairwise correlation inside a CS (size-1 -> 1, empty -> NA)
cs_purity_min_abs <- function(X, idx) {
  k <- length(idx)
  if (k == 0) return(NA_real_)
  if (k == 1) return(1.0)
  C <- suppressWarnings(cor(X[, idx, drop = FALSE], use = "pairwise.complete.obs"))
  min(abs(C[upper.tri(C)]))
}

# Combine per-effect PIPs to single per-variable PIP: 1 - prod_l (1 - α_lj)
combine_pips <- function(alpha_mat) {
  one_minus <- pmax(0, 1 - alpha_mat)
  as.numeric(1 - apply(one_minus, 2, prod))
}

# Overlap rate among a list of index vectors (fraction of CS pairs that overlap)
overlap_rate_from_sets <- function(list_of_sets) {
  m <- length(list_of_sets)
  if (m < 2) return(NA_real_)
  pairs <- utils::combn(m, 2)
  overlaps <- apply(pairs, 2, function(col) {
    i <- list_of_sets[[col[1]]]
    j <- list_of_sets[[col[2]]]
    length(intersect(i, j)) > 0
  })
  mean(overlaps)
}

# Per-effect metrics given one alpha vector
cs_metrics_one_effect <- function(alpha_vec, X, causal_idx, rho = 0.95) {
  cs <- get_credible_set(alpha_vec, rho)
  size <- length(cs)
  purity <- cs_purity_min_abs(X, cs)
  coverage <- as.integer(size > 0 && any(cs %in% causal_idx))
  list(indices = cs, size = size, purity = purity, coverage = coverage)
}

# ---------- classification-style metrics on combined PIPs ----------

# Average Precision (AUPRC) for binary labels using scores in [0,1]
# Implemented as mean precision at each positive's rank.
auprc_average_precision <- function(scores, labels) {
  if (!is.numeric(scores) || !is.numeric(labels))
    stop("scores and labels must be numeric")
  if (length(scores) != length(labels))
    stop("lengths differ")
  labels <- as.integer(labels > 0)
  P <- sum(labels)
  if (P == 0L) return(NA_real_)            # undefined if no positives
  if (P == length(labels)) return(1.0)     # all positive -> AP = 1

  o <- order(scores, decreasing = TRUE)
  y <- labels[o]
  tp <- cumsum(y)
  fp <- cumsum(1L - y)
  pos_idx <- which(y == 1L)
  mean(tp[pos_idx] / (tp[pos_idx] + fp[pos_idx]))
}

# Cross-entropy (log loss) for binary labels and predicted probs
cross_entropy_loss <- function(scores, labels, eps = 1e-12) {
  if (length(scores) != length(labels))
    stop("lengths differ")
  p <- pmin(pmax(scores, eps), 1 - eps)
  y <- as.integer(labels > 0)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

# SNP-heritability / variance explained estimate
# Prefers fitted_y; else uses residual variance; else uses coef.
estimate_hg2 <- function(fit, y, X) {
  if (missing(y) || is.null(y)) return(NA_real_)
  vy <- stats::var(y)
  if (!is.finite(vy) || vy <= 0) return(NA_real_)

  # 1) fitted_y available?
  fy <- tryCatch(fit$model_fit$fitted_y, error = function(e) NULL)
  if (!is.null(fy)) {
    h2 <- stats::var(as.numeric(fy)) / vy
    return(max(0, min(1, h2)))
  }

  # 2) residual variance path available?
  s2 <- tryCatch(fit$model_fit$sigma_2, error = function(e) NULL)
  if (!is.null(s2) && length(s2)) {
    s2_last <- tail(s2[is.finite(s2)], 1)
    if (length(s2_last)) {
      h2 <- 1 - s2_last / vy
      return(max(0, min(1, h2)))
    }
  }

  # 3) posterior mean coefficients available?
  b <- tryCatch(fit$model_fit$coef, error = function(e) NULL)
  if (!is.null(b)) {
    mu <- as.numeric(X %*% b)
    h2 <- stats::var(mu) / vy
    return(max(0, min(1, h2)))
  }

  NA_real_
}

# ---------- main evaluator ----------

evaluate_model <- function(fit, X, y = NULL, causal_idx = integer(0),
                           rho = 0.95, purity_threshold = 0.5,
                           compute_curves = TRUE) {
  stopifnot(is.matrix(X))
  L <- fit$settings$L
  A <- fit$effect_fits$alpha   # L x p

  # Per-effect CS & metrics (UNFILTERED)
  eff <- vector("list", L)
  for (l in seq_len(L)) {
    eff[[l]] <- cs_metrics_one_effect(A[l, ], X, causal_idx, rho = rho)
  }

  effects_unfiltered <- data.frame(
    effect   = seq_len(L),
    size     = sapply(eff, `[[`, "size"),
    purity   = sapply(eff, `[[`, "purity"),
    coverage = sapply(eff, `[[`, "coverage"),
    stringsAsFactors = FALSE
  )
  effects_unfiltered$indices <- I(lapply(eff, `[[`, "indices"))

  # FILTERED view (keep CS with purity >= threshold)
  keep_idx <- which(!is.na(effects_unfiltered$purity) &
                      effects_unfiltered$purity >= purity_threshold)
  effects_filtered <- effects_unfiltered[keep_idx, , drop = FALSE]

  # Unions for power
  union_unfiltered <- sort(unique(unlist(effects_unfiltered$indices)))
  union_filtered   <- sort(unique(unlist(effects_filtered$indices)))

  # Power (fraction of true causals hit at least once)
  power_unfiltered <- if (length(causal_idx) > 0)
    length(intersect(union_unfiltered, causal_idx)) / length(causal_idx) else NA_real_
  power_filtered <- if (length(causal_idx) > 0)
    length(intersect(union_filtered, causal_idx)) / length(causal_idx) else NA_real_

  # Overlap rate
  overlap_unfiltered <- overlap_rate_from_sets(effects_unfiltered$indices)
  overlap_filtered   <- overlap_rate_from_sets(effects_filtered$indices)

  # Effective #effects
  L_eff_unfiltered <- sum(effects_unfiltered$size > 0, na.rm = TRUE)
  L_eff_filtered   <- nrow(effects_filtered)

  # Variable-level PIPs (combined across effects)
  combined_pip <- combine_pips(A)

  # -------- classification-style metrics vs ground truth (if provided) --------
  labels <- rep(0L, ncol(X))
  if (length(causal_idx) > 0) labels[causal_idx] <- 1L

  auprc <- if (sum(labels) > 0) auprc_average_precision(combined_pip, labels) else NA_real_
  xent  <- if (length(labels) == length(combined_pip))
              cross_entropy_loss(combined_pip, labels) else NA_real_

  # -------- h_g^2 --------
  hg2 <- estimate_hg2(fit, y, X)

  # Optional traces if present
  elbo <- tryCatch(fit$model_fit$elbo, error = function(e) NULL)
  sigma2_path <- tryCatch(fit$model_fit$sigma_2, error = function(e) NULL)

  # Model-level summaries
  model_unfiltered <- data.frame(
    L_nominal      = L,
    L_effective    = L_eff_unfiltered,
    mean_size      = mean(effects_unfiltered$size, na.rm = TRUE),
    mean_purity    = mean(effects_unfiltered$purity, na.rm = TRUE),
    mean_coverage  = mean(effects_unfiltered$coverage, na.rm = TRUE),
    power          = power_unfiltered,
    overlap_rate   = overlap_unfiltered,
    AUPRC          = auprc,
    cross_entropy  = xent,
    hg2            = hg2,
    stringsAsFactors = FALSE
  )

  model_filtered <- data.frame(
    L_nominal      = L,
    L_effective    = L_eff_filtered,
    mean_size      = mean(effects_filtered$size, na.rm = TRUE),
    mean_purity    = mean(effects_filtered$purity, na.rm = TRUE),
    mean_coverage  = mean(effects_filtered$coverage, na.rm = TRUE),
    power          = power_filtered,
    overlap_rate   = overlap_filtered,
    AUPRC          = auprc,     # classification metrics don't change with filtering
    cross_entropy  = xent,
    hg2            = hg2,
    stringsAsFactors = FALSE
  )

  out <- list(
    params = list(rho = rho, purity_threshold = purity_threshold),
    # per-effect
    effects_unfiltered = effects_unfiltered,
    effects_filtered   = effects_filtered,
    # unions
    cs_union_indices_unfiltered = union_unfiltered,
    cs_union_indices_filtered   = union_filtered,
    # model summaries
    model_unfiltered = model_unfiltered,
    model_filtered   = model_filtered,
    # extras
    combined_pip = combined_pip,
    traces = list(elbo = elbo, sigma2 = sigma2_path)
  )

  if (!compute_curves) return(out)

  # -------- optional curves (useful for regression tests & calibration) --------
  # ρ-sweep coverage curve (empirical coverage vs rho)
  rho_grid <- seq(0.75, 0.99, by = 0.02)
  coverage_curve <- if (length(causal_idx) > 0) {
    sapply(rho_grid, function(rh) {
      cov_vec <- sapply(seq_len(L), function(l) {
        cs_l <- get_credible_set(A[l, ], rh)
        as.integer(length(cs_l) > 0 && any(cs_l %in% causal_idx))
      })
      mean(cov_vec, na.rm = TRUE)
    })
  } else rep(NA_real_, length(rho_grid))

  # PIP calibration (bin by combined PIP)
  if (length(causal_idx) > 0) {
    brks <- seq(0, 1, length.out = 21)  # 20 bins
    bins <- cut(combined_pip, brks, include.lowest = TRUE)
    pip_calibration <- aggregate(
      data.frame(pip = combined_pip,
                 is_causal = seq_along(combined_pip) %in% causal_idx),
      by = list(bin = bins),
      FUN = mean
    )
  } else {
    pip_calibration <- NULL
  }

  out$curves <- list(
    rho = rho_grid,
    coverage = coverage_curve,
    pip_calibration = pip_calibration
  )
  out
}

# ---------- convenience printer ----------

print_model_summary <- function(res) {
  cat("Credible-set evaluation (rho =", res$params$rho,
      ", purity_threshold =", res$params$purity_threshold, ")\n\n")

  cat("== Unfiltered ==\n")
  print(res$model_unfiltered, row.names = FALSE)

  cat("\n== Purity-filtered ==\n")
  print(res$model_filtered, row.names = FALSE)

  cat("\n#CS kept (filtered):", nrow(res$effects_filtered),
      " | #vars in CS union (filtered):", length(res$cs_union_indices_filtered), "\n")
}

# ---------- example (comment or adapt) ----------
# res <- evaluate_model(fit = fit_a_i, X = X, y = y, causal_idx = causal_idx,
#                       rho = 0.95, purity_threshold = 0.5, compute_curves = TRUE)
# print_model_summary(res)
# head(res$combined_pip)
