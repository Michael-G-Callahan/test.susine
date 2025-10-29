# Task execution and modelling ----------------------------------------------

#' Run all simulation jobs assigned to a SLURM array task.
#'
#' @param job_name Character scalar.
#' @param task_id Integer SLURM array index (1-based).
#' @param job_root Root directory that houses `run_history` and outputs.
#' @param config_path Optional explicit path to `job_config.json`.
#' @param quiet Suppress console output when TRUE.
#'
#' @return Invisible tibble summarising run outcomes.
#' @export
run_task <- function(job_name,
                     task_id,
                     job_root = "output",
                     config_path = NULL,
                     quiet = FALSE) {
  cfg_path <- config_path %||% file.path(job_root, "run_history", job_name, "job_config.json")
  job_config <- load_job_config(cfg_path)
  runs_tbl <- job_config$tables$runs
  task_id <- as.integer(task_id)
  task_runs <- dplyr::filter(runs_tbl, task_id == !!task_id)
  if (!nrow(task_runs)) {
    stop(sprintf("Task %s not found in job config.", task_id))
  }

  if (!quiet) {
    message(sprintf("Running job '%s' task %s with %d run(s).", job_name, task_id, nrow(task_runs)))
  }

  results <- purrr::map_dfr(seq_len(nrow(task_runs)), function(i) {
    run_row <- task_runs[i, , drop = FALSE]
    execute_single_run(run_row, job_config, quiet = quiet)
  })

  invisible(results)
}

#' Load a job configuration JSON file and coerce tibbles.
#' @keywords internal
load_job_config <- function(path) {
  cfg <- jsonlite::read_json(path, simplifyVector = TRUE)
  cfg$tables$runs <- tibble::as_tibble(cfg$tables$runs)
  cfg$tables$scenarios <- tibble::as_tibble(cfg$tables$scenarios)
  cfg$tables$tasks <- tibble::as_tibble(cfg$tables$tasks)
  cfg$tables$use_cases <- tibble::as_tibble(cfg$tables$use_cases)
  cfg
}

#' Retrieve metadata row for a given use case ID.
#' @keywords internal
lookup_use_case <- function(job_config, use_case_id) {
  dplyr::filter(job_config$tables$use_cases, use_case_id == !!use_case_id)
}

#' Execute a single run row: simulate data, fit model, and persist outputs.
#' @keywords internal
execute_single_run <- function(run_row, job_config, quiet = FALSE) {
  run_row <- tibble::as_tibble(run_row)
  use_case <- lookup_use_case(job_config, run_row$use_case_id)
  if (!nrow(use_case)) {
    stop(sprintf("Unknown use_case_id '%s'.", run_row$use_case_id))
  }

  data_bundle <- generate_data_for_run(run_row, job_config)

  model_result <- run_use_case(
    use_case = use_case,
    run_row = run_row,
    data_bundle = data_bundle,
    job_config = job_config
  )

  eval_res <- evaluate_model(
    fit = model_result$fit,
    X = data_bundle$X,
    y = data_bundle$y,
    causal_idx = data_bundle$causal_idx,
    rho = job_config$job$credible_set_rho %||% 0.95,
    purity_threshold = job_config$job$purity_threshold %||% 0.5,
    compute_curves = FALSE
  )

  write_run_outputs(
    run_row = run_row,
    job_config = job_config,
    evaluation = eval_res,
    data_bundle = data_bundle,
    model_result = model_result
  )

  if (!quiet) {
    summary_msg <- sprintf(
      "run_id=%s use_case=%s seed=%s L=%s | power=%.3f size=%.2f purity=%.2f",
      run_row$run_id,
      run_row$use_case_id,
      run_row$seed,
      run_row$L,
      eval_res$model_filtered$power,
      eval_res$model_filtered$mean_size,
      eval_res$model_filtered$mean_purity
    )
    message(summary_msg)
  }

  tibble::tibble(
    run_id = run_row$run_id,
    task_id = run_row$task_id,
    use_case_id = run_row$use_case_id,
    seed = run_row$seed,
    power = eval_res$model_filtered$power,
    mean_size = eval_res$model_filtered$mean_size,
    mean_purity = eval_res$model_filtered$mean_purity
  )
}

#' Generate data bundle for a given run_row based on scenario.
#' @keywords internal
generate_data_for_run <- function(run_row, job_config) {
  scenario <- run_row$data_scenario
  if (scenario %in% c("simulation_n3", "sim_n3")) {
    spec <- list(
      seed = as.integer(run_row$seed),
      p_star = as.integer(run_row$p_star),
      y_noise = as.numeric(run_row$y_noise),
      prior_noise_causal = as.numeric(run_row$prior_noise_causal),
      prior_noise_nonc = as.numeric(run_row$prior_noise_nonc),
      effect_sd = run_row[["effect_sd"]] %||% 1
    )
    return(generate_simulation_data(spec))
  }

  stop(sprintf("Unsupported data_scenario '%s'.", scenario))
}

#' Fit a use case with the susine backend, handling optional annealing/model averaging.
#' @keywords internal
run_use_case <- function(use_case, run_row, data_bundle, job_config) {
  L <- as.integer(run_row$L)
  mu_strategy <- use_case$mu_strategy[[1]]
  sigma_strategy <- use_case$sigma_strategy[[1]]

  mu_0 <- if (mu_strategy %in% c("functional", "eb_mu")) data_bundle$mu_0 else 0
  sigma_0_2 <- if (sigma_strategy %in% c("functional", "eb_sigma")) data_bundle$sigma_0_2 else NULL
  pi_weights <- if (identical(sigma_strategy, "functional")) data_bundle$prior_inclusion_weights else NULL

  extra <- use_case$extra_compute[[1]]
  if (length(extra) && is.na(extra)) {
    extra <- NULL
  }
  anneal <- job_config$job$compute$anneal
  model_avg <- job_config$job$compute$model_average

  base_args <- list(
    L = L,
    X = data_bundle$X,
    y = data_bundle$y,
    mu_0 = mu_0,
    sigma_0_2 = sigma_0_2,
    prior_inclusion_weights = pi_weights,
    prior_update_method = use_case$prior_update_method[[1]] %||% "none",
    verbose = FALSE
  )
  if (resolve_flag(use_case$auto_scale_mu[[1]], FALSE)) {
    base_args$auto_scale_mu_0 <- TRUE
  }
  if (resolve_flag(use_case$auto_scale_sigma[[1]], FALSE)) {
    base_args$auto_scale_sigma_0_2 <- TRUE
  }

  if (identical(extra, "anneal")) {
    base_args$anneal_start_T <- anneal$anneal_start_T %||% 5
    base_args$anneal_schedule_type <- anneal$anneal_schedule_type %||% "geometric"
    base_args$anneal_burn_in <- anneal$anneal_burn_in %||% 5
  }

  if (!identical(extra, "model_avg")) {
    fit <- do.call(susine::susine, base_args)
    fit$settings$L <- fit$settings$L %||% L
    fit$metadata <- list(
      use_case_id = use_case$use_case_id[[1]],
      label = use_case$label[[1]]
    )
    return(list(fit = fit, extra = NULL))
  }

  n_inits <- model_avg$n_inits %||% 5
  init_sd <- model_avg$init_sd %||% 0.05
  subfits <- vector("list", n_inits)
  pips_mat <- matrix(0, nrow = n_inits, ncol = length(data_bundle$beta))
  coef_mat <- matrix(0, nrow = n_inits, ncol = length(data_bundle$beta))

  for (i in seq_len(n_inits)) {
    base_args$init_random <- TRUE
    base_args$init_seed <- as.integer(run_row$seed) + i
    base_args$init_random_sd <- init_sd
    subfits[[i]] <- do.call(susine::susine, base_args)
    pips_mat[i, ] <- subfits[[i]]$model_fit$PIPs
    coef_mat[i, ] <- subfits[[i]]$model_fit$coef
  }

  agg_fit <- subfits[[1]]
  agg_fit$model_fit$PIPs <- colMeans(pips_mat)
  agg_fit$model_fit$coef <- colMeans(coef_mat)
  agg_fit$settings$L <- agg_fit$settings$L %||% L
  agg_fit$model_fit$std_coef <- colMeans(vapply(subfits, function(sf) sf$model_fit$std_coef, numeric(ncol(coef_mat))))
  agg_fit$model_fit$fitted_y <- as.vector(data_bundle$X %*% agg_fit$model_fit$coef)
  agg_fit$metadata <- list(
    use_case_id = use_case$use_case_id[[1]],
    label = use_case$label[[1]],
    model_averaging = list(
      n_inits = n_inits,
      init_sd = init_sd,
      seeds = as.integer(run_row$seed) + seq_len(n_inits)
    )
  )

  list(fit = agg_fit, extra = subfits)
}

#' Persist metrics, truth tables, and model fits for a run.
#' @keywords internal
write_run_outputs <- function(run_row,
                              job_config,
                              evaluation,
                              data_bundle,
                              model_result) {
  run_id <- sprintf("run-%05d", run_row$run_id)
  task_dir <- file.path(job_config$paths$slurm_output_dir, sprintf("task-%03d", run_row$task_id))
  run_dir <- file.path(task_dir, run_id)
  ensure_dir(run_dir)

  model_metrics <- dplyr::mutate(
    dplyr::bind_rows(
      dplyr::mutate(evaluation$model_unfiltered, filtering = "unfiltered"),
      dplyr::mutate(evaluation$model_filtered, filtering = "purity_filtered")
    ),
    run_id = run_row$run_id,
    task_id = run_row$task_id,
    use_case_id = run_row$use_case_id,
    seed = run_row$seed
  )
  readr::write_csv(model_metrics, file.path(run_dir, "model_metrics.csv"))

  format_effects <- function(effects_df, label) {
    if (is.null(effects_df) || !nrow(effects_df)) {
      return(NULL)
    }
    indices <- vapply(effects_df$indices, function(idx) paste(idx, collapse = " "), character(1))
    dplyr::mutate(
      dplyr::select(effects_df, -indices),
      indices = indices,
      filtering = label,
      run_id = run_row$run_id,
      task_id = run_row$task_id,
      use_case_id = run_row$use_case_id,
      seed = run_row$seed
    )
  }
  effect_metrics <- dplyr::bind_rows(
    format_effects(evaluation$effects_unfiltered, "unfiltered"),
    format_effects(evaluation$effects_filtered, "purity_filtered")
  )
  if (!is.null(effect_metrics) && nrow(effect_metrics)) {
    readr::write_csv(effect_metrics, file.path(run_dir, "effect_metrics.csv"))
  }

  pip_tbl <- tibble::tibble(
    snp_index = seq_along(evaluation$combined_pip),
    pip = evaluation$combined_pip,
    run_id = run_row$run_id,
    task_id = run_row$task_id,
    use_case_id = run_row$use_case_id,
    seed = run_row$seed
  )
  readr::write_csv(pip_tbl, file.path(run_dir, "pip.csv"))

  truth_tbl <- tibble::tibble(
    snp_index = seq_along(data_bundle$beta),
    beta = data_bundle$beta,
    mu_0 = data_bundle$mu_0,
    sigma_0_2 = data_bundle$sigma_0_2,
    prior_weight = data_bundle$prior_inclusion_weights,
    causal = as.integer(seq_along(data_bundle$beta) %in% data_bundle$causal_idx)
  )
  readr::write_csv(truth_tbl, file.path(run_dir, "truth.csv"))

  saveRDS(model_result$fit, file.path(run_dir, "fit.rds"))
  if (!is.null(model_result$extra)) {
    saveRDS(model_result$extra, file.path(run_dir, "subfits.rds"))
  }
}
