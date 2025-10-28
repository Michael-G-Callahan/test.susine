# Run-control builders ------------------------------------------------------

#' Build a tidy grid of prior-quality settings.
#'
#' @param causal_levels Numeric vector (either fractions or percentages) for
#'   causal SNP prior noise.
#' @param noncausal_levels Numeric vector (fractions or percentages) for
#'   non-causal SNP prior noise.
#' @return Tibble with columns `prior_noise_causal`, `prior_noise_nonc`, and
#'   `prior_quality_id`.
#' @export
prior_quality_grid <- function(causal_levels = c(0.2, 0.4, 0.6, 0.8),
                               noncausal_levels = causal_levels) {
  normalize <- function(x) {
    x <- unique(x)
    x <- x[!is.na(x)]
    if (any(x > 1)) {
      x <- x / 100
    }
    x
  }
  c_levels <- normalize(causal_levels)
  nc_levels <- normalize(noncausal_levels)
  grid <- tidyr::expand_grid(
    prior_noise_causal = c_levels,
    prior_noise_nonc = nc_levels
  )
  dplyr::mutate(grid, prior_quality_id = dplyr::row_number())
}

#' Construct the full run table for a job.
#'
#' @param use_case_ids Character vector of use-case identifiers.
#' @param L_grid Integer vector of SuSiE/SuSiNE L values.
#' @param y_noise_grid Numeric vector of noise fractions (0-1).
#' @param prior_quality Tibble from [prior_quality_grid()].
#' @param p_star_grid Integer vector for number of causal SNPs.
#' @param seeds Integer vector of RNG seeds.
#' @param data_scenarios Character vector naming the data sources.
#' @return List with elements `scenarios`, `runs`, and `tasks`.
#' @keywords internal
make_run_tables <- function(use_case_ids,
                            L_grid,
                            y_noise_grid,
                            prior_quality,
                            p_star_grid,
                            seeds,
                            data_scenarios,
                            grid_mode = c("full", "minimal")) {
  if (!is.data.frame(prior_quality) ||
      !all(c("prior_noise_causal", "prior_noise_nonc") %in% names(prior_quality))) {
    stop("prior_quality must be a tibble with prior_noise_causal and prior_noise_nonc columns.")
  }
  grid_mode <- match.arg(grid_mode)
  use_cases <- resolve_use_cases(use_case_ids)
  if (!nrow(use_cases)) {
    stop("No valid use cases selected.")
  }

  build_full_grid <- function() {
    tidyr::expand_grid(
      data_scenario = data_scenarios,
      L = unique(L_grid),
      y_noise = unique(y_noise_grid),
      p_star = unique(p_star_grid),
      prior_quality_id = prior_quality$prior_quality_id
    ) %>%
      dplyr::left_join(prior_quality, by = "prior_quality_id") %>%
      dplyr::arrange(data_scenario, L, y_noise, p_star, prior_quality_id) %>%
      dplyr::mutate(scenario_id = dplyr::row_number())
  }

  build_minimal_grid <- function() {
    causal_vals <- unique(prior_quality$prior_noise_causal)
    noncausal_vals <- unique(prior_quality$prior_noise_nonc)
    if (!length(causal_vals)) {
      causal_vals <- 0
    }
    if (!length(noncausal_vals)) {
      noncausal_vals <- 0
    }
    values <- list(
      data_scenario = unique(data_scenarios),
      L = unique(L_grid),
      y_noise = unique(y_noise_grid),
      p_star = unique(p_star_grid),
      prior_noise_causal = causal_vals,
      prior_noise_nonc = noncausal_vals
    )
    lengths <- vapply(values, length, integer(1))
    n_rows <- max(lengths)
    tibble::tibble(
      data_scenario = rep_len(values$data_scenario, n_rows),
      L = rep_len(values$L, n_rows),
      y_noise = rep_len(values$y_noise, n_rows),
      p_star = rep_len(values$p_star, n_rows),
      prior_noise_causal = rep_len(values$prior_noise_causal, n_rows),
      prior_noise_nonc = rep_len(values$prior_noise_nonc, n_rows),
      prior_quality_id = seq_len(n_rows)
    ) %>%
      dplyr::mutate(scenario_id = seq_along(data_scenario))
  }

  scenarios <- if (grid_mode == "minimal") {
    build_minimal_grid()
  } else {
    build_full_grid()
  }

  seed_values <- unique(seeds)
  if (!length(seed_values)) {
    stop("At least one seed must be supplied.")
  }
  if (grid_mode == "minimal") {
    seed_values <- seed_values[1]
  }

  runs <- tidyr::expand_grid(
    scenario_id = scenarios$scenario_id,
    use_case_id = use_cases$use_case_id,
    seed = seed_values
  ) %>%
    dplyr::left_join(scenarios, by = "scenario_id") %>%
    dplyr::left_join(
      dplyr::select(use_cases, use_case_id, requires_prior_quality),
      by = "use_case_id"
    )

  default_scenarios <- scenarios %>%
    dplyr::group_by(data_scenario, L, y_noise, p_star) %>%
    dplyr::summarise(default_scenario_id = dplyr::first(scenario_id), .groups = "drop")

  runs <- runs %>%
    dplyr::left_join(default_scenarios,
      by = c("data_scenario", "L", "y_noise", "p_star")
    ) %>%
    dplyr::mutate(
      scenario_id = dplyr::if_else(
        requires_prior_quality,
        scenario_id,
        default_scenario_id
      ),
      prior_quality_id = dplyr::if_else(
        requires_prior_quality,
        prior_quality_id,
        as.integer(NA)
      ),
      prior_noise_causal = dplyr::if_else(
        requires_prior_quality,
        prior_noise_causal,
        as.numeric(NA)
      ),
      prior_noise_nonc = dplyr::if_else(
        requires_prior_quality,
        prior_noise_nonc,
        as.numeric(NA)
      )
    ) %>%
    dplyr::select(-default_scenario_id) %>%
    dplyr::distinct() %>%
    dplyr::arrange(scenario_id, use_case_id, seed) %>%
    dplyr::mutate(run_id = dplyr::row_number())

  scenarios <- dplyr::semi_join(scenarios, runs, by = "scenario_id")

  list(
    scenarios = scenarios,
    runs = runs,
    use_cases = use_cases
  )
}

#' Attach task identifiers to the run table.
#'
#' @param runs Run tibble from [make_run_tables()].
#' @param seeds_per_task Integer number of seeds executed inside each task.
#' @return Tibble with `task_id` column added and supporting task summary.
#' @keywords internal
assign_task_ids <- function(runs, seeds_per_task) {
  if (seeds_per_task < 1) {
    stop("seeds_per_task must be >= 1")
  }
  runs_aug <- runs %>%
    dplyr::group_by(scenario_id, use_case_id) %>%
    dplyr::mutate(
      seed_index = dplyr::row_number(),
      task_local = ((seed_index - 1L) %/% seeds_per_task) + 1L
    ) %>%
    dplyr::ungroup()

  tasks <- runs_aug %>%
    dplyr::distinct(scenario_id, use_case_id, task_local) %>%
    dplyr::arrange(scenario_id, use_case_id, task_local) %>%
    dplyr::mutate(task_id = dplyr::row_number())

  runs_final <- runs_aug %>%
    dplyr::left_join(tasks, by = c("scenario_id", "use_case_id", "task_local")) %>%
    dplyr::arrange(task_id, scenario_id, use_case_id, seed) %>%
    dplyr::mutate(run_id = dplyr::row_number()) %>%
    dplyr::select(-seed_index, -task_local)

  tasks_summary <- tasks %>%
    dplyr::left_join(
      runs_final %>%
        dplyr::count(task_id, name = "runs_per_task"),
      by = "task_id"
    )

  list(runs = runs_final, tasks = tasks_summary)
}

#' Create an in-memory job configuration.
#'
#' @param job_name Character scalar.
#' @param use_case_ids Character vector of selected use cases.
#' @param L_grid Integer vector of L values.
#' @param y_noise_grid Numeric vector of noise fractions.
#' @param prior_quality Tibble with prior noise settings.
#' @param p_star_grid Integer vector.
#' @param seeds Integer vector of RNG seeds.
#' @param data_scenarios Character vector.
#' @param seeds_per_task Integer seeds handled by each SLURM array task.
#' @param email Notification email address.
#' @param output_root Root directory for outputs (default `output`).
#' @param credible_set_rho Credible set cumulative PIP threshold.
#' @param purity_threshold Minimum purity to keep CS in filtered summary.
#' @param anneal_settings Named list for tempering runs.
#' @param model_average_settings Named list for model averaging runs.
#'
#' @return Job configuration list ready to serialize.
#' @export
make_job_config <- function(job_name,
                            use_case_ids,
                            L_grid,
                            y_noise_grid,
                            prior_quality,
                            p_star_grid,
                            seeds,
                            data_scenarios = "simulation_n3",
                            seeds_per_task = 1,
                            email = "mgc5166@psu.edu",
                            output_root = "output",
                            credible_set_rho = 0.95,
                            purity_threshold = 0.5,
                            grid_mode = c("full", "minimal"),
                            anneal_settings = list(
                              anneal_start_T = 5,
                              anneal_schedule_type = "geometric",
                              anneal_burn_in = 5
                            ),
                            model_average_settings = list(
                              n_inits = 5,
                              init_sd = 0.05
                            )) {
  grid_mode <- match.arg(grid_mode)
  tables <- make_run_tables(
    use_case_ids = use_case_ids,
    L_grid = L_grid,
    y_noise_grid = y_noise_grid,
    prior_quality = prior_quality,
    p_star_grid = p_star_grid,
    seeds = seeds,
    data_scenarios = data_scenarios,
    grid_mode = grid_mode
  )
  runs_tasks <- assign_task_ids(tables$runs, seeds_per_task = seeds_per_task)

  list(
    job = list(
      name = job_name,
      email = email,
      created_at = timestamp_utc(),
      seeds_per_task = seeds_per_task,
      credible_set_rho = credible_set_rho,
      purity_threshold = purity_threshold,
      compute = list(
        anneal = anneal_settings,
        model_average = model_average_settings
      ),
      slurm = list(
        time = "00:30:00",
        mem = "4G",
        cpus_per_task = 1,
        partition = NULL
      )
    ),
    paths = list(
      output_root = output_root,
      run_history_dir = file.path(output_root, "run_history", job_name),
      slurm_output_dir = file.path(output_root, "slurm_output", job_name),
      slurm_prints_dir = file.path(output_root, "slurm_prints"),
      slurm_scripts_dir = file.path(output_root, "slurm_scripts")
    ),
    tables = list(
      scenarios = tables$scenarios,
      runs = runs_tasks$runs,
      tasks = runs_tasks$tasks,
      use_cases = tables$use_cases
    )
  )
}

#' Write job configuration and scripts to disk.
#'
#' @param job_config Output of [make_job_config()].
#' @param run_task_script Path to the general task-runner script.
#' @return List with paths to the serialized artifacts.
#' @export
write_job_artifacts <- function(job_config,
                                run_task_script) {
  paths <- job_config$paths
  ensure_dir(paths$output_root)
  ensure_dir(paths$run_history_dir)
  ensure_dir(paths$slurm_output_dir)
  ensure_dir(paths$slurm_prints_dir)
  ensure_dir(paths$slurm_scripts_dir)

  job_json_path <- file.path(paths$run_history_dir, "job_config.json")
  jsonlite::write_json(
    job_config,
    path = job_json_path,
    auto_unbox = TRUE,
    digits = NA,
    pretty = TRUE
  )

  run_table_path <- file.path(paths$run_history_dir, "run_table.csv")
  readr::write_csv(job_config$tables$runs, run_table_path)

  scenario_path <- file.path(paths$run_history_dir, "scenario_table.csv")
  readr::write_csv(job_config$tables$scenarios, scenario_path)

  use_case_path <- file.path(paths$run_history_dir, "use_cases.csv")
  readr::write_csv(job_config$tables$use_cases, use_case_path)

  task_path <- file.path(paths$run_history_dir, "task_table.csv")
  readr::write_csv(job_config$tables$tasks, task_path)

  slurm_path <- file.path(paths$slurm_scripts_dir, paste0(job_config$job$name, ".slurm"))
  writeLines(
    render_slurm_script(job_config, run_task_script = run_task_script),
    con = slurm_path
  )

  list(
    job_config = job_json_path,
    run_table = run_table_path,
    scenarios = scenario_path,
    use_cases = use_case_path,
    tasks = task_path,
    slurm_script = slurm_path
  )
}

#' Render a SLURM submission script string.
#'
#' @param job_config Job configuration list.
#' @param run_task_script Path to `run_task.R`.
#' @return Character vector of script lines.
#' @keywords internal
render_slurm_script <- function(job_config, run_task_script) {
  job <- job_config$job
  paths <- job_config$paths
  tasks <- job_config$tables$tasks
  n_tasks <- nrow(tasks)
  slurm <- job$slurm

  partition_line <- if (!is.null(slurm$partition)) {
    paste0("#SBATCH --partition=", slurm$partition)
  } else {
    NULL
  }

  script <- c(
    "#!/bin/bash",
    sprintf("#SBATCH --job-name=%s", job$name),
    sprintf("#SBATCH --array=1-%d", n_tasks),
    sprintf("#SBATCH --time=%s", slurm$time),
    sprintf("#SBATCH --mem=%s", slurm$mem),
    sprintf("#SBATCH --cpus-per-task=%s", slurm$cpus_per_task),
    sprintf("#SBATCH --mail-user=%s", job$email),
    "#SBATCH --mail-type=BEGIN,END,FAIL",
    partition_line,
    sprintf("#SBATCH --output=%s/%%x-%%j-%%a.out", paths$slurm_prints_dir),
    sprintf("#SBATCH --error=%s/%%x-%%j-%%a.err", paths$slurm_prints_dir),
    "",
    "set -euo pipefail",
    "",
    "echo \"[$(date -Is)] Starting task ${SLURM_ARRAY_TASK_ID} for job ${SLURM_JOB_NAME}\"",
    sprintf("JOB_ROOT=\"%s\"", normalizePath(paths$output_root, winslash = "/", mustWork = FALSE)),
    sprintf("CONFIG_PATH=\"%s\"", normalizePath(file.path(paths$run_history_dir, "job_config.json"), winslash = "/", mustWork = FALSE)),
    sprintf("RUN_TASK_SCRIPT=\"%s\"", normalizePath(run_task_script, winslash = "/", mustWork = FALSE)),
    "",
    "Rscript \"$RUN_TASK_SCRIPT\" \\",
    "  --job-name \"$SLURM_JOB_NAME\" \\",
    "  --task-id \"$SLURM_ARRAY_TASK_ID\" \\",
    "  --job-root \"$JOB_ROOT\" \\",
    "  --config-path \"$CONFIG_PATH\"",
    "",
    "# Reorganize stdout/stderr into nested folders",
    "RAW_OUT_FILE=\"${SLURM_JOB_NAME}-${SLURM_JOB_ID}-${SLURM_ARRAY_TASK_ID}.out\"",
    "RAW_ERR_FILE=\"${SLURM_JOB_NAME}-${SLURM_JOB_ID}-${SLURM_ARRAY_TASK_ID}.err\"",
    sprintf("RAW_OUT_PATH=\"%s/${RAW_OUT_FILE}\"", paths$slurm_prints_dir),
    sprintf("RAW_ERR_PATH=\"%s/${RAW_ERR_FILE}\"", paths$slurm_prints_dir),
    sprintf("FINAL_DIR=\"%s/${SLURM_JOB_NAME}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}\"", paths$slurm_prints_dir),
    "mkdir -p \"$FINAL_DIR\"",
    "if [ -f \"$RAW_OUT_PATH\" ]; then mv \"$RAW_OUT_PATH\" \"$FINAL_DIR/stdout.out\"; fi",
    "if [ -f \"$RAW_ERR_PATH\" ]; then mv \"$RAW_ERR_PATH\" \"$FINAL_DIR/stderr.err\"; fi",
    "echo \"[$(date -Is)] Completed task ${SLURM_ARRAY_TASK_ID}\""
  )
  script[!is.na(script)]
}

#' Summarise expected run counts per use case and scenario.
#'
#' @param job_config Job configuration list.
#' @return Tibble summarising runs.
#' @export
summarise_job_config <- function(job_config) {
  runs <- job_config$tables$runs
  runs %>%
    dplyr::group_by(use_case_id, L, y_noise, p_star, prior_noise_causal, prior_noise_nonc) %>%
    dplyr::summarise(n_runs = dplyr::n(), .groups = "drop")
}
