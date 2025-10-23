#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (startsWith(arg, "--")) {
      key_val <- sub("^--", "", arg)
      if (grepl("=", key_val, fixed = TRUE)) {
        kv <- strsplit(key_val, "=", fixed = TRUE)[[1]]
        out[[kv[1]]] <- kv[2]
      } else {
        if (i == length(args)) stop("Missing value for argument ", arg)
        out[[key_val]] <- args[[i + 1L]]
        i <- i + 1L
      }
    }
    i <- i + 1L
  }
  out
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  job_name <- args[["job-name"]]
  task_id <- args[["task-id"]]
  config_path <- args[["config-path"]]
  job_root <- args[["job-root"]]
  if (is.null(job_root)) job_root <- "output"
  quiet_flag <- args[["quiet"]]
  if (is.null(quiet_flag)) quiet_flag <- "FALSE"

  if (is.null(job_name) || is.null(task_id)) {
    stop("Arguments --job-name and --task-id are required.")
  }

  task_id <- as.integer(task_id)
  if (is.na(task_id)) {
    stop("task-id must be an integer.")
  }

  suppressPackageStartupMessages(library(test_susine))

  if (is.null(config_path)) {
    config_path <- file.path(job_root, "run_history", job_name, "job_config.json")
  }

  quiet <- tolower(quiet_flag) %in% c("1", "true", "yes")

  tryCatch(
    {
      test_susine::run_task(
        job_name = job_name,
        task_id = task_id,
        job_root = job_root,
        config_path = config_path,
        quiet = quiet
      )
      invisible(NULL)
    },
    error = function(e) {
      message("Error during task execution: ", conditionMessage(e))
      quit(save = "no", status = 1L)
    }
  )
}

main()
