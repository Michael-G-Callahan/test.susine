# Utility helpers -----------------------------------------------------------

#' Coalesce-like helper for NULL/length-0 inputs.
#'
#' Returns the left-hand side if it is non-null and has length > 0,
#' otherwise returns the right-hand side.
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) {
    return(y)
  }
  x
}

#' Ensure a directory exists (creating it recursively when needed).
#'
#' @param path Character scalar path.
#' @return The input path (invisibly), after creating it if necessary.
#' @keywords internal
ensure_dir <- function(path) {
  stopifnot(length(path) == 1L, !is.na(path))
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

#' Format a timestamp in ISO8601 UTC.
#'
#' @return Character scalar timestamp.
#' @keywords internal
timestamp_utc <- function() {
  format(Sys.time(), tz = "UTC", usetz = TRUE)
}

#' Convert a logical flag to scalar TRUE/FALSE, applying default when NA.
#' @param x Logical input.
#' @param default Logical default used when x is NA.
#' @keywords internal
resolve_flag <- function(x, default = FALSE) {
  if (length(x) == 0L || is.na(x)) {
    return(isTRUE(default))
  }
  isTRUE(x)
}
