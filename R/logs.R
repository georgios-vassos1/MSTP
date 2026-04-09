#' Parse SDDP log files for a batch of instances
#'
#' Reads the summary block from each SDDP log file and extracts: total time,
#' numeric issues, best bound, simulation CI (mean ± half-width).
#'
#' @param log_dir     Directory containing log files named
#'                    `SDDPx{iterations}_{idx}.log`.
#' @param iterations  Number of SDDP iterations (used to build file names).
#' @param n_instances Number of instances (default 100).
#' @return Numeric matrix with columns:
#'   `[time, numeric_issues, bound, sim_mean, sim_ci]`
#' @export
parse_logs <- function(log_dir, iterations = 1500L, n_instances = 100L) {
  results <- matrix(NA_real_, nrow = n_instances, ncol = 5L,
                    dimnames = list(NULL, c("time", "numeric_issues",
                                            "bound", "sim_mean", "sim_ci")))

  kv_pat <- "^([^:]+):\\s*(.+)$"
  ci_pat <- "(.+) \\u00b1 (.+)"   # ±

  for (idx in seq_len(n_instances)) {
    log_file <- file.path(log_dir,
      sprintf("SDDPx%d_%03d.log", iterations, idx))
    if (!file.exists(log_file)) next

    lines <- readLines(log_file, warn = FALSE)

    # Find the summary block (lines after "status :")
    start <- grep("^status", lines)
    if (length(start) == 0L) next
    summary_lines <- lines[(start + 1L):(length(lines) - 2L)]

    col <- 1L
    for (line in summary_lines) {
      if (!grepl(kv_pat, line)) next
      m     <- regmatches(line, regexec(kv_pat, line))[[1L]]
      value <- trimws(m[3L])

      if (grepl("\u00b1", value)) {
        ci_m <- regmatches(value, regexec("(.+) \u00b1 (.+)", value))[[1L]]
        results[idx, col]     <- as.numeric(trimws(ci_m[2L]))
        results[idx, col + 1L] <- as.numeric(trimws(ci_m[3L]))
        col <- col + 2L
      } else {
        results[idx, col] <- suppressWarnings(as.numeric(value))
        col <- col + 1L
      }
    }
  }

  results
}

#' Summarise log results
#'
#' @param log_results Matrix from `parse_logs()`.
#' @return Named list with mean/CI for time, bound bias, and CI ratio.
#' @export
summarise_logs <- function(log_results) {
  n    <- nrow(log_results)
  ci95 <- function(x) 1.96 * stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

  bound    <- log_results[, "bound"]
  sim_mean <- log_results[, "sim_mean"]
  sim_ci   <- log_results[, "sim_ci"]
  time     <- log_results[, "time"]

  bias     <- (sim_mean - bound) / bound
  ci_ratio <- sim_ci / bound

  list(
    time      = c(mean = mean(time,     na.rm = TRUE), ci = ci95(time)),
    bias      = c(mean = mean(bias,     na.rm = TRUE), ci = ci95(bias)),
    ci_ratio  = c(mean = mean(ci_ratio, na.rm = TRUE), ci = ci95(ci_ratio)),
    n_issues  = sum(log_results[, "numeric_issues"] > 0L, na.rm = TRUE)
  )
}
