#' Load a batch of instances from JSON files
#'
#' @param dir         Directory containing instance JSON files.
#' @param pattern     Glob pattern for file names (default all `*.json`).
#' @param n_instances Maximum number of instances to load (NULL = all).
#' @return Named list of instance lists, one per file.
#' @export
load_instances <- function(dir, pattern = "*.json", n_instances = NULL) {
  files <- sort(list.files(dir, pattern = glob2rx(pattern), full.names = TRUE))
  if (!is.null(n_instances)) files <- head(files, n_instances)
  lapply(files, jsonlite::fromJSON)
}
