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

#' Generate a correlated covariance matrix (R-side)
#'
#' Mirrors Julia's `gen_cov_mat`. Produces a block-diagonal covariance matrix
#' with uniform intra-block correlations and fixed cross-block correlation.
#'
#' @param n_blocks   Number of blocks.
#' @param block_size Nodes per block.
#' @param cross_corr Off-diagonal cross-block correlation (default 0.4).
#' @return A numeric matrix.
#' @export
gen_corrmat <- function(n_blocks, block_size, cross_corr = 0.4) {
  make_block <- function(n) {
    rho   <- sample(seq(0.6, 0.8, by = 0.1), 1L)
    block <- matrix(rho, nrow = n, ncol = n)
    diag(block) <- 1.0
    block
  }
  blocks <- replicate(n_blocks, make_block(block_size), simplify = FALSE)
  covmat <- as.matrix(Matrix::bdiag(blocks))
  covmat[covmat == 0.0] <- cross_corr
  covmat
}
