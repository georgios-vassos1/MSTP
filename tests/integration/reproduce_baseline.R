# Integration: the reproduction harness (R/reproduce.R) is deterministic and
# produces certified-valid bounds. Runs a small multi-spec table TWICE with the
# same seed and asserts the two tables are identical, every bound is valid, and
# the locked baseline values are reproduced. Run from the package root:
#   R_HOME=... KMP_DUPLICATE_LIB_OK=TRUE Rscript tests/integration/reproduce_baseline.R
#
# This is the fast (tiny-instance) gate. The paper-scale table is the same
# harness called with paper specs (see REPRODUCE.md); it is not run here.

suppressMessages(library(MSTP))
root <- normalizePath(".")
options(warn = 1)

# Tiny specs spanning two axes the paper varies (size and horizon).
specs <- list(
  list(label = "2x2x3", tau = 4L, nOrigins = 2L, nDestinations = 2L,
       nCarriers = 3L, lambda = 30.0, iters = 40L, trials = 120L, n_scenarios = 20L),
  list(label = "2x2x3-t6", tau = 6L, nOrigins = 2L, nDestinations = 2L,
       nCarriers = 3L, lambda = 30.0, iters = 40L, trials = 120L, n_scenarios = 20L)
)

t1 <- mstp_reproduce_bounds(specs, seed = 7L)
t2 <- mstp_reproduce_bounds(specs, seed = 7L)

cat("--- reproduction table (seed 7) ---\n")
print(t1, row.names = FALSE)

ok <- TRUE
check <- function(cond, msg) { cat(if (cond) "PASS: " else "FAIL: ", msg, "\n", sep=""); ok <<- ok && cond }

check(identical(t1, t2),          "table identical across two same-seed runs")
check(all(t1$valid),              "all bounds certified valid (within 3 SE)")
check(nrow(t1) == 2L,             "one row per spec")

# Regression lock: exact baseline values (filled in from the first verified run).
LOCK <- list("2x2x3" = 2263.808966, "2x2x3-t6" = 3930.404253)
for (lbl in names(LOCK)) {
  got <- t1$LB[t1$label == lbl]
  check(isTRUE(all.equal(got, LOCK[[lbl]], tolerance = 1e-5)),
        sprintf("LB baseline locked for %s (%.6f)", lbl, got))
}

cat(if (ok) ">>> M5 REPRODUCTION BASELINE PASSED\n" else ">>> M5 REPRODUCTION BASELINE FAILED\n")
quit(status = if (ok) 0L else 1L)
