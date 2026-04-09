## Reproduction script for arxiv:2505.01813
## "Multistage Stochastic Transportation Problem via SDDP"
##
## Requires: MSTP package installed (or sourced), TLPR, data.table, ggplot2
## Data:     ~/tsproj/logs/   -- SDDP log files (SDDPx{500,1000,1500}_NNN.log)
##           ~/tsproj/output/ -- precomputed RDS files
##           ~/tsproj/instances/ -- instance JSON files

library(MSTP)

LOG_DIR  <- "~/tsproj/logs"
OUT_DIR  <- "~/tsproj/output"
INST_DIR <- "~/tsproj/instances"

# ─────────────────────────────────────────────────────────────────────────────
# Table 4.1 — SDDP performance across iteration counts
# Columns: bias (%), simulation CI ratio (%), wall time (s)
# ─────────────────────────────────────────────────────────────────────────────
cat("=== Table 4.1 ===\n")
table_rows <- lapply(c(500L, 1000L, 1500L), function(iters) {
  res <- parse_logs(LOG_DIR, iterations = iters, n_instances = 100L)
  s   <- summarise_logs(res)
  c(iterations  = iters,
    bias_pct    = round(s$bias["mean"]    * 100, 1),
    bias_ci     = round(s$bias["ci"]      * 100, 1),
    ci_ratio_pct= round(s$ci_ratio["mean"]* 100, 1),
    ci_ratio_ci = round(s$ci_ratio["ci"]  * 100, 1),
    time_s      = round(s$time["mean"],     2),
    time_ci     = round(s$time["ci"],       2),
    n_issues    = s$n_issues)
})
tab <- do.call(rbind, table_rows)
print(tab)
cat("\n")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 4.3 — Sensitivity analysis (inflow and spot rates)
# SensvInflow.RDS / SensvSpot.RDS: flat numeric vectors of length 100*1000
# reshape to (1000 samples × 100 instances) for plot_sensitivity()
# ─────────────────────────────────────────────────────────────────────────────
cat("=== Figure 4.3 — Sensitivity Analysis ===\n")

inflow_vec <- readRDS(file.path(OUT_DIR, "SensvInflow.RDS"))
spot_vec   <- readRDS(file.path(OUT_DIR, "SensvSpot.RDS"))

inflow_mat <- matrix(inflow_vec, nrow = 1000L)   # 1000 samples × 100 instances
spot_mat   <- matrix(spot_vec,   nrow = 1000L)

pdf("fig4.3_sensitivity.pdf", width = 12, height = 5)
par(mfrow = c(1L, 2L))
plot_sensitivity(inflow_mat, title = "100 Instances, 1000 Inflow Samples")
plot_sensitivity(spot_mat,   title = "100 Instances, 1000 Spot Rate Samples")
dev.off()
cat("Written: fig4.3_sensitivity.pdf\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 4.4 — Regret distributions at three iteration counts
# regrets100x1000.RDS: named list (SDDPx500/1000/1500), each 100×1000 matrix
# plot_regret() applies density by row (each row = 1000 trials for one instance)
# ─────────────────────────────────────────────────────────────────────────────
cat("=== Figure 4.4 — Regret Distributions ===\n")

regrets <- readRDS(file.path(OUT_DIR, "regrets100x1000.RDS"))

pdf("fig4.4_regrets.pdf", width = 15, height = 5)
plot_regret(regrets, max_regret = 1.5e-4)
dev.off()
cat("Written: fig4.4_regrets.pdf\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 4.1 / 4.2 — Policy visualization for a 2x2x2x2 instance
# Requires JuliaCall + MSTP Julia engine + TLPR
# ─────────────────────────────────────────────────────────────────────────────
cat("=== Figure 4.1/4.2 — Policy Visualization (requires Julia) ===\n")

tryCatch({
  setup_engine()

  # Small 2×2×2×2 instance for policy visualization
  inst_path <- file.path(INST_DIR, "instance_12x2x2x2_001.json")
  if (!file.exists(inst_path)) {
    cat("instance_12x2x2x2_001.json not found, skipping Figure 4.1/4.2\n")
    stop("skip", call. = FALSE)
  }

  inst   <- jsonlite::fromJSON(inst_path)
  config <- mstp_config(inst)
  model  <- mstp_train(config, iterations = 500L)
  sims   <- mstp_simulate(model, config, trials = 1L)

  cat("Policy simulation complete (1 trial).\n")
  cat("obj =", sims$obj, "\n")
  cat("(Full policy plot requires TLPR::plot2x2instance — see tsproj/R/policy.R)\n\n")
}, error = function(e) {
  if (conditionMessage(e) != "skip")
    message("Julia engine not available: ", conditionMessage(e))
})

cat("Reproduction complete.\n")
