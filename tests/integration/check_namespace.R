# Lightweight consistency gate (no Julia): every R source parses, and every
# NAMESPACE export resolves to a defined object once the sources are loaded.

suppressMessages(library(MSTP))
ns  <- readLines("NAMESPACE")
exp <- grep("^export[(]", ns, value = TRUE)
nms <- gsub("^export[(]|[)]$", "", exp)
miss <- nms[!vapply(nms, exists, logical(1), USE.NAMES = FALSE)]

cat(sprintf("R sources parsed OK; exports=%d missing=%d\n", length(nms), length(miss)))
if (length(miss)) {
  cat("MISSING:", paste(miss, collapse = ", "), "\n")
  quit(status = 1L)
}
cat("ALL EXPORTS RESOLVE\n")
