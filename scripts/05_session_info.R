# 05_session_info.R
# Just captures the R version and package versions used.
# Anyone trying to reproduce this should match these versions
# — mixOmics in particular can behave differently across versions.

sink("results/session_info.txt")
cat("ROSMAP DIABLO Analysis — Session Info\n")
cat("======================================\n")
cat("Date run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
print(sessionInfo())
sink()

cat("Session info saved: results/session_info.txt\n")
