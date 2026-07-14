# ============================================================
# RNA-seq 流程 — R 脚本公共函数
# 被 deseq2.R / clusterprofiler.R / gene2symbol.R source()
# ============================================================

# If R_LIBS_ONLY is set, use only that library path (for conda/container isolation)
if (nzchar(Sys.getenv("R_LIBS_ONLY"))) {
  .libPaths(Sys.getenv("R_LIBS_ONLY"))
}

parse_arg <- function(args, flag) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(NULL)
  if (idx == length(args)) {
    stop(sprintf("Missing value for flag '%s' (cannot be the last argument)", flag))
  }
  args[idx + 1]
}
