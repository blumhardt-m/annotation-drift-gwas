# Adapter around the retained viability script's functions; no source edits.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3L)
suppressMessages(library(fgsea))
env <- new.env(parent=globalenv())
expressions <- parse(args[2])
wanted <- c("read_stats", "read_pathways", "run_version")
for (expr in expressions) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% wanted) {
    eval(expr, envir=env)
  }
}
stopifnot(all(vapply(wanted, exists, logical(1), envir=env, inherits=FALSE)))
env$OUT <- args[1]
env$pathways <- env$read_pathways(file.path(env$OUT, "pathways_ensembl.tsv"))
# The imported runner sorts the preordered input and makes one seed call.
# Bind that call to the literal v1.2 seed and explicit multilevel parameters.
env$set.seed <- function(seed) {
  stopifnot(seed == 42)
  base::set.seed(42L)
}
env$fgsea <- function(pathways, stats, minSize, maxSize, eps) {
  stopifnot(minSize == 15, maxSize == 500, eps == 0)
  result <- fgsea::fgseaMultilevel(pathways=pathways, stats=stats, minSize=15,
                        maxSize=500, eps=0, sampleSize=101, nproc=1,
                        scoreType="pos", nPermSimple=1000)
  # Preserve every returned value before the imported display selects columns.
  saveRDS(result, file.path(env$OUT, sprintf("fgsea_%s_all_values.rds", args[3])))
  exported <- as.data.frame(result)
  exported$leadingEdge <- vapply(exported$leadingEdge, paste, collapse=",", character(1))
  write.table(exported, file=file.path(env$OUT, sprintf("fgsea_%s_all_values.tsv", args[3])),
              sep="\t", quote=TRUE, row.names=FALSE, na="NA")
  result
}
env$run_version(args[3])
