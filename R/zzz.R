# Package-level internal environment for session state (same pattern as
# codeminer's .codeminer_env). The namespace is locked on load but bindings
# inside this environment are not, so helpers can record run-time state here
# -- currently only the most recent allele orientation check, read back by
# last_allele_check().
.mrpipeline_env <- new.env(parent = emptyenv())

# Load ieugwasr now, with the caller's RNG state put back afterwards
# (issue #43). ieugwasr's own .onLoad() calls sample() to build a googleAuthR
# security code, and NAMESPACE imports nothing from ieugwasr, so without
# this the namespace loads inside the first run_mr() of a session --
# after the caller's set.seed() -- and that call's weighted-median bootstrap
# SE and TwoSampleMR id.outcome differ from an identical later call.
# ieugwasr is in Imports, its load is silent, and run_mr() needs it on
# almost every path anyway (clumping, LD).
.onLoad <- function(libname, pkgname) {
  withr::with_preserve_seed(requireNamespace("ieugwasr", quietly = TRUE))
  invisible()
}
