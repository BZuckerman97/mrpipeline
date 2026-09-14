# Package-level internal environment for session state (same pattern as
# codeminer's .codeminer_env). The namespace is locked on load but bindings
# inside this environment are not, so helpers can record run-time state here
# -- currently only the most recent allele orientation check, read back by
# last_allele_check().
.mrpipeline_env <- new.env(parent = emptyenv())
