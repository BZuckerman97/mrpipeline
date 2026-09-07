#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(".data", ".env"))

# Package-level mutable state (the httr2 `the` pattern). The namespace is
# locked on load but bindings inside this environment are not, so helpers can
# record run-time state here -- currently only the most recent allele
# orientation check, read back by last_allele_check().
the <- new.env(parent = emptyenv())
