local({
  repos <- getOption("repos")
  # On CI, r-lib/actions/setup-r has already pointed CRAN at Posit Package
  # Manager's Linux *binary* URL (.../__linux__/<distro>/latest). Replacing it
  # with the generic .../cran/latest, which serves Linux only source, made
  # every Ubuntu job compile every dependency (duckdb alone: ~34 min).
  if (!nzchar(Sys.getenv("CI"))) {
    repos["CRAN"] <- "https://packagemanager.posit.co/cran/latest"
    repos["cloudr"] <- "https://cloud.r-project.org"
  }
  repos["mrcieu"] <- "https://mrcieu.r-universe.dev"
  options(repos = repos)
})
