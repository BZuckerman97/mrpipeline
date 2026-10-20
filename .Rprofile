local({
  if (nzchar(Sys.getenv("CI"))) {
    # R reads only one of ./.Rprofile and ~/.Rprofile, so this file hides the
    # ~/.Rprofile that r-lib/actions/setup-r writes: its RSPM repository
    # (Posit Package Manager's Linux *binary* URL) and the HTTPUserAgent that
    # makes that URL serve binaries. Without them every Ubuntu job compiled
    # every dependency from source (duckdb alone: ~34 min). Source it first.
    user_profile <- path.expand("~/.Rprofile")
    if (file.exists(user_profile)) {
      source(user_profile)
    }
    repos <- getOption("repos")
  } else {
    repos <- c(
      CRAN = "https://packagemanager.posit.co/cran/latest",
      cloudr = "https://cloud.r-project.org"
    )
  }
  repos["mrcieu"] <- "https://mrcieu.r-universe.dev"
  options(repos = repos)
})
