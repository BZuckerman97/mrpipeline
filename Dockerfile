FROM ghcr.io/rocker-org/devcontainer/r-ver:4.4

# Install R packages
RUN R -e 'install.packages(c("renv", "languageserver"), repos="https://cran.rstudio.com/")'
