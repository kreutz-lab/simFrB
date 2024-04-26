
<!-- README.md is generated from README.Rmd. Please edit that file -->

# simulationFromBenchmark

## DIMAR

DIMAR needs to be installed manually from Github

``` r
# Need to install dimar seperatly
install.packages("https://cran.r-project.org/src/contrib/Archive/imputation/imputation_1.3.tar.gz", repos=NULL, type='source')
#> The following package(s) will be installed:
#> - imputation [1.3]
#> These packages will be installed into "C:/Users/mengerj/AppData/Local/R/cache/R/renv/library/simulationFromBenchmark-3ea28dc4/R-4.3/x86_64-w64-mingw32".
#> 
#> # Installing packages --------------------------------------------------------
#> - Installing imputation ...                     OK [linked from cache]
#> Successfully installed 1 package in 17 milliseconds.

devtools::install_github("cran/DMwR")
#> Using GitHub PAT from the git credential store.
#> Skipping install of 'DMwR' from a github remote, the SHA1 (6fd4f0cd) has not changed since last install.
#>   Use `force = TRUE` to force installation

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(pkgs = c("pcaMethods", "impute", "SummarizedExperiment"))
#> 'getOption("repos")' replaces Bioconductor standard repositories, see
#> 'help("repositories", package = "BiocManager")' for details.
#> Replacement repositories:
#>     CRAN: https://packagemanager.posit.co/cran/latest
#> Bioconductor version 3.18 (BiocManager 1.30.22), R 4.3.3 (2024-02-29 ucrt)
#> Warning: package(s) not installed when version(s) same as or greater than current; use
#>   `force = TRUE` to re-install: 'pcaMethods' 'impute' 'SummarizedExperiment'
#> Old packages: 'openssl', 'pkgdown', 'renv', 'boot', 'codetools', 'lattice',
#>   'survival'

install.packages("devtools")
#> The following package(s) will be installed:
#> - devtools [2.4.5]
#> These packages will be installed into "C:/Users/mengerj/AppData/Local/R/cache/R/renv/library/simulationFromBenchmark-3ea28dc4/R-4.3/x86_64-w64-mingw32".
#> 
#> # Installing packages --------------------------------------------------------
#> - Installing devtools ...                       OK [linked from cache]
#> Successfully installed 1 package in 15 milliseconds.
devtools::install_github("kreutz-lab/DIMAR")
#> Using GitHub PAT from the git credential store.
#> Skipping install of 'DIMAR' from a github remote, the SHA1 (b77040ab) has not changed since last install.
#>   Use `force = TRUE` to force installation
```
