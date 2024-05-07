
<!-- README.md is generated from README.Rmd. Please edit that file -->

# simFrB

### Why this simulation Tool?

Simulated omics data often serves as a basis to evaluate or create
analysis tools, plan experiments and get a deeper understanding of the
data. Simulations usually make many assumptions about the distribution
of the data and often times simplify complex relationships observed in
experimental data.

One major challenge in simulating omics data is simulating missing value
patterns. Simulations might take data dependencies into consideration by
including values missing non at random, but we found that such
approaches still fail to reproduce results observed in experimental
data.

While trying to validate results of statistical analysis tools used in
Fröhlich et. al (), we observed a strong influence of the missing value
pattern on statistical results if data was prepossessed with imputation,
as done in the SAM method. While for experimental data, imputation
strongly increased the pAUC of all statistical tests, it decreased the
pAUC for simulated data.

This divergence was found to be caused by **feature dependent
correlations of the missing value pattern**, that got lost in other
simulation approaches. The implemented simulation tool leverages allows
to leverage benchmark data in order to simulate new count matrices with
realistic missing value patterns, by estimating coefficients for both
intensity (linear model) and missingness (logit model) for each feature
and drawing jointly to simulate data that keeps correlations found in
the experimental data.

## Usage

For now this package can be installed from this github repository

``` r
library(remotes, quietly = T)
remotes::install_github("kreutz-lab/simFrB")
#> Using GitHub PAT from the git credential store.
#> Skipping install of 'simFrB' from a github remote, the SHA1 (b5c3eaa4) has not changed since last install.
#>   Use `force = TRUE` to force installation
```

Generally there or two modes of operation:

1.  Use Available Data to estimate coefficients from. This will be most
    accurate for Benchmark data with known ground truth.

``` r
library(simFrB)
#load Experimental Data
exp.df.all <-
  simFrB::loadDataFromGit(
    "https://raw.githubusercontent.com/kreutz-lab/dia-benchmarking/main/data/diaWorkflowResults_allDilutions.rds")

exp.df <- simFrB::subsetDIAWorkflowData(
  data = exp.df.all,
  DIAWorkflow = "DIANN_DIANN_AI_GPF",
  experimentalComparisonGroups = c("1-12","1-25"),
  rowSubset = seq(1,500),
  colSubset = c(seq(1,6),seq(24,29))
)

exp.df <- exp.df[rowSums(exp.df, na.rm = TRUE) > 0, ]
exp.mtx <- as.matrix(exp.df)
DE_idx <- grep("ECOLI",rownames(exp.df)) #This needs to be adapted depending on how you know which feature is differentially expressed

sim <- simFrB::msb.applyFunctionWithSeed(simFrB::msb.simulateDataFromBenchmark,
                                         seed = 123,
                                         mtx = exp.mtx,
                                      DE_idx = DE_idx,
                                      int.mean = mean(exp.mtx,na.rm = T),
                                      nFeatures = 1000,
                                      nSamples = 10)
#> No groupDesign for mtx provided, assuming first half of samples
#>                 are in group 1 and second half in group 2. If this is not the case
#>                 please provide groupDesign_mtx.
#> fitting from feature 1 to 454
#> using input data to estimate feature standard deviations.
#> By default DIMAR is assuming that the first half of the columns belong to group 1 and the second half to group 2.
#>           If this is not the case, please provide the group vector to the dimarConstructDesignMatrix function.
#> [1] "Pattern of MVs is learned by logistic regression."
#> Joining with `by = join_by(rowIDs)`
#> applying dimar coefficients from 1 to 334
#> applying dimar coefficients from 335 to 667
#> applying dimar coefficients from 668 to 1000
```

1.  Use Pre calculated Coefficients stored in data/jointCoefs.rds

``` r
load("data/jointCoefs.rda")

sim <- simFrB::msb.simulateDataFromBenchmark(jointCoefs = jointCoefs,
                                      nFeatures = 6000,
                                      nSamples = 30)
#> No input matrix provided
#> applying dimar coefficients from 1 to 316
#> applying dimar coefficients from 317 to 632
#> applying dimar coefficients from 633 to 948
#> applying dimar coefficients from 949 to 1264
#> applying dimar coefficients from 1265 to 1579
#> applying dimar coefficients from 1580 to 1895
#> applying dimar coefficients from 1896 to 2211
#> applying dimar coefficients from 2212 to 2527
#> applying dimar coefficients from 2528 to 2843
#> applying dimar coefficients from 2844 to 3158
#> applying dimar coefficients from 3159 to 3474
#> applying dimar coefficients from 3475 to 3790
#> applying dimar coefficients from 3791 to 4106
#> applying dimar coefficients from 4107 to 4422
#> applying dimar coefficients from 4423 to 4737
#> applying dimar coefficients from 4738 to 5053
#> applying dimar coefficients from 5054 to 5369
#> applying dimar coefficients from 5370 to 5685
#> applying dimar coefficients from 5686 to 6000
```

The object returned by msb.simulateDataFromBenchmark contains the
original coefficients, as well as the drawn coefficients that were used
the simulate the given data as attributes.

``` r
str(attributes(sim)$expCoefs)
#> List of 3
#>  $ rowCoefs      :'data.frame':  8742 obs. of  7 variables:
#>   ..$ rowIDs           : chr [1:8742] "rowID1" "rowID2" "rowID3" "rowID4" ...
#>   ..$ lmCoefs          : num [1:8742] 4.7 3.87 2.73 3.31 2.97 ...
#>   ..$ lmCoefs_sd       : num [1:8742] 0.288 0.192 0.256 0.271 0.341 ...
#>   ..$ ifDE             : chr [1:8742] "notDE" "notDE" "notDE" "notDE" ...
#>   ..$ FCCoefs          : num [1:8742] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ dimarCoefs_group1: num [1:8742] 4.63 3.63 3.67 3.88 4.09 ...
#>   ..$ dimarCoefs_group2: num [1:8742] 4.63 3.63 3.67 3.88 4.09 ...
#>   ..- attr(*, "dimarXtype")= num [1:1008] 0 1 2 2 2 2 2 2 2 2 ...
#>  $ colCoefs      :'data.frame':  46 obs. of  2 variables:
#>   ..$ lmCoefs   : num [1:46] 0 -0.0141 -0.0445 -0.0584 -0.0901 ...
#>   ..$ dimarCoefs: num [1:46] -7.53 -8.02 -7.2 -7.09 -7.92 ...
#>  $ dimarIntercept: Named num [1:2] 0 -3.86
#>   ..- attr(*, "names")= chr [1:2] "Intercept" "mean"
```

``` r
str(attributes(sim)$usedCoefs)
#> List of 3
#>  $ newLmCoefs   :List of 4
#>   ..$ featureCoefs: Named num [1:6000] -1.6 2.39 -2.63 1.57 -2.64 ...
#>   .. ..- attr(*, "names")= chr [1:6000] "rowID6126" "rowID2993_DE" "rowID727" "rowID134" ...
#>   ..$ sampleCoefs : Named num [1:30] -0.0542 -0.0846 -0.1324 -0.1713 -0.0602 ...
#>   .. ..- attr(*, "names")= chr [1:30] "colID45" "colID44" "colID23...3" "colID38...4" ...
#>   ..$ FCCoefs     : Named num [1:6000] 0 -1.4 0 0 0 ...
#>   .. ..- attr(*, "names")= chr [1:6000] "rowID6126" "rowID2993_DE" "rowID727" "rowID134" ...
#>   ..$ featureSd   : Named num [1:6000] 0.675 0.551 1.054 0.315 0.945 ...
#>   .. ..- attr(*, "names")= chr [1:6000] "rowID6126" "rowID2993_DE" "rowID727" "rowID134" ...
#>  $ newDimarCoefs: Named num [1:12032] 0 -3.86 -7.23 -7.55 -7.16 ...
#>   ..- attr(*, "names")= chr [1:12032] "Intercept" "mean" "colID45" "colID44" ...
#>   ..- attr(*, "xtype")= num [1:12032] 0 1 2 2 2 2 2 2 2 2 ...
#>  $ newDE_idx    : int [1:900] 2010 3137 4736 3844 4771 2912 2476 2063 1115 2094 ...
```
