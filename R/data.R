#' Coefficients estimated with linear model and logistic regression from
#' experimental data submitted in Fröhlich et. al (2022)
#' Fröhlich, K., Brombacher, E., Fahrner, M. et al. 
#' Benchmarking of analysis strategies for data-independent acquisition 
#' proteomics using a large-scale dataset comprising 
#' inter-patient heterogeneity. Nat Commun 13, 2622 (2022). 
#' https://doi.org/10.1038/s41467-022-30094-0
#' 
#' Specifically these coefficients were estimated from data processed with
#' DIANN_DIANN_AI_GPF and comparison c("1-12","1-25"). They were estimated
#' with the simFrB::msb.simulateFromBenchmark() function.
#' 
#'
#' @format ## `jointCoefs`
#' A list with 3 elements:
#' \describe{
#'   \item{rowCoefs}{all feature coefficients}
#'   \item{colCoefs}{all column coefficients}
#'   \item{dimarIntercept}{Intercept of logistic regression (from DIMAR pkg)}
#' }
"jointCoefs"