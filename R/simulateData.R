#' Simulate Data from Benchmark Data or Coefficients
#'
#' This function simulates proteomics data either based on a given count matrix
#' or by previously estimated coefficients.
#' It utilizes a linear model to estimate coefficients for the intensity values
#' and a logistic regression as implemented in DIMAR to estimate coefficients
#' for the amount of missing values in each feature.
#' The function is designed to estimate coefficients from benchmark data where
#' it is known which features are differentially expressed. If no benchmark data
#' is available, the function can still be used to simulate data based on coefficients
#' previously estimated on data presented in
#' Fröhlich, K., Brombacher, E., Fahrner, M. et al.
#' Benchmarking of analysis strategies for data-independent acquisition
#' proteomics using a large-scale dataset
#' comprising inter-patient heterogeneity.
#' Nat Commun 13, 2622 (2022). https://doi.org/10.1038/s41467-022-30094-0
#'
#' @param mtx Optional matrix used to derive LM and DIMAR coefficients if not provided.
#' @param DE_idx Optional indices for differentially expressed (DE) features.
#' @param jointCoefs list of coefficients from linear model and DIMAR. Either
#' estimated from input mtx or found in data/.
#' @param nFeatures Number of features to simulate in the output matrix.
#' @param nSamples Number of samples to simulate in the output matrix.
#' @param nDE Number of DE features, computed as a proportion of `nFeatures` by default.
#' @param npat Number of patterns for missing value simulations, defaults to 10.
#' @param groupDesign_mtx Group design for the input matrix `mtx` if it exists.
#' @param groupDesign_new Group design for the new simulated matrix. By default
#' first half of samples are assigned to group 1 and the second half to group 2.
#' @param int.mean The mean intensity value added to simulated data.
#' @param drawReplace Logical indicating whether sampling of coefficients is with replacement.
#'
#' @return An array of simulated data with dimensions nFeatures x nSamples x npat,
#'         with attributes detailing the used coefficients.
#'
#' @details If `lmCoefs` or `dimarCoefs` are not provided, and an input matrix `mtx` is available,
#'          the function will compute these coefficients based on `mtx` and `DE_idx`.
#'          Make sure to install and load the DIMAR package as it is crucial for the function's operation.
#'
#' @examples
#'
#' \dontrun{
#' url <- "https://raw.githubusercontent.com/kreutz-lab/dia-benchmarking/main/data/diaWorkflowResults_allDilutions.rds"
#' data <-
#' loadDataFromGit(url)
#' exp.df <- subsetDIAWorkflowData(data = data, DIAWorkflow = "DIANN_DIANN_AI",
#' experimentalComparisonGroups = c("1-12","1-25"),
#' rowSubset = seq(1,1000), colSubset = c(seq(1,5),seq(24,28)))
#'
#' exp.df <- exp.df[rowSums(exp.df, na.rm = TRUE) > 0, ]
#' DE_idx <- grep("ECOLI", rownames(exp.df))
#'   sim_data <- msb.simulateDataFromBenchmark(
#'    mtx = as.matrix(exp.df),
#'    DE_idx = DE_idx,
#'    nFeatures = 800,
#'    nSamples = 20,
#'    npat = 2)
#' }
#'
#'
#' @importFrom DIMAR dimarLearnPattern dimarAssignPattern
#' @export
msb.simulateDataFromBenchmark <- function(mtx = NULL,
                                          groupDesign_mtx = NULL,
                                          DE_idx = NULL,
                                          jointCoefs = NULL,
                                          int.mean = 19,
                                          nFeatures = 2000,
                                          nSamples = 10,
                                          nDE = NULL,
                                          npat = 10,
                                          groupDesign_new = rep(c(1, 2), each = nSamples / 2),
                                          drawReplace = T) {



  if (is.null(mtx)) {
    message("No input matrix provided")
    if (is.null(jointCoefs))
      stop("Coefficients need to be provided if no matrix is given,
           Coefficients calculated from Froehlich et. al 2022 Benchmark Data
           are available in data (pkg: simFrB).")
  } else {
    if (!"matrix" %in% class(mtx))
      mtx <- as.matrix(mtx)
      mtxAssertions(mtx)
      if (is.null(groupDesign_mtx)){
        groupDesign_mtx <- rep(c(1, 2), each = ncol(mtx) / 2)
        message("No groupDesign for mtx provided, assuming first half of samples
                are in group 1 and second half in group 2. If this is not the case
                please provide groupDesign_mtx.")
      }
      assertthat::assert_that(is.vector(groupDesign_mtx),
                              msg = "groupDesign needs to be a vector")
      assertthat::assert_that(all(groupDesign_mtx %in% c(1,2)),
                              msg = "currently groupDesign can only handle two
                          groups, labelled with integers 1 and 2")
  }


  assertthat::assert_that(all(groupDesign_new %in% c(1,2)),
                          msg = "currently groupDesign can only handle two
                          groups, labelled with integers 1 and 2")



  DEidxAssertions(DE_idx, mtx)

  assertthat::assert_that(nSamples %% 2 == 0,
                          msg = "nSamples needs to be an even number")
  assertthat::assert_that(length(groupDesign_new) == nSamples,
                          msg = "groupDesign_new needs to have one entry for each
                          sample, specifing the group assigment of that sample")

  assertthat::assert_that(nFeatures > 0,
                          msg = "nFeatures needs to be a positive number")

  if(is.null(nDE)){
    nDE <- round(nFeatures * 0.15)
    message("no nDE given, using default value of 15% of nFeatures")
  }
  assertthat::assert_that(nDE > 0,
                          msg = "nDE needs to be a positive number")

  assertthat::assert_that(is.numeric(int.mean) & length(int.mean) == 1,
                          msg = "int.mean needs to be a numeric value")



  #Estimate Linear Model coefficients from experimental data
  if (is.null(jointCoefs))
    jointCoefs <- estimateAndJoinCoefficients(mtx = mtx,
                                                     DE_idx = DE_idx,
                                                     groupDesign_mtx =
                                                       groupDesign_mtx)


  #draw randomly from dataframes for new data of custom size
  newCoefs <- drawCoefs(jointCoefficients = jointCoefs,
                        nFeatures = nFeatures,
                        nSamples = nSamples,
                        nDE = nDE,
                        drawReplace = drawReplace)
  #Use drawn linear model coefficients to simulate intensity values
  #you can also include custom values for mean int or add a value to
  #multiply to the coefficients
  full.mtx <- msb.predictFromLmCoefs(coefs = newCoefs$newLmCoefs,
                                     nFeatures = nFeatures,
                                     groupDesign = groupDesign_new,
                                     DE_idx = newCoefs$newDE_idx,
                                     int.mean = int.mean)
  #Use drawn logit coefficients to include missing values
  sim.df <- dimarAssignInChunks(ref = full.mtx,
                                coef = newCoefs$newDimarCoefs,
                                groupDesign = groupDesign_new,
                                npat = npat)


  attr(sim.df, 'full') <- full.mtx
  attr(sim.df,"usedCoefs") <- newCoefs
  attr(sim.df,"expCoefs") <- jointCoefs
  return(sim.df)
}

#' Estimate and Join LM and DIMAR Coefficients
#'
#' This internal function estimates linear model (LM) coefficients for intensity values
#' and Logit Coefficients (with DIMAR package) for the probability of missing values based on
#' experimental data. It then combines these coefficients into a joint dataframe,
#' used for simulation.
#'
#' @param mtx Matrix of experimental data used to compute coefficients.
#' @param DE_idx Numeric vector of indices for differentially expressed (DE) features.
#' @param groupDesign_mtx Group design vector explaining the group assignment of samples in `mtx`.
#'
#' @return A list containing the combined LM and DIMAR coefficients.
#'         The list includes separate elements for row coefficients (with DE labels)
#'         and column coefficients.
#'
#' @details The function performs two main tasks:
#'          - It calls `msb.omicsGetLmCoefs` to estimate LM coefficients based on
#'            provided matrix data (`mtx`), DE indices (`DE_idx`), and the group design (`groupDesign_mtx`).
#'          - It estimates DIMAR coefficients using `dimarLearnPattern` for modeling
#'            the chance of missing data, also based on `mtx` and `DE_idx`.
#'          The estimated coefficients are then combined using `joinCoefs`.
#'
#' @noRd
#'
estimateAndJoinCoefficients <- function(mtx, DE_idx, groupDesign_mtx){

  #estimate linear model coefficients for intensity values
  #from experimental data
  lmCoefs <- msb.omicsGetLmCoefs(mtx = mtx, DE_idx = DE_idx,
                                 groupDesign = groupDesign_mtx)

  #Estimate Logit Coefficients (DIMAR) for chance of missing value from
  #experimental data
  dimarCoefs <- DIMAR::dimarLearnPattern(mtx = mtx, DE_idx = DE_idx,
                                         group = groupDesign_mtx,
                                  orderCoefByName = T)

  #Create a joint dataframe matching coefficients obtained for intensity values
  #and missing values. One data frame for row coefficients with DE label, one
  #for column coefficients
  jointCoefficients <- joinCoefs(lmCoefs = lmCoefs,
                               dimarCoefs = dimarCoefs)
  return(jointCoefficients)
}

#' Simulate Matrix from Linear Model Coefficients
#'
#' This function simulates an intensity matrix based on specified linear model coefficients
#' for a given experimental design. It optionally adds noise to the simulated values.
#'
#' @param groupDesign A vector indicating the group design for each sample.
#' @param DE_idx Optional vector of indices for differentially expressed features.
#' @param nFeatures Integer specifying the number of features (rows) in the simulated matrix.
#' @param percentDE Proportion of features that are differentially expressed.
#' @param coefs List of coefficients including `featureCoefs`, `sampleCoefs`, and `FCCoefs`.
#' @param int.mean Mean intensity value to be added to each simulated data point.
#' @param addNoise Logical indicating whether to add simulated noise to the matrix.
#' @param rowScale Scaling factor for row coefficients.
#' @param rowSdScale Scaling factor for the standard deviation of rows.
#' @param colScale Scaling factor for column coefficients.
#' @param FCScale Scaling factor for fold change coefficients.
#' @return A matrix representing simulated intensity values.
#' @export
#' @examples
#'   groupDesign <- rep(1:2, each = 5)
#'   DE_idx <- sample(1:500, 100)
#'   coefs <- list(featureCoefs = stats::rnorm(500),
#'   sampleCoefs = stats::rnorm(10), FCCoefs = stats::rnorm(100))
#'   sim_matrix <- msb.predictFromLmCoefs(groupDesign, DE_idx = DE_idx, coefs = coefs)
msb.predictFromLmCoefs <- function(groupDesign, DE_idx = NULL,
                                   nFeatures = 500,
                                   percentDE = 0.2,
                                   coefs, int.mean = 19,
                                   addNoise = T, rowScale = 1,
                                   rowSdScale = 1, colScale = 1, FCScale = 1){



  assertthat::assert_that(is.list(coefs),
                          msg = "coefs needs to be a list")

  rowNames <- paste0(rep("feature_", nFeatures),1:nFeatures)
  #include DE label in rownames
  rowNames[DE_idx] <- paste0(rowNames[DE_idx], "_DE")

  nSamples <- length(groupDesign)
  colNames <- paste(paste0("group", groupDesign),
                    1:length(groupDesign), sep = "_")

  full.mtx <- mtxSimulate(nFeatures, nSamples, coefs, int.mean, groupDesign,
                          rowScale, colScale, FCScale, rowNames, colNames)


  if (addNoise){
    #draw noise from normal distribution with
    #previously calculated standard deviations
    #optionally scale the noise by rowSdScale
    for (rowID in 1:nrow(full.mtx)){
      full.mtx[rowID,] <- stats::rnorm(ncol(full.mtx),
                                       mean = full.mtx[rowID,],
                                       sd = coefs$featureSd[rowID] * rowSdScale)
    }
  }
  return(full.mtx)
}

#' Internal Helper Function to Simulate Matrix
#'
#' This function constructs a matrix based on the provided coefficients and scales.
#' It is used internally by `msb.predictFromLmCoefs` to generate the base matrix before
#' noise addition.
#'
#' @param nFeatures Integer, number of features (rows) in the matrix.
#' @param nSamples Integer, number of samples (columns) in the matrix.
#' @param coefs List of coefficients including `featureCoefs`, `sampleCoefs`, and `FCCoefs`.
#' @param int.mean Mean intensity value to be added to each simulated data point.
#' @param rowScale Scaling factor for row coefficients.
#' @param colScale Scaling factor for column coefficients.
#' @param FCScale Scaling factor for fold change coefficients.
#' @param rowNames Vector of row names for the matrix.
#' @param colNames Vector of column names for the matrix.
#' @return A matrix with dimensions [nFeatures x nSamples] filled with simulated values.
#' @noRd
mtxSimulate <- function(nFeatures, nSamples, coefs, int.mean, groupDesign,
                        rowScale, colScale, FCScale, rowNames, colNames){
  #initialize matrix
  full.mtx <- matrix(data = NA, nrow = nFeatures, ncol = nSamples,
                     dimnames = list(rowNames, colNames))
  for (rowID in 1:nFeatures) {
    for (colID in 1:nSamples) {
      full.mtx[rowID,colID] <- coefs$featureCoefs[rowID] * rowScale +
        coefs$sampleCoefs[colID] * colScale + int.mean
      #for samples in group 2, add fold change
      if (groupDesign[colID] == 2) {
        full.mtx[rowID,colID] <- coefs$featureCoefs[rowID] +
          coefs$sampleCoefs[colID] + coefs$FCCoefs[rowID] + int.mean
      }
    }
  }
  return(full.mtx)
}

#' Assign DIMAR Coefficients in Chunks
#'
#' This function applies DIMAR coefficients to a reference matrix (`ref`) in chunks,
#' optionally using an existing matrix (`mtx`) for additional computations. It is designed
#' for large datasets where applying coefficients in chunks helps manage memory usage.
#' The function processes data by group as defined in the coefficients.
#'
#' @param ref A matrix or data frame used as the reference for assigning coefficients.
#' @param coef A named vector of coefficients including both row and other coefficients.
#'             Row coefficients should have names marked with "#" to denote groups.
#' @param chunksize The number of rows to process per chunk; defaults to the minimum of
#'                  400 or the number of rows in `ref`.
#' @param npat The number of patterns to simulate, defaults to 1.
#'
#' @return A data frame combining all chunks with the DIMAR pattern assignment applied.
#'         Each row in the returned data frame corresponds to a row in the `ref` matrix,
#'         with additional dimensions if `npat` > 1.
#'
#' @details The function segregates coefficients into groups based on their names, applies
#'          these coefficients to chunks of the reference matrix. It outputs a single data frame
#'          that aggregates the results from all chunks.
#' @noRd
dimarAssignInChunks <- function(ref,
                                coef,
                                chunksize = min(nrow(ref), 300),
                                groupDesign,
                                npat = 1){

  dimarRow_idx <- which(attributes(coef)$xtype == 3)
  rowCoefs <- coef[dimarRow_idx]
  otherCoefs <- coef[-dimarRow_idx]

  #seperate rowCoefs by group
  rowCoefs.lst <- list()
  rowGroupVector <- gsub(".*#", "", names(rowCoefs))
  for (g in unique(rowGroupVector)) {
    rowCoefs.group <- rowCoefs[rowGroupVector == g]
    rowCoefs.lst[[paste0("dimarCoefs_group",g)]] <- rowCoefs.group
  }
  # Get chunks indices to apply linear model chunk by chunk
  chunks <- ceiling(seq(0, nrow(ref),
                        length.out = ceiling(nrow(ref) / chunksize)))
  if (length(chunks) == 1){
    chunks <- c(0, nrow(ref))}

  sim <- list()

  for (n in 1:(length(chunks) - 1)) {
    idx1 <- chunks[n] + 1
    idx2 <- chunks[n + 1]
    message(paste("applying dimar coefficients from", idx1, "to", idx2))
    ref.chunk <- ref[idx1:idx2, ]
    rowCoefs.chunk <- c()
    for (g in unique(rowGroupVector)) {
      rowCoefs.chunk.group <- rowCoefs.lst[[paste0("dimarCoefs_group",g)]]
      rowCoefs.chunk <- c(rowCoefs.chunk, rowCoefs.chunk.group[idx1:idx2])
    }
    coefs.chunk <- c(otherCoefs, rowCoefs.chunk)


    sim.chunk <- suppressMessages(DIMAR::dimarAssignPattern(ref = ref.chunk,
                                           coef = coefs.chunk,
                                           group = groupDesign,
                                           npat = npat))

    sim.chunk <- array(sim.chunk,c(dim(ref.chunk),npat), dimnames = list(rownames(ref.chunk),
                                                                         colnames(ref.chunk),
                                                                         1:npat))

    sim[[n]] <- sim.chunk
  }
  sim.array <- abind::abind(sim, along = 1)
  return(sim.array)
}
