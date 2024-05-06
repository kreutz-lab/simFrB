#' Extract Linear Model Coefficients for Omics Data
#'
#' This function extracts linear model coefficients from a provided matrix
#' representing omics data, such as proteomics or genomics,
#' by fitting linear models
#' in chunks.
#'
#'
#' @param mtx A matrix of omics data where rows represent features and
#' columns represent samples.
#' @param DE_idx A vector indicating the indices of differentially
#' expressed features.
#' @param groupDesign A vector indicating group assignments
#' for each column in the matrix.
#' @param chunksize The number of rows to process in each chunk.
#' Defaults to the minimum of 400 or the number of rows in the matrix.
#' @return A list containing the coefficients of the linear models.
#' @export
#' @examples
#' nFeatures <- 100
#' nCols <- 10
#' data <- stats::rnorm(nFeatures * nCols)
#' mtx <- matrix(data, nrow = nFeatures, ncol = nCols)
#' DE_idx <- sample(1:nFeatures, 0.2*nFeatures)
#' msb.omicsGetLmCoefs(mtx, DE_idx)

msb.omicsGetLmCoefs <- function(mtx = NULL, DE_idx = NULL,
                                groupDesign = rep(c(1, 2),
                                                  each = (ncol(mtx) / 2)),
                                chunksize = min(nrow(mtx), 400)) {

  mtxAssertions(mtx = mtx)
  DEidxAssertions(mtx = mtx, DE_idx = DE_idx)
  groupDesignAssertions(groupDesign = groupDesign, mtx = mtx)


  if (is.null(DE_idx))
    warning("no index of DE features (DE_idx) provided to omicsGetLmCoefs.
            No FC coefficient can be learned.")


  predictors <- makePredictors(mtx = mtx, DE_idx = DE_idx, groupDesign = groupDesign)

  #run fits chunk by chunk
  fit <- fitLmModel(mtx = mtx, chunksize = chunksize, predictors = predictors)
  coefs <- stats::coef(fit)

  assertthat::assert_that(length(coefs) == (nrow(predictors$rowIDs.mtx) + ncol(predictors$colIDs.mtx) + nrow(predictors$FC.mtx)- 1),
                        msg = "number of coefficients does not match number of predictors")

  coefs <- aggregateLmCoefs(coefs = coefs, DE_idx = DE_idx)

  lmCoefs <- addSdToCoefs(mtx = mtx,
                          coefs = coefs,
                          nFeatures = nrow(mtx),
                          groupDesign = groupDesign,
                          DE_idx = DE_idx)

  return(lmCoefs)
}
#' Calculate Feature Standard Deviations and Add to Coefficients
#'
#' This function estimates the standard deviations of features from a given
#' matrix and adds them to a provided list of coefficients.
#' If no matrix is provided, it uses
#' default values to generate standard deviations.
#'
#' @param mtx Optional matrix from which to calculate standard deviations.
#'            If NULL, standard deviations are generated from default parameters.
#' @param coefs List containing coefficient data,
#' to which standard deviations will be added.
#' @param DE_idx Optional vector of indices for
#' differentially expressed features.
#' If NULL, a warning is issued and calculations may be less accurate.
#' @param groupDesign Optional vector specifying the group design for the matrix columns.
#'                    If NULL, it assumes two groups of equal size.
#' @param sd.mean.default Mean value used for generating standard deviations
#' if no matrix is provided.
#' @param sd.sd.default Standard deviation value used for
#' generating standard deviations if no matrix is provided.
#'
#' @return The `coefs` list, augmented with `featureSd`,
#' which contains the calculated standard deviations.
#' @details If `mtx` is not provided, standard deviations are generated using a normal distribution with
#'          specified mean and standard deviation. The function first tries to calculate standard deviations
#'          from the residuals of a simulated full matrix (calculated using provided coefficients) subtracted
#'          from `mtx`. If the matrix `mtx` is not available, it falls back to using default values.
#' @examples
#'   mtx <- matrix(stats::rnorm(100), nrow = 10, ncol = 10)
#'   coefs <- list(featureCoefs = runif(10), sampleCoefs = runif(10))
#'   DE_idx <- c(1, 5, 10)
#'   updated_coefs <- addSdToCoefs(mtx, coefs, DE_idx)
#' @noRd

###### Helper functions

addSdToCoefs <- function(mtx = NULL, coefs, DE_idx = NULL, nFeatures,
                         groupDesign = NULL, sd.mean.default = 0.6,
                         sd.sd.default = 0.4){

  if (is.null(mtx)){
    message("no input data provide to estimate protein
            standard deviations from, using default values.")
    coefs$featureSd <- abs(stats::rnorm(nFeatures, mean = sd.mean.default,
                                 sd = sd.sd.default))
    names(coefs$featureSd) <- paste0("rowID",nFeatures)
    return(coefs)
  }
  else
    message("using input data to estimate feature standard deviations.")

  assertthat::assert_that(is.list(coefs),
                          msg = "coefs needs to be a list")
  if (is.null(groupDesign)) {
    message("no group design provided, assuming group 1 and 2 of equal size")
    groupDesign <- rep(c(1,2), each = ncol(mtx)/2)
  }
  if (is.null(DE_idx))
    warning("no DE_idx provided, calculation of sd will be inacurate")


  #simulate full matrix as close to exp as possible. No noise included.
  full.exact <- msb.predictFromLmCoefs(nFeatures = nrow(mtx),
                                       groupDesign = groupDesign,
                                       int.mean = mean(mtx,na.rm = T),
                                       DE_idx = DE_idx,
                                       coefs = coefs,
                                       addNoise = F)

  #get residuals by subtracting from original matrix
  res.matrix <- mtx - full.exact

  # Compute standard deviations, handling NAs and infinite values
  featureSd <- apply(res.matrix, 1, function(x) {
    valid_data <- x[!is.na(x)]  # Exclude NA values for computation
    if (length(valid_data) > 2) {
      sqrt(sum(valid_data^2) / (length(valid_data) - 2))
    } else {
      NA  # Return NA if not enough data to compute SD
    }
  })

  # Replace NA and infinite values with the mean of the absolute residuals
  # Calculate mean of absolute residuals ignoring NA and infinite
  mean_abs_res <- mean(abs(res.matrix[is.finite(res.matrix)]), na.rm = TRUE)
  # Replace NA and Inf in protein.sd with mean_abs_res
  featureSd[is.na(featureSd) | is.infinite(featureSd)] <- mean_abs_res

  names(featureSd) <- paste0("rowID",1:nFeatures)
  names(featureSd)[DE_idx] <- paste0(names(featureSd)[DE_idx],"_DE")

  coefs$featureSd <- featureSd
  return(coefs)
}

#' Generate predictors for the regression model
#'
#' This internal function generates predictors for regression modeling from the matrix,
#' considering differentially expressed indices.
#'
#' @param mtx The data matrix.
#' @param DE_idx Indices of differentially expressed features.
#' @return A dataframe of predictors.
#' @noRd
makePredictors <- function(mtx, DE_idx, groupDesign){

  rowIDs <- rep(1:dim(mtx)[1])
  colIDs <- rep(1:dim(mtx)[2])

  rowIDs.mtx <- matrix(rep(rowIDs, length(colIDs)),
                       nrow = length(rowIDs),
                       ncol = length(colIDs))
  colIDs.mtx <- matrix(rep(colIDs, each = length(rowIDs)),
                       nrow = length(rowIDs), ncol = length(colIDs))
  FC.mtx <- matrix(data = rep(as.factor(groupDesign), each =  length(rowIDs)),
                   nrow = length(rowIDs), ncol = length(colIDs))

  #Only proteins of group 2 which are differntially expressed are also
  #assigned to group 2 here. Others are set to 1.
  FC.mtx[!rowIDs %in% DE_idx,] <- 1

  #combine in list
  predictors <- list(rowIDs.mtx = rowIDs.mtx,
                     colIDs.mtx = colIDs.mtx,
                     FC.mtx = FC.mtx)
  return(predictors)
}

#' Prepare a data chunk for regression
#'
#' This internal function prepares a subset of the data matrix for fitting
#' a regression model, adjusting based on the provided indices.
#'
#' @param mtx The full data matrix.
#' @param idx1 The starting index of the chunk.
#' @param idx2 The ending index of the chunk.
#' @param predictors A list containing the predictors (rowIDs, colIDs, FC).
#' @return A dataframe ready for regression fitting.
#' @noRd
prepareDataChunk <- function(mtx, idx1, idx2, predictors) {
  assertthat::assert_that(!(is.na(idx1) | is.na(idx2)),
                          msg = "idx1 or idx2 are NA")

  #assertthat::assert_that(idx2 > idx1,
  #                        msag = "idx2 is not larger than idx1")
  data.chunk <- data.frame(
    intensity = as.vector(mtx[idx1:idx2, ] - mean(mtx, na.rm = TRUE))
  )
  data.chunk$rowID <- factor(predictors$rowIDs.mtx[idx1:idx2,],
                             levels = levels(as.factor(predictors$rowIDs.mtx)))
  data.chunk$colID <- factor(predictors$colIDs.mtx[idx1:idx2,],
                             levels = levels(as.factor(predictors$colIDs.mtx)))
  data.chunk$FC <- factor(predictors$FC.mtx[idx1:idx2,],
                          levels = levels(as.factor(predictors$FC.mtx)))
  return(data.chunk)
}


#' Fit the linear model in chunks
#'
#' This internal function applies a linear regression model to chunks of the data matrix.
#'
#' @param mtx The data matrix.
#' @param chunksize The size of each chunk.
#' @param predictors Predictors prepared for the regression model.
#' @return A linear model object.
#' @noRd
fitLmModel <- function(mtx, chunksize, predictors) {
  # Get chunks indices to apply linear model chunk by chunk
  chunks <- ceiling(seq(0, nrow(mtx),
                        length.out = ceiling(nrow(mtx) / chunksize)))
  if (length(chunks) == 1)
    chunks <- c(0, nrow(mtx))
  for (n in 1:(length(chunks) - 1)) {
    idx1 <- chunks[n] + 1
    idx2 <- chunks[n + 1]
    message(paste("fitting from feature", idx1, "to", idx2))

    data.chunk <- prepareDataChunk(mtx, idx1, idx2, predictors)

    if (n == 1) {
      fit <- biglm::biglm(
        formula = stats::formula(intensity ~ rowID + colID + rowID:FC - 1),
        data = data.chunk)
    } else {
      fit <- stats::update(fit, moredata = data.chunk)
    }

  }
  return(fit)
}

#' Aggregate and organize linear model coefficients
#'
#' This internal function aggregates coefficients from a fitted linear model and organizes them.
#'
#' @param coefs Coefficients from the fitted model.
#' @param DE_idx Indices of differentially expressed features.
#' @return A list of structured coefficients.
#' @noRd
aggregateLmCoefs <- function(coefs, DE_idx){
  lmCoefs <- list(
    featureCoefs =
      coefs[grepl("rowID", names(coefs)) &
                  !grepl("FC", names(coefs))],
    sampleCoefs =
      c(colID1 = 0, coefs[grep("colID", names(coefs))]),
    FCCoefs =
     coefs[grep("FC", names(coefs))]
  )
  # set NA coefficients to 0
  lmCoefs$FCCoefs[is.na(lmCoefs$FCCoefs)] <- 0
  names(lmCoefs$FCCoefs) <- gsub(":FC2","",names(lmCoefs$FCCoefs))
  #add DE_idx to names
  names(lmCoefs$FCCoefs)[DE_idx] <- paste0(names(lmCoefs$FCCoefs)[DE_idx],"_DE")
  names(lmCoefs$featureCoefs)[DE_idx] <- paste0(names(lmCoefs$featureCoefs)[DE_idx],"_DE")
  return(lmCoefs)
}

#' Join LM and DIMAR Coefficients
#'
#' This function combines linear model coefficients and DIMAR coefficients into a single
#' structured list, categorizing and labeling DE and non-DE features based on specified tags and keywords.
#'
#' @param lmCoefs List of linear model coefficients including `featureCoefs`, 'sampleCoefs', FCCoefs, `featureSd`.
#' @param dimarCoefs List of DIMAR coefficients.
#' @return A list containing combined and structured coefficients.
#' @examples
#' # Assuming lmCoefs and dimarCoefs are already defined
#' joined_coefs <- joinCoefs(lmCoefs, dimarCoefs, c("1-25", "1-12"), "ECOLI")
#' @noRd
joinCoefs <- function(lmCoefs,
                      dimarCoefs) {

  assertthat::assert_that(!is.null(attributes(dimarCoefs)$xtype),
                          msg = "dimarCoefs must have xtype attribute.")
  #Join column coefficients
  dimarCol_idx <- which(attributes(dimarCoefs)$xtype == 2)
  dimarColCoef <- dimarCoefs[dimarCol_idx]
  jointColCoefs <- data.frame(lmCoefs = lmCoefs$sampleCoefs,
                              dimarCoefs = dimarColCoef)

  jointRowCoefs <- joinRowCoefs(lmCoefs = lmCoefs,
                                dimarCoefs = dimarCoefs,
                                nSamples = length(dimarColCoef))

  #combine everything
  jointCoefs <- list(rowCoefs = jointRowCoefs,
                     colCoefs = jointColCoefs,
                     dimarIntercept = dimarCoefs[1:2])
  return(jointCoefs)
}



#' Join Row Coefficients
#'
#' Internal function to combine LM coefficients with DIMAR row coefficients.
#'
#' @param lmCoefs List containing LM coefficients.
#' @param dimarCoefs List containing DIMAR coefficients.
#' @param nSamples number of samples in dimarCoefs
#' @return Data frame of combined row coefficients.
#' @noRd
joinRowCoefs <- function(lmCoefs, dimarCoefs, nSamples) {
  #initialise data frame
  jointRowCoefs <- data.frame(rowIDs = names(lmCoefs$featureCoefs),
                              lmCoefs = lmCoefs$featureCoefs,
                              lmCoefs_sd = lmCoefs$featureSd)

  #DE idx is computed during dimar row coef extraction
  #and then used to label DE proteins
  jointRowCoefs[["ifDE"]] <- "notDE"


  #join FC coefficients to row Coefficients
  FCCoefs.df <- data.frame(rowIDs = names(lmCoefs$FCCoefs),
                           FCCoefs = lmCoefs$FCCoefs, row.names = NULL)
  FCCoefs.df$rowIDs <- sub(":FC*[0-9]","",FCCoefs.df$rowIDs)
  jointRowCoefs <- dplyr::inner_join(jointRowCoefs, FCCoefs.df)
  #set NA values to 0 to enable adding to feature coefs
  #jointRowCoefs$FCCoefs[is.na(jointRowCoefs$FCCoefs)] <- 0

  DE_idx_dimar <- attributes(dimarCoefs)$DE_idx
  DE_idx <- grep("DE",names(lmCoefs$FCCoefs))
  assertthat::are_equal(DE_idx, DE_idx_dimar,
                        msg = "DE_idx from lmCoefs and dimarCoefs do not match!")

  #include DIMAR coefs
  dimarRowCoef <- dimarCoefs[seq(from = nSamples + 3,
                                 to = length(dimarCoefs))]

  rowGroupVector <- gsub(".*#", "", names(dimarRowCoef))
  for (g in unique(rowGroupVector)) {
    dimarRowCoef.group <- dimarRowCoef[rowGroupVector == g]
    jointRowCoefs[[paste0("dimarCoefs_group",g)]] <- dimarRowCoef.group
  }
  jointRowCoefs$ifDE[DE_idx] <- "DE"

  attr(jointRowCoefs, "dimarXtype") <- attributes(dimarCoefs)$xtype
  return(jointRowCoefs)

}

#' Draw Coefficients
#'
#' This function handles the sampling of coefficients from provided joint coefficients,
#' including handling of DE and non-DE differentiations.
#'
#' @param jointCoefficients List containing both LM and DIMAR coefficients.
#' @param nFeatures Total number of features to include.
#' @param nSamples Total number of samples to include.
#' @param nDE Number of differentially expressed features.
#' @param drawReplace Logical indicating whether sampling should be with replacement.
#' @return A list containing new LM coefficients and DIMAR coefficients.
#' @noRd
drawCoefs <- function(jointCoefficients,
                      nFeatures, nSamples, nDE,
                      drawReplace = TRUE) {

  # Draw row coefficients
  rowResults <- drawRowCoefficients(jointCoefficients, nFeatures, nDE, drawReplace)
  newRowCoefs <- rowResults$newRowCoefs
  newDE_idx <- rowResults$newDE_idx

  # Draw column coefficients
  newColCoefs <- drawColumnCoefficients(jointCoefficients, nSamples, drawReplace)

  ##combine lmCoefs
  newLmCoefs <- list(featureCoefs = stats::setNames(newRowCoefs[, "lmCoefs"],
                                             newRowCoefs$rowIDs),
                     sampleCoefs = stats::setNames(newColCoefs[, "lmCoefs"],
                                            rownames(newColCoefs)),
                     FCCoefs = stats::setNames(newRowCoefs[,"FCCoefs"], newRowCoefs$rowIDs),
                     featureSd =  stats::setNames(newRowCoefs[, "lmCoefs_sd"],
                                           newRowCoefs$rowIDs))

  ## combine dimar coefs
  nGroups <- length(newRowCoefs) - 5
  newDimarCoefs <- combineNewDimarCoefs(jointCoefficients, newRowCoefs, newColCoefs, nGroups = nGroups)
  newXtype <- c(0,1,rep(2,nSamples),rep(3,nFeatures*nGroups))
  attr(newDimarCoefs, "xtype") <- newXtype

  # Return final results
  return(list(newLmCoefs = newLmCoefs, newDimarCoefs = newDimarCoefs, newDE_idx = rowResults$newDE_idx))
}

#' Combine DIMAR Coefficients
#'
#' This function integrates DIMAR intercept coefficients with newly sampled
#' column and row coefficients. It handles the integration across multiple
#' specified groups, adding appropriate naming conventions to distinguish
#' between differentially expressed (DE) and non-DE features within these groups.
#'
#' @param jointCoefficients A list containing the original joint coefficients
#'        from which the DIMAR intercepts are extracted.
#' @param newRowCoefs A data frame of new row coefficients that have been sampled.
#'        It must include `rowIDs` and
#'        a column that specifies if a row is DE or not (`ifDE`).
#' @param newColCoefs A data frame of new column coefficients that have been sampled
#'        from the original set. This should include a column named `dimarCoefs` which
#'        will be used in the combination process.
#' @param nGroups Integer specifying the number of groups to process. This is used to
#'        determine how many times to loop over the group-specific coefficients in
#'        `newRowCoefs` for integration.
#'
#' @return A vector of combined DIMAR coefficients, starting with the intercepts,
#'         followed by the integrated column and row coefficients for each group.
#'
#' @noRd
combineNewDimarCoefs <- function(jointCoefficients, newRowCoefs, newColCoefs, nGroups){
  #first intercept and mean (unchanged)
  newDimarCoefs <- jointCoefficients$dimarIntercept
  #then add column coefficients
  newDimarCoefs <- c(newDimarCoefs,
                     stats::setNames(newColCoefs[, "dimarCoefs"],
                              rownames(newColCoefs)))
  #loop over groups to add row coefficients
  for (g in 1:nGroups) {
    dimarCoefs.group <- newRowCoefs[,ncol(newRowCoefs) - (g - 1)]
    names(dimarCoefs.group) <- newRowCoefs$rowIDs
    groupDE_idx <- which(newRowCoefs$ifDE == "DE")
    #names(dimarCoefs.group)[groupDE_idx] <-
    #  paste0(names(dimarCoefs.group)[groupDE_idx],"_","DE")
    names(dimarCoefs.group) <- paste0(names(dimarCoefs.group),"#",g)
    newDimarCoefs <- c(newDimarCoefs, dimarCoefs.group)
  }
  return(newDimarCoefs)
}

#' Draw Row Coefficients
#'
#' Sample row coefficients for DE and non-DE features from joint coefficients.
#'
#' @param jointCoefficients List containing joint coefficients.
#' @param nFeatures Integer, total number of features to sample.
#' @param nDE Integer, number of DE features.
#' @param drawReplace Logical, whether to sample with replacement.
#' @return A data frame of sampled row coefficients.
#' @noRd
drawRowCoefficients <- function(jointCoefficients, nFeatures, nDE, drawReplace) {
  newRowCoefs <- as.data.frame(matrix(NA,
                                      nrow = nFeatures,
                                      ncol = length(jointCoefficients$rowCoefs),
                                      ))
  colnames(newRowCoefs) <- colnames(jointCoefficients$rowCoefs)

  newDE_idx <- sample(1:nFeatures, nDE)
  newnotDE_idx <- which(!(1:nrow(newRowCoefs) %in% newDE_idx))

  newRowCoefs[newnotDE_idx,] <- dplyr::filter(jointCoefficients$rowCoefs, jointCoefficients$rowCoefs$ifDE == "notDE") %>%
    dplyr::sample_n(nFeatures - nDE, replace = drawReplace)
  newRowCoefs[newDE_idx,] <- dplyr::filter(jointCoefficients$rowCoefs, jointCoefficients$rowCoefs$ifDE == "DE") %>%
    dplyr::sample_n(nDE, replace = drawReplace)

  return(list(newRowCoefs = newRowCoefs, newDE_idx = newDE_idx))
}

#' Draw Column Coefficients
#'
#' Sample column coefficients from joint coefficients.
#'
#' @param jointCoefficients List containing joint coefficients.
#' @param nSamples Integer, number of samples to draw.
#' @param drawReplace Logical, whether to sample with replacement.
#' @return A vector of sampled column coefficients.
#' @noRd
drawColumnCoefficients <- function(jointCoefficients, nSamples, drawReplace) {
  newColCoefs <- dplyr::sample_n(jointCoefficients$colCoefs, nSamples, replace = drawReplace)
  return(newColCoefs)
}
