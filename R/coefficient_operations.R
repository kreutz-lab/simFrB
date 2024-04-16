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
#' data <- rnorm(nFeatures * nCols)
#' mtx <- matrix(data, nrow = nFeatures, ncol = nCols)
#' DE_idx <- sample(1:nFeatures, 0.2*nFeatures)
#' msb.omicsGetLmCoefs(mtx, DE_idx)

msb.omicsGetLmCoefs <- function(mtx, DE_idx = NULL,
                                groupDesign = rep(c(1, 2),
                                                  each = (dim(mtx)[2] / 2)),
                                chunksize = min(nrow(mtx), 400)) {

  assertthat::assert_that(is.matrix(mtx),
                          msg = "Expected a matrix for 'mtx',
                          got something else.")

  if (!is.null(DE_idx))
    assertthat::assert_that(all(DE_idx %in% 1:nrow(mtx)),
                            msg = "indices contained in DE_idx are outside
                           the range of the provided input matrix")
  else
    warning("no index of DE features (DE_idx) provided to omicsGetLmCoefs.
            No FC coefficient can be learned.")

  assertthat::assert_that(is.vector(groupDesign),
                          msg = "groupDesign needs to be a vector")


  assertthat::assert_that(length(groupDesign) == ncol(mtx),
                          msg = "groupDesign needs to have one entry for each
                          column of mtx,
                          specifing the group assigment of that sample")

  assertthat::assert_that(all(groupDesign %in% c(1,2)),
                          msg = "currently groupDesign can only handle two
                          groups, labelled with integers 1 and 2")

  predictors <- makePredictors(mtx = mtx, DE_idx = DE_idx, groupDesign = groupDesign)

  #run fits chunk by chunk
  fit <- fitLmModel(mtx = mtx, design = design, chunksize = chunksize, predictors = predictors)

  lmCoefs <- aggregateLmCoefs(coefs = stats::coef(fit), DE_idx = DE_idx)

  return(lmCoefs)
}

###### Helper functions

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
fitLmModel <- function(mtx, design, chunksize, predictors) {
  # Get chunks indices to apply linear model chunk by chunk
  chunks <- ceiling(seq(0, nrow(mtx),
                        length.out = ceiling(nrow(mtx) / chunksize)))
  if (length(chunks == 1))
    chunks <- c(0, nrow(mtx))
  for (n in 1:(length(chunks) - 1)) {
    idx1 <- chunks[n] + 1
    idx2 <- chunks[n + 1]
    message(paste("fitting from feature", idx1, "to", idx2))

    data.chunk <- prepareDataChunk(mtx, idx1, idx2, predictors)

    #The design can only include FC coefficient if DE proteins present in chunk
    if (2 %in% data.chunk$FC)
      design <- stats::as.formula(intensity ~ rowID + colID + rowID:FC - 1)
    else
      design <- stats::as.formula(intensity ~ rowID + colID - 1)

    if (n == 1) {
      fit <- biglm::biglm(design, data.chunk)
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
      stats::setNames(coefs[grep("FC", names(coefs))],
               paste0("feature_", DE_idx))
  )
}

