#test_that("msb.omicsGetLmCoefs return proper error messages", {
#  expect_error(msb.omicsGetLmCoefs(mtx = data.frame(1:5, 5:9)), "mtx needs to be a matrix")
#})


test_that("Error is thrown for 'DE_idx' out of bounds", {
  correct_matrix_input <- matrix(rnorm(100), 10, 10)
  out_of_bounds_DE_idx <- c(1, 2, 11)  # 11 is out of bounds for 10 features

  expect_error(
    msb.omicsGetLmCoefs(mtx = correct_matrix_input,
                        DE_idx = out_of_bounds_DE_idx,
                        groupDesign = rep(1:2, 5)),
    regexp = "indices contained in DE_idx are outside the range of the provided input matrix"
  )
})

test_that("Error is thrown for not supported groupDesign", {
  correct_matrix_input <- matrix(rnorm(100), 10, 10)
  incorrect_group_design <- rep(c(2,3), each = 5)
  correct_DEidx <- sample(1:nrow(correct_matrix_input),
                          0.2*nrow(correct_matrix_input))
  expect_error(
    msb.omicsGetLmCoefs(correct_matrix_input,
                        DE_idx = correct_DEidx,
                        groupDesign = incorrect_group_design),
    "currently groupDesign can only handle two
                          groups, labelled with integers 1 and 2"
  )
})


# Test for incorrect groupDesign length
test_that("Error is thrown for incorrect 'groupDesign' length", {
  correct_matrix_input <- matrix(rnorm(100), 10, 10)
  incorrect_group_design <- rep(1:2, 4)  # Only 8 elements, needs 10
  correct_DEidx <- sample(1:nrow(correct_matrix_input), 0.2*nrow(correct_matrix_input))
  expect_error(
    msb.omicsGetLmCoefs(correct_matrix_input, DE_idx = correct_DEidx, groupDesign = incorrect_group_design),
    "groupDesign needs to have one entry for each
                          column of mtx,
                          specifing the group assigment of that sample"
  )
})

test_that("makePredictors outputs correct predictors", {
  # Create a small example matrix
  mtx <- matrix(1:100, nrow = 10, ncol = 10)

  # Define DE indices and group design
  DE_idx <- c(1, 5, 10)  # Assume features 1, 5, and 10 are differentially expressed
  groupDesign <- c(rep(1, 5), rep(2, 5))  # Simple group design for illustration

  # Call the function
  predictors <- makePredictors(mtx, DE_idx, groupDesign)

  # Check the dimensions and content of the predictors
  expect_length(predictors, 3)  # Expect 100 rows
  #dimension of each matrix in predictors should be same as of mtx
  predictor_dims <- lapply(predictors, dim)
  expect_true(all(sapply(predictor_dims, function(x) all(x == dim(mtx)))))

  #FC predictors for samples of group 1 should always be 1
  expect_true(all(predictors$FC[predictors$colIDs.mtx %in%
                                  seq(1,ncol(mtx))[groupDesign == 1]] == 1))
  #FC predictors for proteins not DE should always be 1
  expect_true(all(predictors$FC[!predictors$rowIDs %in% DE_idx] == 1))
})

test_that("msb.omicsGetLmCoefs runs as expected chunk does not contain any DE", {
  # Create a small example matrix
  nFeatures = 250
  nSamples = 10
  mtx <- matrix(rnorm(20), nrow = nFeatures, ncol = nSamples)
  DE_idx <- sample(1:100,20)

  # Define group design
  groupDesign <- c(rep(1, 5), rep(2, 5))  # Simple group design for illustration

  # Call the function
  lmCoefs <- msb.omicsGetLmCoefs(mtx, DE_idx = DE_idx, groupDesign, chunksize = 100)

  expect_length(lmCoefs, 4)  # Expect 4 elements in the output list
  expect_length(lmCoefs$featureCoefs, nFeatures)
  expect_length(lmCoefs$sampleCoefs, nSamples)
  expect_length(lmCoefs$FCCoefs, nFeatures)
  expect_length(lmCoefs$featureSd, nFeatures)

})
