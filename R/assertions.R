' Assert Validity of Group Design
#'
#' Ensures that the group design vector is valid according to specified conditions,
#' including length and permissible values.
#'
#' @param groupDesign The group design vector to validate.
#' @param mtx The matrix whose columns should match the length of `groupDesign`.
#' @return Invisibly returns TRUE if all assertions pass.
#' @noRd
groupDesignAssertions <- function(groupDesign, mtx){
  assertthat::assert_that(is.vector(groupDesign),
                          msg = "groupDesign needs to be a vector")


  assertthat::assert_that(length(groupDesign) == ncol(mtx),
                          msg = "groupDesign needs to have one entry for each
                          column of mtx,
                          specifing the group assigment of that sample")

  assertthat::assert_that(all(groupDesign %in% c(1,2)),
                          msg = "currently groupDesign can only handle two
                          groups, labelled with integers 1 and 2")
}
#' Assert Validity of Differential Expression Indices
#'
#' Checks validity of DE indices against a matrix. Ensures indices are numeric,
#' within the correct range, and properly formatted.
#'
#' @param DE_idx Vector of indices for differentially expressed features or NULL.
#' @param mtx The matrix to check DE indices against.
#' @return Invisibly returns TRUE if all assertions pass.
#' @noRd
DEidxAssertions <- function(DE_idx, mtx){
  assertthat::assert_that(is.null(DE_idx) | is.vector(DE_idx),
                          msg = "DE_idx needs to be a vector or NULL")
  #needs to be numeric
  assertthat::assert_that(is.null(DE_idx) | all(is.numeric(DE_idx)),
                          msg = "DE_idx needs to be numeric")
  #needs to be within bounds
  if (!is.null(mtx))
    assertthat::assert_that(is.null(DE_idx) | all(DE_idx %in% 1:nrow(mtx)),
                          msg = "indices contained in DE_idx are outside the range of the provided input matrix")

}
#' Assert Validity of a Matrix
#'
#' Verifies that the provided object is a matrix with numeric entries and no rows
#' containing only NA values.
#'
#' @param mtx The matrix to validate.
#' @return Invisibly returns TRUE if all assertions pass.
#' @noRd
mtxAssertions <- function(mtx){
  assertthat::assert_that(!is.null(mtx),
                          msg = "specify a matrix as input")
  assertthat::assert_that(is.matrix(mtx),
                          msg = "mtx needs to be a matrix")
  assertthat::assert_that(all(is.numeric(mtx)),
                          msg = "mtx needs to be a numeric matrix")
  #no row in mtx should contain only NAs
  assertthat::assert_that(all(apply(mtx, 1, function(x) !all(is.na(x)))),
                          msg = "no row in mtx should contain only NAs")
}
