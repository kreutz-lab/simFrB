#' @title msb.applyFunctionWithSeed
#'
#' @description Function to apply seed to a function which has a
#' random component, and
#' reset to previous seed afterwards
#'
#'
#' @param functionName Function name as string
#' @param seed Seed defined outside this function, default = 123
#' @param ... further parameters passed to function
#'
#' @return Result of the function the seed is applied to (res) and seed (seed)
#'
#' @examples
#' mtx <- matrix(runif(100000),nrow=1000)
#' randomCols <- msb.applyFunctionWithSeed(sample, seed = 83434,
#' x = 1:ncol(mtx), size= min(50, ncol(mtx)))
#'
#'
#' @keywords seed
#' @export

msb.applyFunctionWithSeed <- function(functionName, seed = 123,  ...){
  oldseed <- .Random.seed

  #oldseed <- .Random.seed
  oldseed <- mlr3misc::get_seed()

  if (is.null(seed)) {
    # seed <- .Random.seed
    seed <- mlr3misc::get_seed()
  }

  set.seed(seed)
  res <- functionName(...)

  .Random.seed <- oldseed

  return(list(res = res, seed = seed))
}
