#' Prepare Data for Logistic Regression
#'
#' This function prepares a dataset for logistic regression analysis by transforming
#' a given matrix into a long format. It creates a binary response variable indicating
#' the presence (1) or absence (0) of data (NA). Predictors include row ID, column ID,
#' group assignments based on a provided group design, and group-specific row means.
#'
#' @param mtx A matrix with features as rows and samples as columns.
#' @param groupDesign A vector indicating the group assignment for each column of the matrix.
#'        The default splits columns evenly into two groups.
#'
#' @return A list containing two elements: a data frame suitable for logistic regression
#'         and a formula object. The data frame includes binary indicators for data presence,
#'         row IDs, column IDs, group assignments, and group-specific means as predictors.
#'         The formula defines the model to be fitted, including all predictors.
#'
#' @examples
#' set.seed(123)
#' # Create a matrix with random data and some NAs
#' data_matrix <- matrix(rnorm(100), nrow = 10)
#' data_matrix[sample(1:100, 20)] <- NA
#'
#' # Prepare data for logistic regression
#' result <- prepareLogisticRegressionData(data_matrix)
#'
#' # View the prepared data
#' head(result$data)
#'
#' # View the regression formula
#' print(result$formula)
#'
#' @importFrom reshape2 melt
#' @importFrom tidyr pivot_wider
#' @export
prepareLogisticRegressionData <- function(
    mtx,
    groupDesign = rep(c(1, 2),each = (ncol(mtx) / 2))) {
  if (length(groupDesign) != ncol(mtx)) {
    stop("The length of 'groupDesign' must be equal to the number of columns in 'mtx'.")
  }

  rownames(mtx) <- seq_len(nrow(mtx))
  colnames(mtx) <- seq_len(ncol(mtx))

  # Create a data frame from the matrix in a long format
  data_long <- as.data.frame(as.table(mtx))

  # Rename columns for clarity
  names(data_long) <- c("RowID", "ColID", "Value")

  # Create a response variable where 1 indicates data is present (not NA)
  data_long$Presence <- as.integer(!is.na(data_long$Value))

  data_long <- data_long %>% dplyr::group_by(RowID) %>% dplyr::mutate(Mean = mean(Value, na.rm = TRUE))

  data_long$Group <- as.factor(groupDesign[data_long$ColID])


  #include group assignment in rowID
  data_long$RowID <- paste0(data_long$RowID,"#",data_long$Group)
  data_long$RowID <- as.factor(data_long$RowID)
  # Create the model formula
  # Remove the original 'Value' column as it's not needed for the model
  data_long$Value <- NULL
  dsgn <- formula(Presence ~ Mean + ColID + RowID)
  # Return the data frame and the formula
  return(list(data = data_long, design = dsgn))
}
#' Perform Logistic Regression Using speedglm
#'
#' This function performs logistic regression on a provided dataset using the speedglm package.
#' For datasets larger than 500 rows, the function fits the model in chunks to manage memory
#' and computation efficiency.
#'
#' @param data_long The data frame prepared for logistic regression.
#' @param formula The logistic regression formula.
#'
#' @return The logistic regression model object.
#'
#' @examples
#' # Assuming 'data_long' and 'formula' are defined as outputs from previous function:
#' # fit_model <- perform_logistic_regression(data_long, formula)
#' # summary(fit_model)
#'
#' @export
perform_logistic_regression <- function(data_long, design) {

  n_rows <- nrow(data_long)
  chunk_size <- 500000  # Define the chunk size

  tmp_data1 <- tempfile("data1",fileext=".rds")
  saveRDS(data_long,tmp_data1)

  # If data is less than or equal to 500 rows, fit the model on the entire dataset
  if (n_rows <= chunk_size) {
    fit <- speedglm::speedglm(formula = design, data = data_long, family = binomial())
  } else {
    # If data is larger, perform logistic regression in chunks
    fit <- NULL
    for (start_row in seq(1, n_rows, by = chunk_size)) {
      end_row <- min(start_row + chunk_size - 1, n_rows)
      chunk_data <- data_long[start_row:end_row, ]

      if (is.null(fit)) {
        # Initial model fitting
        fit <- speedglm::speedglm(formula = design, data = chunk_data, family = binomial())
      } else {
        # Update the model with new data chunk
        fit <- speedglm::updateWithMoreData(fit, formula = design, data = chunk_data)
      }
    }
  }

  return(fit)
}
#' create data for chunks for shglm
#'
make.data <- function(filename, chunksize,...){
  conn <- NULL
  function(reset = FALSE){
    if (reset){
      if (!is.null(conn)) close(conn)
      conn <<- gzfile(filename)
    } else{
      all <- readRDS(conn)
      rval <- all[1:chunksize,]
      if ((nrow(rval) == 0)) {
        close(conn)
        conn <<- NULL
        rval <- NULL
      }
      return(rval)
    }
  }
}


