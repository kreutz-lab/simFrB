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

  oldseed <- mlr3misc::get_seed()

  if (is.null(seed)) {
    seed <- mlr3misc::get_seed()
  }

  set.seed(seed)
  res <- functionName(...)

  .Random.seed <- oldseed

  return(list(res = res, seed = seed))
}

#' Load Data from a GitHub Raw URL
#'
#' This function downloads and loads data stored in an RDS file hosted on GitHub,
#' specifically from a raw content URL. It is designed to ensure the URL is correct,
#' downloads the file, loads it into R, and then cleans up by deleting the downloaded file.
#'
#' @param githubURL The full URL to the raw RDS file on GitHub.
#' This URL should start with "https://raw.githubusercontent.com/".
#' The correct URL can be obtained by navigating to the file on GitHub,
#' clicking "Raw", and then copying the URL from the browser's download manager.
#'
#' @return Returns the data loaded from the RDS file.
#'
#' @examples
#' \dontrun{
#'   # Example URL (this is a placeholder, replace with a real URL for actual use)
#'   url <-
#'   "https://raw.githubusercontent.com/kreutz-lab/
#'   dia-benchmarking/main/data/diaWorkflowResults_allDilutions.rds"
#'   data <- loadDataFromGit(url)
#' }
#'
#' @export
loadDataFromGit <- function(githubURL){
  #check that url starts with https://raw.githubusercontent.com/
  if(grep("https://raw.githubusercontent.com/", githubURL) == 0){
    stop("URL must start with https://raw.githubusercontent.com/.
         Get this download link by clicking the 'Download' button in github
         and left clicking the file in the download manager.")
  }
  utils::download.file(url = githubURL,
                destfile = "data.rds", method = "curl")

  gitData <- readRDS("data.rds")

  #delete downloaded file again
  file.remove("data.rds")
  return(gitData)
}


#'
#' This function subsets proteomics data based on specified Data-Independent Acquisition (DIA) workflow and
#' experimental comparison groups. It also allows further subsetting by specific rows and columns.
#'
#' @param data List of datasets keyed by DIA workflow names. Can be loaded with
#' loadDatarFromGit()
#' @param DIAWorkflow Character string specifying the DIA workflow to subset from `data`.
#' @param experimentalComparisonGroups Character vector of groups to include based on their
#'        tags in column names which might include specific concentrations of proteins or control setups.
#'        This is highly adapted to the data from https://github.com/kreutz-lab/dia-benchmarking/tree/main
#' @param rowSubset Optional numeric vector indicating the rows to be retained.
#' @param colSubset Optional numeric vector indicating the columns to be retained.
#'
#' @return A subset of the data frame corresponding to the specified DIA workflow and experimental groups.
#'         Further subsets by rows and columns if specified.
#'
#' @details The function requires a valid `diaWorkflow` which should be one of the keys in `data`.
#'          The `experimentalComparisonGroups` should match parts of the column names in the data frame
#'          corresponding to the groups of interest. This function stops and throws an error
#'          if the required workflow or comparison groups are not found or if multiple workflows are specified.
#'
#' @examples
#' \dontrun{
#' rownames <- paste0("feature_",1:20)
#' colnames <- c(sample(1:5,5), paste0("group1_",sample(6:10,5)), paste0("group2_",sample(6:10,5)))
#'   data_list <- list(
#'     DIANN_DIANN_AI = data.frame(matrix(ncol = 15, nrow = 20, dimnames = list(rownames, colnames))),
#'     DIANN_MaxQuant = data.frame(matrix(ncol = 15, nrow = 20, dimnames = list(rownames, colnames)))
#'   )
#'   DIAWorkflow <- "DIANN_DIANN_AI"
#'   groups <- c("control", "group2")
#'   subsetted_data <- subsetDIAWorkflowData(data = data_list, DIAWorkflow, groups)
#' }
#'
#' @export
subsetDIAWorkflowData <- function(data,
                                  DIAWorkflow,
                                  experimentalComparisonGroups,
                                  rowSubset = NULL,
                                  colSubset = NULL){


  if (is.null(DIAWorkflow)) {
    stop("Please provide a DIAWorkflow, options are: DIANN_DIANN_AI,
    DIANN_DIANN_AI_GPF, DIANN_MaxQuant, DIANN_MSFragger,
         DIANN_PROSIT_EDIA_GPF, OSW_DIANN_AI_GPF, OSW_MaxQuant, OSW_MSFragger,
         SkylineDIANN_DIANN_AI_GPF, SkylineDIANN_MaxQuant,
         SkylineDIANN_MSFragger, SkylineDIANN_PROSIT_EDIA_GPF,
         SkylineMSstats_DIANN_AI_GPF, SkylineMSstats_MaxQuant,
         SkylineMSstats_MSFragger, SkylineMSstats_PROSIT_EDIA_GPF,
         Spectronaut_DIANN_AI_GPF, Spectronaut_DirectDIA,
         Spectronaut_MaxQuant, Spectronaut_MSFragger,
         Spectronaut_PROSIT_EDIA_GPF")
  }
  if (is.null(experimentalComparisonGroups)) {
    stop("Please provide experimentalComparisonGroups, options are:
         c(control, 1-25, 1-12, 1-6), these refer to the groups, and their
         respective concentration of ecoli proteins in the dilution series.
         Control contains no Ecoli proteins.
         1-25 contains 25 parts human proteins and 1 part Ecoli proteins.")
  }

  allWorkflows <- names(data)

  if (length(DIAWorkflow) != 1)
    stop("Only one DIAWorkflow is allowed.")
  else if (!DIAWorkflow %in% allWorkflows)
    stop(paste("DIAWorkflow not found in data. Please provide a valid DIAWorkflow. Options are:",
               paste(allWorkflows, collapse = ", ")))
  else
    exp.df.all <- data[[DIAWorkflow]]

  colIdx <- c()
  for (grouptag in experimentalComparisonGroups) {
    colIdx.group <- grep(grouptag, colnames(exp.df.all))
    if (grouptag == "control") {
      #control groups should not contain "_" in colnames
      colIdx.group <- which(!grepl("_", colnames(exp.df.all)))
      if (is.null(colIdx.group)){
        stop("No control group found in data" )
      }
    }else if (length(colIdx.group) == 0) {
      stop(paste("No group found in data for", grouptag, "."))
    }
    colIdx <- c(colIdx, colIdx.group)
  }
  exp.df <- exp.df.all[,colIdx]

  exp.df <- subsetRowsCols(exp.df, rowSubset, colSubset)

  return(exp.df)
}


#' Subset Rows and Columns of a Data Object
#'
#' This function subsets rows and/or columns of a data object based on the provided indices.
#' It allows for flexible subsetting by either rows, columns, or both.
#'
#' @param data The data object to be subsetted, typically a data frame or matrix.
#' @param rowSubset Numeric or logical vector specifying the rows to be retained.
#'                  If empty, no subsetting is done on rows.
#' @param colSubset Numeric or logical vector specifying the columns to be retained.
#'                  If empty, no subsetting is done on columns.
#'
#' @return The subsetted data object with the specified rows and/or columns retained.
#'         If no subsetting parameters are provided, the original data is returned unchanged.
#'
#' @examples
#' data_matrix <- matrix(1:9, nrow = 3, ncol = 3)
#' rows_to_keep <- c(1, 3)
#' cols_to_keep <- c(2, 3)
#' subsetted_data <- subsetRowsCols(data_matrix, rows_to_keep, cols_to_keep)
#'
#' @noRd
subsetRowsCols <- function(data, rowSubset, colSubset) {
  if (length(rowSubset) > 0 & length(colSubset) == 0)
    data <- data[rowSubset,]
  if (length(rowSubset) == 0 & length(colSubset) > 0)
    data <- data[,colSubset]
  if (length(rowSubset) > 0 & length(colSubset) > 0)
    data <- data[rowSubset,colSubset]
  return(data)
}

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL
