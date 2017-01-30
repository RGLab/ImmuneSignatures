#' Function loads hipc eset objects from rds files and performs GE adjustment and subsetting, if directed via params.
#'
#' @param rds_data_dir path to the directory containing the dataset .rds files
#' @param sdy ImmuneSpace study id
#' @param group sample subset to load, either young or old
#' @param response HAI discrete variable to load into variable "Response".
#' Responders are the samples in the last class, while non-responders are
#' the samples in the first class. Samples from the intermediate class (if any)
#' get an NA value.
#' @param adjusted if TRUE, baseline expression data is adjusted for age (but NOT gender) with a linear regression model.  (intercept + residual) are given as the adjusted values for the baseline data only and any other time points will have missing values (NA).  Default is TRUE for Old subjects but FALSE for young or all subjects. NOTE: Behavior is different than in previous version!
#' @param baselineOnly When `TRUE`, only baseline (day 0 or pre-vaccination) time points will be returned. Default is `FALSE` for compatibility with previous versions.
#' @return an ExpressionSet object with only one matrix in assayData. Which matrix depends on the combination of the group (age) and whether adjusted is set to TRUE or FALSE.
#' @export

# EH notes: edited function to have clearer arguments
getSDY <- function(rds_data_dir, sdy, group, response, adjusted = FALSE, baselineOnly = FALSE){

  # Problem: The way the original .getSDY() was written, adjusted was supposed to be
  # set in the arguments, but because the group has not yet been set, then it remains FALSE
  # when the intention was to have it be TRUE for "old" cohort.
  # Original below as it was in the argument.
  # adjusted <- ifelse(group == "old", TRUE, FALSE)

  eset <- readRDS(file.path(rds_data_dir, paste0(gsub("\\.rds$", "", sdy), '.rds')))

  ## Limit to basline if baselineOnly flag is TRUE
  if(baselineOnly) {
    eset <- eset[,grep("^((d0)|(dneg))", eset$Condition)]
  }

  ## change from "lockedEnvironment" to "list" to modify
  storageMode(assayData(eset)) <- "list"

  ## massage phenotypic data and select which assayData to use (adjusted or not)
  hai_var <- grep("^((fc)|(d0))_", varLabels(eset), value = TRUE)
  if(!adjusted) {
    assayData(eset) <- list(exprs=assayData(eset)[["exprs"]])
  }
  if( group == 'all' ){
    # to keep all different end metrics, comment out next line or use RDS file directly
    phenoData(eset) <- phenoData(eset)[, grep("^((young)|(old))_", varLabels(eset), invert = TRUE)]
    if(adjusted) {
      assayData(eset) <- list(exprs=assayData(eset)[["adjustedExprs-All"]])
    }
  }else if( group == 'young' ){
    if(any(eset$Age.class %in% 'Young') && adjusted) {
      assayData(eset) <- list(exprs=assayData(eset)[["adjustedExprs-Young"]])
    }
    eset <- eset[, eset$Age.class %in% 'Young']
    if( !ncol(eset) ) warning("No samples found in group 'Young'")
    phenoData(eset) <- phenoData(eset)[, setdiff(varLabels(eset),
                                                 c(hai_var, paste0("old_", hai_var)))]
    varLabels(eset) <- gsub("^young_", '', varLabels(eset))
  } else if( group == 'old' ){
    if(any(eset$Age.class %in% 'Old') && adjusted) {
      assayData(eset) <- list(exprs=assayData(eset)[["adjustedExprs-Old"]])
    }
    eset <- eset[, eset$Age.class %in% 'Old']
    if( !ncol(eset) ) warning("No samples found in group 'Old'")
    phenoData(eset) <- phenoData(eset)[, setdiff(varLabels(eset),
                                                 c(hai_var, paste0("young_", hai_var)))]
    varLabels(eset) <- gsub("^old_", '', varLabels(eset))
  }
  ## change back to "lockedEnvironment" to prevent chagnes
  storageMode(assayData(eset)) <- "lockedEnvironment"

  ## define response variable: use extreme classes
  resp <- eset[[response]]
  lev <- levels(resp)
  lev <- c(head(lev, 1L), tail(lev, 1L))
  eset$Response <- factor(c('Non-responder', 'Responder')[match(resp, lev)])

  ## result
  return(eset)
}

