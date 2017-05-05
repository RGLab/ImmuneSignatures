#' Loading HIPC Datasets
#' 
#' @param dir path to the directory containing the dataset .rds files
#' @param sdy SDY study id
#' @param group sample subset to load
#' @param response HAI discrete variable to load into variable "Response".
#' Responders are the samples in the last class, while non-responders are 
#' the samples in the first class. Samples from the intermediate class (if any)
#' get an NA value.
#' @param adjusted if TRUE, baseline expression data is adjusted for age (but NOT gender) with a linear regression model.  (intercept + residual) are given as the adjusted values for the baseline data only and any other time points will have missing values (NA).  Default is TRUE for Old subjects but FALSE for young or all subjects. NOTE: Behavior is different than in previous version!
#' @param baselineOnly When `TRUE`, only baseline (day 0 or pre-vaccination) time points will be returned. Default is `FALSE` for compatibility with previous versions.
#' @param keep.all logical that indicates if all the response measures should be kept (only used when \code{group = 'all'}).
#' @param feature.id type of feature identifier. 
#' If \code{'symbol'} then only probes mapped to symbols are returned. 
#' Multiple probes mapping to the same genes are averaged.
#' @return an ExpressionSet object with only one matrix in assayData. Which matrix depends on the combination of the group (age) and whether adjusted is set to TRUE or FALSE.
.getSDY <- function(dir = '.'){
  function(sdy, group = c('all', 'young', 'old'),
           response = c("fc_res_max_d20", "fc_res_max_d30", "fc_max_4fc",
               "fc_norm_max_d20", "fc_norm_max_d30"),
           swapSDY404 = TRUE,
           adjusted=ifelse(group == "old", TRUE, FALSE), baselineOnly=FALSE
           , keep.all = FALSE, feature.id = c('probes', 'symbol')) {
    ## parameters
    group <- match.arg(group)
    response <- match.arg(response)
    
    ## load data
    suppressMessages(library(Biobase))
    suppressMessages(library(stringr))  # required for str_match function
    
    # process SDY id
    sdy <- gsub("\\.rds$", "", sdy)
    sdy0 <- sdy
    # check if corrected data is requested
    sdy <- gsub("_c0$", '', sdy)
    
    if( sdy != sdy0 ){
      dir <- file.path(dir, 'corrected')
    }
    
   # look for group specific data file
    sdy_file_prefix <- sdy
    if( file.exists(file.path(dir, gprefix <- sprintf("%s-%s.rds", sdy_file_prefix, group))) ){
      sdy_file_prefix <- paste0(sdy_file_prefix, "-", group)
    }
    data_file <- file.path(dir, paste0(sdy_file_prefix, '.rds'))
    message(sprintf("Loading data from file '%s' [%s]", data_file, group))
    eset <- readRDS(data_file)

    ## change from "lockedEnvironment" to "list" to modify
    storageMode(assayData(eset)) <- "list" 

    ## Perform sample swaps to correct errors
    if(swapSDY404) {
      print("I am swapping 404!")
      if(sdy != "SDY404") {
        warning("Setting swapSDY404 to TRUE has no effect when sdy != 'SDY404'")
      }
      expr <- assayData(eset)$exprs
      ## Swap 2 with incorrect sex
      inds1 <- grep("SUB120473_d[027]", colnames(expr))
      inds2 <- grep("SUB120474_d[027]", colnames(expr))
      print(c(inds1, inds2))
      tmp1 <- expr[,inds1]
      expr[,inds1] <- expr[,inds2]
      expr[,inds2] <- tmp1
      ## Swap 2 additional samples with strong evidence
      inds1 <- grep("SUB120460_d0", colnames(expr))
      inds2 <- grep("SUB120485_d28", colnames(expr))
      print(c(inds1,inds2))
      tmp1 <- expr[,inds1]
      expr[,inds1] <- expr[,inds2]
      expr[,inds2] <- tmp1
      assayData(eset)$exprs <- expr
    }
    
    ## Limit to basline if baselineOnly flag is TRUE
    if(baselineOnly) {
      #eset <- eset[,grep("^((d0)|(dneg))", eset$Condition)]
      eset <- eset[,grep("^d0", eset$Condition)]
    }
    
    ## massage phenotypic data and select which assayData to use (adjusted or not)
    hai_var <- grep("^((fc)|(d0))_", varLabels(eset), value = TRUE)
    if(!adjusted) {
      assayData(eset) <- list(exprs=assayData(eset)[["exprs"]])
    }
    if( group == 'all' ){
      # to keep all different end metrics, comment out next line or use RDS file directly
      if( !keep.all ) phenoData(eset) <- phenoData(eset)[, grep("^((young)|(old))_", varLabels(eset), invert = TRUE)]
      if(adjusted) {
        assayData(eset) <- list(exprs=assayData(eset)[["adjustedExprs-All"]])
      } 
    }else if( group == 'young' ){
      if(any(eset$Age.class %in% 'Young') && adjusted) {
        assayData(eset) <- list(exprs=assayData(eset)[["adjustedExprs-Young"]])
      } 
      eset <- eset[, eset$Age.class %in% 'Young']
      if( !ncol(eset) ) warning("No samples found in group 'Young'") 
      # remove HAI variables computed on other groups (only if necessary)
      if( any(grepl("^young_", varLabels(eset))) ){
        phenoData(eset) <- phenoData(eset)[, setdiff(varLabels(eset),
                                                     c(hai_var, paste0("old_", hai_var)))]
        varLabels(eset) <- gsub("^young_", '', varLabels(eset))
      }
      
    } else if( group == 'old' ){
      if(any(eset$Age.class %in% 'Old') && adjusted) {
        assayData(eset) <- list(exprs=assayData(eset)[["adjustedExprs-Old"]])
      }
      eset <- eset[, eset$Age.class %in% 'Old']
      if( !ncol(eset) ) warning("No samples found in group 'Old'")
      
      # remove HAI variables computed on other groups (only if necessary)
      if( any(grepl("^old_", varLabels(eset))) ){
        phenoData(eset) <- phenoData(eset)[, setdiff(varLabels(eset),
                                                     c(hai_var, paste0("young_", hai_var)))]
        varLabels(eset) <- gsub("^old_", '', varLabels(eset))
      }
    }

    ## change back to "lockedEnvironment" to prevent chagnes    
    storageMode(assayData(eset)) <- "lockedEnvironment"

    ## define response variable: use extreme classes
    resp <- eset[[response]]
    lev <- levels(resp)
    lev <- c(head(lev, 1L), tail(lev, 1L))
    eset$Response <- factor(c('Non-responder', 'Responder')[match(resp, lev)])
    
    # convert to SYMBOL if requested
    feature.id <- match.arg(feature.id)
    if( feature.id == 'symbol' ){
      symb <- 'geneSymbol'
      if( is.null(fData(eset)[[symb]]) ) symb <- 'gs' 
      eset <- eset[!is.na(fData(eset)[[symb]]) & nzchar(as.character(fData(eset)[[symb]])), ]
      fData(eset) <- droplevels(fData(eset))
      # convert to symbol
      library(applyBy)
      if( xbioc::is_logscale(eset) ){
        eset <- logb(colMeansBy(expb(eset, 2), fData(eset)[[symb]]), 2)
      }else eset <- colMeansBy(eset, fData(eset)[[symb]])
    }
    
    ## result
    eset
  }
}
