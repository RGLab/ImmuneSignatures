# Original information
#'---
#' title: "HIPC Datasets - Processing"
#' author: "Renaud Gaujoux & Stefan Avey"
#' date: "06/11/2014"
#' output:
#'  pdf_document:
#'    toc: true
#'    toc_depth: 1
#'
#' ---
#+setup, include = FALSE

# Evan Henrich Notes: This script was edited to work as a set of functions that could
# be called by the main hipc_full_pipeline() function.  Original work that was not needed
# for the core functionality was moved to the bottom of the script and commented out.  This
# included a large amount of plotting work.  Also, the make_rds() function was changed
# so that it would take the study name as the input instead of looking through a directory
# to build a list of study names.  Variable names were changed for clarity in places too.

#**********TESTING********************************************************************
#-----------DEPENDENCIES------------------------------------------------
# library(knitr)
# library(plyr)
# library(stringr)
# library(Biobase)
# ## library(CLIR)
# library(ggplot2)
# library(reshape2)
#
# # #EH NOTE: These packages are necessary too, but only need sourcing once
# source("https://bioconductor.org/biocLite.R")
# biocLite("hugene10sttranscriptcluster.db")
# biocLite("hgu133plus2.db")
# biocLite("illuminaHumanv4.db")
# biocLite("DESeq")
# # biocLite("illuminaHumanv1.db") # rm b/c changed only search for v4 for yale
# # biocLite("illuminaHumanv2.db") # rm b/c changed only search for v4 for yale
# # biocLite("illuminaHumanv3.db") # rm b/c changed only search for v4 for yale


#****************************************************************************************

#-------------HELPER METHODS----------------------------------------
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)), {s <- substring(s, 2); if(strict) tolower(s) else s}, sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

# process pids and infer annotation
.match_db <- function(x, db){
  n <- sapply(db, function(pkg){
    suppressMessages(library(pkg, character.only = TRUE))
    length(intersect(x, AnnotationDbi::keys(get(pkg, envir = asNamespace(pkg), inherits = FALSE))))
  })
  if( all(n == 0 ) ) setNames(0, '')
  else tail(sort(n), 1L)
}


#----------------MAIN METHODS-----------------------------------------------------------------
#' Function creates combined hai table and rds file from all six studies
#'
#' @param hai_dir Directory holding HAI data tables
#' @param output_dir Output Directory
#' @export
combine_hai_data <- function(hai_dir, output_dir){
  file_list <- unlist(lapply(hai_dir,
                            list.files,
                            recursive = TRUE,
                            full.names = TRUE))

  message("************* MAKING COMBINED HAI RDS ************************")

  sdy_list <- str_match(basename(file_list), "^((CHI)|(SDY[0-9]+))-?[^_]*")
  hai <- ldply(split(file_list, sdy_list[, 1L]), .id = 'Study', function(d){

    # extract info from filename: SDY, cohort, mfc titer base
    m <- str_match(basename(d), "(((CHI)|(SDY[0-9]+))-?([^_]*))_((combined)|(young)|(old))?_?hai")
    sdy <- unique(m[, 3L])
    cohort <- unique(m[, 6L])
    mfc_base <- m[, 7L]
    if( nzchar(cohort) ){
      cohort <- sprintf("%s-%s", sdy, str_trim(gsub("[-_)(]", "", cohort)))
    }else{
      cohort <- sdy
    }

    message(sprintf("  - %s [%s]", sdy, cohort))

    f <- d
    # load each file
    subjects <- character()
    res <- lapply(seq_along(f), function(i){
      f <- f[i]
      message("    -> Processing ", mfc_base[i], ": ", basename(f))

      hai <- read.table(f, header = TRUE, sep = "\t", row.names = 1L)

      subjects <<- unique(c(subjects, rownames(hai)))
      col <- grep("_max", colnames(hai))
      res <- hai[ , col, drop = FALSE]
      # flag with mfc base if necessary -- and rename/flag
      if( nzchar(mfc_base[i]) && mfc_base[i] != 'combined' ){
        colnames(res) <- paste0(mfc_base[i], "_", colnames(res))
      }
      res <- res[subjects,, drop = FALSE]
      rownames(res) <- subjects
      as.matrix(res)
    })

    # combine
    stopifnot( all(sapply(res, function(x) identical(rownames(x), rownames(res[[1]])))) )
    rn <- rownames(res[[1L]])
    if( !grepl("[^0-9]", rn[1L]) ){
      rn <- paste0("sub_", rn)
    }
    res <- do.call(cbind, res)
    data.frame(Cohort = cohort, SubjectID = rn, res, stringsAsFactors = FALSE)
  })

  # add age class
  age <- rep(NA, nrow(hai))
  if( length(i_age <- grep("^old_", names(hai))) ){
    age[!is.na(hai[[i_age[1]]])] <- 'Old'
  }
  if( length(i_age <- grep("^young_", names(hai))) ){
    age[!is.na(hai[[i_age[1]]])] <- 'Young'
  }
  hai$Age.class <- factor(age)

  # make factors and save data.frame object
  l_ply(grep("((4fc)|(d[0-9]+))$", names(hai)), function(i){
    hai[[i]] <<- factor(ifelse(is.nan(hai[[i]]), NA, hai[[i]]))
  })

  saveRDS(hai, file = paste0(output_dir,"/HAI.rds"))
  return(hai)
}

#' Function to generate expressionSet object and save as Rds file
#'
#' @param sdy Study
#' @param ge_dir Directory holding gene expression data tables
#' @param combined_hai Combined hai data frame
#' @param output_dir Output directory
#' @export
make_rds <- function(sdy, ge_dir, combined_hai, output_dir){

  message(paste0("**************** Processing ", sdy, " *******************"))

  if(!(sdy %in% c("SDY80","SDY67"))){
    ge_file <- paste0(ge_dir, "/", sdy, ".GEMatrix.txt")
  }else if(sdy == "SDY80"){
    ge_file <- paste0(ge_dir, "/", "CHI-nih.GEMatrix.txt")
  }else if(sdy == "SDY67"){
    ge_file <- paste0(ge_dir, "/", sdy, "-batch2.GEMatrix.txt")
  }

  # expression data
  if(sdy == "SDY400"){
    e <- as.data.frame(orig_400_GE)
  }else{
    e <- read.table(ge_file, sep = "\t", header = TRUE)
  }

  fdata <- e[, 1:2]
  e <- as.matrix(e[, -(1:2)])
  rownames(e) <- as.character(fdata[, 1L])

  ## Normalize RNA-seq Data (only for SDY67)
  if(sdy == "SDY67") {
    if(require(DESeq)) {
      ##  countTable is a data frame and condition is the colnames of countTable
      countTable <- e
      condition <- colnames(e)
      cds <-newCountDataSet(countTable, condition)
      cds <- estimateSizeFactors(cds) ## estimate size factor
      cdsBlind <- estimateDispersions(cds, method="blind" )
      vsd <- varianceStabilizingTransformation(cdsBlind)
      e <- exprs(vsd)                   # normalized data into `e`
    }
    message("    - ", sdy, " counts normalized using DESeq (version ", packageVersion("DESeq"), ")")
  }

  annot <- setNames(0, '')
  pids <- rownames(e)
  if( any(grepl("_PM_", pids, fixed = TRUE)) ){ pids <- gsub("_PM_", "_", pids, fixed = TRUE)}
  if( any(grepl("_at$", pids)) ) { annot <- .match_db(pids, c('hgu133plus2.db', 'hgu133plus2.db'))}
  if( all(grepl("^[0-9]+$", pids)) ) { annot <- .match_db(pids, c('hugene10sttranscriptcluster.db'))
  } else if( any(grepl("^ILMN", pids)) ) {annot <- .match_db(pids, sprintf('illuminaHumanv4.db'))}

  if( !annot ) warning("Could not infer annotation package for ", sdy)
  message(sprintf("    - Annotation: %s [%s/%s probes matched]", names(annot), annot, nrow(e)))
  annot <- names(annot)
  rownames(e) <- rownames(fdata) <- pids
  eset <- ExpressionSet(e, featureData = AnnotatedDataFrame(fdata), annotation = annot)
  eset$Study <- sdy

  # pheno data
  m <- str_match(sampleNames(eset), "(((SUB)|(sub))_?[0-9]+)(\\.[0-9])?_(.*)") # modified for SDY212
  subjectID <- m[, 2L]
  eset$Condition <- factor(m[, 7L])

  if( file.exists(pfile <- gsub("GEMatrix", "demographics", ge_file, fixed = TRUE)) ){
    p <- read.table(pfile, sep = "\t", header = TRUE)
  }else{
    warning("No phenotypic data for ", sdy)
  }

  p_subjectID <- as.character(p[, 1L])
  with_annot <- subjectID %in% p_subjectID
  message(sprintf("    - %s samples | %s/%s with annotation [%s available subjects]",
                  ncol(eset), sum(with_annot), ncol(eset), nrow(p)))
  if( !all(with_annot) ){
    stop("Missing sample annotation for ", paste0(sampleNames(eset)[!with_annot], collapse = ", "))
  }
  p <- p[match(subjectID, p_subjectID), ]
  # uniformise
  names(p) <- capwords(names(p))
  p <- cbind(pData(eset), p)
  rownames(p) <- sampleNames(eset)
  pData(eset) <- droplevels(p)

  # add HAI to the phenotypic data
  if(sdy != "SDY80"){
    cohort <- combined_hai[combined_hai$Study == sdy, ]
  }else{
    cohort <- combined_hai[combined_hai$Study == "CHI-nih", ]
  }

  if( !nrow(cohort) ){
    warning("No HAI data found for ", cohortID)
  }else{
    rownames(cohort) <- as.character(cohort$SubjectID)
    with_hai <- subjectID %in% rownames(cohort)
    message(sprintf("    - %s samples | %s/%s with HAI [%s - %s available subjects]",
                    ncol(eset), sum(with_hai), ncol(eset), sdy, nrow(cohort)))
    # do not include columns Study and SubjectID which are already there
    p <- cohort[subjectID, setdiff(colnames(cohort), c('Study', "SubjectID")), drop = FALSE]
    p <- cbind(pData(eset), p)
    rownames(p) <- sampleNames(eset)
    pData(eset) <- droplevels(p)
    varMetadata(phenoData(eset))$HAI <- isHAI <- grepl("_max", varLabels(eset))
    HAI_type <- rep('', ncol(pData(eset)))
    HAI_type[isHAI] <- ifelse(sapply(pData(eset)[isHAI], is.factor), "class", 'continuous')
    varMetadata(phenoData(eset))$HAI_type <- factor(HAI_type)
  }

   storageMode(assayData(eset)) <- "list" # change from "lockedEnvironment" to "list" to modify

  ## Store raw count data for SDY67
  if(sdy == "SDY67") {
    ## countTable should only exist if study is SDY67
    if(all.equal(colnames(countTable), colnames(e))) { # make sure samples are in same order
      assayData(eset)[["counts"]] <- countTable
    }
  }

  ## Adjust for age (not gender) and add these matrices to the eset
  ageClass <- pData(eset)$Age.class

  ## Fit 3 models, one for all, and one for old and young separately if data is available
  classes <- "All"
  if(any(!is.na(ageClass))) {           # If any are not missing (NA)
    if(any(tolower(ageClass) == "young", na.rm=TRUE)) {
      classes <- c(classes, "Young")
    }
    if(any(tolower(ageClass) == "old", na.rm=TRUE)) {
      classes <- c(classes, "Old")
    }
  }
  for(cl in classes) {
    message(sprintf("    - Calculating age-adjusted expression for %s subjects...", cl))

    ## Select only day 0 and the age class specified by cl
    if(cl == "All") {
      sel <- grep("^((d0)|(dneg))", eset$Condition)
    } else {
      sel <- which( grepl("^((d0)|(dneg))", eset$Condition) & (ageClass == cl) )
    }
    eMat <- exprs(eset)[,sel]
    age <- as.numeric(pData(eset)$Age[sel])

    ## Fit a linear model to each probe and use the intercept + residuals as
    ## the age-corrected values. Was previously correcting for gender effects as well but
    ## not any longer. If any values are NA (as currently in SDY212), then do not
    ## even attempt to correct them
    lm1 <- apply(eMat, 1, function(probeExpression) {
      if(!any(is.na(probeExpression) | is.infinite(probeExpression)))
        ## lm(probeExpression ~ age + gender, na.action="na.fail")
        lm(probeExpression ~ age , na.action="na.fail")
      else
        # return missing value for residuals because no model is fit
        list(residuals=rep(NA, length(sel)),
             coefficients=c(`(Intercept)`=NA, NA, NA))
    })
    resList <- lapply(lm1, '[[', "residuals")
    intercepts <- sapply(lm1, function(model) model$coefficients["(Intercept)"])
    rm(list=c("lm1", "eMat"))
    omit <- which(sapply(resList, function(x) any(is.na(x)))) # only include if all values are present
    if(length(omit) > 0) {
      message("    + The following probes have missing or infinite values and were not adjusted: ",
              paste(names(omit), collapse=', '))
    }
    combMat <- t(data.frame(resList, check.names=FALSE)) + intercepts
    finalMat <- matrix(nrow=nrow(eset), ncol=ncol(eset), dimnames=list(rownames(eset), colnames(eset)))
    finalMat[,sel] <- combMat
    assayData(eset)[[paste("adjustedExprs", cl, sep='-')]] <- finalMat
    rm(list=c("resList", "combMat", "finalMat"))
  }

  storageMode(assayData(eset)) <- "lockedEnvironment" # prevent changes

  if(!(sdy %in% c("SDY80","SDY67"))){
    file_out <- paste0(output_dir, "/", sdy, ".rds")
  }else if(sdy == "SDY80"){
    file_out <- paste0(output_dir, "/CHI-nih.rds")
  }else if(sdy == "SDY67"){
    file_out <- paste0(output_dir, "/", sdy, "-batch2.rds")
  }

  saveRDS(eset, file = file_out)
  message("    - Saved in file '", file_out, "'")
}



#*************************ORIGINAL CODE********************************
#' # Objective
#' Provide uniform access to datasets and end points, until all the data is available in ImmuneSpace
#'
#' # Data processing
#' ## HAI data
#+make_hai, results='asis', cache=TRUE

# if( !use_hai_combined ){
#   sdy_data <- unlist(lapply(HAI_dir, list.files, pattern = "^SDY.*_hai_.*\\.txt$", recursive = TRUE, full.names = TRUE))
#   sdy_data <- grep("hai_titer", sdy_data, invert = TRUE, value = TRUE)
# }else{
#   sdy_data <- unlist(lapply(HAI_dir, list.files, pattern = "^((SDY)|(CHI)).*_hai_titer_table\\.txt$", recursive = TRUE, full.names = TRUE))
# }
#----------------PARAMETERS-------------------------------------
# opts_chunk$set(echo = FALSE, size="scriptsize", dev='pdf', error = TRUE)
# options(width = 120)
# use_hai_combined <- TRUE
# force <- TRUE  # default is FALSE, whether to force recreating objects if they already exist
# save_data <- TRUE

# str(hai)

#-------------gene expression-------------------------------------
# #' ## Gene expression data
# #+make_eset, results='asis', cache=TRUE
# dfile <- list.files(DATA_dir, pattern = "GEMatrix", recursive = TRUE, full.names = TRUE)
# l_ply(dfile, function(f){
#   m <- str_match(basename(f), "^([^- _]+)(.*)\\.GEMatrix")
#   sdy <- m[, 2L]
#   sdy_extra <- str_trim(gsub("[-_)(]", "", m[, 3L]))
#   message(sprintf("  * %s [%s] - Processing %s", sdy, sdy_extra, basename(f)), appendLF = FALSE)
#   eset_file <- sprintf('%s%s.rds', sdy, if( nzchar(sdy_extra) ) paste0("-", sdy_extra) else '')
#   eset_file <- file.path(getwd(), "data", eset_file) # EH edit to save directly to data folder for analysis
#   if( file.exists(eset_file) && !force ){
#     message(" ... SKIP")
#     return()
#   }
#   message()

#---------------------Testing Code------------------------------------
# #' # Data overview
# #' The number of samples at baseline (day 0 or -7) are shown below.
# #' \scriptsize
# #+ desc, fig.width=10
# ## rds <- list.files(".", pattern = "^SDY.*\\.rds$")
# rds <- list.files(".", pattern = "^((SDY)|(CHI)).*\\.rds$")
# ## pd <- ldply(rds, function(x) pData(readRDS(x)))
# ## Limit to d0 only
# pd <- ldply(rds, function(x)
#             {
#               eset <- readRDS(x)
#               eset <- eset[,grep("^d0", eset$Condition)]
#               pData(eset)
#             })
# # Age dist
# library(RColorBrewer)
# ageCut <- 65
# pd$Age.group <- ifelse(pd$Age >= ageCut, paste("Over", ageCut), paste("Under", ageCut))
# p <- ggplot(data = pd, aes(x = Study, y=Age.group)) + geom_dotplot(binaxis="y", fill="#0838FF", binwidth=1, stackdir="center", aes(y = Age)) + ggtitle("Age distribution") + ylab("Age") + theme_bw() +
# theme(plot.background = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.border = element_blank(),
#       axis.line = element_line(colour = "black"))
# plot(p)
# p <- ggplot(data = pd, aes(x = Study, fill = paste(Gender, Age.group))) + scale_fill_manual('Groups', values = brewer.pal(4, "Paired")[1:4])
# p + geom_boxplot(aes(y = Age)) + ggtitle("Age distribution") + geom_jitter(aes(y = Age), position = position_jitterdodge(), show_guide = FALSE)
# p + geom_bar(stat = 'bin', position = position_dodge()) + geom_bar(stat = 'bin', position = position_dodge(), colour="black", show_guide=FALSE) + ggtitle("Sample size")
#
#
# #' # Loading datasets
# #' ## Formatted
# #' The file `getSDY.R` defines a function `.getSDY` that enables loading datasets, specifying which subset (e.g., young) and
# #' response variable to use, available in phenotypic variable `"Response"`.
# #' The returned `ExpressionSet` object contains all other HAI measures as well, but computed based on the selected subset.
# #' You can specify whether you want the expression to be adjusted for age by setting the `adjusted` argument in the `getSDY` function to `TRUE` or `FALSE` accordingly.  The adjustment is only done on the baseline data, so if there are other time points in the original data, they will be missing in the adjusted data.  NOTE: Adjustment is only correcting for age effects and not gender effects as these were found to not be important and to average out the signal.
# #+ getSDY, echo = TRUE
# # split data into young/old
# source('getSDY.R')
# # configure rds directory
# getSDY <- .getSDY('.')
# # all
# eset <- getSDY('SDY212')
# varLabels(eset)
# summary(pData(eset)[c('Age.class', 'Response')])
# # young only
# eset <- getSDY('SDY212', 'young')
# varLabels(eset)
# summary(pData(eset)[c('Age.class', 'Response')])
# # young only using another response definition
# eset <- getSDY('SDY212', 'young', 'fc_res_max_d30')
# varLabels(eset)
# summary(pData(eset)[c('Age.class', 'Response')])
# # What expression looks like without age adjustment
# head(exprs(eset))
# # What expression looks like with age adjustment
# eset <- getSDY('SDY212', 'young', 'fc_res_max_d30', adjusted=TRUE)
# head(exprs(eset))
#
#
# #' ## Using raw object
# #+ example, echo = TRUE
# # load data
# eset <- readRDS('SDY212.rds')
# eset
# varLabels(eset)
#
# # summarize HAI variables only
# hai_type <- varMetadata(phenoData(eset))$HAI_type
# summary(pData(eset)[hai_type == 'continuous'])
#
# # To get the different matrices out of assayData
# head(exprs(eset))
# head(assayData(eset)[["exprs"]]) # same as above using exprs(eset)
# head(assayData(eset)[["adjustedExprs-All"]]) # all subjects
# head(assayData(eset)[["adjustedExprs-Young"]]) # only those with Age.class == "Young"
# head(assayData(eset)[["adjustedExprs-Old"]]) # only those with Age.class == "Old"
#
# #' # All datasets
# #+ demographics, results='asis'
# # all demographics
# l_ply(rds, function(x){
#   eset <- readRDS(x)
#   cat("## ", as.character(eset$Study[1L]), "\n")
#   l_ply(capture.output(print(summary(pData(eset)[!varMetadata(phenoData(eset))$HAI]))), message)
# })
