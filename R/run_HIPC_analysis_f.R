# Evan Henrich
# ehenrich@fredhutch.org
# Gottardo Lab - Fred Hutchinson Cancer Research Center
# January 2017

#----------------------HELPER METHODS--------------------------------------
#' Generates TSV text tables for use in the hipc pipeline from ImmuneSpace data
#'
#' @param studies a vector of characters holding the study names for which to generate eSets
#' @param hai_dir a path for the preprocessed HAI data
#' @param ge_dir a path for the preprocessed GE and demographic data
#' @param rawdata_dir a path for directory to download rawdata into
#' @export
hipc_preprocess <- function(studies, hai_dir, ge_dir, rawdata_dir, orig_params = T){
  # get user input
  user <- readline(prompt = "Username for ImmuneSpace: ")
  pwd <- readline(prompt = "Password for ImmuneSpace: ")

  if(orig_params == F){
    yale.anno <- readline(prompt =
                            "Yale Studies' Gene Annotation: [H = HIPC manuscript, l = R library, m = Illumina manifest]  ")
    yale.anno <- ifelse(yale.anno %in% c("H", "h", ""), "o", yale.anno)
    sdy80.process <- readline(prompt =
                                "Process SDY80 rawdata or use pre-processed tables from HIPC manuscript? [p = process / H = HIPC manuscript] ")
    sdy80.process <- ifelse(sdy80.process %in% c("H", "h", ""), FALSE, TRUE)

    if(sdy80.process == TRUE){
      sdy80.anno <- readline(prompt = "SDY80 Gene Annotation: [H = HIPC manuscript, l = R library]  ")
      sdy80.anno <- ifelse(sdy80.anno %in% c("H", "h", ""), "o", sdy80.anno)
      sdy80.norm <- readline(prompt = "Re-normalize SDY80 GE data? [T / f]  ")
      sdy80.norm <- ifelse(sdy80.norm %in% c(T, "", "t", "T"), TRUE, FALSE)
    }
  }else{
    yale.anno <- "o"
    sdy80.process <- F
    sdy80.anno <- "o"
    sdy80.norm <- F
  }

  input_map <- list()
  input_map$o <- "original"
  input_map$l <- "library"
  input_map$m <- "manifest"

  yale.anno <- input_map[[yale.anno]]
  sdy80.anno <- input_map[[sdy80.anno]]

  for(sdy in studies){
    if( !(sdy %in% c("SDY400","SDY80")) | (sdy == "SDY80" & sdy80.process == TRUE) ){
        print(paste0("Generating GE for ", sdy))
        makeGE(sdy, user, pwd,
               yale.anno = yale.anno,
               sdy80.anno = sdy80.anno,
               sdy80.norm = sdy80.norm,
               output_dir = ge_dir,
               rawdata_dir = rawdata_dir)
        print(paste0("Generating HAI for ", sdy))
        makeHAI(sdy, output_dir = hai_dir)
        print(paste0("Generating Demo for ", sdy))
        makeDemo(sdy, output_dir = ge_dir)
    }else if(sdy == "SDY400"){
      print(paste0("Generating HAI for ", sdy))
      makeHAI(sdy, output_dir = hai_dir)
      print(paste0("Generating Demo for ", sdy))
      makeDemo(sdy, output_dir = ge_dir)
    }
  }
}

#' Generates rds objects from studies given the study names and directories holding the rawdata.
#'
#' @param studies a vector of characters holding the study names for which to generate eSets
#' @param hai_dir a path for the preprocessed HAI data
#' @param ge_dir a path for the preprocessed GE and demographic data
#' @param rds_dir a path for the output files (rds)
#' @export
hipc_make_rds <- function(studies, hai_dir, ge_dir, rds_dir){
  message("Now combining files for each study into rds / eset file")
  message("\n")
  combined_hai <- combine_hai_data(hai_dir, output_dir = rds_dir)
  for(sdy in studies){
    make_rds(sdy, ge_dir, combined_hai, output_dir = rds_dir)
  }
}

#' Perform HIPC Meta Analysis using given parameters
#'
#' Takes in various parameters related to the cutoff points for gene pathway selection, HAI endpoint selection,
#' adjustments to gene expression data, and output types.
#' @param rds_dir the path for the directory holding the RDS files with esets of each study
#' @param cohort the age cohort, either 'young' or 'old', to use in analysis
#' @param orig_params a boolean that defaults to TRUE, but when set to false allows the user to specify non-original parameters for the analysis
#' @param output_dir the path for the pdf and text outputs of the analysis
#' @export
hipc_meta_analysis <- function(rds_dir, cohort, orig_params = T, output_dir){

  if(orig_params == F){
    FDR.cutoff <- as.float(readline(prompt = "Input FDR.cutoff (orig = 0.5):  "))
    pvalue.cutoff <- as.float(readline(prompt = "Input pvalue.cutoff (orig = 0.01):  "))
    endPoint <- as.integer(readline(prompt = "Input HAI discretization low % cutoff [20 or 30]:  "))
    endPoint <- paste0("fc_res_max_d", endPoint)
    adjusted <- readline(prompt = "Adjusted = [T / f] ")
    adjusted <- ifelse(adjusted %in% c(T, TRUE, "t", t), TRUE, FALSE)
    baselineOnly <- readline(prompt = "baselineOnly = [T / f] ")
    baselineOnly <- ifelse(baselineOnly %in% c(T, TRUE, "t", t), TRUE, FALSE)
    indiv_rds <- readline(prompt = "output individual Rds of gene module analysis? [T / F]  ")
    indiv_rds <- ifelse(indiv_rds %in% c(T, TRUE, "t", t), TRUE, FALSE)

    message(paste0("Running meta analysis for ", cohort, " cohort"))
    message("\n")
    meta_analysis(geneSetDB = geneSetDB,
                  rds_data_dir = rds_dir,
                  cohort = cohort,
                  FDR.cutoff = FDR.cutoff,
                  pvalue.cutoff = pvalue.cutoff,
                  endPoint = endPoint,
                  adjusted = adjusted,
                  baselineOnly = baselineOnly,
                  indiv_rds = indiv_rds,
                  output_dir = output_dir)
  }else{
    message(paste0("Running meta analysis for ", cohort, " cohort"))
    meta_analysis(geneSetDB = geneSetDB,
                  rds_data_dir = rds_dir,
                  cohort = cohort,
                  output_dir = output_dir)
  }
}

#----------------MAIN METHOD-----------------------------------------
#' Runs full pipeline from start to finish with user input required at certain points
#' @importFrom xml2 read_html
#' @importFrom XML xmlToList xmlParse
#' @importFrom plyr ldply l_ply
#' @import dplyr
#' @importFrom httr GET authenticate write_disk
#' @importFrom stringr str_sub str_match str_trim
#' @importFrom data.table fread setnames
#' @importFrom hash hash has.key
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom GEOquery gunzip getGEO Table
#' @importFrom tibble as_tibble
#' @import qusage
#' @import knitr
#' @import Biobase
#' @import hgu133plus2.db
#' @import DESeq
#' @import illuminaHumanv4.db
#' @import hugene10sttranscriptcluster.db
#'
#' @return text and pdf files wtih significant gene pathways
#' @export
hipc_full_pipeline <- function(){

  studies <- c("SDY212", "SDY63", "SDY404", "SDY400", "SDY80", "SDY67")
  directory <- readline(prompt = "Save output files to new directory called 'ImmSig_Analysis'? [ T / f ] ")
  if(directory %in% c(T, "T", "t", "")){
    wkdir <- file.path(getwd(),"ImmSig_Analysis")
    dir.create(path = wkdir)
  }else if(directory %in% c(F, "F", "f")){
    wkdir <- readline(prompt = "Please specify full path for directory where you want to save files: ")
    if(!dir.exists(wkdir)){
      stop("Directory does not exist. Exiting analysis")
    }
  }else{
    stop("Terminating: direction not understood")
  }


  # 1. Extract and pre-process data from ImmuneSpace, ImmPort, or GEO databases
  gen_files <- readline(prompt = "Generate and preprocess rawdata? [T / f]  ")
  if(gen_files %in% c(T, "T", "t", "")){
    message("Creating PreProc_Data directory and sub-directories ... ")
    # create preproc dir and subdirectories
    pp_dir <- file.path(wkdir,"PreProc_Data")
    dir.create(path = pp_dir)

    hai_dir <- file.path(pp_dir,"HAI")
    dir.create(path = hai_dir)

    ge_dir <- file.path(pp_dir,"GE")
    dir.create(path = ge_dir)

    rawdata_dir <- file.path(pp_dir,"rawdata")
    dir.create(path = rawdata_dir)

    orig_params <- readline(prompt = "Use original parameters from HIPC manuscript? [ T / f ] ")
    orig_params <- ifelse(orig_params %in% c(T, "", "t", "T"), TRUE, FALSE)

    hipc_preprocess(studies, hai_dir, ge_dir, rawdata_dir, orig_params)
  }

  cat("Moving on to RDS Generation")

  # Step 2: Combine pre-processed files into Rds (BioConductor eset)
  run_rds <- readline(prompt = "Are you ready to run rds generation file? [T / f]  ")
  rds_dir <- file.path(wkdir,"Rds_data")
  dir.create(path = rds_dir)

  if(run_rds %in% c(T, "T", "t", "")){
    hipc_make_rds(studies, hai_dir, ge_dir, rds_dir)
  }

  cat("Moving on to Meta Analysis")

  # Step 3: Run meta analysis script
  run_meta <- readline(prompt = "Are you ready to run meta analysis? [T / f]  ")
  if(run_meta %in% c(T, "T", "t", "")){
    output_dir <- file.path(wkdir,"meta_analysis_output")
    dir.create(path = output_dir)
    hipc_meta_analysis(rds_dir, cohort = "young", orig_params = T, output_dir = output_dir)
    hipc_meta_analysis(rds_dir, cohort = "old", orig_params = T, output_dir = output_dir)
  }else{
    stop("process stopped")
  }
  cat(paste0("Output saved to ", output_dir))
}

#-----------------------------------------------------------------------




