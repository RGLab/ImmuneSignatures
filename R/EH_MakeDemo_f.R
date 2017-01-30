# Evan Henrich
# January 2017
# Fred Hutchinson Cancer Research Center
# Gottardo Lab
# ehenrich@fredhutch.org

# PURPOSE: generate DEMO tables for downstream HIPC Gene Module Analysis using raw data from Immunespace

# NOTES: The original code to perform these operations was developed by collaborators of 
# the HIPC Immune Signatures project.  This code was developed to provide the same functionality
# as the original in terms of statistical operation, but use data pulled directly from the 
# ImmuneSpace portal at www.immunespace.org instead of data shared among the collaborating
# labs via a google drive and also to handle all the studiesused in the meta-analysis in an 
# automated format.

#**************TESTING ONLY*****************************************
# setwd("/home/ehenrich/R/ImmSig_Testing/")
# library(ImmuneSpaceR)
# library(hash)
#*******************************************************************

#------MAIN METHOD---------------
makeDemo <- function(sdy, output_dir){
  
  # Setup directory vars
  wk_dir <- getwd()
  
  # Get rawdata from ImmuneSpace
  con <- CreateConnection(sdy)
  data <- data.frame(con$getDataset("demographics"))
  
  cols_to_keep <- c("participant_id","age_reported","gender")
  if(sdy %in% c("SDY212","SDY63","SDY404","SDY400")){
    cols_to_keep <- c(cols_to_keep,"phenotype","race","ethnicity")
  }else if(sdy == "SDY67"){
    cols_to_keep <- c(cols_to_keep,"race")
  }else if(sdy == "SDY80"){
    cols_to_keep <- c(cols_to_keep,"race")
  }
  
  # General data cleaning
  data <- data[ , names(data) %in% cols_to_keep]
  
  data$participant_id <- sapply(data$participant_id, FUN = function(x){
    res <- strsplit(x, ".", fixed = T)
    return(res[[1]][1])
  })
  
  data$age_reported <- sapply(data$age_reported, FUN = function(x){
    return(round(x, digits = 0))
  })
  
  # Study specific data cleaning, including removal of subjects to mimic original files
  if(sdy %in% c("SDY404","SDY400")){
    data <- data[ , c(1,3,4,5,6,2)]
    names(data) <- c("SubjectID","Gender","Age","Ethnicity","Race","frailty.Pheno")
    
    if(sdy == "SDY404"){
      data$frailty.Pheno <- sapply(data$frailty.Pheno, FUN = function(x){
        if(x == "Older adult, not-frail"){
          return("Old, not-frail")
        }else if(x == "Older adult, frail"){
          return("Old, frail")
        }else if(x == "Older adult, pre-frail"){
          return("Old, pre-frail")
        }else if(x == "Younger adult"){
          return("Young")
        }
      })
      subs_rm <- c("SUB120432","SUB120448","SUB120451")
      data <- subset(data, !(data$SubjectID %in% subs_rm))
    
      }else{
      data$frailty.Pheno <- sapply(data$frailty.Pheno, FUN = function(x){
        if(x == "Old, Non-frailty"){
          return("Old, not-frail")
        }else if(x == "Old Frailty"){
          return("Old, frail")
        }else if(x == "Old Pre-frailty" | x == "Old, Pre-frailty"){
          return("Old, pre-frail")
        }else{ 
          return(x)
        }
      })
      subs_rm <- c("SUB123271","SUB123319","SUB123334","SUB123338","SUB123343")
      data <- subset(data, !(data$SubjectID %in% subs_rm))
    }
  
  }else if(sdy == "SDY63"){
    # remove observations to mimic original
    subs_rm <- c("SUB112972", "SUB112976", "SUB112983", "SUB112988",
                 "SUB113013", "SUB113017", "SUB113018")
    data <- subset(data, !(data$participant_id %in% subs_rm))
    
    data <- data[, -2] # remove phenotype column to mimic original
    names(data) <- c("SubjectID","Gender","Age","Ethnicity","Race")
  
  }else if(sdy == "SDY212"){
    names(data) <- c("subjectID", "phenotype","gender","age","ethnicity","race")
  
  }else if(sdy == "SDY67"){
    # add exp samples
    subs_rm <- c("SUB113458","SUB113463","SUB113470","SUB113473","SUB113474",
                 "SUB113476","SUB113483","SUB113487","SUB113490","SUB113494",
                 "SUB113495","SUB113496","SUB113498","SUB113504","SUB113505",
                 "SUB113513","SUB113514","SUB113524","SUB113526","SUB113527",
                 "SUB113532","SUB113535","SUB113537","SUB113545","SUB113548",
                 "SUB113555","SUB113558","SUB113559","SUB113561","SUB113566",
                 "SUB113567","SUB113568","SUB113571","SUB113572","SUB113582",
                 "SUB113583","SUB113588","SUB113595","SUB113610")
    data <- subset(data, !(data$participant_id %in% subs_rm))
    data <- data[ , c(1,3,2,4)]
    names(data) <- c("SubjectID","Age","Gender","Race")
    exp_hash <- hash(SDY67_exp_map$old_res.SubjectID, SDY67_exp_map$old_res.exp_sample_acc)
    data$exp_sample_acc <- unlist(unname(sapply(data$SubjectID, FUN = function(x){
      return(exp_hash[[x]])
    })))
    sdy <- "SDY67-batch2"
  
  }else if(sdy == "SDY80"){
    names(data) <- c("SubjectID","Gender","Age","Race")
    sdy <- "CHI-nih"
  }
  
  write.table(data, file = paste0(output_dir, "/", sdy, ".demographics.txt"), 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
}
  
