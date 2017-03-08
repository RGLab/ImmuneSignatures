# Evan Henrich
# January 2017
# Fred Hutchinson Cancer Research Center
# Gottardo Lab
# ehenrich@fredhutch.org

# PURPOSE: generate HAI tables for downstream HIPC Gene Module Analysis using raw data from Immunespace

# NOTES: The original code to perform these operations was developed by Yuri Kotliarov at
# NIH (yuri.kotliarov@nih.gov).  This version was developed to provide the same functionality as
# the original in terms of statistical operation, but use data pulled directly from the ImmuneSpace
# portal at www.immunespace.org instead of data shared among the collaborating labs via a google
# drive and also to handle all the studies used in the meta-analysis in an automated format.


#------Helper Functions---------

nm_builder <- function(prefix, strains){
  tmp <- unlist(lapply(strains, FUN = function(strain){ paste0(prefix, strain)}))
}

# Calc median and SD, return as part of list of lists
med_sd_calc <- function(prefix, strains, glob_vals, titer_data){
  suf <- c("med","sd")
  for(virus in strains){
    for(s in suf){
      tar_list <- paste0(prefix,s)
      tar_col <- paste0(prefix,virus)
      if(s == "med"){
        glob_vals[[tar_list]][[virus]] <- median(titer_data[[tar_col]], na.rm = TRUE)
      }else{
        glob_vals[[tar_list]][[virus]] <- sd(titer_data[[tar_col]], na.rm = TRUE)
      }
    }
  }
  return(glob_vals)
}

# to mimic matlab median(x,1) to bypass NA and infinite vals
robust_med <- function(vec){
  return(median(vec[!is.na(vec) & !is.infinite(vec)]))
}

# to mimic matlab median absolute deviation func - mad(x,1)
r_mad <- function(vec){
  robmed <- robust_med(vec)
  mad <- median(abs(vec - robmed))
  return(mad)
}

# Seems computationally costly, but enough edge cases merit checking for valid data before extraction
# Valid means at least two titers per strain, with one being zero and the other greater than zero.
sub_check <- function(sub_df, strains){
  result <- list()
  for(vir in strains){
    d_zero <- sub_df[which(sub_df$virus == vir & sub_df$study_time_collected == 0),]
    d_other <- sub_df[which(sub_df$virus == vir & sub_df$study_time_collected > 0),]
    if(nrow(d_zero) > 0 & nrow(d_other) > 0){
      result[[vir]] <- TRUE
    }else{
      result[[vir]] <- FALSE
    }
  }
  final <- !(FALSE %in% result)
  return(final)
}

max_select <- function(subid, trg_col){
  tmp_ls <- list()
  for(virus in gl_strains){
    tmp_ls[virus] <- gl_tdata[gl_tdata$subject == subid, paste0(trg_col, virus)]
  }
  return(max(unlist(tmp_ls), na.rm = TRUE))
}

# Method by Yuri Kotliarov to categorize an observation based on low and high percentiles
# Changed slightly to round fc_res_max values to 7 digits prior to comparison with quantile
# values which are interpolated and therefore can throw off assignment if intention is to
# use them as if they are nearest order statistic (similar to type = 3 in quantiles args).
discretize <- function(df, input_col, low_perc, sdy, name){
  xq <- quantile(df[[input_col]],
                 c( (low_perc/100), 1 - (low_perc/100) ),
                 na.rm = T,
                 type = 7)
  xd <- ""
  if(sdy == "SDY67" && name == "combined"){
    xd <- sapply(df[[input_col]], FUN = function(x){
      if(is.na(x)){
        return(NaN)
      }else if(round(x, digits = 7) < round(xq[[1]], digits = 7)){
        return(0)
      }else if(round(x, digits = 7) >= round(xq[[2]], digits = 7)){
        return(2)
      }else{
        return(1)
      }
    })
  }else if(sdy == "SDY404" && name == "young"){
    xd <- sapply(df[[input_col]], FUN = function(x){
      if(is.na(x)){
        return(NaN)
      }else if(round(x, digits = 7) < round(xq[[1]], digits = 7)){
        return(0)
      }else if(round(x, digits = 7) > round(xq[[2]], digits = 7)){
        return(2)
      }else{
        return(1)
      }
    })
  }else{
    xd <- sapply(df[[input_col]], FUN = function(x){
      if(is.na(x)){
        return(NaN)
      }else if(round(x, digits = 7) <= round(xq[[1]], digits = 7)){
        return(0)
      }else if(round(x, digits = 7) >= round(xq[[2]], digits = 7)){
        return(2)
      }else{
        return(1)
      }
    })
  }
  return(xd)
}

# Original discretize function from Yuri Kotliarov
# NOTE: Does not return original results for $fc_res_max_d30
# discretize <- function(df, input_col, low_perc) {
#   x <- df[[input_col]]
#   xq = quantile(x, c((low_perc/100), 1 - (low_perc/100)), na.rm=T)
#   xd = ifelse(is.na(x),NA,1)
#   xd[x<=xq[1]] = 0
#   xd[x>=xq[2]] = 2
#   return(xd)
# }

# for splitting participant_id string to remove sdy info
sub_split <- function(x){
  tmp <- (strsplit(x, split = "[.]"))
  return(tmp[[1]][1])
}

# Drop columns by name
drop_cols <- function(df, cols_to_drop){
  df <- df[ , !(names(df) %in% cols_to_drop)]
  return(df)
}

#-----Main method--------------

#' Function to generate HAI data table from ImmuneSpace Connection
#'
#' @param rawdata dataframe of raw HAI data
#' @param sdy ImmuneSpace Study Name
#' @export
adjust_hai <- function(sdy, rawdata){

  res_list <- list() # holder for output of 3 df for combined, young, old hai data

  # Get strain names, which are unique variations for each study unfortunately
  strains <- list()
  if(sdy == "SDY80"){
    strains <- c("H1N1", "A_Brisbane", "A_Uruguay","B_Brisbane")
  }else{
    strains <- unique(rawdata$virus)
    strains <- sapply(strains, FUN = function(x){
      x <- gsub("\\.|\\/| |-|\\(|\\)", "_", x)
      return(x)
    })
  }

  # Setup all colnames
  std_names <- nm_builder("d0_std_norm_", strains)
  std_names <- c(std_names, nm_builder("fc_std_norm_", strains))

  d28_names <- nm_builder("d28_", strains)
  d0_names <- nm_builder("d0_", strains)
  fc_names <- nm_builder("fc_", strains)

  keep_names <- c("d0_norm_max","fc_norm_max","fc_norm_max_ivt", "d0_max",
                  "fc_max","fc_max_4fc", "fc_norm_max_d20","fc_norm_max_d30",
                  "fc_res_max","fc_res_max_d20","fc_res_max_d30")

  if(sdy == "SDY80"){
    allnms <- c("cohort", std_names, keep_names)
    tmpmat <- matrix(NA, ncol = 20, nrow = 64)
    colnames(tmpmat) <- allnms

    # This removes all subjects without original demo information, 223 has no IS sub id
    # and is removed post-calculations
    rawdata <- rawdata [ which(
      (rawdata$subject >=200 & rawdata$subject <= 284)
      ), ]

    titer_data <- data.frame(rawdata, tmpmat)

    # Adding cohort definition from Yuri comments re: code
    old_subs <- c(212, 229, 232, 233, 244, 245, 250, 251, 260, 261, 273, 277, 280)
    titer_data$cohort <- ifelse(titer_data$subject %in% old_subs, "Old", "Young")

    subids <- titer_data$subject
    cohorts <- c("Young", "Old")

  }else{
    # Generate vectors from rawdata or instantiate variables based on
    # nomenclature of original column headers
    subids <- unique(rawdata$participant_id)
    cohorts <- unique(rawdata$cohort)
    days_collected <- c(0,28)

    allnms <- c("subject","cohort", d0_names, d28_names, fc_names, std_names, keep_names)

    titer_data <- data.frame(matrix(vector(),
                                    nrow = 0,
                                    ncol = length(allnms)),
                             stringsAsFactors = F)
    colnames(titer_data) <- allnms

    # Parse day 0 (initial) and day 28 (follow-up) titer data into df as basis for all future calculations.
    # NOTE 1: Because the follow-up titer measurement is not always collected on day 28 exactly,
    # it is assumed that any value greater than 0 represents the 28th day value if there is
    # not a day 28 value present.
    iterator <- 1
    for(id in subids){
      sub_data <- rawdata[which(rawdata$participant_id == id),]
      valid <- sub_check(sub_data, names(strains))
      if(valid){
        titer_data[iterator,1] <- id
        cohort_df <- rawdata[which(rawdata$participant_id == id),]
        cohort_val <- unique(cohort_df$cohort)
        titer_data[iterator,2] <- cohort_val
        for(vir_name in names(strains)){
          for(day in days_collected){
            col_to_find <- paste0("d", day, "_", strains[[vir_name]])
            rowid <- which(rawdata$participant_id == id &
                             rawdata$study_time_collected == day &
                             rawdata$virus == vir_name)
            if(length(rowid) == 0){
              rowid <- which(rawdata$participant_id == id &
                               rawdata$study_time_collected > 0 &
                               rawdata$virus == vir_name)
            }
            if(length(rowid) > 1){
              rowid <- rowid[[1]]
            }
            target_row <- rawdata[rowid, ]
            titer_data[iterator, col_to_find] <- as.integer(target_row$value_reported)
          }
        }
        iterator <- iterator + 1
      }
    }
  }

  # calc fold change
  for(i in 1:length(strains)){
    titer_data[, fc_names[[i]] ] <- titer_data[ , d28_names[[i]] ] / titer_data[ , d0_names[[i]] ]
  }

  glob_vals <- list()
  if(sdy != "SDY80"){
    # setup list to hold median and sd values for later use in calculations
    glob_names <- c("d0_med", "d0_sd","fc_med", "fc_sd")
    li <- vector("list",length = length(strains))
    names(li) <- strains

    for(l in glob_names){
      glob_vals[[l]] <- li
    }

    # calc median and sd for d0 cols
    glob_vals <- med_sd_calc("d0_", strains , glob_vals, titer_data)

    # calc fold change med and sd
    glob_vals <- med_sd_calc("fc_", strains , glob_vals, titer_data)
  }

  # calc standardized and normalized value for each possibility of (d0,fc) x (strains)
  # sd used in all non-SDY80 studies because more than 50% of subjects were defined
  # as non-responders which causes denominator to be zero
  opts <- c("d0_","fc_")
  for(ver in opts){
    for(virus in strains){
      std_norm_col <- paste0(ver, "std_norm_", virus)
      tcol <- titer_data[[paste0(ver, virus)]]
      if(sdy == "SDY80"){
        robmed <- robust_med(tcol)
        titer_data[std_norm_col] <- (tcol - robmed) / r_mad(tcol)
        rep_ls <- titer_data[std_norm_col] == Inf
        titer_data[std_norm_col][rep_ls] <- NaN
      }else{
        colmed <- glob_vals[[paste0(ver,"med")]][[virus]]
        colsd <- glob_vals[[paste0(ver,"sd")]][[virus]]
        titer_data[std_norm_col] <- (tcol - colmed) / colsd
      }
    }
  }

  # Assign snapshots of variables to global_env for use with mapply.
  assign("gl_tdata", titer_data, envir = .GlobalEnv)
  assign("gl_strains", strains, envir = .GlobalEnv)

  # Select maxima for (d0,fc) x ("","std_norm") columns
  for(ver in opts){
    titer_data[[paste0(ver,"max")]] <- mapply(max_select, titer_data$subject, ver)
    titer_data[[paste0(ver,"norm_max")]] <- mapply(max_select, titer_data$subject, paste0(ver,"std_norm_"))
  }

  # determine fc_max_4fc, which is categorization based on fold change > 4
  titer_data$fc_max_4fc <- unlist(lapply(titer_data$fc_max, FUN = function(x){
    if(x > 4){return(1)}else{return(0)}
    }))

  # Inverse normal transformation of standardized/normalized max fold change column
  # Done by quantile normalization on a modified ranking of observations
  # fc_norm_max_ivt << provided from Yuri Kotliarov
  ranked <- rank(titer_data$fc_norm_max, na.last = "keep")
  P <- ranked / (sum(!is.na(ranked)) + 1) # +1 needed to avoid 0,1 values that generate inf, -inf
  titer_data$fc_norm_max_ivt <- qnorm(P)

  # Need to remove SDY80 subjects that have low flow results, NA for B-Brisbane, or not in original results
  if(sdy == "SDY80"){
    low_flow <- c(206, 226, 243, 247, 249, 252, 254, 263, 270, 275, 281, 282)
    titer_data <- titer_data[ which(!(titer_data$subject %in% low_flow)), ]
  }

  # setup meta-list for possible titer tables based on cohorts: young, old, and combined are possible
  submxs <- list()

  # Generate subset matrices based on age and perform statistical work on each separately
  # SDY212 and the other studies use different nomenclature for categorization, therefore
  # need to check against lists
  yng_ls <- c("Cohort_1", "Young adults 21-30 years old", "Young")
  old_ls <- c("Cohort_2", "Cohort2", "Older adults >= 65 years old", "healthy adults, 50-74 yo", "Old")

  titer_data$Age.class <- sapply(titer_data$cohort, FUN = function(x){
    if(x %in% yng_ls){
      return("young")
    }else if(x %in% old_ls){
      return("old")
    }else{
      return(NA)
    }
  })

  # SDY67-old is only subjects over 60 yrs old.  These subject IDs are from the demographic file.
  if(sdy == "SDY67"){
    subs_keep <-  c("SUB113453", "SUB113458", "SUB113460", "SUB113461", "SUB113463", "SUB113464",
                    "SUB113466", "SUB113467", "SUB113468", "SUB113470", "SUB113471", "SUB113486",
                    "SUB113492", "SUB113493", "SUB113494", "SUB113495", "SUB113496", "SUB113497",
                    "SUB113499", "SUB113500", "SUB113501", "SUB113502", "SUB113503", "SUB113504",
                    "SUB113505", "SUB113508", "SUB113509", "SUB113510", "SUB113512", "SUB113516",
                    "SUB113517", "SUB113525", "SUB113526", "SUB113527", "SUB113528", "SUB113529",
                    "SUB113532", "SUB113533", "SUB113535", "SUB113536", "SUB113540", "SUB113541",
                    "SUB113543", "SUB113545", "SUB113546", "SUB113547", "SUB113549", "SUB113553",
                    "SUB113555", "SUB113556", "SUB113564", "SUB113571", "SUB113572", "SUB113574",
                    "SUB113575", "SUB113577", "SUB113578", "SUB113580", "SUB113586", "SUB113593",
                    "SUB113594", "SUB113595", "SUB113596", "SUB113599", "SUB113602", "SUB113603",
                    "SUB113604", "SUB113605", "SUB113606")
    subs_keep <- sapply(subs_keep, FUN = function(x){x <- paste0(x,".67")})
    titer_data <- titer_data[(titer_data$subject %in% subs_keep),]
  }

  submxs[["combined"]] <- titer_data

  for(coh in cohorts){
    if(coh %in% yng_ls){
      submxs[["young"]]  <- titer_data[which(titer_data$cohort == coh),]
    }else if(coh %in% old_ls){
      submxs[["old"]]  <- titer_data[which(titer_data$cohort == coh),]
    }
  }

  # fc_res_max is generated by first binning the subjects' d0_norm_max data
  # according to a manual selection and then normalizing/standardizing
  # the inverse normal transformation values within those bins.
  # This code is based on Yuri Kotliarov's work.
  SDY404_bins <- list(c(1,5),c(0.25,1,5),c())
  SDY212_bins <- list(c(0.01),c(0.05),c(0.01))
  # Although SDY400 only shows bins 1,5 on google drive png, it has two NaN vals indicative of
  # single member bins, therefore third bin was added (6) in young / combined. Old appeared to
  # need very different bins.
  SDY400_bins <- list(c(1,5,6),c(1,5,6),c(0.1,3))
  SDY63_bins <- list(c(2,6),c(2),c(2))
  # SDY67 / 80 not derived from pngs in HIPC google drive, but rather determined by following method:
  # file <- fread("original_hai_tbl")
  # df <- data.frame(file$fc_res_max,file$fc_norm_max_ivt,file$d0_norm_max)
  # df <- df[order($fc_norm_max_ivt,$d0_norm_max),]
  # By looking at where ivt is the same but fc_res_max was different, one can guess that they
  # were in different bins and estimate the cutoff points by looking at the d0_norm_max value.
  SDY67_bins <- list(c(0.1,0.5,1.5,4,6),c(),c(0.1,0.5,1.5,6))
  SDY80_bins <- list(c(-10,1,50),c(-10,1,50),c()) # From Yuri's code

  bname <- paste0(sdy,"_bins")
  bins <- get(bname)
  names(bins) <- c("combined", "young", "old")

  for(name in names(submxs)){
    df <- submxs[[name]]
    tmp_mad <- r_mad(df$fc_norm_max_ivt)
    df$bin <- cut(df$d0_norm_max,
                  breaks = c(-Inf, bins[[name]], Inf),
                  labels = 1:( length(bins[[name]]) + 1) )

    # need count of bin members to force NA value for bins with < 3 members.
    # Original code was in Matlab and may not have generated median / sd vals for < 3 groups.
    df = df %>%
      group_by(bin) %>%
      dplyr::mutate(count = n()) %>% # MUST name 'dplyr' in front of mutate to get dplyr::n()
      ungroup()

    if(sdy == "SDY80"){
      df = df %>%
        group_by(bin) %>%
        mutate(fc_res_max = (fc_norm_max_ivt - median(fc_norm_max_ivt)) / r_mad(fc_norm_max_ivt)) %>%
        ungroup()

      df$fc_res_max <- ifelse(df$fc_res_max %in% c(Inf, NA), NaN, df$fc_res_max)

    }else{
      df = df %>%
        group_by(bin) %>%
        mutate(fc_res_max =
                 ifelse( count > 2 | sdy == "SDY63", # SDY63 Old allows a 2 count group to be processed
                         (fc_norm_max_ivt - median(fc_norm_max_ivt, na.rm=T)) / sd(fc_norm_max_ivt, na.rm=T),
                         NaN)) %>%
        ungroup()
    }


    # discretize for all combinations of ("fc_norm_max","fc_res_max") x ("d20","d30")
    in_cols <- c("fc_norm_max", "fc_res_max")
    in_percs <- c(20, 30)
    for(cl in in_cols){
      for(perc in in_percs){
        targ_col <- paste0(cl,"_d", perc)
        df[[targ_col]] <- discretize(df, cl, perc, sdy, name)
      }
    }

    # Remove d28, bin, and cohort columns b/c not present in results of original manual versions.
    df <- drop_cols(df, c(d0_names, d28_names, fc_names, std_names, "bin", "cohort", "count"))

    if(sdy == "SDY80"){
      id_hsh <- hash(SDY80_IDmap$bioSampleID, SDY80_IDmap$participantID)
      df$subject <- sapply(df$subject, FUN = function(x){
        val <- id_hsh[[as.character(x)]]
        return(ifelse(is.null(val), 223, val))
      })
      ifelse(name != "old", df <- df[-which(df$subject == 223), ] , df <- df ) # 223 cannot be mapped to IS subject
    }

    df$subject <- unlist(lapply(df$subject, sub_split))

    if(name %in% c("young","old")){
      tmp_nms <- paste0(name,"_",colnames(df))
      tmp_nms[1] <- "subject"
      colnames(df) <- tmp_nms
    }
    res_list[[name]] <- df
  }

  # combine all three into one df (should just be subject and keep names columns here)
  merged_df <- merge(res_list$combined, res_list$old, by = "subject", all.y = T, all.x = T)
  if(sdy != "SDY67"){
    merged_df <- merge(merged_df, res_list$young, by = "subject", all.y = T, all.x = T)
    merged_df <- drop_cols(merged_df, c("young_Age.class", "old_Age.class"))
  }else{
    merged_df <- drop_cols(merged_df, c("old_Age.class"))
  }

  return(merged_df)
}



