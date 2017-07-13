

# Dependencies
library(ImmuneSpaceR)
library(data.table)
library(tidyr)

# Get ImmuneSpace data
con <- CreateConnection("SDY80")
isData <- con$getDataset("neut_ab_titer")
isData <- data.frame(isSubject = isData$participant_id,
                     day = isData$study_time_collected,
                     isStrain = isData$virus,
                     value = isData$value_reported)

# Load rawdata from Yuri
pathToData <- "/home/ehenrich/R/CHI_titer_calc/chi_data/"
d0 <- "day0_d0_d70match_08_2013_raw.txt"
d70 <- "day70_d0_d70match_08_2013_raw.txt"
d0Data <- fread(file.path(pathToData, d0))
d70Data <- fread(file.path(pathToData, d70))



# load SDY80 id map
idMap <- read.table("/home/ehenrich/R/ImmSig2/OrigCode/SDY80_IDmap.tsv", 
                    header = T,
                    stringsAsFactors = F,
                    sep = "\t")
idMap$bioSampleID <- as.character(idMap$bioSampleID)

# map strain names
chiNames <- c("swine",
              "a_brisbane", 
              "uruguay", 
              "b_brisbane")
isNames <- c("A_Ca_07_swine", 
             "A/Brisbane/59/2007" , 
             "A_Uruguay_716_2007", 
             "B_Brisbane_60_2001")
strainNms <- data.frame(chiNames, isNames)

# Helper for updating things
updateHai <- function(df, day, strainNms){
  df$sub <- as.character(df$sub)
  newDf <- data.table(gather(df, strain, value, 2:5))
  newDf[idMap, isSubject := participantID, on=c(sub = "bioSampleID")]
  newDf$day <- day
  newDf[strainNms, isStrain := isNames, on=c(strain = "chiNames")]
  return(newDf)
}

# Make CHI data match IS
chiD0 <- updateHai(d0Data, 0)
chiD70 <- updateHai(d70Data, 70)
chiData <- rbind(chiD0, chiD70)

# Subjects removed prior to calculations
older <- c(212, 229, 232, 233, 244, 245, 250, 251, 260, 261, 273, 277, 280)
noDemoLow <- 200
noDemoHigh <- 284
chiData$sub <- as.numeric(chiData$sub)
chiData <- chiData[!(chiData$sub %in% older), ]
chiData <- chiData[(chiData$sub < noDemoHigh & chiData$sub > noDemoLow ), ]

# Compare data - observations not found at all in ImmuneSpace
notSubMap <- chiData[(is.na(chiData$isSubject)), ] # just sub 223

# rm chiData cols and 223 row to prep for comparison
chiData <- chiData[ !(is.na(chiData$isSubject)), ]
chiData <- chiData[,-(1:2)]

chiData$isVal <- apply(chiData, 1, FUN = function(row){
  isRow <- isData[(isData$isSubject == row[2] &
                   isData$day == row[3] &
                   isData$isStrain == row[4]), ]
  if(dim(isRow)[1] != 0){
    return(isRow$value)
  }else{
    return(NA)
  }
})

# 

