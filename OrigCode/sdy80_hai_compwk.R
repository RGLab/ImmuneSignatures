

# Creating the output table from the RAW IS data
library(ImmuneSpaceR)
library(dplyr)
library(hash)

con <- CreateConnection("SDY80")
rawdata <- con$getDataset("neut_ab_titer")
rawdata <- as_tibble(rawdata)
rawdata <- select(rawdata, participant_id, study_time_collected, virus, value_reported)

ugh <- matrix(NA, ncol=4, nrow=688)
ugh <- data.frame(ugh)
names(ugh) <- c("H1N1", "A_Brisbane", "B_Brisbane", "A_Uruguay")
just2 <- select(rawdata, participant_id, study_time_collected)
output <- cbind(just2, ugh)
output <- output[ !duplicated(output), ]

rawkeys <- c(  "A_Uruguay_716_2007", 
              "A_Ca_07_swine",
              "A/Brisbane/59/2007",
              "B_Brisbane_60_2001")
rawvals <- c("A_Uruguay","H1N1","A_Brisbane","B_Brisbane")
rawhash <- hash(rawkeys,rawvals)

for(i in 1:nrow(rawdata) ){
  x <- rawdata[i, ]
  sub <- x[1][[1]]
  day <- x[2][[1]]
  virus <- x[3][[1]]
  value <- x[4][[1]]
  targ_row <- which(output$participant_id == sub & output$study_time_collected == day)
  col_nm <- which(colnames(output) == rawhash[[ virus ]])
  output[ targ_row, col_nm ] <- value 
}

# cleanup, only use day 0 and day 7
output <- output[ !(output$study_time_collected %in% c(-7,70) ), ]
chk <- output$participant_id[ duplicated(output$participant_id) ]
output <- output[ output$participant_id %in% chk, ]

cleaner <- function(df, day){
  df <- df[df$study_time_collected == day, ]
  df <- df[ , -(which(colnames(df) == "study_time_collected")) ]
  df <- df[ order(df$participant_id), ]
  tmp <- colnames(df)[2:5]
  if(day == 7){ day <- 28 }
  tmp <- paste0("d", day, "_", tmp)
  colnames(df) <- c("subject", tmp)
  return(df)
}

day0 <- cleaner(output, 0)
day7 <- cleaner(output, 7)
haiIS <- merge(day0, day7, by = "subject")

# Create similar table from Yuri Raw data
library(ImmSig2)
haiOr <- SDY80_rawtiterdata_v2
id_tbl <- data.frame(read.table("/home/ehenrich/R/ImmSig2/origCode/SDY80_IDmap.tsv",
                     stringsAsFactors = F,
                     header = T,
                     sep = "\t"))
haiOr$ISsub <- id_tbl$participantID[ match(haiOr$subject, id_tbl$bioSampleID) ]

# Note: this causes sub 223 to be dropped, which is used later
haiMatch <- haiOr[ !is.na(haiOr$ISsub), -(which(colnames(haiOr) == "subject")) ]
colnames(haiMatch)[9] <- "subject"

# compare original to IS data
x <- colnames(haiIS)
y <- colnames(haiMatch)
haiMatch <- haiMatch[order(haiMatch$subject), order(match(y,x)) ]
haiIS <- haiIS[order(haiIS$subject), ]

