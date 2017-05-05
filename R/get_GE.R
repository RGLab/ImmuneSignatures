
# The following code was used to generate the expression matrix (EM) for each study that
# is included as data in the ImmSig2 package.  The EMs are pre-packaged as data to save time
# in running the report on ImmuneSpace.  However, a user may run this code by hand to see how
# the information was extracted from various online sources according to the methods of the
# the original work by the HIPC ImmuneSignatures collaborators.  It is important to note that
# the original annotation is kept even though it has different geneSymbol mapping than the
# official probe-to-geneSymbol mapping in the various bioconductor packages (e.g. illuminaHumav4).
# The reason for this is that the use of more current mappings alters the pathways selected
# in the meta-analysis.

# # Expression Matrix Creation ---------------------------------------------------------------
# library(ImmuneSpaceR)
# library(httr)
# library(RCurl)
# library(data.table)
# library(hash)
# library(preprocessCore)
# library(GEOquery)
# library(plyr)
# library(ImmSig2) # to get the biosmpl2sub table for mapping colnames
# library(DESeq) # for SDY67 normalizaion
#
# #*************************************************************************************
#
# #-----------Setup directory vars------------------------------------------
# rawdata_dir <- "/home/ehenrich/R/Gen_Test/ImmSig2_work/raw/"
#
# #------HELPER METHODS ------------------------------------------------------
# # Get gene expression file names from ImmuneSpace and use only baseline data
# get_gef <- function(con, sdy){
#   gef <- con$getDataset("gene_expression_files")
#   if(sdy != "SDY80"){
#     gef <- gef[ which(gef$file_info_name != "NA" & gef$study_time_collected == 0),  ]
#   }
#   return(gef)
# }
#
# # For SDY212 / SDY67 pull from ImmuneSpace via RCurl to use netrc instead of Auth with httr::GET
# make_handle <- function(con){
#   opts <- con$config$curlOptions
#   opts$netrc <- 1L
#   handle <- getCurlHandle(.opts = opts)
#   return(handle)
# }
#
# # download microarray or gene expression files from ImmuneSpace and return list of filenames
# dl_IS_GE <- function(handle, ge_flnms, sdy, rawdata_dir){
#   fl_list <- llply(ge_flnms, .fun = function(x){
#     func <- basicTextGatherer()
#     link <- paste0("https://www.immunespace.org/_webdav/Studies/",
#                    sdy,
#                    "/%40files/rawdata/gene_expression/",
#                    x)
#
#     curlPerform(url = link, curl = handle, writefunction = func$update)
#     fl <- file.path(rawdata_dir, x)
#     write(func$value(), file = fl)
#     return(fl)
#   })
#   return(fl_list)
# }
#
# get_GE_vals <- function(fl_list, subids){
#   vals_list <- lapply(fl_list, FUN = function(x){
#     vals <- read.table(x, header = T, stringsAsFactors = F, sep = "\t", fill = T)
#     vals <- vals[order(vals[,1]),]
#     vals <- vals[,2]
#     return(vals)
#   })
#   names(vals_list) <- subids
#   return(vals_list)
# }
#
# # download gene expression matrix file from GEO - faster than GEOquery::getGSEtable()
# get_GSE_files <- function(sdy, rawdata_dir){
#   link <- ""
#   if(sdy == "SDY63"){
#     link <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59635/suppl/GSE59635_non-normalized.txt.gz"
#   }else if(sdy == "SDY404"){
#     link <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59654/suppl/GSE59654_non-normalized.txt.gz"
#   }else if(sdy == "SDY400"){
#     link <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59743/suppl/GSE59743_series_matrix.txt.gz"
#   }else if(sdy == "SDY80"){
#     link <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE47nnn/GSE47353/matrix/GSE47353_series_matrix.txt.gz"
#   }
#
#   inputFile <- paste0(rawdata_dir, "/", sdy, ".txt.gz")
#   suppressMessages(GET(url = link, write_disk(inputFile, overwrite = TRUE)))
#   GEOquery::gunzip(inputFile, overwrite = TRUE)
#   inputFile <- paste0(rawdata_dir, "/", sdy, ".txt")
#   return(inputFile)
# }
#
# # Significantly faster than using data.table or which statements for searching paired-values
# # NOTE 1: The key must be presented as a character; it does not auto-cast.
# # NOTE 2: dflt_ret needed for for Yale Studies.
# hashmap <- function(keys, values, targ_vec, dflt_ret = NA){
#   tmp_hash <- hash(keys, values)
#   mapped_res <- sapply(targ_vec, FUN = function(x){
#     val <- tmp_hash[[as.character(x)]]
#     val <- ifelse(is.null(val), dflt_ret, val)
#     return(val)
#   })
#   return(mapped_res)
# }
#
# # Generate subject ID values that can be used downstream and mapped to expr values
# map_headers <- function(sdy, exprs_headers, is_col, smpl_col){
#   if(sdy == "SDY212" | sdy == "SDY67"){
#     id_map <- biosmpl2sub
#   } else {
#     id_map <- get(paste0(sdy, "_IDmap"))
#   }
#
#   if(sdy %in% c("SDY212","SDY67", "SDY80")){
#     mapped_headers <- hashmap(id_map[[smpl_col]], id_map[[is_col]], exprs_headers)
#   }else if(sdy %in% c("SDY63","SDY400","SDY404")){
#     mapped_headers <- sapply(exprs_headers, FUN = function(x){
#       spl_x <- strsplit(x, split = "_")
#       if(sdy == "SDY63"){
#         # sample ids were off according to original code, therefore correct
#         spl_x[[1]][1] = as.integer(spl_x[[1]][1]) + 91000
#       }
#       sub_no_date <- id_map[[is_col]][which(id_map[[smpl_col]] == spl_x[[1]][1])]
#       res <- paste0(sub_no_date,"_d",spl_x[[1]][3])
#       return(res)
#     })
#   }
#
#   # for the case where two instances / columns of one participant "SUB134307"
#   if(sdy == "SDY212"){
#     dbl_subid <- which(mapped_headers == "SUB134307_d0")
#     mapped_headers[dbl_subid[1]] <- "SUB134307.1_d0"
#     mapped_headers[dbl_subid[2]] <- "SUB134307.2_d0"
#   }
#   return(mapped_headers)
# }
#
# # log2 transform, quantile normalize raw expression values and map sample to IS ids
# norm_map_mx <- function(sdy, exprs, is_col, smpl_col){
#   cnames <- colnames(exprs)
#   norm_exprs <- preprocessCore::normalize.quantiles(log2(as.matrix(exprs)))
#   colnames(norm_exprs) <- map_headers(sdy, cnames, is_col, smpl_col)
#   return(norm_exprs)
# }
#
# # Remove subjects not found in original files.  Unclear why subjects were removed.
# remove_subs <- function(expr_matrix, subs_to_rm){
#   expr_matrix <- expr_matrix[ , !(colnames(expr_matrix) %in% subs_to_rm)]
#   return(expr_matrix)
# }
#
# # original main method
# makeGE <- function(sdy, rawdata_dir){
#
#   # Get raw data and process it uniquely for each study
#   if(sdy == "SDY212" | sdy == "SDY67"){
#     # build links from Immunespace and then download
#     con <- CreateConnection(sdy)
#     gef <- get_gef(con, sdy)
#     handle <- make_handle(con)
#     fl_list <- dl_IS_GE(handle, unique(gef$file_info_name), sdy, rawdata_dir)
#
#     if(sdy == "SDY212"){
#       # Clean and Prep
#       rawdata <- read.table(fl_list[[1]], header = T, stringsAsFactors = F, sep = "\t", fill = T)
#
#       # remove gene with NA values. These NAs are found in ImmPort as well.
#       # ImmPort staff (Patrick.Dunn@nih.gov) said this was an unresolved issue
#       # with ticket "HDB-13" last discussed in October 2015 with study authors.
#       rawdata <- rawdata[ -c(which(rawdata[,"PROBE_ID"] == "ILMN_2137536")), ]
#       probe_ids <- rawdata[, "PROBE_ID"]
#       gene_syms <- rawdata[ , "SYMBOL"]
#       sigcols <- grep("Signal", colnames(rawdata), value = TRUE)
#       rawdata <- rawdata[ , sigcols ]
#       setnames(rawdata, gsub(".AVG.*$", "", colnames(rawdata)))
#
#       # normalize, then remove subjects not found in original file
#       final_expr_vals <- norm_map_mx(sdy, rawdata, "GE_ID", "Biosample")
#       final_expr_vals <- remove_subs(final_expr_vals, c("SUB134242_d0","SUB134267_d0"))
#
#     }else if(sdy == "SDY67"){
#       # Build correct subids for colnames to have both ID and day value
#       subids <- unname(sapply(gef$file_info_name, FUN = function(x){
#         fname <- basename(x)
#         targ_row <- gef[which(gef$file_info_name == fname),]
#         subid <- substr(targ_row$participant_id[[1]],1,9)
#         day_val <- as.integer(targ_row$study_time_collected)
#         subid <- paste0(subid, "_d", day_val)
#         return(subid)
#       }))
#
#       raw_smpl <- read.table(fl_list[[1]], header = T, stringsAsFactors = F, sep = "\t", fill = T)
#       gene_syms <- sort(toupper(raw_smpl[,1]))
#       probe_ids <- seq(from = 1, to = length(gene_syms), by = 1) # Like HIPC manuscript, probe_ids are row number
#
#       val_list <- get_GE_vals(fl_list, subids)
#       rawdata <- quickdf(val_list)
#
#       # these subjects were removed from original file because they were not processed
#       # at the same time as the others and had a significant batch effect.
      # subs_rm <- c("SUB113458_d0","SUB113463_d0","SUB113470_d0","SUB113473_d0","SUB113474_d0",
      #              "SUB113476_d0","SUB113483_d0","SUB113487_d0","SUB113490_d0","SUB113494_d0",
      #              "SUB113495_d0","SUB113496_d0","SUB113498_d0","SUB113504_d0","SUB113505_d0",
      #              "SUB113513_d0","SUB113514_d0","SUB113524_d0","SUB113526_d0","SUB113527_d0",
      #              "SUB113532_d0","SUB113535_d0","SUB113537_d0","SUB113545_d0","SUB113548_d0",
      #              "SUB113555_d0","SUB113558_d0","SUB113559_d0","SUB113561_d0","SUB113566_d0",
      #              "SUB113567_d0","SUB113568_d0","SUB113571_d0","SUB113572_d0","SUB113582_d0",
      #              "SUB113583_d0","SUB113588_d0","SUB113595_d0","SUB113610_d0")
      #
      # # Normalized using DESeq b/c RNAseq data not microarray like other studies
      # countTable <- exmx <- remove_subs(rawdata, subs_rm)
      # condition <- colnames(exmx)
      # cds <- newCountDataSet(countTable, condition)
      # cds <- estimateSizeFactors(cds) ## estimate size factor
      # cdsBlind <- estimateDispersions(cds, method="blind" )
      # vsd <- varianceStabilizingTransformation(cdsBlind)
      # final_expr_vals <- exprs(vsd)                   # normalized data into `e`
#     }
#
#   }else if(sdy %in% c("SDY63","SDY404","SDY400","SDY80")){
#
#     inputFiles <- get_GSE_files(sdy, rawdata_dir)
#
#     if(sdy != "SDY80"){
#
#       if(sdy == "SDY400"){
#         rawdata <- read.table(inputFiles, skip = 64, header = T, fill = T)
#         rawdata <- rawdata[ 1:(dim(rawdata)[1] - 1), ] # drop the last row with !series_matrix_table_end and NA vals
#         colnames(rawdata) <- read.table(inputFiles,
#                                         sep = "\t",
#                                         header = F,
#                                         stringsAsFactors = F,
#                                         fill = T,
#                                         nrows = 1,
#                                         skip = 31)
#         colnames(rawdata)[1] <- "ID_REF"
#
#       }else{
#         rawdata <- as.data.frame(fread(inputFiles))
#       }
#
#       colnames(SDY63_anno_tbl) <- c("probeIDs","geneSymbol")
#       keys <- SDY63_anno_tbl$probeIDs
#       values <- SDY63_anno_tbl$geneSymbol
#
#       # Clean and Prep
#       probe_ids <- rawdata$ID_REF
#       gene_syms <- hashmap(keys,
#                            values,
#                            probe_ids,
#                            dflt_ret = "")
#
#       # SDY400 has different headers
#       rawdata <- rawdata[ , grepl("PBMC", names(rawdata))]
#
#       if(sdy == "SDY400"){
#         final_expr_vals <- rawdata
#         colnames(final_expr_vals) <- unname(map_headers(sdy,
#                                                         colnames(rawdata),
#                                                         "Sub.Org.Accession",
#                                                         "User.Defined.ID"))
#
#       }else{
#         final_expr_vals <- norm_map_mx(sdy, rawdata, "Sub.Org.Accession", "User.Defined.ID")
#       }
#
#       # SDY404 ids from Immport have actual day the sample was taken (e.g. 4 instead of 2) according to
#       # personal communication from Stefan Avey at Yale. However this creates conflicts with the original
#       # file that had standard days (e.g. 2 and 28). Therefore, changing here to ensure reproducibility
#       # of original information.
#       if(sdy == "SDY404"){
#         colnms <- colnames(final_expr_vals)
#         bad_ids <-  c("SUB120470_d24", "SUB120472_d32", "SUB120480_d4", "SUB120481_d4", "SUB120483_d4",
#                       "SUB120484_d4", "SUB120485_d4", "SUB120487_d4", "SUB120488_d4")
#         good_ids <- c("SUB120470_d28", "SUB120472_d28", "SUB120480_d2", "SUB120481_d2", "SUB120483_d2",
#                       "SUB120484_d2", "SUB120485_d2", "SUB120487_d2", "SUB120488_d2")
#         id_hash <- hash(bad_ids, good_ids)
#         colnms <- sapply(colnms, FUN = function(x){
#           if(has.key(x, id_hash)){
#             return(id_hash[[x]])
#           }else{
#             return(x)
#           }
#         })
#         colnames(final_expr_vals) <- colnms
#       }
#
#     }else if(sdy == "SDY80"){
#       rawdata <- as.data.frame(fread(inputFiles, skip = 96, header = T, fill = T))
#       probe_ids <- rawdata$ID_REF
#       gene_syms <- hashmap(SDY80_orig_anno$probeID, SDY80_orig_anno$gs, probe_ids)
#
#       # To mimic original file, I remove all rows that were not
#       # successfully mapped, which is approximately 45%.  Unsure why this was done.
#       gene_syms <- gene_syms[which(sapply(gene_syms, FUN = function(x){!is.na(x)}))]
#       probe_ids <- probe_ids[which(probe_ids %in% names(gene_syms))]
#       rawdata <- rawdata[ which(rawdata$ID_REF %in% names(gene_syms)), ]
#       rawdata <- rawdata[,-1]
#       colnames(rawdata) <- unname(map_headers(sdy,
#                                               colnames(rawdata),
#                                               is_col = "GE_name",
#                                               smpl_col = "GSM"))
#       final_expr_vals <- rawdata
#     }
#   }
#
#   em <- as.data.frame(final_expr_vals)
#   em$geneSymbol <- gene_syms
#   rownames(em) <- probe_ids
#   return(em)
# }

# MAIN METHOD -----------------------------------------------------------------------
#' Function to return the prepackaged expression matrix held in library. Kept
#' outside of Rmd to allow for original code to be shown above.
#'
#' @import data.table Rlabkey
#' @param sdy name of study
#' @export
#'
getGE <- function(sdy){
  em <- get(paste0(tolower(sdy), "_exprMx"))
}
