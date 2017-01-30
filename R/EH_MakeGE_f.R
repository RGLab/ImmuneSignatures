# Evan Henrich
# January 2017
# Fred Hutchinson Cancer Research Center
# Gottardo Lab
# ehenrich@fredhutch.org

#--------GEN NOTES-------------------------------------------------------
# This script preprocesses the raw Gene Expression data from ImmuneSpace.org, ImmPort 
# or GEO into the format of GEMatrix.txt file used in the HIPC ImmuneSignatures meta 
# analysis. Much of the code is based on work done by Damian Fermin [dfermin.yale@gmail.com] 
# from Yale in 2014 that was used for the meta analysis paper, particularly the 
# statistical operations. 

#**************TESTING ONLY************************************************
# Comment out when running as part of full analysis
# setwd("/home/ehenrich/R/ImmSig_Testing//")
# 
# library(ImmuneSpaceR)
# library(httr)
# library(XML)
# library(xml2)
# library(R.utils)
# library(Rlabkey)
# library(tools)
# library(data.table)
# library(hash)
# source("https://bioconductor.org/biocLite.R")
# biocLite("preprocessCore")
# library(GEOquery)
# library(illuminaHumanv4.db)
# library(hugene10sttranscriptcluster.db)

#*************************************************************************************

#-----------Setup directory vars-------------------------------------------
# wk_dir <- getwd()
# pp_dir <- file.path(wk_dir,"PreProc")
# rawdata_dir <- file.path(pp_dir,"rawdata")

#------HELPER METHODS ------------------------------------------------------
# Get gene expression file names from ImmuneSpace and return as data frame
get_gef <- function(sdy){
  con <- CreateConnection(sdy)
  gef <- con$getDataset("gene_expression_files")
  if(sdy != "SDY80"){
    gef <- gef[ which(gef$file_info_name != "NA" & gef$study_time_collected == 0),  ]
  }
  return(gef)
}

# download microarray or gene expression files from ImmuneSpace and return list of file paths
get_is_files <- function(gef, sdy, rawdata_dir, user, pwd){
  inputFiles <- unique(gef$file_info_name)
  
  # download files from IS to local machine
  links <- paste0("https://www.immunespace.org/_webdav/Studies/",
                  sdy,
                  "/%40files/rawdata/gene_expression/", 
                  inputFiles)
  
  # only setting to an object to limit output to console
  dump <- sapply(links, FUN = function(x){
    GET(url = x, 
        write_disk(file.path(rawdata_dir,basename(x)), 
        overwrite = T), 
        authenticate(user,pwd))
    })
  inputFiles <- file.path(rawdata_dir, inputFiles)
  return(inputFiles)
}

# for SDY80, generate a list of tables by scraping GEO websites
get_geo_files <- function(gef){
  
  # list of accessions to pull
  geo_acc <- gef$geo_accession
  
  # create links with following format and pull data:
  # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GSM1147759&id=13600&db=GeoDb_blob98
  base <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc="
  joint <- "&id="
  tail <- "&db=GeoDb_blob98"
  
  tbl_list <- lapply(geo_acc, FUN = function(x){
    # generate id and link
    digits <- gsub("GSM", "", x)
    id <- as.integer(digits) - 1134159
    link <- paste0(base, x, joint, id, tail)
    
    # scrape link for content and parse to table, return with probe ID
    dump <- GET(link)
    dump_contents <- read_html(dump$content)
    dump_contents_list <- xmlToList(xmlParse(dump_contents))
    raw_table <- dump_contents_list$body$font$pre[[4]]
    tbl_con<- textConnection(raw_table)
    data <- read.table(tbl_con, sep = "\t", stringsAsFactors = F, header = F)
    colnames(data) <- c("probeID", x)
    return(data)
  })
  
  # return list of tables
  return(tbl_list)
}

# download gene expression files from immport for yale studies
get_immport_files <- function(sdy, rawdata_dir){
  link <- ""
  if(sdy == "SDY63"){
    link <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59635/suppl/GSE59635_non-normalized.txt.gz"
  }else if(sdy == "SDY404"){
    link <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59654/suppl/GSE59654_non-normalized.txt.gz"
  }else if(sdy == "SDY400"){
    print("data not available yet")
  }
  inputFile <- paste0(rawdata_dir, "/", sdy, ".txt.gz")
  dump <- GET(url = link, write_disk(inputFile, overwrite = TRUE))
  gunzip(inputFile, overwrite = TRUE)
  inputFile <- paste0(rawdata_dir, "/", sdy, ".txt")
  return(inputFile)
}

# Significantly faster than using data.table or which statements for searching paired-values
# NOTE 1: The key must be presented as a character; it does not auto-cast.
# NOTE 2: keep_names needed for SDY80 when matching probes and dflt_ret for Yale Studies.
hashmap <- function(keys, values, targ_vec, keep_names = FALSE, dflt_ret = NA){
  tmp_hash <- hash(keys, values)
  mapped_res <- sapply(targ_vec, FUN = function(x){
    val <- tmp_hash[[as.character(x)]]
    val <- ifelse(is.null(val), dflt_ret, val)
    if(keep_names){names(val) <- x}
    return(val)
  })
  return(mapped_res)
}

# Important to maintain has-many relationship of gene symbol to probes
inv_alias2probe <- function(a2p_list){
  a2p_list <- a2p_list[!is.na(a2p_list)]
  p2a_list <- setNames(unlist(a2p_list, use.names=F), rep(names(a2p_list), lengths(a2p_list)))
  return(p2a_list)
}

# Generate subject ID values that can be used downstream and mapped to expr values
map_headers <- function(sdy, exprs_headers, is_col, smpl_col){
  id_map <- get(paste0(sdy,"_IDmap"))
  
  if(sdy %in% c("SDY212","SDY67", "SDY80")){
    mapped_headers <- hashmap(id_map[[smpl_col]], id_map[[is_col]], exprs_headers)
  }else if(sdy %in% c("SDY63","SDY400","SDY404")){
    mapped_headers <- sapply(exprs_headers, FUN = function(x){
      spl_x <- strsplit(x, split = "_")
      if(sdy == "SDY63"){
        # sample ids were off according to original code, therefore correct
        spl_x[[1]][1] = as.integer(spl_x[[1]][1]) + 91000
      }
      sub_no_date <- id_map[[is_col]][which(id_map[[smpl_col]] == spl_x[[1]][1])]
      res <- paste0(sub_no_date,"_d",spl_x[[1]][3])
      return(res)
    })
  }
  
  # for the case where two instances / columns of one participant "SUB134307"
  if(sdy == "SDY212"){
    dbl_subid <- which(mapped_headers == "SUB134307_d0")
    mapped_headers[dbl_subid[1]] <- "SUB134307.1_d0"
    mapped_headers[dbl_subid[2]] <- "SUB134307.2_d0"
  }
  return(mapped_headers)
}

# log2 transform, quantile normalize raw expression values and map sample to IS ids
norm_map_mx <- function(sdy, exprs, is_col, smpl_col){
  cnames <- colnames(exprs)
  norm_exprs <- preprocessCore::normalize.quantiles(log2(as.matrix(exprs)))
  colnames(norm_exprs) <- map_headers(sdy, cnames, is_col, smpl_col)
  return(norm_exprs)
}

# Remove subjects not found in original files.  Unclear why subjects were removed.
remove_subs <- function(expr_matrix, subs_to_rm){
  expr_matrix <- expr_matrix[ , !(colnames(expr_matrix) %in% subs_to_rm)]
  return(expr_matrix)
}

# write out final file containing prove ids, gene symbols, and expression values
write_out <- function(probe_ids, gene_syms, norm_map_exprs, sdy, output_dir){
  out <- data.frame(
    probeID=probe_ids,
    geneSymbol=gene_syms,
    round(norm_map_exprs,6),
    stringsAsFactors=FALSE,
    row.names = NULL)
  
  write.table(out, 
              file = paste0(output_dir, "/" , sdy, ".GEMatrix.txt"), 
              sep = "\t", 
              quote=FALSE, 
              row.names=FALSE)
}


#---------MAIN METHOD--------------------------------------------------------
makeGE <- function(sdy, 
                   user, 
                   pwd, 
                   yale.anno = "original", 
                   sdy80.anno = "original", 
                   sdy80.norm = F,
                   output_dir,
                   rawdata_dir){
  
  # instantiate variables to hold data for output at end
  probe_ids <- ""
  gene_syms <- ""
  final_expr_vals <- ""
  
  # Get raw data and process it uniquely for each study
  if(sdy == "SDY212" | sdy == "SDY67"){
    # Get filenames from Immunespace and then download
    gef <- get_gef(sdy)
    inputFiles <- get_is_files(gef, sdy, rawdata_dir, user, pwd)
    
    if(sdy == "SDY212"){
      # Clean and Prep
      rawdata <- fread(inputFiles, header = TRUE)
      
      # remove gene with NA values. These NAs are found in ImmPort as well.
      # ImmPort staff (Patrick.Dunn@nih.gov) said this was an unresolved issue
      # with ticket "HDB-13" last discussed in October 2015 with study authors.
      rawdata <- rawdata[ -c(which(rawdata$PROBE_ID == "ILMN_2137536")),]
      probe_ids <- rawdata[, PROBE_ID]
      gene_syms <- rawdata[ , SYMBOL]
      sigcols <- grep("Signal", colnames(rawdata), value = TRUE)
      rawdata <- rawdata[, sigcols, with = FALSE]
      setnames(rawdata, gsub(".AVG.*$", "", colnames(rawdata)))
      
      # normalize, then remove subjects not found in original file
      final_expr_vals <- norm_map_mx(sdy, rawdata, "final_id", "BioSampleID")
      final_expr_vals <- remove_subs(final_expr_vals, c("SUB134242_d0","SUB134267_d0"))
      
    }else if(sdy == "SDY67"){
      # inputFiles here only contain 1 subject worth of data per file, therefore
      # need to iterate through and build a master file with all subjects to return.
      # NOTE: no normalization was done on the original file and that is mimicked her.
      gs_tbl <- fread(inputFiles[[1]])
      gs_tbl$GENE_SYMBOL <- toupper(gs_tbl$GENE_SYMBOL)
      gs_tbl <- gs_tbl[ order(GENE_SYMBOL),]
      gene_syms <- gs_tbl$GENE_SYMBOL
      
      final_expr_vals <- as.data.frame(do.call(cbind,lapply(inputFiles, FUN = function(x){
        fname <- basename(x)
        targ_row <- gef[which(gef$file_info_name == fname),]
        subid <- substr(targ_row$participant_id[[1]],1,9)
        day_coll <- as.integer(targ_row$study_time_collected)
        subid <- paste0(subid,"_d",day_coll)
        
        #read in the table and relabel colnames
        tmp <- fread(x)
        colnames(tmp) <- c("SYMBOL",subid)
        tmp$SYMBOL <- toupper(tmp$SYMBOL)
        tmp <- tmp[ order(SYMBOL),]
        
        # make sure the gene symbol vec for the curr subject matches original
        if(!all.equal(tmp$SYMBOL,gs_tbl$GENE_SYMBOL)){
          stop(paste0("Gene Symbols in ",fname,
                      " do not match those in first file. Please check before re-running."))
        }
        
        # pull out expr values and return as vec
        return(tmp[,2])
      })))
      
      # these subjects were removed from original file because they were not processed
      # at the same time as the others and had a significant batch effect.
      subs_rm <- c("SUB113458_d0","SUB113463_d0","SUB113470_d0","SUB113473_d0","SUB113474_d0",
                    "SUB113476_d0","SUB113483_d0","SUB113487_d0","SUB113490_d0","SUB113494_d0",
                    "SUB113495_d0","SUB113496_d0","SUB113498_d0","SUB113504_d0","SUB113505_d0",
                    "SUB113513_d0","SUB113514_d0","SUB113524_d0","SUB113526_d0","SUB113527_d0",
                    "SUB113532_d0","SUB113535_d0","SUB113537_d0","SUB113545_d0","SUB113548_d0",
                    "SUB113555_d0","SUB113558_d0","SUB113559_d0","SUB113561_d0","SUB113566_d0",
                    "SUB113567_d0","SUB113568_d0","SUB113571_d0","SUB113572_d0","SUB113582_d0",
                    "SUB113583_d0","SUB113588_d0","SUB113595_d0","SUB113610_d0")
      final_expr_vals <- remove_subs(final_expr_vals, subs_rm)
      
      # to replicate original file, probe_ids are simply the row number
      probe_ids <- seq(from = 1, to = length(gene_syms), by = 1)
      
      # to be consistent with naming to have Datasets.R match them correctly
      sdy <- "SDY67-batch2"
    }
    
  }else if(sdy %in% c("SDY63","SDY404","SDY400")){
    # download files from ImmPort and read in data
    inputFiles <- get_immport_files(sdy, rawdata_dir)
    rawdata <- as.data.frame(fread(inputFiles))
    
    if(yale.anno == "original"){
      # "SDY63_anno_tbl.tsb" was generated from original file (SDY63.GEMatrix.txt) 
      anno_tbl <- SDY63_anno_tbl
      colnames(anno_tbl) <- c("probeIDs","geneSymbol")
      keys <- anno_tbl$probeIDs
      values <- anno_tbl$geneSymbol
      
    }else if(yale.anno == "library"){
      id_ls <- as.list(illuminaHumanv4ALIAS2PROBE)
      keys <- inv_alias2probe(id_ls)
      values <- names(keys) 
      # gene_syms <- hashmap(exp_ids, names(exp_ids), probe_ids, dflt_ret = "")
      
    }else if(yale.anno == "manifest"){
      anno_tbl <- ilv4_anno_tbl
      keys <- anno_tbl$ID
      values <- anno_tbl$ILMN_Gene
    }
    
    # Clean and Prep
    probe_ids <- rawdata$ID_REF
    gene_syms <- hashmap(keys, 
                         values, 
                         probe_ids, 
                         dflt_ret = "")
    rawdata <- rawdata[ , grepl("PBMC", names(rawdata))]
    
    # Norm / Map
    final_expr_vals <- norm_map_mx(sdy, rawdata, "Sub.Org.Accession", "User.Defined.ID")
    
    # SDY404 ids from Immport have actual day the sample was taken (e.g. 4 instead of 2) according to
    # personal communication from Stefan Avey at Yale. However this creates conflicts with the original 
    # file that had standard days (e.g. 2 and 28). Therefore, changing here to ensure reproducibility 
    # of original information.
    if(sdy == "SDY404"){
      colnms <- colnames(final_expr_vals)
      bad_ids <-  c("SUB120470_d24", "SUB120472_d32", "SUB120480_d4", "SUB120481_d4", "SUB120483_d4",
                    "SUB120484_d4", "SUB120485_d4", "SUB120487_d4", "SUB120488_d4")
      good_ids <- c("SUB120470_d28", "SUB120472_d28", "SUB120480_d2", "SUB120481_d2", "SUB120483_d2",
                    "SUB120484_d2", "SUB120485_d2", "SUB120487_d2", "SUB120488_d2")
      id_hash <- hash(bad_ids, good_ids)
      colnms <- sapply(colnms, FUN = function(x){
        if(has.key(x, id_hash)){
          return(id_hash[[x]])
        }else{
          return(x)
        }
      })
      colnames(final_expr_vals) <- colnms
    }
    
  }else if(sdy == "SDY80"){
    # get data via scraping geo site and parsing into list of tables
    gef <- get_gef(sdy)
    inputFiles <- get_geo_files(gef)
    probe_ids <- unname(inputFiles[[1]][,1])
    expr_vals <- as.data.frame(do.call(cbind,
                                       lapply(inputFiles, FUN = function(x){
                                         return(x[,2])
                                        })
                                ))
    colnames(expr_vals) <- sapply(inputFiles, FUN = function(x){colnames(x)[[2]]})
    rawdata <- cbind(probe_ids, expr_vals)
    
    if(sdy80.anno == "original"){
      gene_syms <- hashmap(CHI_nih_gene_map$probeID, CHI_nih_gene_map$gs, probe_ids, keep_names = T)
      
    }else if(sdy80.anno == "library"){
      id_ls <- as.list(hugene10sttranscriptclusterALIAS2PROBE)
      keys <- inv_alias2probe(id_ls)
      gene_syms <- hashmap(keys, names(keys), probe_ids)
      
    }
    
    # To mimic original file, I remove all rows that were not 
    # successfully mapped, which is approximately 45%.  Unsure why this was done.
    gene_syms <- gene_syms[which(sapply(gene_syms, FUN = function(x){!is.na(x)}))]
    probe_ids <- probe_ids[which(probe_ids %in% names(gene_syms))]
    rawdata <- rawdata[ which(rawdata$probe_ids %in% names(gene_syms)), ]
    rawdata <- rawdata[,-1]
    
    if(sdy80.norm == F){
      # to mimic original file data, but use ImmuneSpace subjectIDs instead of biosampleIDs
      colnames(rawdata) <- unname(map_headers(sdy, 
                                           colnames(expr_vals), 
                                           is_col = "GE_name", 
                                           smpl_col = "GSM"))
      final_expr_vals <- rawdata
    }else{
      # log2 untransform rawdata and send down regular pipeline
      rawdata <- 2 ^ rawdata
      final_expr_vals <- norm_map_mx(sdy, rawdata, is_col = "GE_name", smpl_col = "GSM")
    }
    sdy <- "CHI-nih"
  }
  
  write_out(probe_ids, gene_syms, final_expr_vals, sdy, output_dir)
  
}


#-----------DEPRECATED METHODS--------------------------------------------------------
# FOR SDY67
# Using this manifest file from Illumina because bioconductor library for 450k methylation is defunct
# and no longer mainted.Available for download: 
# http://support.illumina.com/array/array_kits/infinium_humanmethylation450_beadchip_kit/downloads.html
# hmn450k_ann_tbl <- read.table(file.path(pp_dir,"IlluminaHuman_450kMethylation_anno_table.txt"),
#                               stringsAsFactors = F, sep = "\t", header = T)