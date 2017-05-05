#######################################################################################
# HIPC gene-level multi-cohort analysis script
# Original Developer: Francesco Vallania, Stanford University (2014)
# Refactoring Developer: Evan Henrich, Fred Hutchinson Cancer Research Center (2017)
# purpose: perform gene-level analysis for HIPC
#######################################################################################

#######################################################################################
###                       HELPER FUNCTIONS                                          ###
#######################################################################################
# ensembl_gene_connector_hs76
#~fetch list of all genes in Homo sapiens
# loading as part of pkg
# ensembl_gene_connector_hs76 <- function(){
#
#   #Establish connection to EnsEMBL
#   con <- dbConnect(dbDriver("MySQL"),
#                    username= 'anonymous',
#                    host    = 'ensembldb.ensembl.org',
#                    dbname  = 'homo_sapiens_core_76_38',
#                    password= '',
#                    port    = 3306);
#
#   #build the parts for the query string
#   query  <- "select gene.stable_id,
#   xref.display_label
#   from xref
#   inner join gene
#   where xref.xref_id = gene.display_xref_id";
#
#   #run query and get table
#   ens_table <- dbGetQuery(con,query);
#
#   #Close connection
#   dbDisconnect(con);
#
#   #return table
#   return(ens_table);
# }

# Convert eset into appropriate list type for MetaIntegrator --------------
makeGEM <- function(sdy, eset_list, ens_table, valProc){
  eset <- eset_list[[sdy]]
  if( !(sdy %in% c("SDY67", "SDY80")) | valProc == TRUE ){
    eset <- eset[ , eset$fc_res_max_d30 != 1 ] # rm moderates
  }

  #create key vector
  keys        <- fData(eset)$gene_symbol
  names(keys) <- rownames(fData(eset))

  #remove probes that do not map back to genes
  keys[ !(keys %in% ens_table$display_label) ] <- NA

  # create class vector and change high-responders to "1" from "2" for MetaIntegrator
  # Able to do this because all moderate responders have already been removed
  endPoint <- pData(eset)[["fc_res_max_d30"]]
  endPoint <- if( !(sdy %in% c("SDY67","SDY80")) ){ gsub(2, 1, endPoint) }else{ endPoint }
  classV <- as.numeric(endPoint)
  names(classV) <- row.names(pData(eset))

  # Ensure that exprs is numeric
  em <- exprs(eset)
  class(em) <- "numeric"

  #create list
  gem <- list(expr          = em,
              pheno         = pData(eset),
              name          = sdy,
              formattedName = sdy,
              keys          = keys,
              class         = classV,
              comment       = "HIPC gem object created from SDY GEM matrix + demographics file")

  #return gem object
  return(gem)
}

# Performs any pairwise ROC comparison in one dataset -----------------------------
# input-> data.frame containint class and score as columns
multi_roc_cmp <- function(input_df,
                          title_name=input_df$name[1]){

  #Get input classes and form pairs (combinations of n elements taken 2 at a time)
  input_classes <- unique(input_df$class);
  input_pairs   <- t(combn(input_classes,2));

  test_function <- function(sets){

    # select data for test
    roc_obj <- roc(subset(input_df,class %in% sets)$class,
                   subset(input_df,class %in% sets)$score,
                   lty=5,
                   lwd=5,
                   col='royalblue',
                   percent=T
                )

    return(list(roc=roc_obj,
                name=paste(title_name,"\n",
                           subset(input_df,class==sets[1])$label[1],
                           "VS",
                           subset(input_df,class==sets[2])$label[1],sep=" ")
    ))
  }

  #apply function
  res <- lapply(1:nrow(input_pairs),
                function(i) test_function(input_pairs[i,]))
}

multi_stat_cmp <- function(input_df,test.type = 't.test'){
  #Get input classes and form pairs (combinations of n elements taken 2 at a time)
  input_classes <- unique(input_df$class);
  input_pairs   <- t(combn(input_classes,2));

  test_function <- function(sets){
    #select data for test
    set_a <- input_df$score[input_df$class==sets[1]];
    set_b <- input_df$score[input_df$class==sets[2]];

    #perform stats
    two.sided <- NA;
    if(test.type=='wilcox.test'){
      two.sided <- wilcox.test(set_a,set_b)$p.value;
    }
    else if(test.type=='t.test'){
      two.sided  <- t.test(set_a,set_b)$p.value;
    }
    else{
      warning('Unsupported test type.');
    }
    return(as.numeric(two.sided));
  }

  #apply function
  res <- sapply(1:nrow(input_pairs),
                function(i) test_function(input_pairs[i,]));

  #return data in the form of a data.frame
  return(data.frame(Var1   = input_pairs[,1],
                    Var2   = input_pairs[,2],
                    pvalue = res,
                    bonf   = res*nrow(input_pairs)));
}

# computes geometric mean in real space [arithmetic mean in log space] ---------------
arithMean_comb_validate <- function(GEM,
                                    genes_p = NULL,
                                    genes_n = NULL,
                                    scale   = 1){

  #compute score with mean
  score <- rep(0, length(GEM$class) )
  if( !is.null(genes_p) ){
    currGenes_p <- getSampleLevelGeneData(GEM,genes_p)
    score <- score + colMeans(currGenes_p, na.rm=T)
  }
  if( !is.null(genes_n) ){
    currGenes_n <- getSampleLevelGeneData(GEM,genes_n)
    score <- score - colMeans(currGenes_n, na.rm=T)
  }
  if( !is.null(genes_p) || !is.null(genes_n) ){
    if( scale == 1 ){
      score <- scale(score)[,1];
    }
  }

  #return data.frame
  return(data.frame(class = GEM$class,
                    score = score,
                    name  = GEM$name));
}

# Creates violin plot using GGPLOT2 -------------------------------------------------
violinPlotter <- function(inputDF,
                          labels,
                          colorVect,
                          shapeVect){

  #format data into data.table
  inputDT       <- as.data.table(inputDF)
  inputDT$label <- labels[inputDT$class+1]
  inputDT$label <- factor(inputDT$label,
                          levels  = labels,
                          ordered = TRUE)

  #function for standard error
  sterr <- function(vector){sd(vector)/sqrt(length(vector))}

  #add column with mean and standard error
  inputDT[, mean:=mean(score), by=class]
  inputDT[,sterr:=sterr(score),by=class]

  #create plot object with ggplot
  finalPlot <- ggplot(inputDT,aes(y=score,x=label)) +
    geom_violin(trim = F) +
    geom_jitter(width = 0.2,aes(col=label,shape=label),size=4) +
    scale_color_manual(values=colorVect)                       +
    scale_shape_manual(values=shapeVect)                       +
    geom_point(aes(y=mean),size=6,col='grey35') +
    geom_errorbar(aes(ymax=mean+sterr,ymin=mean-sterr),
                  col='grey35',width = 0.2) + xlab("") +
    ylab("Response score") +
    theme(axis.text.x=element_blank())

}

# This function computes ROC curves from an input list of pROC objects --------------
rocGGMultiRoc <- function(validation_roc_list,
                          rocStyle      = 1,
                          rocSize       = 1.4,
                          newLineLegend = TRUE){

  roc_data_table <- rocDT(validation_roc_list)

  # generate plot
  p <- ggplot(roc_data_table,aes(x= Specificity,y = Sensitivity,colour = name)) +
    geom_line(size = rocSize,linetype = rocStyle)                               +
    scale_x_reverse()                                                           +
    geom_segment(x=-100, y=0, xend=0, yend=100, linetype=3, colour="grey75")    +
    theme(legend.key.height=unit(2,"line"),
          legend.title    = element_blank(),
          legend.text     = element_text(size = 10),
          legend.position ="bottom")
}

# reformat in a data.table -----------------------------------------------------------
rocDT <- function(validation_roc_list){

  roc_data_table <- rbindlist(lapply(validation_roc_list, function(i){
                          final_name <- paste(i$name,
                                              "\n",
                                              paste("AUC = ",
                                                    round(i$roc$auc[[1]],2),
                                                    "%",
                                                    collapse=""),
                                              paste(" CI = ",
                                                    paste(c(round(ci(i$roc)[1],2),
                                                            round(ci(i$roc)[3],2)),
                                                          collapse = "-")))

                             x <- data.table(name        = final_name,
                                             AUC         = paste("AUC = ",
                                                                 round(i$roc$auc[[1]],2),
                                                                 "%",
                                                                 collapse=""),
                                             Specificity = i$roc$specificities,
                                             Sensitivity = i$roc$sensitivities,
                                             FPR         = 100 - i$roc$specificities)
                             x <- x[ order(Sensitivity), ]
                           }))

}

# get metaIntegratorResults ----------------------------------------------------------
metaRes <- function(gems){
  res <- list()

  # Structure it as a list to correctly input into MetaIntegrator
  mRes <- MetaIntegrator::runMetaAnalysis( list("originalData" = gems) )

  # filter meta-analysis results by 10% FDR
  mResFilt <- MetaIntegrator::filterGenes(mRes,
                                                 isLeaveOneOut           = FALSE,
                                                 FDRThresh               = 0.10,
                                                 heterogeneityPvalThresh = 0)

  # store up/down regulated genes as results, filtering FDR thresholds as above
  tmp <- mRes$metaAnalysis$pooledResults
  tmp <- tmp[ tmp$effectSizeFDR < 0.1, ]
  res$upReg <- tmp[ tmp$fisherFDRUp < 0.1, ]
  res$dnReg <- tmp[ tmp$fisherFDRDown < 0.1, ]

  # select set of positive genes to be used from this point on
  res$posGenes <- mResFilt$filterResults$FDR0.1_es0_nStudies1_looaFALSE_hetero0$posGeneNames

  return(res)
}


#######################################################################################
###                       EXPORTED FUNCTIONS                                        ###
#######################################################################################
#' helper function SDY404 column swap
#' 
#' @param em expression matrix with colnames including sub1 and sub2
#' @param sub1 subject 1
#' @param sub2 subject 2
#' @export
swapCols <- function(em, sub1, sub2){
  idx1 <- grep(sub1, colnames(em))
  idx2 <- grep(sub2, colnames(em))
  tmp1 <- em[, idx1]
  em[, idx1] <- em[ , idx2]
  em[ , idx2] <- tmp1
  return(em)
}

#***********MAIN METHOD***************************************************************
#' Function to perform single gene meta analysis for HIPC ImmuneSignatures Project
#'
#' @import ggplot2 pROC clinfun
#' @param adjEsetList list of expressionSet objects holding expr, pheno data adjusted for cohort
#' @param cohort Study cohort, young or old
#' @export
singleGeneAnalysis <- function(adjEsetList, cohort){

  # Create holder list for results
  finalRes <- list("disc","val")

  # reformat esets for use in MetaIntegrator functions
  # Need to take out moderate responders here for discovery studies,
  # because runMetaAnalysis only allows 1 or 0 in class aka endPoint,
  # but leave for validation studies for plotting work.
  gems <- lapply(names(adjEsetList), function(sdy){
              makeGEM(sdy = sdy,
                      eset_list = adjEsetList,
                      ens_table = ens_table,
                      valProc = FALSE)
  })
  names(gems) <- names(adjEsetList)

  # Assign discovery set and validation
  discGems <- gems[grep("212|400|404|63", names(gems))]
  finalRes$disc <- metaRes(discGems)

  sdyId <- ifelse(cohort == "young", "SDY80", "SDY67")
  valGemUnProc <- gems[[ sdyId  ]]

  # compute validation object using validation study gem
  valObj <- arithMean_comb_validate(valGemUnProc, finalRes$disc$posGenes, NULL)

  # ViolinPlot for Val - fig 6A in manuscript
  finalRes$val$ViolinPlot <- violinPlotter(valObj,
                                 c("Low responders","Moderate responders","High responders"),
                                 c("deepskyblue","magenta","red"),
                                 c(15,19,17))

  # ROCPlot for Val - fig 6B in manuscript
  valRoc <- multi_roc_cmp(valObj)
  valRocLs <- list(valRoc[[1]], valRoc[[2]])
  valRocLs[[1]]$name <- "High Responders Vs Low Responders"
  valRocLs[[2]]$name <- "Moderate Responders Vs Low Responders"
  finalRes$val$RocPlot <- rocGGMultiRoc(valRocLs) + theme_bw() +
    scale_colour_manual(values=c("springgreen1","darkgreen"))

  return(finalRes)
}
