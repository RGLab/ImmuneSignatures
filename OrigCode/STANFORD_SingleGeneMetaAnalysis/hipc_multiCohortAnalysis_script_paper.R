#######################################################################################
#HIPC gene-level multi-cohort analysis script
#by Francesco Vallania
#@Stanford, CA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#the purpose of this code is to perform all the gene-level analysis for HIPC
#######################################################################################

#######################################################################################
#HIPC Analysis functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~Cleaning workspace
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~Load all the necessary libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(pROC)
library(MADAM)
library(RMySQL)
library(ggplot2)
library(Biobase)
library(clinfun)
library(data.table)
library(MetaIntegrator)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ensembl_gene_connector_hs76
#~fetch list of all genes in Homo sapiens
#######################################################################################
ensembl_gene_connector_hs76 <- function(){

  #Establish connection to EnsEMBL
  con <- dbConnect(dbDriver("MySQL"),
                   username= 'anonymous',
                   host    = 'ensembldb.ensembl.org',
                   dbname  = 'homo_sapiens_core_76_38',
                   password= '',
                   port    = 3306);

  #build the parts for the query string
  query  <- "select gene.stable_id,
  xref.display_label
  from xref
  inner join gene
  where xref.xref_id = gene.display_xref_id";

  #run query and get table
  ens_table <- dbGetQuery(con,query);

  #Close connection
  dbDisconnect(con);

  #return table
  return(ens_table);
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#hipc_renaud2gem_corr
#~convert Reanud's HIPC objects into MetaIntegrator dataset inpt format
#######################################################################################
hipc_renaud2gem_corr <- function(renaud_file,name=NULL,force_log = FALSE){
  #read in the RDS
  object <- readRDS(renaud_file)

  #create key vector
  keys        <- gsub(";| /// ",",",as.character(fData(object)[,2]))
  names(keys) <- fData(object)[,1]

  #remove probes that do not map back to genes
  genes      <- ensembl_gene_connector_hs76()
  fake_genes <- which(!(keys %in% genes$display_label))
  keys[fake_genes] <- NA

  #create class vector
  classV <- rep(0,nrow(pData(object)))
  names(classV) <- row.names(pData(object))

  #create list
  gem <- list(expr          = exprs(object),
              pheno         = pData(object),
              name          = name,
              formattedName = gsub(".rds","",name),
              keys          = keys,
              class         = classV,
              comment       = "HIPC gem object created from SDY GEM matrix+demographics file")

  if(force_log == TRUE){
    #make sure there are not 0s [set to pseudocounts]
    gem$expr[which(gem$expr == 0)] <- 0.0001

    gem$expr    <- log2(gem$expr)
    gem$comment <- "Data was forced into log scale"
  }

  #final check
  if(!all(colnames(gem$expr)==row.names(gem$pheno))){
    stop("Warning! Error with the input file!")
  }

  #update here-> if SDY404
  if(length(grep("SUB120473_d[027]", colnames(gem$expr)))>0){

    #Swap 2 with incorrect sex
    inds1 <- grep("SUB120473_d[027]", colnames(gem$expr))
    inds2 <- grep("SUB120474_d[027]", colnames(gem$expr))
    tmp1 <- gem$expr[,inds1]
    gem$expr[,inds1] <- gem$expr[,inds2]
    gem$expr[,inds2] <- tmp1

    #Swap 2 additional samples with strong evidence
    inds1 <- grep("SUB120460_d0", colnames(gem$expr))
    inds2 <- grep("SUB120485_d28", colnames(gem$expr))
    tmp1 <- gem$expr[,inds1]
    gem$expr[,inds1] <- gem$expr[,inds2]
    gem$expr[,inds2] <- tmp1
  }

  #return gem object
  return(gem)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#GEM_formatter
#~formats GEM to be used in input with MetaIntegrator
#######################################################################################
GEM_formatter <- function(GEM, field = "fc_res_max_d20",baseline_only=T){
  GEM$class <- as.integer(as.character(GEM$pheno[[field]]))

  GEM$class[which(is.na(GEM$class))]          <- -1

  if(baseline_only==TRUE){
    GEM$class[which(GEM$pheno$Condition!='d0')] <- -1
  }

  GEM <- geo_subset_cutter(GEM)

  return(GEM)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#geo_subset_cutter [by Francesco Vallania 04/25/2014]
#~this function removes GSM entries in a GSE object
#######################################################################################
geo_subset_cutter <- function(GSE){

  #re-cut everything
  survivors <- which(GSE$class!=-1);
  GSE$class <- GSE$class[survivors];
  GSE$pheno <- GSE$pheno[survivors,];
  GSE$expr  <- GSE$expr[,survivors];

  #return GSE
  return(GSE);
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#multi_stat_cmp
#this function performs any pairwise comparison in one dataset
#input-> data.frame containing class and score as columns
#######################################################################################
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#multi_roc_cmp
#~this function performs any pairwise ROC comparison in one dataset
#~input-> data.frame containint class and score as columns
#######################################################################################
multi_roc_cmp <- function(input_df,
                          title_name=input_df$name[1]){

  #Get input classes and form pairs (combinations of n elements taken 2 at a time)
  input_classes <- unique(input_df$class);
  input_pairs   <- t(combn(input_classes,2));

  test_function <- function(sets){
    #select data for test

    roc_obj <- roc(subset(input_df,class %in% sets)$class,
                   subset(input_df,class %in% sets)$score,
                   lty=5,
                   lwd=5,
                   col='royalblue',
                   percent=T
    );

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

  #return data in the form of a data.frame
  return(res)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#arithMean_comb_validate
#~computes geometric mean in real space [arithmetic mean in log space]
#######################################################################################
arithMean_comb_validate <- function(GEM,
                                    genes_p = NULL,
                                    genes_n = NULL,
                                    scale   = 1){

  #compute score with mean
  score <- rep(0,length(GEM$class));
  if(!is.null(genes_p)){
    currGenes_p <- getSampleLevelGeneData(GEM,genes_p)
    score <- score + colMeans(currGenes_p,na.rm=T);
  }
  if(!is.null(genes_n)){
    currGenes_n <- getSampleLevelGeneData(GEM,genes_n)
    score <- score - colMeans(currGenes_n,na.rm=T);
  }
  if(!is.null(genes_p) || !is.null(genes_n)){
    if(scale==1){
      score <- scale(score)[,1];
    }
  }

  #return data.frame
  return(data.frame(class=GEM$class,
                    score=score,
                    name =GEM$name));
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#validator
#~Computes a validation object from a signature [returns data.frame,ROCs, and t-tests]
#######################################################################################
validator <- function(GEM,posG,negG){

  #compute validation
  validation <- arithMean_comb_validate(GEM,posG,negG)

  #compute stats and ROC curves for each comparison
  rocs  <- multi_roc_cmp(validation)
  stats <- multi_stat_cmp(validation)

  #return output
  return(list('v'=validation,'r'=rocs,'s'=stats))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#violinPlotter
#~Creates violin plot using GGPLOT2
#######################################################################################
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

  #return plot object
  return(finalPlot)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#negCorPlot
#~Creates correlation plots between old and young effect size vectors
#######################################################################################
negCorPlot <- function(youngES,oldES){

  #prepare the data
  youngESdT <- as.data.table(melt(youngES))
  youngESdT$age <- "young"

  oldESdT <- as.data.table(melt(oldES))
  oldESdT$age <- "older"

  finalDT <- rbind(oldESdT,youngESdT)
  finalDT <- dcast.data.table(finalDT,Var1+Var2~age,value.var='value')
  finalDT <- finalDT[!is.na(young)]
  finalDT <- finalDT[!is.na(older)]

  inputD <- unique(finalDT$Var2)
  indP   <- lapply(inputD,
                   function(i){
                     newCOrrPlotBW <- ggplot(finalDT[Var2==i],aes(y = older,x = young)) +
                       geom_vline(xintercept=0,size=0.1, colour="black")   +
                       geom_hline(yintercept=0,size=0.1, colour="black")   +
                       geom_point(color = 'gray65')                        +
                       theme_bw() +
                       theme(axis.title        = element_text(size=16,face="bold"),
                             legend.title      = element_blank(),
                             legend.key        = element_rect(colour = "white"),
                             legend.position   = c(0.25, 0.9),
                             legend.text       = element_text(size = 15,face="bold"),
                             legend.background = element_rect(colour = "lightgray")) +
                       geom_point(data = finalDT[Var2==i][grep("^(RAB24|GRB2|DPP3|ACTB|MVP)$",Var1),],
                                  aes(x = young,y = older),
                                  color = 'gray20',
                                  size  = 3,
                                  shape = 15) +
                       geom_point(data = finalDT[Var2==i][grep("^(PTPN22|PURA|SP4)$",Var1),],
                                  aes(x = young,y = older),
                                  color = 'gray20',
                                  size  = 3,
                                  shape = 17) +
                       xlab("Gene effect size (young)") + ylab("Gene effect size (older)") +
                       ggtitle(i)

                     return(newCOrrPlotBW)
                   })
  names(indP) <- inputD

  finalDT[,cor(young,older),by=Var2]

  #return the data
  return(list('cor' = finalDT[,cor(young,older),by=Var2],
              'plot'= indP))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#rocGGMultiRoc
#~This function computes ROC curves from an input list of pROC objects
#######################################################################################
rocGGMultiRoc <- function(validation_roc_list,
                          rocStyle      = 1,
                          rocSize       = 1.4,
                          newLineLegend = TRUE){

  #reformat in a data.table
  roc_data_table <- lapply(validation_roc_list,
                           function(i){
                             #
                             final_name <- NULL#paste(i$name,"\n",paste("AUC = ",round(i$roc$auc[[1]],2),"%",collapse=""))
                             if(newLineLegend == FALSE){
                               final_name <- paste(i$name," ",paste("AUC = ",round(i$roc$auc[[1]],2),"%",collapse=""))
                             }else{
                               final_name <- paste(i$name,"\n",paste("AUC = ",round(i$roc$auc[[1]],2),"%",collapse=""),
                                                   paste(" CI = ",paste(c(round(ci(i$roc)[1],2),round(ci(i$roc)[3],2)),collapse = "-")))
                             }

                             #
                             x <- data.table(name        = final_name,
                                             AUC         = paste("AUC = ",round(i$roc$auc[[1]],2),"%",collapse=""),
                                             Specificity = i$roc$specificities,
                                             Sensitivity = i$roc$sensitivities,
                                             FPR         = 100 - i$roc$specificities)
                             x <- x[order(Sensitivity),]
                             return(x)
                           })
  roc_data_table <- rbindlist(roc_data_table)

  #generate plot
  p <- ggplot(roc_data_table,aes(x= Specificity,y = Sensitivity,colour = name)) +
    geom_line(size = rocSize,linetype = rocStyle)                               +
    scale_x_reverse()                                                           +
    geom_segment(x=-100, y=0, xend=0, yend=100, linetype=3, colour="grey75")    +
    theme(legend.key.height=unit(2,"line"),
          legend.title    = element_blank(),
          legend.text     = element_text(size = 10),
          legend.position ="bottom")

  #return plot
  return(p)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#SDY80_time_course_plot_input
#~this funciton create timeplot for SDY80
#######################################################################################
SDY80_time_course_plot_input <- function(GSE47353_gem_object,
                                         ...){

  GSE47353_gem_object$pheno$timepoint.days  <- as.numeric(gsub("neg","-",gsub("d","",GSE47353_gem_object$pheno$Condition)))

  #select gene of interest and compute validation
  GSE47353_gem_object$pheno$score <- arithMean_comb_validate(GSE47353_gem_object,
                                                             ...)$score

  GSE47353_gem_object$pheno$resp <- c("Low responder",
                                      "Moderate responder",
                                      "High responder")[GSE47353_gem_object$class+1]

  return(GSE47353_gem_object$pheno[,43:45])
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#rocDT
#~this function turns a ROC list object into a ROC data.table
#######################################################################################
rocDT <- function(validation_roc_list){

  newLineLegend <- TRUE

  #reformat in a data.table
  roc_data_table <- lapply(validation_roc_list,
                           function(i){
                             #
                             final_name <- NULL#paste(i$name,"\n",paste("AUC = ",round(i$roc$auc[[1]],2),"%",collapse=""))
                             if(newLineLegend == FALSE){
                               final_name <- paste(i$name," ",paste("AUC = ",round(i$roc$auc[[1]],2),"%",collapse=""))
                             }else{
                               final_name <- paste(i$name,"\n",paste("AUC = ",round(i$roc$auc[[1]],2),"%",collapse=""),
                                                   paste(" CI = ",paste(c(round(ci(i$roc)[1],2),round(ci(i$roc)[3],2)),collapse = "-")))
                             }

                             #
                             x <- data.table(name        = final_name,
                                             AUC         = paste("AUC = ",round(i$roc$auc[[1]],2),"%",collapse=""),
                                             Specificity = i$roc$specificities,
                                             Sensitivity = i$roc$sensitivities,
                                             FPR         = 100 - i$roc$specificities)
                             x <- x[order(Sensitivity),]
                             return(x)
                           })
  roc_data_table <- rbindlist(roc_data_table)

  return(roc_data_table)
}

#######################################################################################
#######################################################################################
#R Main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################################################################################
#######################################################################################

#######################################################################################
#######################################################################################
#Prepare/Process Input
#######################################################################################
#######################################################################################

#######################################################################################
#Data Loading and Pre-processing
#######################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Load new Objects
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#read and process new files
renaud_path                <- '/input_data/'
hipc_gems_new              <- lapply(list.files(renaud_path,full.names = FALSE),
                                     function(i)
                                       hipc_renaud2gem_corr(paste(renaud_path,i,sep=""),i))
names(hipc_gems_new)       <- list.files(renaud_path,full.names = FALSE)

#train only on SDY212,SDY63,SDY400, and SDY404, leave SDY80 for validation
hipc_gems                  <- hipc_gems_new[grep("212|400|404|63",names(hipc_gems_new))]
names(hipc_gems)           <- gsub(".rds","",names(hipc_gems))

hipc_gems_valid            <- hipc_gems_new[grep("CHI",names(hipc_gems_new))]
names(hipc_gems_valid)     <- "SDY80"

hipc_gems_valid_old        <- hipc_gems_new[grep("SDY67.rds",names(hipc_gems_new))]
names(hipc_gems_valid_old) <- "SDY67"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Prepare young training cohort for MetaIntegrator: High Vs Low Responders at Day0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#remove anything but day 0
hipc_youngInput <- lapply(hipc_gems,
                          function(j){
                            #select pre-vaccination
                            excluded <- which(j$pheno$Condition != 'd0')

                            #label it for removal
                            j$class[excluded] <- -1
                            j <- geo_subset_cutter(j)

                            #return it
                            return(j)
                          })

#tag based on young_fc_res_max_d30 (2 is cases, 0 is controls, remove NAs)
hipc_young <- lapply(hipc_youngInput,
                     function(j){

                       #select pre-vaccination
                       excluded  <- which(j$pheno$young_fc_res_max_d30 == 1)
                       excluded2 <- which(is.na(j$pheno$young_fc_res_max_d30))
                       positive  <- which(j$pheno$young_fc_res_max_d30 == 2)

                       #label it for removal
                       j$class[excluded]  <- -1
                       j$class[excluded2] <- -1
                       j$class[positive]  <-  1
                       j <- geo_subset_cutter(j)

                       #return things that have been classified
                       if(!all(j$class==0)){

                         #return if there are more than one obs per group
                         if(length(which(table(j$class)>1))==2){
                           #return it
                           return(j)
                         }
                       }
                     })
hipc_young <- hipc_young[which(sapply(hipc_young,function(i) !is.null(i)))]

#Structure it as a list to correctly input into MetaIntegrator
hipc_young_new_Meta <- list("originalData" = hipc_young)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Prepare older training cohort for MetaIntegrator: High Vs Low Responders at Day0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#remove anything but day 0
hipc_OldInput <- lapply(hipc_gems,
                        function(j){
                          #select pre-vaccination
                          excluded <- which(j$pheno$Condition != 'd0')

                          #label it for removal
                          j$class[excluded] <- -1
                          j <- geo_subset_cutter(j)

                          #return it
                          return(j)
                        })
#tag based on old_fc_res_max_d30 (2 is cases, 0 is controls, remove NAs)
hipc_Old <- lapply(hipc_OldInput,
                   function(j){

                     #select pre-vaccination
                     excluded  <- which(j$pheno$old_fc_res_max_d30 == 1)
                     excluded2 <- which(is.na(j$pheno$old_fc_res_max_d30))
                     positive  <- which(j$pheno$old_fc_res_max_d30 == 2)

                     #label it for removal
                     j$class[excluded]  <- -1
                     j$class[excluded2] <- -1
                     j$class[positive]  <-  1
                     j <- geo_subset_cutter(j)

                     #return things that have been classified
                     if(!all(j$class==0)){

                       #return if there are more than one obs per group
                       if(length(which(table(j$class)>1))==2){
                         #return it
                         return(j)
                       }
                     }
                   })
hipc_Old <- hipc_Old[which(sapply(hipc_Old,function(i) !is.null(i)))]

#Structure it as a list to correctly input into MetaIntegrator
hipc_old_new_Meta <- list("originalData" = hipc_Old)


#######################################################################################
#Run MetaIntegrator on Young and Older individuals
#######################################################################################

#compute meta-analysis on young and older
hipc_young_new_Meta <- MetaIntegrator::runMetaAnalysis(hipc_young_new_Meta)
hipc_old_new_Meta   <- MetaIntegrator::runMetaAnalysis(hipc_old_new_Meta)

#filter meta-analysis results by 10% FDR
hipc_young_new_Meta <- MetaIntegrator::filterGenes(hipc_young_new_Meta,
                                                   isLeaveOneOut           = FALSE,
                                                   FDRThresh               = 0.10,
                                                   heterogeneityPvalThresh = 0)
hipc_old_new_Meta   <- MetaIntegrator::filterGenes(hipc_old_new_Meta,
                                                   isLeaveOneOut           = FALSE,
                                                   FDRThresh               = 0.10,
                                                   heterogeneityPvalThresh = 0)

#######################################################################################
#######################################################################################
#Analysis
#######################################################################################
#######################################################################################

#######################################################################################
#Performance on discovery cohorts
#######################################################################################

#run on positive genes only
youngSignPosPerf <- lapply(hipc_young_new_Meta$originalData,
                           function(i)
                             validator(i,
                                       hipc_young_new_Meta$filterResults$FDR0.1_es0_nStudies1_looaFALSE_hetero0$posGeneNames,
                                       NULL))

#run on positive genes only
oldSignPosPerf <- lapply(hipc_old_new_Meta$originalData,
                         function(i)
                           validator(i,
                                     hipc_young_new_Meta$filterResults$FDR0.1_es0_nStudies1_looaFALSE_hetero0$posGeneNames,
                                     NULL))

#run on negative genes only
youngSignNegPerf <- lapply(hipc_young_new_Meta$originalData,
                           function(i)
                             validator(i,
                                       NULL,
                                       hipc_young_new_Meta$filterResults$FDR0.1_es0_nStudies1_looaFALSE_hetero0$negGeneNames))
#create violin plots for both sets
youngSignPosViolinPlot <- lapply(youngSignPosPerf,
                                 function(i)
                                   violinPlotter(i$v,
                                                 c("Low responders","High responders"),
                                                 c("deepskyblue","red"),
                                                 c(19,19)))
youngSignNegViolinPlot <- lapply(youngSignNegPerf,
                                 function(i)
                                   violinPlotter(i$v,
                                                 c("Low responders","High responders"),
                                                 c("deepskyblue","red"),
                                                 c(19,19)))

oldSignNegViolinPlot <- lapply(oldSignPosPerf,
                               function(i)
                                 violinPlotter(i$v,
                                               c("Low responders","High responders"),
                                               c("deepskyblue","red"),
                                               c(19,19)))

#######################################################################################
#Performance on validation cohort
#######################################################################################

#select set of positive genes to be used from this point on
inputPosGenes10p <- hipc_young_new_Meta$filterResults$FDR0.1_es0_nStudies1_looaFALSE_hetero0$posGeneNames

#select day0 samples from SDY80 with adjMFC labeling [use fc_res_max_d30]
SDY80day0 <- GEM_formatter(hipc_gems_valid$SDY80,field = "fc_res_max_d30")
SDY80day0$class[which(is.na(SDY80day0$pheno$Age.class))] <- -1
SDY80day0 <- geo_subset_cutter(SDY80day0)

#compute validation object for SDY80
SDY80day0Validation10p <- validator(SDY80day0,
                                    inputPosGenes10p,
                                    NULL)

#create violin plots
SDY80day0ViolinPlot <- violinPlotter(SDY80day0Validation10p$v,
                                     c("Low responders","Moderate responders","High responders"),
                                     c("deepskyblue","magenta","red"),
                                     c(15,19,17))

#test on SDY67
#select day0 samples from SDY67 with adjMFC labeling [use fc_res_max_d30]
SDY67day0 <- GEM_formatter(hipc_gems_valid_old$SDY67,field = "fc_res_max_d30")
SDY67day0$class[which(is.na(SDY67day0$pheno$Age.class))] <- -1
SDY67day0$class[which(SDY67day0$class==1)] <- -1
SDY67day0 <- geo_subset_cutter(SDY67day0)
SDY67day0$class[which(SDY67day0$class==2)] <-  1
SDY67day0Validation <- validator(SDY67day0,
                                 inputPosGenes10p,
                                 NULL)
SDY67day0ViolinPlot <- violinPlotter(SDY67day0Validation$v,
                                     c("Low responders","High responders"),
                                     c("deepskyblue","red"),
                                     c(19,19))

#######################################################################################
#Create ROC curve plots
#######################################################################################

#training rocs using positive signature
youngSignPosRocs <- list(youngSignPosPerf$SDY212$r[[1]],
                         youngSignPosPerf$SDY400$r[[1]],
                         youngSignPosPerf$SDY404$r[[1]],
                         youngSignPosPerf$SDY63$r[[1]])
youngSignPosRocs[[1]]$name <- names(youngSignPosPerf)[1]
youngSignPosRocs[[2]]$name <- names(youngSignPosPerf)[2]
youngSignPosRocs[[3]]$name <- names(youngSignPosPerf)[3]
youngSignPosRocs[[4]]$name <- names(youngSignPosPerf)[4]
youngSignPosRocP           <- rocGGMultiRoc(youngSignPosRocs) + theme_bw()

#training rocs using negative signature
youngSignNegRocs <- list(youngSignNegPerf$SDY212$r[[1]],
                         youngSignNegPerf$SDY400$r[[1]],
                         youngSignNegPerf$SDY404$r[[1]],
                         youngSignNegPerf$SDY63$r[[1]])
youngSignNegRocs[[1]]$name <- names(youngSignNegPerf)[1]
youngSignNegRocs[[2]]$name <- names(youngSignNegPerf)[2]
youngSignNegRocs[[3]]$name <- names(youngSignNegPerf)[3]
youngSignNegRocs[[4]]$name <- names(youngSignNegPerf)[4]
youngSignNegRocP           <- rocGGMultiRoc(youngSignNegRocs) + theme_bw()

#validation roc on SDY80
SDY80day0ValidationRocs <- list(SDY80day0Validation10p$r[[1]],
                                SDY80day0Validation10p$r[[2]])
SDY80day0ValidationRocs[[1]]$name <- "High responders Vs Low responders"
SDY80day0ValidationRocs[[2]]$name <- "Moderate Responders Vs Low responders"
SDY80day0ValidationRocP           <- rocGGMultiRoc(SDY80day0ValidationRocs) +
  theme_bw() + scale_colour_manual(values=c("springgreen1","darkgreen"))

#validation roc on SDY80
SDY67day0ValidationRocs <- list(SDY67day0Validation$r[[1]])#,
#SDY67day0Validation$r[[2]])
SDY67day0ValidationRocs[[1]]$name <- "High responders Vs Low responders"
#SDY67day0ValidationRocs[[2]]$name <- "Moderate Responders Vs Low responders"
SDY67day0ValidationRocP           <- rocGGMultiRoc(SDY67day0ValidationRocs) +
  theme_bw() + scale_colour_manual(values=c("red"))

#######################################################################################
#Time-course plots
#######################################################################################

#create input data table with all time points from SDY80
SDY80allDays             <- GEM_formatter(hipc_gems_valid$SDY80,
                                          field = "fc_res_max_d30",
                                          baseline_only = F)
SDY80allDays$class[which(is.na(SDY80allDays$pheno$Age.class))] <- -1
SDY80allDays             <- geo_subset_cutter(SDY80allDays)
SDY80_TC_input_data      <- SDY80_time_course_plot_input(SDY80allDays,
                                                         inputPosGenes10p,
                                                         NULL)
SDY80_TC_input_data$resp <- factor(SDY80_TC_input_data$resp,
                                   c("Low responder",
                                     "Moderate responder",
                                     "High responder"),
                                   ordered = TRUE)
#compute plot object using GGPLOT2
SDY80_TC_plot  <- ggplot(SDY80_TC_input_data,
                         aes(x    = factor(timepoint.days),
                             y    = score,
                             fill = resp))                     +
  geom_boxplot(outlier.size = 0)                               +
  geom_point(position=position_jitterdodge(dodge.width=0.8),
             aes(color = resp,shape = resp),size = 3)          +
  scale_fill_manual(values=c("white","white","white"))         +
  scale_shape_manual(values=c(15,19,17))                       +
  scale_color_manual(values=c("deepskyblue","magenta","red"))  +
  xlab("Time (Days)") + ylab("Response Score")

#compute tests and ROC curves
syd80TCdT <- as.data.table(SDY80_TC_input_data,keep.rownames = T)
syd80TCdT$rn <- gsub("_d.*$","",syd80TCdT$rn)

#compute day0 vs day1 paired t-test
t01  <- dcast(syd80TCdT[timepoint.days %in% c(0,1)],rn~timepoint.days,value.var='score')
t01p <- t.test(t01$`0`,t01$`1`,paired = T,alternative = 'less')

#compute comparison of low vs high responders at each time point
timeSDY80 <- lapply(unique(syd80TCdT$timepoint.days),
                    function(t)
                      t.test(syd80TCdT[timepoint.days==t & resp=='Low responder']$score,
                             syd80TCdT[timepoint.days==t & resp=='High responder']$score))
names(timeSDY80) <- unique(syd80TCdT$timepoint.days)

#compute AUCs for each time point
timeSDY80auc <- lapply(unique(syd80TCdT$timepoint.days),
                       function(t){
                         tempDT <- syd80TCdT[timepoint.days==t][resp!='Moderate responder']
                         tempDT$class <- 0
                         tempDT[resp=='High responder']$class <- 1
                         return(roc(tempDT$class,tempDT$score))
                       })
names(timeSDY80auc) <- unique(syd80TCdT$timepoint.days)


#######################################################################################
#Sex analysis
#######################################################################################

#create data.table
ageSexDT <- rbindlist(lapply(hipc_youngInput,
                             function(i){

                               #remove NAs
                               excluded1 <- which(is.na(i$pheno$young_fc_res_max_d30))
                               i$class[excluded1] <- -1
                               i$class[which(i$pheno$young_fc_res_max_d30 == 1)]  <-  1
                               i$class[which(i$pheno$young_fc_res_max_d30 == 2)]  <-  2
                               i                  <- geo_subset_cutter(i)

                               #compute score
                               mfcScore   <- data.table('sub' = row.names(i$pheno),
                                                        'Age' = i$pheno$Age,
                                                        "Sex" = i$pheno$Gender)

                               validation <- as.data.table(arithMean_comb_validate(i,
                                                                                   inputPosGenes10p,
                                                                                   NULL),keep.rownames = T)
                               outDT <- merge(validation,mfcScore,by.y='sub',by.x='rn')
                               return(outDT)
                             }))
ageSexDT$name <- gsub(".rds","",ageSexDT$name)

#Test for M/F differences [do not use SDY63 because they have only 1 male]
SexTtest <- lapply(setdiff(unique(ageSexDT$name),"SDY63"),
                   function(i)
                     t.test(ageSexDT[name==i & Sex=='Male']$score,ageSexDT[name==i & Sex=='Female']$score))
SexFisherP <- fisher.method(t(as.data.frame(sapply(SexTtest,function(i) i$p.value))))


#######################################################################################
#Cell-Proportion correction analysis
#######################################################################################

#load and log data
SDY80propCorr  <- readRDS('corr_data/SDY80_corrected.rds')
SDY80propCorrL <- log2(sweep(SDY80propCorr,2,apply(SDY80propCorr,2,min)-1))

#Compute geometric mean [using log scale]
SYD80propCorrPS  <- SDY80propCorrL[which(row.names(SDY80propCorrL) %in% inputPosGenes10p),]

if(is.null(nrow(SYD80propCorrPS))){
  SYD80propCorrPS  <- as.data.table(scale(SYD80propCorrPS),keep.rownames=T)
}else{
  SYD80propCorrPS  <- as.data.table(scale(colMeans(SYD80propCorrPS)),keep.rownames=T)
}

#create annotation
SDY80_PC_ann  <- data.table('sample' = row.names(SDY80day0$pheno),
                            'resp'   = SDY80day0$pheno$fc_res_max_d30)

#merge annotation with geometric mean
SDY80_PC_final                  <- merge(SDY80_PC_ann,SYD80propCorrPS,by.x='sample',by.y='rn')
setnames(SDY80_PC_final,c('sample','class','score'))
SDY80_PC_final$name             <- "SDY80"
SDY80_PC_final$labels           <- "High responder"
SDY80_PC_final[class==0]$labels <- "Low responder"
SDY80_PC_final[class==1]$labels <- "Moderate responder"
SDY80_PC_final$class <- as.numeric(SDY80_PC_final$class)

#compute ROCs of High Vs Low and Moderate Vs Low
SDY80_PC_ROCs <- multi_roc_cmp(SDY80_PC_final)

for(i in 1:length(SDY80_PC_ROCs)){
  SDY80_PC_ROCs[[i]]$name <- gsub("SDY80 \n ","",SDY80_PC_ROCs[[i]]$name)
}
SDY80_PC_ROCs_correctedP   <- rocGGMultiRoc(SDY80_PC_ROCs[1:2]) + theme_bw() +
  scale_colour_manual(values=c("springgreen1","darkgreen"))

#combine them for a final plot
finalSYD80 <- c(SDY80_PC_ROCs[1:2],SDY80day0ValidationRocs[1:2])
finalSYD80[[1]]$name <- "High VS Low Responders: Corrected"
finalSYD80[[2]]$name <- "Moderate VS Low Responders: Corrected"
finalSYD80[[3]]$name <- "High VS Low Responders: Uncorrected"
finalSYD80[[4]]$name <- "Moderate VS Low Responders: Uncorrected"

#add set label to spearate original data from proportion corrected data
finalSYD80pOne      <-  rocDT(finalSYD80[1:2])
finalSYD80pOne$set  <- "Corrected"
finalSYD80pTwo      <-  rocDT(finalSYD80[3:4])
finalSYD80pTwo$set  <- "Uncorrected"

#create final combined plot using GGPLOT2
SDY80_PC_ROCsCombinedPLot <- ggplot(rbind(finalSYD80pOne,finalSYD80pTwo),
                                    aes(x        = Specificity,
                                        y        = Sensitivity,
                                        colour   = name,
                                        linetype = set))                      +
  geom_line(size = 1)                                                         +
  scale_x_reverse()                                                           +
  geom_segment(x=-100, y=0, xend=0, yend=100, linetype=3, colour="grey75")    +
  theme(legend.key.height=unit(2,"line"),
        legend.title    = element_blank(),
        legend.text     = element_text(size = 10),
        legend.position ="bottom") +
  theme_bw() +
  scale_colour_manual(values=c("springgreen1",
                               "springgreen1",
                               "darkgreen",
                               "darkgreen")) +
  theme(legend.key.size = unit(2, 'lines'))

#add violin and correlation plot
SDY80_PC_final$class <- SDY80_PC_final$class-1
SDY80_PC_violin <- violinPlotter(SDY80_PC_final[,c('class','score'),with=F],
                                 c("Low responders","Moderate responders","High responders"),
                                 c("deepskyblue","magenta","red"),
                                 c(15,19,17))

SDY80day0Validation10p$v$sample <- row.names(SDY80day0Validation10p$v)
combinedPCandOriDT <- merge(SDY80_PC_final,SDY80day0Validation10p$v,by='sample')



#######################################################################################
#Old Vs Young Correlation Plots
#######################################################################################

#create input dataset
young_vs_old <- as.data.table(merge(hipc_young_new_Meta$metaAnalysis$pooledResults,
                                    hipc_old_new_Meta$metaAnalysis$pooledResults,by=0))

#create figure
newCOrrPlotBW <- ggplot(young_vs_old,aes(y = effectSize.y,x = effectSize.x)) +
  geom_vline(xintercept=0,size=0.1, colour="black")   +
  geom_hline(yintercept=0,size=0.1, colour="black")   +
  geom_point(color = 'gray65')                        +
  theme_bw() +
  theme(axis.title        = element_text(size=16,face="bold"),
        legend.title      = element_blank(),
        legend.key        = element_rect(colour = "white"),
        legend.position   = c(0.25, 0.9),
        legend.text       = element_text(size = 15,face="bold"),
        legend.background = element_rect(colour = "lightgray")) +
  geom_point(data = young_vs_old[grep("^(RAB24|GRB2|DPP3|ACTB|MVP|DPP7|ARCP4)$",young_vs_old$Row.names),],
             aes(x = effectSize.x,y = effectSize.y),
             color = 'gray20',
             size  = 3,
             shape = 15) +
  geom_point(data = young_vs_old[grep("^(PTPN22|PURA|SP4)$",young_vs_old$Row.names),],
             aes(x = effectSize.x,y = effectSize.y),
             color = 'gray20',
             size  = 3,
             shape = 17) +
  xlab("Gene effect size (young)") + ylab("Gene effect size (older)")

#compute plots @ a single dataset level
newCOrrPlotBWpD <- negCorPlot(hipc_young_new_Meta$metaAnalysis$datasetEffectSizes,
                              hipc_old_new_Meta$metaAnalysis$datasetEffectSizes)


#######################################################################################
#######################################################################################
#Save Output
#######################################################################################
#######################################################################################

save(file = 'hipc_analysis_output_multicohort_analysis.RData',
     #results
     hipc_gems_new,
     hipc_young_new_Meta,
     hipc_old_new_Meta,
     #validation on SDY80
     SDY80day0,
     SDY80day0Validation10p,
     SDY80day0ViolinPlot,
     SDY80day0ValidationRocs,
     SDY80day0ValidationRocP,
     SDY80_TC_plot,
     timeSDY80,
     #SDY67
     SDY67day0ValidationRocP,
     SDY67day0ViolinPlot,
     SDY67day0Validation,
     #SDY80 proportion corrected
     SDY80_PC_ROCsCombinedPLot,
     SDY80_PC_violin,
     combinedPCandOriDT,
     #correlation plots
     young_vs_old,
     newCOrrPlotBW,
     newCOrrPlotBWpD)
