####################
## Set up R session for co-deletions and amplifications project
####################

## This scripts installs packages, data and functions that I might need for the R session.

#################
###Load packages
library(tidyr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(org.Hs.eg.db)
library(biomaRt)
library(parallel)
library(ggplot2)

###################
###Set working directory
setwd("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/")


#################
### Load data
#load("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/R workspaces/co-deletions workspace.RData")
load("./R workspaces/co-deletions workspace.RData")
co.deletions.per.cytoband.circle.plots.table <- readRDS("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/total.co-deletion.events.per.cytoband.rds")
deletions.per.cytoband.circle.plots.table<- readRDS( "/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/total.deletion.events.per.cytoband.rds")
short.cnv.list<- readRDS("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/target.cancer.list.rds")
co_deletions_removed_zeros_plot_table2<- readRDS("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/BRCA_co_deletion_distance_plot_table.rds")
gene_information_list<- readRDS("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/target_gene_information_list.rds")
## Table of distance from target gene and proportion of deletions
co_deletions_distance_from_target_gene_plot_table<- readRDS("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/BRCA_co_deletion_distance_from_target_gene_plot_table.rds")





################
### Load Functions.

###F1: Function to import multiple files from multiple folders

import.files.from.directories<-function(x,file.to.import){
  currentwd<- getwd()
  setwd(x)
  directory.names<-dir()
  my.list <- vector("list", length(directory.names))
  for (i in 1: length(directory.names)){
    move.to.directory<-paste0(x,"/",directory.names[i])
    setwd(move.to.directory)
    my.list[[i]]<-read.delim(file.to.import, stringsAsFactors = FALSE, header = TRUE)
    print(directory.names[i])
  }
  names(my.list)<- directory.names
  setwd(currentwd)
  return(my.list)
}

##########
##F2: Function to select some or all CNV datasets from CNV list object 
#and combining data into one large dataframe.
join.cnv.datasets<- function(x, column, data.sets = "all data sets"){
  
  if(data.sets == "all data sets"){
    
    index<- seq(1:length(x))
    
  } else {
    names.of.tables<-names(x)
    names.of.tables
    index<-which(names.of.tables %in% data.sets)
    print("DO NOT WORRY ABOUT THE FOLLOWING WARNING MESSAGE:")
    
  }
  
  df<- x[[1]] %>% dplyr::select(c(Gene.Symbol,Locus.ID,Cytoband))
  
  for (i in index){
    
    df2<-x[[i]][,c(1,column:ncol(x[[i]]))]
    
    df<- full_join(df, df2, by = "Gene.Symbol")
    
  }
  return(df)
}

##########
### F3: Obtain chromosomal locations of genes

chromosomal_location<- function(object_name){
  
  ##Variables I want:
  keys<- as.character(object_name$Locus.ID)
  columns<- c("CHR", "CHRLOC", "CHRLOCEND")
  
  ## Search for gene chromosomal locations:
  genes.of.interest<- AnnotationDbi::select(org.Hs.eg.db, keys, columns, keytype = "ENTREZID")
  
  ## Remove duplicated entries:
  genes.of.interest<- na.omit(genes.of.interest[!duplicated(genes.of.interest$ENTREZID), ])
  
  ## Add additional columns to dataframe including strand, gene start and end:
  genes.of.interest$strand<- ifelse(genes.of.interest$CHRLOC <0, "-", "+")
  genes.of.interest$start<- abs(genes.of.interest$CHRLOC)
  genes.of.interest$end<- abs(genes.of.interest$CHRLOCEND)
  
  ##########
  ###Join gene location dataframe to original CNV data table
  
  ##Rename ENTREZID column so it can be joined to oringinal dataframe:
  genes.of.interest<- dplyr::rename(genes.of.interest, Locus.ID = ENTREZID)
  
  genes.of.interest$Locus.ID<- as.integer(genes.of.interest$Locus.ID)
  
  genes.of.interest<- full_join(genes.of.interest,object_name, by="Locus.ID")
  
  ##########
  ### Order genes by chromosome and location:
  
  ## Rename X and Y chromosomes to integer so they can be properly sorted.
  genes.of.interest$CHR<-sub("X", "23", genes.of.interest$CHR)
  genes.of.interest$CHR<-sub("Y", "24", genes.of.interest$CHR)
  genes.of.interest$CHR<- as.integer(genes.of.interest$CHR)
  genes.of.interest<- dplyr::arrange(genes.of.interest, CHR, start)
  genes.of.interest$CHR<-sub("23","X",  genes.of.interest$CHR)
  genes.of.interest$CHR<-sub("24", "Y", genes.of.interest$CHR)
  
  return(genes.of.interest)
}


#################
### Events per cytoband or chromosomal interval
events.per.cytoband<- function(object_name, threshold = -1, cytoband_column = 10,
                               column_data_start = 11, chromosome_interval = 0, 
                               select_chromosome = 1, deletion = TRUE){
  
  ## Code to get proportion of deletions per cytoband
  cnv.matrix<- as.matrix(object_name[,column_data_start:ncol(object_name)])
  
  if (deletion == TRUE) {
    
    cnv.matrix<- ifelse(cnv.matrix <= threshold, 1, 0)
    
  } else {
    
    cnv.matrix<- ifelse(cnv.matrix >= threshold, 1, 0)
  }
  
  ##Create list to store data:
  if (chromosome_interval == 0){
    results <- vector("list", 2)
  } else {
    results <- vector("list", 4)
  }
  
  ##Create matrix with raw deletions/amplification values and number of 
  #deletions/amplifications per gene:
  results[[1]]<- cnv.matrix %>% 
    as.data.frame() %>% 
    dplyr::mutate(sum.of.deletions = rowSums(.)) %>%
    dplyr::mutate(number.of.tumours = ncol(.)-1) %>%
    cbind(cytoband=object_name[,cytoband_column], 
          chromosome = object_name$CHR, 
          start = object_name$start, .) 
  
  ##Create dataframe with number of deletions/amplifications per gene, number of genes and number of 
  #potential deletion/amplification events that could occur:
  results[[2]]<- results[[1]]%>%
    group_by(cytoband) %>%
    dplyr::summarise(sum.of.genes.deleted = sum(sum.of.deletions),
                     total.number.of.genes=n(),
                     total.number.of.events = sum(number.of.tumours)) %>%
    dplyr::mutate(proportion.of.deletions = sum.of.genes.deleted/total.number.of.events) %>%
    tidyr::separate(cytoband, c("chromosome", "band"), sep = "[p:q]", remove = FALSE, convert = TRUE)
  
  results[[2]]$chromosome<-sub("X", "23", results[[2]]$chromosome)
  results[[2]]$chromosome<-sub("Y", "24", results[[2]]$chromosome)
  results[[2]]$chromosome<- as.integer(results[[2]]$chromosome)
  results[[2]]<- dplyr::arrange(results[[2]], chromosome, cytoband)
  results[[2]]$chromosome<-sub("23", "X", results[[2]]$chromosome)
  results[[2]]$chromosome<-sub("24", "Y", results[[2]]$chromosome)
  
  #############
  #### Code to get proportion of deletions per unit size of chromosome:
  
  if (chromosome_interval > 0){
    ## Very rough estimate of lengths of chromosomes:
    results[[3]]<- results[[1]] %>%
      group_by(chromosome) %>%
      summarise(chromosome_start = min(start),
                chromosome_end = max(start)) %>%
      mutate(estimated_chromosome_length = chromosome_end - chromosome_start,
             intervals_for_kb = estimated_chromosome_length/1000,
             intervals_for_10kb = estimated_chromosome_length/10000,
             intervals_for_100kb = estimated_chromosome_length/100000,
             intervals_for_Mb = estimated_chromosome_length/1000000,
             intervals_for_10Mb = estimated_chromosome_length/10000000)
    
    ## Select chromosome to investigate:
    filter.table<- results[[1]] %>%
      dplyr::filter(chromosome == select_chromosome)
    
    ## Create data table with proportion of deletion/amplification events that have occured:
    results[[4]]<- filter.table$start %>%
      cut(chromosome_interval, labels = seq(1:chromosome_interval)) %>%
      cbind(Intervals = ., filter.table) %>%
      dplyr::group_by(Intervals) %>%
      dplyr::summarise(sum.of.genes.deleted = sum(sum.of.deletions),
                       total.number.of.genes=n(),
                       total.number.of.potential.events = sum(number.of.tumours)) %>%
      dplyr::mutate(proportion.of.deletions = sum.of.genes.deleted/total.number.of.potential.events)
  }
  
  return(results)
}


####################
###F1: Function to create matrix containing proportions of co-deletions or amplifications 
#normalised by the total number of tumours compared

co.deletion_co.amplification_matrix<- function(cnv.table, column_start = 11, threshold = -1, 
                                               selection_criteria, Gene.Symbol = FALSE, start = FALSE, 
                                               Chromosome = 0, Cytoband = FALSE, remove_NA = TRUE, deletion = TRUE, normalisation = "total.tumour.number"){
  if (remove_NA == TRUE){
    
    cnv.table<- cnv.table %>%
      dplyr::filter(!is.na(start))
  }
  
  if (Gene.Symbol == TRUE){
    
    ##Select Genes of interest by Gene Symbol and convert CNV data to matrix:
    matrix<- cnv.table %>% dplyr::filter(Gene.Symbol %in% selection_criteria) 
    cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
    rownames(cnv.matrix)<- matrix$Gene.Symbol
    
  } else if (Chromosome[1] > 0 & start == FALSE){
    
    ##Select Chromosome of interest and convert CNV data to matrix:
    matrix<- cnv.table %>% dplyr::filter(CHR %in% Chromosome) 
    cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
    rownames(cnv.matrix)<- matrix$Gene.Symbol
    
  } else if (start == TRUE & Chromosome[1] == 0){
    
    ##Select Chromosomal region of interest and convert CNV data to matrix:
    matrix<- cnv.table$start %>%
      dplyr::between(selection_criteria[1], selection_criteria[2]) %>%
      cnv.table[.,]
    
    ## Remove NAs which are maintained by dplyr::between function
    matrix <- matrix %>% filter(!is.na(start))
    
    cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
    rownames(cnv.matrix)<- matrix$Gene.Symbol
    
  }  else if (start == TRUE & Chromosome[1] > 0){
    
    ##Select Chromosome of interest and convert CNV data to matrix:
    matrix<- cnv.table %>% dplyr::filter(CHR %in% Chromosome) 
    
    ##Select Chromosomal region of interest and convert CNV data to matrix:
    matrix<- matrix %>%
      dplyr::filter(start >= selection_criteria[1], end <= selection_criteria[2])
    
    ##Convert Chromosome and region of interest into matrix:
    cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
    rownames(cnv.matrix)<- matrix$Gene.Symbol
    
  } else if (Cytoband == TRUE){
    ##Select Cytobands of interest and convert CNV data to matrix:
    matrix<- cnv.table %>% dplyr::filter(Cytoband %in% selection_criteria) 
    cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
    rownames(cnv.matrix)<- matrix$Gene.Symbol
    
  } else {
    ##Convert ALL CNV data to matrix:
    cnv.matrix<- as.matrix(cnv.table[,column_start:ncol(cnv.table)])
    rownames(cnv.matrix)<- cnv.table$Gene.Symbol
  }
  
  ##Create a binary matrix of CNV data such that deletions (deletion = TRUE) or 
  #amplifications (deletion = FALSE) below (for deletions) or above (for amplifications) a threshold = 1 
  if (deletion == TRUE) {
    
    cnv.matrix<- ifelse(cnv.matrix <= threshold, 1, 0)
    
  } else {
    
    cnv.matrix<- ifelse(cnv.matrix >= threshold, 1, 0)
  }
  
  
  ## Calculate total number of co-deletions or amplifications
  heatmap.matrix<- cnv.matrix%*%t(cnv.matrix)
  
  if (normalisation == "total.tumour.number") {
    
    heatmap.matrix<- heatmap.matrix/ncol(cnv.matrix)
    
  } else if (normalisation == "tumours.with.event") {
    tumours.with.del.or.amp<- colSums(cnv.matrix)
    tumours.with.del.or.amp<-sum(tumours.with.del.or.amp >0)
    heatmap.matrix<- heatmap.matrix/tumours.with.del.or.amp
    
  } else if (normalisation == "none") {
    
  } else if (normalisation == "frequency.for.whole.sample"){
    heatmap.matrix<- sum(heatmap.matrix)/(ncol(cnv.matrix)*(nrow(cnv.matrix)^2))
    
  }
  
  ## Remove NAs caused by divideing by 0.
  heatmap.matrix[is.nan(heatmap.matrix)] = 0
  
  return(heatmap.matrix)
}
