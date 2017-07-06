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
## List containing dataframes with each dataframe corresponding to CNVs for one cancer using gistic thresholded values 
threshold_cnv_list<- readRDS("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/threshold_cnv_list.rds")
## Dataframe containing gistic thresholded CNV values for all cancer types.
threshold_CNV_all_table<- readRDS("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/threshold_CNV_all_table.rds")
## Dataframe containing gistic thresholded CNV values for all cancer types and gene chromosomal location information
threshold_CNV_all_table_loc<- readRDS("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/threshold_CNV_all_table_loc.rds")
## List containing dataframes with each dataframe corresponding to CNVs for one cancer using gistic thresholded values and containing chromosomal location
threshold_cnv_list_loc<- readRDS("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/threshold_cnv_list_loc.rds")
## List of CNV Dataframes of gistic thresholded CNV scores for specific cancers of interest 
threshold_short_cnv_list<- readRDS("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/threshold_short_cnv_list.rds")
## List of CNV Dataframes of gistic thresholded CNV scores with chromosomal location information for specific cancers of interest 
threshold_short_cnv_list_loc<- readRDS("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/threshold_short_cnv_list_loc.rds")
## List of target genes 1:
target_genes<- readRDS("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/target_genes1.rds")
## Comprehensive clinical data from firehose using fbget
clinical_fbget_list<- readRDS(file = "./R workspaces/clinical_fbget_list.rds")
## Concise clinical data from cbioportal
clinical_cbioportal_list<- readRDS(file = "./R workspaces/clinical_cbioportal_list")
## List of Overall survival and disease free survival
clinical_survival_list<- readRDS(file = "./R workspaces/clinical_survival_list")



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


###################
##F1: Function to create list of target genes, their chromosome, cytoband, start and end
###################

create_target_gene_information_list<- function(cnv_table, target_genes){
  
  ## Create an empty list to store gene information
  gene_information_list<- vector("list", length(target_genes)) 
  
  ##loop to create list of genes and their start and stop locations
  for (i in 1: length(target_genes)){
    
    gene<- target_genes[i]
    
    gene_information<- cnv.table %>% 
      dplyr::filter(Gene.Symbol == gene) %>%
      dplyr::select(Gene.Symbol, CHR, Cytoband, start, end)
    
    gene_information<- as.list(gene_information)
    gene_information_list[[i]]<- gene_information
    
  }
  
  return(gene_information_list)
  
}





#################
###F2: Function to calculate the distance from each gene in list to a gene of interest (target_gene):
#################


distance_from_target_gene_function<- function(cnv.table, gene_information, distance){
  
  ##Get end site of genes 2.5MB away of 5' start site of target gene
  
  end_sites_5prime_genes<- cnv.table %>%
    dplyr::filter(CHR == gene_information[[2]]) %>%
    dplyr::filter(start <= gene_information[[4]], start >=  gene_information[[4]] - distance) %>%
    dplyr::select(end)
  
  ##Calculate the distances between genes
  distance_5prime_genes<- gene_information[[4]] - end_sites_5prime_genes
  
  ##Any value <0 = 0 i.e. the gene of interest and any overlapping genes
  distance_5prime_genes[distance_5prime_genes < 0]<- 0
  
  ##Get start site of genes 2.5MB away of 3' end of end of target gene
  start_sites_3prime_genes<-cnv.table %>%
    dplyr::filter(CHR == gene_information[[2]]) %>%
    dplyr::filter(start > gene_information[[4]], end <= gene_information[[5]] + distance ) %>%
    dplyr::select(start)
  
  ##calculate the distance to the end of the genes 5' of the start of the target gene.
  distance_3prime_genes<- start_sites_3prime_genes - gene_information[[5]]
  
  ##Any value <0 = 0 i.e. the gene of interest and any overlapping genes
  if (nrow(distance_3prime_genes) == 0){
    
  }else {
    distance_3prime_genes[distance_3prime_genes < 0]<- 0
  }
  
  ##Make sure both tables have same colnames
  colnames(distance_5prime_genes)<- "start"
  
  ## Join data together
  
  distance_from_target_gene<- rbind(distance_5prime_genes, distance_3prime_genes)
  
  return(distance_from_target_gene)
  
}



#############
###F3: Function that takes a dataframe including Gene.Symbol, chromosome(CHR) and start site (start) and 
#creates a matrix of intergene distances using the start site per chromosome
#############

#Modify function to accept target genes and intervals?

intergene_distance_matrix_function<- function(cnv.table, chromosome){
  
  ##Get start location for chromosome of choice
  selected.genes.start<- cnv.table %>%
    dplyr::filter(CHR == chromosome) %>%
    dplyr::select(start)
  
  ##Calculate pairwise distances between gene start sites
  intergene.distance.table <- data.frame(matrix(NA, ncol = nrow(selected.genes.start), nrow = nrow(selected.genes.start)))
  
  for (j in 1:nrow(selected.genes.start)) {
    
    intergene.distance.table[,j]<- abs(selected.genes.start[,1] - selected.genes.start[j,1])
    
  }
  
  ##Get gene names and use to name matrix rows and columns 
  Gene_Symbol<- cnv.table %>%
    dplyr::filter(CHR == chromosome) %>%
    dplyr::select(Gene.Symbol) %>%
    t() %>% #required to convert output from data.frame to vector
    as.character()
  
  colnames(intergene.distance.table)<- Gene_Symbol
  rownames(intergene.distance.table)<- Gene_Symbol
  
  return(intergene.distance.table)
  
}





###################
###F4: Function for gene v's target gene: Creates table with gene distances and 
#proportion of co-deletion/amplifications:
###################

distance_from_target_gene_co_deletion_co_amplification_function<- function(cnv.table, gene_information_list, 
                                                                           distance = 2.5e+06, 
                                                                           deletion = TRUE, 
                                                                           threshold = -2, 
                                                                           compare_all_genes = FALSE, 
                                                                           normalisation = "tumours.with.event", 
                                                                           column_start = 11){
  
  ############
  ### Create a long datafame of co-deletions 2.5MB upstream and downstream of gene of interest.
  
  ##Create co-deletion matricies for each target gene
  co.deletion.per.target.gene<- lapply(gene_information_list, function(x) co.deletion_co.amplification_matrix(cnv.table = cnv.table, column_start = column_start, threshold = threshold, start = TRUE, Chromosome = x[[2]], selection_criteria = c(x[[4]] - distance, x[[5]] + distance), deletion = deletion, normalisation = normalisation))
  
  ##Add gene name to each column to be used with gather function later
  co.deletion.per.target.gene<- lapply(co.deletion.per.target.gene, function(x) as.data.frame(cbind(Gene.Symbol.row = rownames(x), x)))
  
  ##Create a long 3 column wide table with pair-wise proportion of pair wise deletions
  gathered<- lapply(co.deletion.per.target.gene, function(x) tidyr::gather(x, Gene.Symbol.col,proportion, 2:ncol(x)))
  
  ##Keep rows relating to MET v's all genes only and not all genes v's all genes
  if (compare_all_genes == FALSE){
    gathered_target_genes<- vector("list", length(target_genes))
    
    for(i in 1: length(target_genes)){
      
      gene<- target_genes[[i]][1]
      
      gathered_target_genes[[i]]<- gathered[[i]] %>% 
        dplyr::filter(Gene.Symbol.col == gene)
    }
    
  } else if (compare_all_genes == TRUE){
    
    gathered_target_genes<- gathered
    
  }
  
  ##Bind all dataframes in list together
  gathered.co.deletion.per.target.gene<- do.call(rbind, gathered_target_genes)
  
  ############
  ### Create a dataframe of gene distances from gene of interest.
  if (compare_all_genes == FALSE){
    
    distance_from_target_gene_table<- lapply(gene_information_list, function(x) distance_from_target_gene_function(cnv.table = cnv.table, gene_information = x, distance = distance))
    
    ##Bind all dataframes in list together
    distance_between_genes_table<- do.call(rbind, distance_from_target_gene_table)
    
  } else if (compare_all_genes == TRUE){
    
    intergene_distance_table_function<- function(cnv.table, gene_information, distance){
      selected_genes_start<- cnv.table %>%
        dplyr::filter(CHR == gene_information[[2]]) %>%
        dplyr::filter(start >=  gene_information[[4]] - distance, end <= gene_information[[5]] + distance ) %>%
        dplyr::select(start)
      
      intergene_distance_table <- data.frame(matrix(NA, ncol = nrow(selected_genes_start), nrow = nrow(selected_genes_start)))
      
      for (j in 1:nrow(selected_genes_start)) {
        
        intergene_distance_table[,j]<- abs(selected_genes_start[,1] - selected_genes_start[j,1])
        
      }
      
      Gene_Symbol<- cnv.table %>%
        dplyr::filter(CHR == gene_information[[2]]) %>%
        dplyr::filter(start >=  gene_information[[4]] - distance, end <= gene_information[[5]] + distance ) %>%
        dplyr::select(Gene.Symbol) %>%
        t() %>% #required to convert output from data.frame to vector
        as.character()
      
      colnames(intergene_distance_table)<- Gene_Symbol
      rownames(intergene_distance_table)<- Gene_Symbol
      
      return(intergene_distance_table)
      
    }
    
    intergene_distance_table<- lapply(gene_information_list, function(x) intergene_distance_table_function(cnv.table = cnv.table, gene_information = x, distance = distance))
    
    ##Add gene name to each column to be used with gather function later
    intergene_distance_table<- lapply(intergene_distance_table, function(x) as.data.frame(cbind(Gene.Symbol.row = rownames(x), x)))
    
    ##Add target gene to each column to be used with gather function later
    for (i in 1: length(gene_information_list)){
      
      intergene_distance_table[[i]]<- cbind(gene_information_list[[i]][[1]], intergene_distance_table[[i]])
      
    }
    #intergene_distance_table<- lapply(intergene_distance_table, function(x) as.data.frame(cbind(Target_gene = (x), x)))
    
    ##Create a long 3 column wide table with pair-wise proportion of pair wise deletions
    gathered<- lapply(intergene_distance_table, function(x) tidyr::gather(x, Gene.Symbol.col,proportion, 3:ncol(x)))
    
    ##Bind all dataframes in list together
    distance_between_genes_table<- do.call(rbind, gathered)
    
  }
  
  ###########
  ### join pair-wise distance table to pair-wise proportion of co-deletions table:
  if (compare_all_genes == FALSE){
    
    co_deletions_distance_from_target_gene_plot_table<- cbind(proportion_of_co_deletion = gathered.co.deletion.per.target.gene, distance_from_target_gene = distance_between_genes_table)
    
    colnames(co_deletions_distance_from_target_gene_plot_table)<- c("Comparison_gene", "Target_gene", "proportion_co_del_amp", "distance_from_target_genes")
    
  } else if (compare_all_genes == TRUE) {
    
    co_deletions_distance_from_target_gene_plot_table<- cbind(distance_from_target_gene = distance_between_genes_table, proportion_of_co_deletion = gathered.co.deletion.per.target.gene)
    
    colnames(co_deletions_distance_from_target_gene_plot_table)<- c("Target_gene", "Comparison_gene1", "Comparison_gene2", "distance", "Comparison_gene1_again", "Comparison_gene2_again", "proportion_of_co_deletion")
    
    
  }
  
  return(co_deletions_distance_from_target_gene_plot_table)
  
}


##############
###F5: Function to get mean co-deletion of genes +/- n genes away from current gene:
##############

mean_co_deletion_co_amplification_values_around_gene<- function(co_deletion_table, 
                                                                distance_from_gene_to_calculate_mean){
  
  ##Calculate co-del/co-amp for n genes around current gene(sixe of sliding window)
  n<- distance_from_gene_to_calculate_mean
  
  result<-rep(NA, nrow(co_deletion_table))
  
  ##Loop to get co-del/co-amp mean in sliding window
  for(i in 1: nrow(co_deletion_table)){
    
    if(i < n){
      start<- i
      end<- i+n
      result[i]<- mean(co_deletion_table[i, start:end])
      
    } else if((i >= n) & (i+n-1 < ncol(co_deletion_table))) {
      start<- i-n
      end<- i+n
      result[i]<- mean(co_deletion_table[i, start:end])
      
    } else if(i+n > ncol(co_deletion_table)) {
      
      start<- i-n
      end<- ncol(co_deletion_table)
      result[i]<- mean(co_deletion_table[i, start:end])
      
    }
  }
  
  return(result)
  
}


