####################
### Function to create table of intergene distances and proportion of co-deletions/co-amplifications
####################




######################
## Initial function parameters:
target_genes<- c("MET", "CDKN2A", "RB1", "WWOX", 
                 "LRP1B", "PDE4D", "CCNE1", "TP53",
                 "FGFR1", "MYC", "EGFR","WHSC1L1",
                 "ERBB2", "MCL1", "MDM2", "CCND1", "ATM",
                 "NOTCH1", "PPP2R2A", "BRD4", "ARID1A",
                 "STK11", "PARK2")

cnv.table<- threshold_short_cnv_list_loc[[1]]



###################
##Function to create list of target genes, their chromosome, cytoband, start and end
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

###############
##Test function
gene_information_list<- create_target_gene_information_list(cnv_table = cnv.table, target_genes = target_genes)
gene_information_list[[4]]
gene_information_list





#################
### Function to calculate the distance from each gene in list to a gene of interest (target_gene):
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


############
###Test function
distance_from_target_gene_function(cnv.table = cnv.table, gene_information = gene_information_list[[1]], distance = 2.5e+06)


#############
### Function that takes a dataframe including Gene.Symbol, chromosome(CHR) and start site (start) and 
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

#######
###Test function:

chromosome_9_intergene_distance<- intergene_distance_matrix_function(cnv.table = cnv.table, chromosome = 9)
dim(chromosome_9_intergene_distance)
chromosome_9_intergene_distance


chromosome_all_intergene_distance<-lapply(c((1:22), "X"), function(x) intergene_distance_matrix_function(cnv.table, chromosome = x))
length(chromosome_all_intergene_distance)
dim(chromosome_all_intergene_distance[[9]])










################
### Function for all v's all genes??????
###############

intergene_distance_co_deletion_co_amplification_function<- function(cnv.table, Chromosome, deletion = TRUE, threshold = -2, compare_all_genes = FALSE, normalisation = "total.tumour.number", column_start = 11){}









###################
### Function for gene v's target gene:
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

###########
### Test function:
test<- distance_from_target_gene_co_deletion_co_amplification_function(cnv.table = cnv.table, gene_information_list = gene_information_list, distance = 2.5e+06, deletion = TRUE, threshold = -2, compare_all_genes = FALSE, normalisation = "tumours.with.event")
dim(test)
head(test)
tail(test)

test2<- distance_from_target_gene_co_deletion_co_amplification_function(cnv.table = cnv.table, gene_information_list = gene_information_list, distance = 2.5e+06, deletion = TRUE, threshold = -2, compare_all_genes = TRUE, normalisation = "tumours.with.event")
dim(test2)
head(test2)
tail(test2)


