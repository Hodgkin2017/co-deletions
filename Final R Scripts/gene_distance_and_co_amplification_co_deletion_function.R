####################
### Function to create table of intergene distances and proportion of co-deletions/co-amplifications
####################

###Contents:

##Functions:
#F1. create_target_gene_information_list    #Function to create a list containing target gene name, cytoband, chromosome and
                                            #start and end locations
                                            #Function arguments:
                                            #cnv.table            # CNV file containing gene location information
                                            #target_genes         # Vector of gene names used to create list
#F2. distance_from_target_gene_function     #Function to calculate the distances between the target gene and genes within a 
                                            #certain 5' and 3' distance of the genes start and end site.
                                            #Function arguments:
                                            #cnv.table            # CNV file containing gene location information
                                            #gene_information     # List containing genes of interest and location information 
                                                                  #created using create_target_gene_information_list function 
                                            #distance             # distance 5' and 3' of target gene to specify which 
                                                                  #genes to include in analysis.
#F3. intergene_distance_matrix_function     #Function to create intergene distances for all genes on a chromosome.
                                            #Function arguments:
                                            #cnv.table            # CNV file containing gene location information
                                            #chromosome           # Chromosome to get intergene distances for.

#F4. distance_from_target_gene_co_deletion_co_amplification_function    #Function to obtain the distance between each gene and 
                                                                        #the target gene and the proportion of 
                                                                        #co-deletions/co-amplifications
                                            #Function arguments:
                                            #cnv.table                # CNV file
                                            #column_start             # First column in file containing tumour CNV data
                                            #gene_information_list    # List containing genes of interest and location information 
                                                                      #created using create_target_gene_information_list function 
                                            #distance                 # distance 5' and 3' of target gene to specify which 
                                                                      #genes to include.
                                            #threshold                # Value above (deletion = FALSE) or below 
                                                                      #(deletion = TRUE) which CNVs will be included in analysis
                                            #deletion = TRUE          # If deletion = TRUE: Count number of deletion events
                                                                      # If deletion = FALSE: Count number of amplification events
                                            #normalisation            #If normalisation = "total.tumour.number" then normalise by total number of tumours
                                                                      #If normalisation = "tumours.with.event" then normalise by 
                                                                      # number of tumours with deletion or amplification event. e.g. if only 7 out of 90 
                                                                      # individuals had a deletion in a cytoband then divide number of co-deletions by 7
                                                                      #If normalisation = "none" then raw number of co-deletions or co-amplifications returned.
                                                                      #If normalisation = "frequency.for.whole.sample" then normalise number of events 
                                                                      #by total number of possible events. Returns a single value and not a matrix
#F5. mean_co_deletion_co_amplification_values_around_gene             #Function to create a list containing target gene name, cytoband, chromosome and
                                                                      #start and end locations
                                            #Function arguments:
                                            #cnv.table                            # CNV file containing gene location information
                                            #distance_from_gene_to_calculate_mean # distance 5' and 3' of target gene to specify which               
                                                                                  #genes to include in analysis.
 
##Actions and or Loops outside of functions:
#None

##Selected Objects:
#O1. target_genes           #List of gene names used as input for functions
#O2. cnv.table              #A table from TCGA containing CNV data and gene information (generated by chromosomal_location function)
#O3. gene_information_list  #List of genes and their location information required input for functions
#O4. chromosome_9_intergene_distance                  #Output co-deletions for chromosome 9 only
#O5. chromosome_all_intergene_distance                  #Output co-deletions for Cytobands "9p21.2" and "9p21.3" only
#O6. test1                  #Output co-deletions between residues 100,000 and 250,000 only
#O7. test2                  #Output co-deletions on chromosome 9 between 100,000 and 250,000 only






######################
## O1:Initial function parameters:
target_genes<- c("MET", "CDKN2A", "RB1", "WWOX", 
                 "LRP1B", "PDE4D", "CCNE1", "TP53",
                 "FGFR1", "MYC", "EGFR","WHSC1L1",
                 "ERBB2", "MCL1", "MDM2", "CCND1", "ATM",
                 "NOTCH1", "PPP2R2A", "BRD4", "ARID1A",
                 "STK11", "PARK2")
##O2:
cnv.table<- threshold_short_cnv_list_loc[[1]]



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

###############
##O3: Test function
gene_information_list<- create_target_gene_information_list(cnv_table = cnv.table, target_genes = target_genes)
gene_information_list[[4]]
gene_information_list





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


############
###Test function
distance_from_target_gene_function(cnv.table = cnv.table, gene_information = gene_information_list[[1]], distance = 2.5e+06)


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

#######
###Test function:

##O4: One chromosome:
chromosome_9_intergene_distance<- intergene_distance_matrix_function(cnv.table = cnv.table, chromosome = 9)
dim(chromosome_9_intergene_distance)
chromosome_9_intergene_distance

##O5: All chromosomes
chromosome_all_intergene_distance<-lapply(c((1:22), "X"), function(x) intergene_distance_matrix_function(cnv.table, chromosome = x))
length(chromosome_all_intergene_distance)
dim(chromosome_all_intergene_distance[[9]])



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

###########
### Test function:

##O6: Compare each each to target gene
test1<- distance_from_target_gene_co_deletion_co_amplification_function(cnv.table = cnv.table, gene_information_list = gene_information_list, distance = 2.5e+06, deletion = TRUE, threshold = -2, compare_all_genes = FALSE, normalisation = "tumours.with.event")
dim(test1)
head(test1)
tail(test1)

##O7: Compare all genes to each other
test2<- distance_from_target_gene_co_deletion_co_amplification_function(cnv.table = cnv.table, gene_information_list = gene_information_list, distance = 2.5e+06, deletion = TRUE, threshold = -2, compare_all_genes = TRUE, normalisation = "tumours.with.event")
dim(test2)
head(test2)
tail(test2)

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

