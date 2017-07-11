##########
###Parameters:

target_gene_list<- gene_information_list[[2]]#CDKN2A
survival_time_list<- clinical_survival_list[[1]]#BRCA
cnv.table<- threshold_short_cnv_list_loc[[1]]#BRCA
dim(cnv.table)
distance = 2.5e+06 
threshold = -2
deletion = TRUE
time_of_death_column 
death_event_column
column_start = 11
start = TRUE
remove_NA = TRUE
Cytoband = FALSE
print_to_screen = FALSE 
plot_graph = FALSE

#########
### Main Function:

survival_analysis_of_gene_list<- function(target_gene_list, survival_time_list, cnv.table, distance = 2.5e+06, 
                                          threshold = -2, deletion = TRUE, time_of_death_column, 
                                          death_event_column, column_start = 11, start = TRUE, 
                                          remove_NA = TRUE, Cytoband = FALSE, print_to_screen = FALSE, 
                                          plot_graph = FALSE){
  
  ##########################################
  ###Create a list of tables were each table contains deletion categories(1-4) for each tumour(columns) per 
  #gene(rows) surrounding an interval around a target gene.
  
  
  
  
  
  
  #########################################
  
  
  
  
  
  
  
  ##########
  ###Create a list of tables were each table contains deletion categories(1-4) for each tumour(columns) per 
  #gene(rows) surrounding an interval around a target gene.
  
  deletion_category_target_gene_list <- lapply(target_gene_list, function(x) categorise_deletion_type_around_target_gene(cnv.table = cnv.table, target_gene = x[[1]], Chromosome = x[[2]], 
                                                                                                                         selection_criteria = c(x[[4]] - distance, x[[5]]+ distance), threshold = threshold, deletion = deletion, column_start = column_start, start = start, remove_NA = remove_NA, Cytoband = Cytoband))
  
  # deletion_category_target_gene_list[[1]][,1:3]
  # deletion_category_target_gene_list[[2]][,1:3]
  
  
  
  
  
  
  ########
  ### Take a deletion category table for a target gene and join to survival table
  ## Add extra column with patient ID:
  # deletion_category_target_gene<-deletion_category_target_gene_list[[2]]
  # survival_time_list<- clinical_survival_list[[1]]
  # 
  # join_clinical_deletion_category_table(deletion_category_target_gene = deletion_category_target_gene, survival_time_list = survival_time_list)
  
  deletion_category_survival_target_gene_list<- lapply(deletion_category_target_gene_list, function(x) join_clinical_deletion_category_table(x, survival_time_list, time_of_death_column = time_of_death_column, death_event_column = death_event_column))
  # deletion_category_survival_target_gene_list[[1]][1:3, 1:12]
  
  #########################From here!
  return(deletion_category_survival_target_gene_list)
  
}

test_function<- survival_analysis_of_gene_list(target_gene_list,
                                               survival_time_list,
                                               cnv.table,
                                               distance = 2.5e+06,
                                               threshold = -2,
                                               deletion = TRUE,
                                               column_start = 11,
                                               time_of_death_column = 5,
                                               death_event_column = 6,
                                               start = TRUE,
                                               remove_NA = TRUE,
                                               Cytoband = FALSE)
length(test_function)
test_function[[1]]
dim(test_function[[1]])
dim(test_function[[2]])