####################
### Survival analysis of target genes in all cancers and Fishers exact test:
######################

#####################
###Survival analysis of ALL Tumour suppressors in all cancers using a for loop:
###################

################
### Create a list of lists for ALL Tumour suppressors:

gene_information_list

##Import Maria's Tumour suppressor list



##Which Tumour supressors in the list do I not have CNV information for?



##Comment: Do the tumour suppressors not in my CNV list have synonyms

##Check I have the start and end positions for all Tumour suppressors


###############
### Attach the table of all the tumours combined into the per cancer CNV list
length(threshold_short_cnv_list_loc)
length(threshold_cnv_list_loc)
length(threshold_CNV_all_table_loc)
length(list(threshold_CNV_all_table_loc))
threshold_CNV_all_table_loc[1:2, 1:12]

threshold_cnv_list_plus_all_loc<- c(threshold_cnv_list_loc, ALL=list(threshold_CNV_all_table_loc))
length(threshold_cnv_list_plus_all_loc)
dim(threshold_cnv_list_plus_all_loc[[38]])

names(threshold_cnv_list_plus_all_loc)
threshold_selected_cnv_list_plus_all_loc<- threshold_cnv_list_plus_all_loc[c(1,2,3,4,5,7,8,9,10,12,13,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,32,33,34,35,36,37,38)]
names(threshold_selected_cnv_list_plus_all_loc)
#############
### Make a table of all the tumours survival analysis combined and to the end of clinical_survival_list
#COAD and READ combined into COADREAD

length(clinical_survival_list)



################
### Create a for loop to go through all cancers and perform survival analysis on target genes and co-deletions:

##Create an empty list
survival_stats_cancer_list<- vector("list", length(threshold_short_cnv_list_loc))

for (i in 1:2){
  
  ##Get CNV table
  cnv.table<- threshold_short_cnv_list_loc[[i]]
  
  
  ##Get survival data
  survival_time_list<- clinical_survival_list[[i]]
  
  
  ## Perform survival analysis of co-deletions:
  print(paste(names(threshold_short_cnv_list_loc[i]), "Cancer type:"))
  survival_stats_list<- lapply(gene_information_list, function(x) survival_analysis_of_gene_list(target_gene_list = x, survival_time_list = survival_time_list, cnv.table = cnv.table, distance = 2.5e+06,
                                                                                                 threshold = -2, deletion = TRUE, time_of_death_column = 5, 
                                                                                                 death_event_column = 6, column_start = 11, start = TRUE, 
                                                                                                 remove_NA = TRUE, Cytoband = FALSE, print_to_screen = FALSE, 
                                                                                                 plot_graph = FALSE))
  
  ##Combine survival stats together and order by log-rank test p-value
  survival_stats_list_table<- do.call(rbind, survival_stats_list)
  survival_stats_list_table<- survival_stats_list_table[order(survival_stats_list_table$`p-value_logrank_test`, decreasing = FALSE),]
  #survival_stats_list_table[1:100, c(1,2,5,6,8,9)]
  
  ##Save file as .csv
  write.csv(survival_stats_list_table, file = paste0(names(threshold_short_cnv_list_loc[i]),"_co-deletion_survival_stats.csv"), quote = FALSE)
  
  ##Add to list
  survival_stats_cancer_list[[i]]<- survival_stats_list_table
  
  
  
}

## Save object:
saveRDS(survival_stats_cancer_list, file = "./R workspaces/survival_stats_cancer_list")




# target_gene_list<- gene_information_list[[1]]
# survival_time_list<- clinical_survival_list[[1]]
# cnv.table<- threshold_short_cnv_list_loc[[1]]
# 
# 
# test_apply<- lapply(gene_information_list, function(x) survival_analysis_of_gene_list(target_gene_list = x, survival_time_list = survival_time_list, cnv.table = cnv.table, distance = 2.5e+06,
#                                                                                       threshold = -2, deletion = TRUE, time_of_death_column = 5, 
#                                                                                       death_event_column = 6, column_start = 11, start = TRUE, 
#                                                                                       remove_NA = TRUE, Cytoband = FALSE, print_to_screen = FALSE, 
#                                                                                       plot_graph = FALSE))
# 












