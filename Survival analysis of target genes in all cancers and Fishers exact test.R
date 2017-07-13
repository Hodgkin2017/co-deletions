# ####################
# ### Survival analysis of target genes in all cancers and Fishers exact test:
# ######################
# 
# #####################
# ###Survival analysis of ALL Tumour suppressors in all cancers using a for loop:
# ###################
# 
# ################
# ### Create a list of lists for ALL Tumour suppressors:
# 
# gene_information_list
# 
# ##############################
# ##Import Maria's and Christophe's Tumour suppressor list:
# target_genes_Maria<- read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/From Maria/Tumour Suppressors list/tumorSuppressors/TSGlist_compiledChristopheVogelsteinCOSMIC_MS_20160914.txt", stringsAsFactors = FALSE, header = TRUE)
# dim(target_genes_Maria)
# target_genes_Maria
# 
# ##How many of Maria's gene are in my table:
# sum(target_genes_Maria$Gene %in% cnv.table$Gene.Symbol)
# which(!(target_genes_Maria$Gene %in% cnv.table$Gene.Symbol))
# ## Comment: Missing genes: 109 133 134
# 
# target_genes_Maria$Gene[c(109,133,134)]
# ##[1] "FAM123B" "MLL2"    "MLL3" 
# 
# target_genes_Maria$Gene
# grep("FAM123B",target_genes_Maria$Gene)
# grep("FAM123B",cnv.table$Gene.Symbol)
# grep("MLL2",cnv.table$Gene.Symbol)
# grep("MLL3",cnv.table$Gene.Symbol)
# 
# ## FAM123B also known as WTX and AMER1
# grep("WTX",cnv.table$Gene.Symbol)
# grep("AMER1",cnv.table$Gene.Symbol) #name in my list
# 
# ##MLL2 also known as KMT2B
# grep("KMT2B",cnv.table$Gene.Symbol) #name in my list
# 
# ##MLL3 also known a KMT2C
# grep("KMT2C",cnv.table$Gene.Symbol) #name in my list
# 
# ###Replace gene names in target_genes_Maria
# target_genes_Maria$Gene[c(109,133,134)]<- c("AMER1","KMT2B","KMT2C")
# 
# ##How many of Maria's gene are in my table:
# sum(target_genes_Maria$Gene %in% cnv.table$Gene.Symbol)
# which(!(target_genes_Maria$Gene %in% cnv.table$Gene.Symbol))
# ## Comment: Now no missing genes!
# 
# ## Create list of target genes and their information:
# gene_information_long_list<- create_target_gene_information_list(cnv_table = cnv.table, target_genes = target_genes_Maria$Gene)
# length(gene_information_long_list)
# gene_information_long_list[[4]]
# gene_information_long_list
# 
# ##Identify genes without a start:
# genes_without_start<- lapply(gene_information_long_list, function(x) !is.na(x$start))
# genes_without_start
# genes_without_start<- unlist(genes_without_start)
# genes_without_start
# 
# ##REmove genes with no start:
# gene_information_long_list<- gene_information_long_list[genes_without_start]
# length(gene_information_long_list)
# 
# ##Identify genes without a end:
# sapply(gene_information_long_list, function(x) !is.na(x$end))
# 
# ##save object
# saveRDS(gene_information_long_list, file = "./R workspaces/gene_information_long_list")
# 
# 
# 
# 
# ###############
# ### Attach the table of all the tumours combined into the per cancer CNV list
# # length(threshold_short_cnv_list_loc)
# # length(threshold_cnv_list_loc)
# # length(threshold_CNV_all_table_loc)
# # length(list(threshold_CNV_all_table_loc))
# # threshold_CNV_all_table_loc[1:2, 1:12]
# # 
# # threshold_cnv_list_plus_all_loc<- c(threshold_cnv_list_loc, ALL=list(threshold_CNV_all_table_loc))
# # length(threshold_cnv_list_plus_all_loc)
# # dim(threshold_cnv_list_plus_all_loc[[38]])
# # 
# # names(threshold_cnv_list_plus_all_loc)
# # threshold_selected_cnv_list_plus_all_loc<- threshold_cnv_list_plus_all_loc[c(1,2,3,4,5,7,8,9,10,12,13,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,32,33,34,35,36,37,38)]
# # names(threshold_selected_cnv_list_plus_all_loc)
# # names(clinical_survival_long_list)
# # length(threshold_selected_cnv_list_plus_all_loc)
# 
# ##Create one large table with all cancer information CNV information in it:
# names(clinical_survival_long_list)
# CNV_all_table<-join.cnv.datasets(threshold_cnv_list, column = 4, data.sets = c("ACC", "BLCA", "BRCA", "CESC",
#                                                                                "CHOL", "COADREAD", "DLBC", "ESCA",
#                                                                                "GBM", "HNSC", "KICH", "KIRC",
#                                                                                "KIRP", "LAML", "LGG", "LIHC",
#                                                                                "LUAD", "LUSC", "MESO", "OV",
#                                                                                "PAAD", "PCPG", "PRAD", "SARC",
#                                                                                "SKCM", "STAD", "TGCT", "THCA",
#                                                                                "THYM", "UCEC", "UCS", "UVM"))
# dim(CNV_all_table)
# CNV_all_table[1:2, 1:10]
# 
# ##Append chromosomal location information to it:
# CNV_all_table_selected_loc<- chromosomal_location(CNV_all_table)
# dim(CNV_all_table_selected_loc)
# CNV_all_table_selected_loc[1:2, 1:12]
# 
# ## Create list of CNV data matching survival data list (clinical_survival_long_list)
# threshold_selected_cnv_list_loc<- threshold_cnv_list_loc[c(1,2,3,4,5,7,8,9,10,12,13,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,32,33,34,35,36,37)]
# length(threshold_selected_cnv_list_loc)
# names(threshold_selected_cnv_list_loc)
# 
# ## Append ALL Cancer CNV table onto list of cancer CNVs
# threshold_selected_cnv_list_plus_all_loc<- c(threshold_selected_cnv_list_loc, ALL=list(CNV_all_table_selected_loc))
# length(threshold_selected_cnv_list_plus_all_loc)
# names(threshold_selected_cnv_list_plus_all_loc)
# 
# #############
# ### Make a table of all the tumours survival analysis combined and add to the end of clinical_survival_list
# #COAD and READ combined into COADREAD
# 
# length(clinical_survival_long_list)
# names(clinical_survival_long_list)
# 
# clinical_survival_long_all_table<-do.call(rbind, clinical_survival_long_list)
# lapply(clinical_survival_long_list, function(x) ncol(x))
# 
# ##Comment: UCS, PCPG, MESO, ACC do not have disease free survival?!
# ##Fix by adding empty new_tumor_days column to tables with missing data
# clinical_survival_long_list$ACC_clinical.tsv[1:2,]
# clinical_survival_long_list$BRCA_clinical.tsv[1:2,]
# 
# ##Add blank columns to these dataframes
# clinical_survival_long_list$ACC_clinical.tsv<- cbind(patient_IDs = clinical_survival_long_list$ACC_clinical.tsv[,1],
#       new_tumor_days = rep(NA,nrow(clinical_survival_long_list$ACC_clinical.tsv)),
#       clinical_survival_long_list$ACC_clinical.tsv[,2:6])
# 
# clinical_survival_long_list$MESO_clinical.tsv<- cbind(patient_IDs = clinical_survival_long_list$MESO_clinical.tsv[,1],
#                                                      new_tumor_days = rep(NA,nrow(clinical_survival_long_list$MESO_clinical.tsv)),
#                                                      clinical_survival_long_list$MESO_clinical.tsv[,2:6])
# 
# clinical_survival_long_list$PCPG_clinical.tsv<- cbind(patient_IDs = clinical_survival_long_list$PCPG_clinical.tsv[,1],
#                                                      new_tumor_days = rep(NA,nrow(clinical_survival_long_list$PCPG_clinical.tsv)),
#                                                      clinical_survival_long_list$PCPG_clinical.tsv[,2:6])
# clinical_survival_long_list$UCS_clinical.tsv<- cbind(patient_IDs = clinical_survival_long_list$UCS_clinical.tsv[,1],
#                                                      new_tumor_days = rep(NA,nrow(clinical_survival_long_list$UCS_clinical.tsv)),
#                                                      clinical_survival_long_list$UCS_clinical.tsv[,2:6])
# 
# 
# sapply(clinical_survival_long_list, function(x) ncol(x))
# 
# ##Make one large survival table:
# clinical_survival_long_all_table<-do.call(rbind, clinical_survival_long_list)
# dim(clinical_survival_long_all_table)
# length(unique(clinical_survival_long_all_table$patient_IDs))
# 
# ##Remove non-unique patient IDS:
# duplicated_patient_IDs<- !duplicated(clinical_survival_long_all_table$patient_IDs)
# clinical_survival_long_all_table<- clinical_survival_long_all_table[duplicated_patient_IDs,]
# dim(clinical_survival_long_all_table)
# 
# ##Append clinical data for all tumours to list:
# clinical_long_plus_all_survival<- c(clinical_survival_long_list, ALL=list(clinical_survival_long_all_table))
# length(clinical_long_plus_all_survival)
# names(clinical_long_plus_all_survival)

##Save objects
# saveRDS(threshold_selected_cnv_list_plus_all_loc, file = "./R workspaces/threshold_selected_cnv_list_plus_all_loc")
# saveRDS(clinical_long_plus_all_survival, file = "./R workspaces/clinical_long_plus_all_survival")



################
### Create a for loop to go through all cancers and perform survival analysis on target genes and co-deletions:

##CNV list = threshold_selected_cnv_list_plus_all_loc
## clinical list = clinical_long_plus_all_survival
#short_gene_information_list<- gene_information_list[1]

##Create an empty list
survival_stats_cancer_list<- vector("list", length(threshold_selected_cnv_list_plus_all_loc))

for (i in 1:length(threshold_selected_cnv_list_plus_all_loc)){
  
  ##Get CNV table
  cnv.table<- threshold_selected_cnv_list_plus_all_loc[[i]]
  
  
  ##Get survival data
  survival_time_list<- clinical_long_plus_all_survival[[i]]
  
  
  ## Perform survival analysis of co-deletions:
  print(paste(names(threshold_selected_cnv_list_plus_all_loc[i]), "Cancer type:"))
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
  write.csv(survival_stats_list_table, file = paste0(names(threshold_selected_cnv_list_plus_all_loc[i]),"_co-deletion_overall_survival_stats.csv"), quote = FALSE)
  
  ##Add to list
  survival_stats_cancer_list[[i]]<- survival_stats_list_table
  
  
  
}

#sapply(survival_stats_cancer_list, function(x) nrow(x))
#length(survival_stats_cancer_list)


## Save object:
saveRDS(survival_stats_cancer_list, file = "./R workspaces/survival_stats_overall_survival_cancer_list")


#########################
### Repeat loop but perform disease free survival analysis:

##Create an empty list
survival_stats_cancer_list2<- vector("list", length(threshold_selected_cnv_list_plus_all_loc))

for (i in 1:2){
  
  ##Get CNV table
  cnv.table<- threshold_selected_cnv_list_plus_all_loc[[i]]
  
  
  ##Get survival data
  survival_time_list<- clinical_long_plus_all_survival[[i]]
  
  
  ## Perform survival analysis of co-deletions:
  ##Disease free survival = column 7!!!
  print(paste(names(threshold_selected_cnv_list_plus_all_loc[i]), "Cancer type:"))
  survival_stats_list<- lapply(short_gene_information_list, function(x) survival_analysis_of_gene_list(target_gene_list = x, survival_time_list = survival_time_list, cnv.table = cnv.table, distance = 2.5e+06,
                                                                                                       threshold = -2, deletion = TRUE, time_of_death_column = 7, 
                                                                                                       death_event_column = 6, column_start = 11, start = TRUE, 
                                                                                                       remove_NA = TRUE, Cytoband = FALSE, print_to_screen = FALSE, 
                                                                                                       plot_graph = FALSE))
  
  ##Combine survival stats together and order by log-rank test p-value
  survival_stats_list_table<- do.call(rbind, survival_stats_list)
  survival_stats_list_table<- survival_stats_list_table[order(survival_stats_list_table$`p-value_logrank_test`, decreasing = FALSE),]
  #survival_stats_list_table[1:100, c(1,2,5,6,8,9)]
  
  ##Save file as .csv
  write.csv(survival_stats_list_table, file = paste0(names(threshold_selected_cnv_list_plus_all_loc[i]),"_co-deletion_disease_free_survival_stats.csv"), quote = FALSE)
  
  ##Add to list
  survival_stats_cancer_list2[[i]]<- survival_stats_list_table
  
  
  
}

## Save object:
saveRDS(survival_stats_cancer_list2, file = "./R workspaces/survival_stats_disease_free_cancer_list")


###########
###Delete cancers from list that dont have correct disease free survival data:
# UCS, PCPG, MESO, ACC



###################
###Perform fishers exact test per co-deletion (between deletion categories:1, 2, 3 and 4).

## Create a vector to store the results of the test:

fishers_co_deletion_test_results_list<- vector("list", length(survival_stats_cancer_list))

for (i in 1: length(survival_stats_cancer_list)){
  
  survival_stats_table<- survival_stats_cancer_list[[i]]
  
  ##Create an empty table for fishers exact test results:
  fishers_test_results<- data.frame(matrix(NA, ncol = 4, nrow = nrow(survival_stats_table)))
  colnames(fishers_test_results)<- c("target_gene", "proximal_gene", "co-deletion_p-value", "co-deletion_Odds")
  
  ##Perform Fishers exact test between target gene and co-deleted gene (deletion categories: 1,2,3,and 4).
  for(j in 1: nrow(survival_stats_table)) {
    
    ##Create an empty table for fisher exact test results:
    fishers_test_input_table<- data.frame(matrix(NA, ncol = 2, nrow = 2))
    
    ##Fill table:
    fishers_test_input_table[2,1]<- survival_stats_table[j,7]
    fishers_test_input_table[1,1]<- survival_stats_table[j,8]
    fishers_test_input_table[1,2]<- survival_stats_table[j,9]
    fishers_test_input_table[2,2]<- survival_stats_table[j,10]
    
    ##Perform Fishers exact test:
    fishers_test_co_deletions<- fisher.test(fishers_test_input_table)
    
    ##Store results of fishers exact test:
    fishers_test_results[j,1]<- survival_stats_table[j,1]
    fishers_test_results[j,2]<- survival_stats_table[j,2]
    fishers_test_results[j,3]<- fishers_test_co_deletions$p.value
    fishers_test_results[j,4]<- fishers_test_co_deletions$estimate
  }
  
  fishers_co_deletion_test_results_list[[i]]<- fishers_test_results
 
}


##Save object
saveRDS(fishers_co_deletion_test_results_list, file = "./R workspaces/fishers_co_deletion_test_results_list")

##Create .csv
write.csv(fishers_co_deletion_test_results_list, file = "fishers_co_deletion_test_results", quote = FALSE)


###################
###Perform fishers exact test per co-deletion between different cancers:.

## Create a vector to store the results of the test:

# fishers_test_per_cancer_results_list<- vector("list", length(survival_stats_cancer_list))
# 
# for (i in 1: length(survival_stats_cancer_list)-1){
#   
#   survival_stats_table<- survival_stats_cancer_list[[i]]
#   
#   ##Create an empty table for fishers exact test results:
#   fishers_test_results<- data.frame(matrix(NA, ncol = 16, nrow = nrow(survival_stats_table)))
#   colnames(fishers_test_results)<- c("target_gene", "proximal_gene", "cat1", "cat1", "cat2", "cat2",
#                                      "cat3", "cat3", "cat4", "cat4","co-deletion_p-value", "co-deletion_Odds", "co-deletion_p-value", "co-deletion_Odds", "co-deletion_p-value", "co-deletion_Odds")
#   
#   ##Perform Fishers exact test between target gene and co-deleted gene (deletion categories: 1,2,3,and 4).
#   for(j in 1: nrow(survival_stats_table)) {
#     
#     ##Create an empty table for fisher exact test results:
#     fishers_test_input_table<- data.frame(matrix(NA, ncol = 2, nrow = 2))
#     
#     ##Fill table:
#     fishers_test_input_table[1,1]<- survival_stats_table[j,8]
#     fishers_test_input_table[1,2]<- survival_stats_table[j,9]
#     fishers_test_input_table[2,1]<- survival_stats_table[j,7]
#     
#     fishers_test_input_table[2,2]<- survival_stats_table[j,10]
#     
#     ##Perform Fishers exact test:
#     fishers_test_co_deletions<- fisher.test(fishers_test_input_table)
#     
#     ##Store results of fishers exact test:
#     fishers_test_results[j,1]<- survival_stats_table[j,1]
#     fishers_test_results[j,2]<- survival_stats_table[j,2]
#     fishers_test_results[j,3]<- fishers_test_co_deletions$p.value
#     fishers_test_results[j,4]<- fishers_test_co_deletions$estimate
#   }
#   
#   fishers_test_per_cancer_results_list[[i]]<- fishers_test_results
#   
# }
# 
# 
# ##Save object











