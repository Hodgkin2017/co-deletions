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
# 
# 
# 
# ################
# ### Create a for loop to go through all cancers and perform survival analysis on target genes and co-deletions:
# 
# ##CNV list = threshold_selected_cnv_list_plus_all_loc
# ## clinical list = clinical_long_plus_all_survival
# #short_gene_information_list<- gene_information_list[1]
# 
# ##Create an empty list
# survival_stats_cancer_list<- vector("list", length(threshold_selected_cnv_list_plus_all_loc))
# 
# for (i in 1:length(threshold_selected_cnv_list_plus_all_loc)){
#   
#   ##Get CNV table
#   cnv.table<- threshold_selected_cnv_list_plus_all_loc[[i]]
#   
#   
#   ##Get survival data
#   survival_time_list<- clinical_long_plus_all_survival[[i]]
#   
#   
#   ## Perform survival analysis of co-deletions:
#   print(paste(names(threshold_selected_cnv_list_plus_all_loc[i]), "Cancer type:"))
#   survival_stats_list<- lapply(gene_information_long_list, function(x) survival_analysis_of_gene_list(target_gene_list = x, survival_time_list = survival_time_list, cnv.table = cnv.table, distance = 2.5e+06,
#                                                                                                  threshold = -2, deletion = TRUE, time_of_death_column = 5, 
#                                                                                                  death_event_column = 6, column_start = 11, start = TRUE, 
#                                                                                                  remove_NA = TRUE, Cytoband = FALSE, print_to_screen = FALSE, 
#                                                                                                  plot_graph = FALSE))
#   
#   ##Combine survival stats together and order by log-rank test p-value
#   survival_stats_list_table<- do.call(rbind, survival_stats_list)
#   survival_stats_list_table<- survival_stats_list_table[order(survival_stats_list_table$`p-value_logrank_test`, decreasing = FALSE),]
#   #survival_stats_list_table[1:100, c(1,2,5,6,8,9)]
#   
#   ##Save file as .csv
#   write.csv(survival_stats_list_table, file = paste0(names(threshold_selected_cnv_list_plus_all_loc[i]),"_co-deletion_overall_survival_stats_long.csv"), quote = FALSE)
#   
#   ##Add to list
#   survival_stats_cancer_list[[i]]<- survival_stats_list_table
#   
#   
#   
# }
# 
# #sapply(survival_stats_cancer_list, function(x) nrow(x))
# #length(survival_stats_cancer_list)
# 
# 
# ## Save object:
# saveRDS(survival_stats_cancer_list, file = "./R workspaces/survival_stats_overall_survival_cancer_list")
# 
# 
# #########################
# ### Repeat loop but perform disease free survival analysis:
# 
# ##Create an empty list
# survival_stats_cancer_list2<- vector("list", length(threshold_selected_cnv_list_plus_all_loc))
# 
# for (i in 1:length(survival_stats_cancer_list2)){
#   
#   ##Get CNV table
#   cnv.table<- threshold_selected_cnv_list_plus_all_loc[[i]]
#   
#   
#   ##Get survival data
#   survival_time_list<- clinical_long_plus_all_survival[[i]]
#   
#   
#   ## Perform survival analysis of co-deletions:
#   ##Disease free survival = column 7!!!
#   print(paste(names(threshold_selected_cnv_list_plus_all_loc[i]), "Cancer type:"))
#   survival_stats_list<- lapply(gene_information_long_list, function(x) survival_analysis_of_gene_list(target_gene_list = x, survival_time_list = survival_time_list, cnv.table = cnv.table, distance = 2.5e+06,
#                                                                                                        threshold = -2, deletion = TRUE, time_of_death_column = 7, 
#                                                                                                        death_event_column = 6, column_start = 11, start = TRUE, 
#                                                                                                        remove_NA = TRUE, Cytoband = FALSE, print_to_screen = FALSE, 
#                                                                                                        plot_graph = FALSE))
#   
#   ##Combine survival stats together and order by log-rank test p-value
#   survival_stats_list_table<- do.call(rbind, survival_stats_list)
#   survival_stats_list_table<- survival_stats_list_table[order(survival_stats_list_table$`p-value_logrank_test`, decreasing = FALSE),]
#   #survival_stats_list_table[1:100, c(1,2,5,6,8,9)]
#   
#   ##Save file as .csv
#   write.csv(survival_stats_list_table, file = paste0(names(threshold_selected_cnv_list_plus_all_loc[i]),"_co-deletion_disease_free_survival_stats_long.csv"), quote = FALSE)
#   
#   ##Add to list
#   survival_stats_cancer_list2[[i]]<- survival_stats_list_table
#   
#   
#   
# }
# 
# ## Save object:
# saveRDS(survival_stats_cancer_list2, file = "./R workspaces/survival_stats_disease_free_cancer_list")
# 
# 
# ###########
# ###Delete cancers from list that dont have correct disease free survival data:
# # UCS, PCPG, MESO, ACC
# 
# 
# 
# ###################
# ###Perform fishers exact test per co-deletion (between deletion categories:1, 2, 3 and 4).
# 
# ## Create a vector to store the results of the test:
# 
# fishers_co_deletion_test_results_list<- vector("list", length(survival_stats_cancer_list))
# 
# for (i in 1: length(survival_stats_cancer_list)){
#   
#   survival_stats_table<- survival_stats_cancer_list[[i]]
#   
#   ##Create an empty table for fishers exact test results:
#   fishers_test_results<- data.frame(matrix(NA, ncol = 4, nrow = nrow(survival_stats_table)))
#   colnames(fishers_test_results)<- c("target_gene", "proximal_gene", "co-deletion_p-value", "co-deletion_Odds")
#   
#   ##Perform Fishers exact test between target gene and co-deleted gene (deletion categories: 1,2,3,and 4).
#   for(j in 1: nrow(survival_stats_table)) {
#     
#     ##Create an empty table for fisher exact test results:
#     fishers_test_input_table<- data.frame(matrix(NA, ncol = 2, nrow = 2))
#     
#     ##Fill table:
#     fishers_test_input_table[2,1]<- survival_stats_table[j,7]
#     fishers_test_input_table[1,1]<- survival_stats_table[j,8]
#     fishers_test_input_table[1,2]<- survival_stats_table[j,9]
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
#   fishers_co_deletion_test_results_list[[i]]<- fishers_test_results
#  
# }
# 
# 
# ##Save object
# saveRDS(fishers_co_deletion_test_results_list, file = "./R workspaces/fishers_co_deletion_test_results_list")
# 
# 
# ###################
# ###Perform fishers exact test per co-deletion between different cancers:.
# 
# ##Data: 
# # survival_stats_overall_survival_cancer_list<- readRDS(file = "./R workspaces/survival_stats_overall_survival_cancer_list")
# # length(survival_stats_overall_survival_cancer_list)
# # dim(survival_stats_overall_survival_cancer_list[[33]])
# # survival_stats_overall_survival_cancer_list[[33]][1:2,]
# # survival_stats_cancer_list<- survival_stats_overall_survival_cancer_list
# 
# 
# ## Create a vector to store the results of the test:
# 
# fishers_test_per_cancer_results_list<- vector("list", length(survival_stats_cancer_list))
# 
# for (i in 1: (length(survival_stats_cancer_list)-1)){
#   
#   survival_stats_table<- survival_stats_cancer_list[[i]]
#   
#   ##Create an empty table for fishers exact test results:
#   fishers_test_results<- data.frame(matrix(NA, ncol = 16, nrow = nrow(survival_stats_table)))
#   colnames(fishers_test_results)<- c("target_gene", "proximal_gene", "n_cat1_per_cancer", "n_cat1_all_cancers", 
#                                      "n_cat2_per_cancer", "n_cat2_all_cancers",
#                                      "n_cat3_per_cancer", "n_cat3_all_cancers",
#                                      "n_cat4_per_cancer", "n_cat4_all_cancers",
#                                      "co-deletion_p-value_cat2vs1", "co-deletion_Odds_cat2vs1",
#                                      "co-deletion_p-value_cat2vs3", "co-deletion_Odds_cat2vs3", 
#                                      "co-deletion_p-value_cat2vs4", "co-deletion_Odds_cat2vs4")
#   
#   ##Perform Fishers exact test between target gene and co-deleted gene (deletion categories: 1,2,3,and 4).
#   for(j in 1: nrow(survival_stats_table)) {
#     
#     ##Create an empty table for fisher exact test results:
#     fishers_test_input_table<- data.frame(matrix(NA, ncol = 2, nrow = 2))
#     
#     ##Get category data from other cancers
#     ##Cat1:
#     Cat1_sum<- sapply(survival_stats_cancer_list, function(x) x$number_of_samples_cat1[j]) %>%
#       .[-c(i,33)] %>%
#       sum()
#     ##Cat2:
#     Cat2_sum<- sapply(survival_stats_cancer_list, function(x) x$number_of_samples_cat2[j]) %>%
#       .[-c(i,33)] %>%
#       sum()
#     ##Cat3:
#     Cat3_sum<- sapply(survival_stats_cancer_list, function(x) x$number_of_samples_cat3[j]) %>%
#       .[-c(i,33)] %>%
#       sum()
#     ##Cat4:
#     Cat4_sum<- sapply(survival_stats_cancer_list, function(x) x$number_of_samples_cat4[j]) %>%
#       .[-c(i,33)] %>%
#       sum()
#     
#     ################
#     ###Cat 2 v's Cat1
#     ##Fill table for fishers exact test:
#     fishers_test_input_table[1,1]<- survival_stats_table[j,8]
#     fishers_test_input_table[1,2]<- Cat2_sum
#     
#     fishers_test_input_table[2,1]<- survival_stats_table[j,7]
#     fishers_test_input_table[2,2]<- Cat1_sum
#     
#     ##Perform Fishers exact test:
#     fishers_test_co_deletions<- fisher.test(fishers_test_input_table)
#     
#     ##Store results of fishers exact test:
#     fishers_test_results[j,1]<- survival_stats_table[j,1]
#     fishers_test_results[j,2]<- survival_stats_table[j,2]
#     fishers_test_results[j,3]<- survival_stats_table[j,7]
#     fishers_test_results[j,4]<- Cat1_sum
#     fishers_test_results[j,5]<- survival_stats_table[j,8]
#     fishers_test_results[j,6]<- Cat2_sum
#     fishers_test_results[j,7]<- survival_stats_table[j,9]
#     fishers_test_results[j,8]<- Cat3_sum
#     fishers_test_results[j,9]<- survival_stats_table[j,10]
#     fishers_test_results[j,10]<- Cat4_sum
#     fishers_test_results[j,11]<- fishers_test_co_deletions$p.value
#     fishers_test_results[j,12]<- fishers_test_co_deletions$estimate
#     
#     ###############
#     ###Cat 2 v's Cat 3
#     ##Fill table for fishers exact test:
#     fishers_test_input_table[2,1]<- survival_stats_table[j,9]
#     fishers_test_input_table[2,2]<- Cat3_sum
#     
#     ##Perform Fishers exact test:
#     fishers_test_co_deletions<- fisher.test(fishers_test_input_table)
#     
#     ##Store results of fishers exact test:
#     fishers_test_results[j,13]<- fishers_test_co_deletions$p.value
#     fishers_test_results[j,14]<- fishers_test_co_deletions$estimate
#     
#     ###############
#     ###Cat 2 v's Cat 4
#     ##Fill table for fishers exact test:
#     fishers_test_input_table[2,1]<- survival_stats_table[j,10]
#     fishers_test_input_table[2,2]<- Cat4_sum
#     
#     ##Perform Fishers exact test:
#     fishers_test_co_deletions<- fisher.test(fishers_test_input_table)
#     
#     ##Store results of fishers exact test:
#     fishers_test_results[j,15]<- fishers_test_co_deletions$p.value
#     fishers_test_results[j,16]<- fishers_test_co_deletions$estimate
#     
#     
#     
#   }
#   
#   fishers_test_per_cancer_results_list[[i]]<- fishers_test_results
#   
# }
# 
# names(fishers_test_per_cancer_results_list)<- names(threshold_selected_cnv_list_plus_all_loc)
# ##Save object
# saveRDS(fishers_test_per_cancer_results_list, file = "./R workspaces/fishers_test_per_cancer_results_list")

###################
### Perform survival analysis between catagotories 1 and 2 only:
##################


#########
###F2: Function to Create a table of survival p-value statistics for co-deletions or co-amplifications

# survival_analysis_of_gene_list_cat1_and_2<- function(target_gene_list,
#                                                      survival_time_list,
#                                                      cnv.table,
#                                                      distance = 2.5e+06,
#                                                      threshold = -2,
#                                                      deletion = TRUE,
#                                                      time_of_death_column, 
#                                                      death_event_column,
#                                                      column_start = 11,
#                                                      start = TRUE, 
#                                                      remove_NA = TRUE,
#                                                      Cytoband = FALSE,
#                                                      print_to_screen = FALSE, 
#                                                      plot_graph = FALSE){
#   
#   #############
#   ## Get target gene name, chromosome start and end locations and set interval upstream and downstream 
#   #of target gene used to identify genes in close proximety to target gene. 
#   target_gene<- target_gene_list[[1]]
#   Chromosome<- target_gene_list[[2]]
#   selection_criteria<- c(target_gene_list[[4]] - distance, target_gene_list[[5]] + distance)
#   print(target_gene)
#   
#   #############
#   ### Get genes surrounding target gene
#   
#   ##Remove genes with no start site
#   if (remove_NA == TRUE){
#     
#     cnv.table<- cnv.table %>%
#       dplyr::filter(!is.na(start))
#   } 
#   ## Select chromosome
#   if (Chromosome[1] > 0 & start == FALSE){
#     
#     ##Select Chromosome of interest and convert CNV data to matrix:
#     matrix<- cnv.table %>% dplyr::filter(CHR %in% Chromosome) 
#     cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
#     rownames(cnv.matrix)<- matrix$Gene.Symbol
#     
#   }
#   
#   ## Select chromosome and region of interest
#   if (start == TRUE & Chromosome[1] > 0){
#     
#     ##Select Chromosome of interest and convert CNV data to matrix:
#     matrix<- cnv.table %>% dplyr::filter(CHR %in% Chromosome) 
#     
#     ##Select Chromosomal region of interest and convert CNV data to matrix:
#     matrix<- matrix %>%
#       dplyr::filter(start >= selection_criteria[1], end <= selection_criteria[2])
#     
#     ##Convert Chromosome and region of interest into matrix:
#     cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
#     rownames(cnv.matrix)<- matrix$Gene.Symbol
#     
#   } else if (Cytoband == TRUE){
#     ##Select Cytobands of interest and convert CNV data to matrix:
#     matrix<- cnv.table %>% dplyr::filter(Cytoband %in% selection_criteria) 
#     cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
#     rownames(cnv.matrix)<- matrix$Gene.Symbol
#     
#   } else {
#     ##Convert ALL CNV data to matrix:
#     cnv.matrix<- as.matrix(cnv.table[,column_start:ncol(cnv.table)])
#     rownames(cnv.matrix)<- cnv.table$Gene.Symbol
#   }
#   
#   ##Create a binary matrix of CNV data such that deletions (deletion = TRUE) or 
#   #amplifications (deletion = FALSE) below (for deletions) or above (for amplifications) a threshold = 1 
#   if (deletion == TRUE) {
#     
#     cnv.matrix<- ifelse(cnv.matrix <= threshold, 1, 0)
#     
#   } else {
#     
#     cnv.matrix<- ifelse(cnv.matrix >= threshold, 1, 0)
#   }
#   
#   
#   ###########
#   ###Create a table with genes in rows and tumours on columns and each entry has a number from 1 to 4 depending on which 
#   #deletion/amplification catagory it belongs to: 
#   #1 = Deletion of target gene only. 
#   #2 = Deletion of both target gene and proximal gene (co-deleted gene or interest)
#   #3 = Deletion of proximal gene only.
#   #4 = Deletion in neither target gene or proximal gene.
#   
#   ##Add gene names as new column
#   gene_names_cnv.matrix<- cbind.data.frame(gene = row.names(cnv.matrix) ,cnv.matrix)
#   gene_names_cnv.matrix$gene<-as.character(gene_names_cnv.matrix$gene)
#   gene_names<- row.names(cnv.matrix)
#   
#   ##########
#   ###Function to categorise gene deletions:
#   ## If I remove it from within this function it no longer works for some strange reason!
#   categorise_deletion_type_function<- function(gene_names, gene_names_cnv.matrix){
#     
#     ##Determine if target gene is deleted
#     target_gene_deleted<- gene_names_cnv.matrix %>%
#       dplyr::filter(gene == target_gene) %>%
#       dplyr::select(-gene)
#     
#     target_gene_deleted<- as.vector(target_gene_deleted == 1)
#     
#     ##Determine if proximal gene is deleted
#     proximal_gene_deleted<- gene_names_cnv.matrix %>%
#       dplyr::filter(gene == gene_names) %>%
#       dplyr::select(-gene)
#     
#     proximal_gene_deleted<- as.vector(proximal_gene_deleted == 1)
#     
#     ##Determine if target gene and not proximal gene is deleted(category 1)
#     deletion_category1<- proximal_gene_deleted - target_gene_deleted < 0
#     
#     ##Determine if target gene and proximal gene are deleted together (category 2)
#     deletion_category2<- target_gene_deleted + proximal_gene_deleted > 1
#     
#     ##Determine if proximal gene is deleted without target gene (category 3)
#     deletion_category3<- target_gene_deleted - proximal_gene_deleted < 0
#     
#     ##Determine if neither proximal gene or target gene is deleted (category 4)
#     deletion_category4<- target_gene_deleted + proximal_gene_deleted == 0
#     
#     ## Join deletion_category vectors into table then take the non-NA value for each tumour:
#     all_deletion_category_table<- cbind(deletion_category1, deletion_category2, deletion_category3, deletion_category4)
#     
#     final_deletion_category<- max.col(all_deletion_category_table,ties.method="first")
#     
#     return(final_deletion_category)
#     
#   }
#   
#   ## Use categorise_deletion_type_function to loop over all genes around target gene and 
#   #categorise each tumour into group 1,2,3 or 4.
#   deletion_category_table<- lapply(gene_names, function(x) categorise_deletion_type_function(x, gene_names_cnv.matrix))
#   
#   
#   deletion_category_table<- do.call(rbind, deletion_category_table)
#   colnames(deletion_category_table)<- colnames(cnv.matrix)
#   rownames(deletion_category_table)<- rownames(cnv.matrix)
#   
#   
#   #########################################
#   ### Join tumour category table above to tumour survival data:
#   
#   deletion_category<-t(deletion_category_table)
#   
#   ##Convert patient IDs in deletion_category so they match with the patient IDS in the clinical table
#   deletion_category_patient_ID<- rownames(deletion_category) %>%
#     substr(0, 12) %>%
#     gsub("\\.", "-", .) %>%
#     cbind.data.frame(patient_IDs = ., deletion_category)
#   
#   deletion_category_patient_ID$patient_IDs<- as.character(deletion_category_patient_ID$patient_IDs)
#   
#   ##Extract patient ID, time of death and if death occured columns from survival table:
#   short_survival_time_list<- survival_time_list[,c(1, time_of_death_column, death_event_column)]
#   
#   ##Join clinical table with tumour deletion category table:
#   clinical_survival_deletion_category<- dplyr::full_join(short_survival_time_list, deletion_category_patient_ID, by = "patient_IDs") %>%
#     dplyr::select(-patient_IDs)
#   
#   ## Remove entries with NA deletion categories due to missing clinical data
#   vars<- colnames(deletion_category_patient_ID[-1])
#   deletion_category_proximal_genes<- clinical_survival_deletion_category %>%
#     dplyr::select(one_of(vars))
#   
#   clinical_survival_deletion_category<- clinical_survival_deletion_category[complete.cases(deletion_category_proximal_genes), ]
#   
#   
#   ################################################
#   ### Survival analysis for categories 1 and 2 only:
#   
#   # ##Create Surv object and remove NA values
#   # death_time<- clinical_survival_deletion_category[,1]
#   # death_event<- clinical_survival_deletion_category[,2]
#   # surv_object<- Surv(death_time, death_event==1)
#   # ##Remove entries with NA:
#   # NA_object<- !is.na(surv_object)
#   # death_time<- death_time[NA_object]
#   # death_event<- death_event[NA_object]
#   # surv_object<- surv_object[NA_object]
#   
#   ##Create empty table to store stats:
#   survival_stats_table<- data.frame(matrix(NA, ncol = 14, nrow = ncol(clinical_survival_deletion_category) -2))
#   names(survival_stats_table)<-c("target_gene", "proximal_gene", "p-value_Likelihood_ratio_test",
#                                  "p-value_Wald_test", "p-value_logrank_test", "Hazard_ratio",
#                                  "number_of_samples_cat1", "number_of_samples_cat2", "number_of_samples_cat3",
#                                  "number_of_samples_cat4","mean_survival_cat_1","mean_survival_cat_2","mean_survival_cat_3"
#                                  ,"mean_survival_cat_4")
#   
#   ## Loop through genes around target gene and calculate survival stats:
#   for (i in 1:nrow(survival_stats_table)){
#     
#     ##Get covariable object for survfit:
#     covariable_object<- clinical_survival_deletion_category[,i+2]
#     ##Keep only entries with 1 and 2
#     covariable_object_cat1_and_2<- which(covariable_object %in% c(1,2))
#     covariable_object<- covariable_object[covariable_object_cat1_and_2]
#     
#     ##Create Surv object if there are any tumours in catagory 1 and 2
#     if(length(covariable_object) > 1){
#       death_time<- clinical_survival_deletion_category[covariable_object_cat1_and_2,1]
#       death_event<- clinical_survival_deletion_category[covariable_object_cat1_and_2,2]
#       surv_object<- Surv(death_time, death_event==1)
#       
#       ##Fit Kaplain meier graph to one co-variable to compare data with :
#       fittedSurv <- survfit(surv_object~covariable_object, na.action = na.exclude)
#       fittedSurv_mean<- survival:::survmean(fittedSurv, rmean="individual") 
#       #survival:::survmean(fittedSurv, rmean="common")
#       
#       ##Get names of categories:
#       df.categ <- cbind(sapply(names(fittedSurv$strata), function(x) strsplit(x,"=")[[1]][2]),
#                         fittedSurv$n)
#       df.categ <- data.frame(df.categ)
#       
#       ##Get category names for one variable
#       categNames <- apply(df.categ, 1, function(x) paste0(x[1]," (",x[2],")"))
#       
#       if (plot_graph == TRUE) {
#         #save as tiff in folder....
#         plot(fittedSurv, main=plotTitle,
#              xlab="Time (days)", ylab=ylabel,
#              col=brewer.pal(9,"Set1"), mark.time=T)
#         legend("topright", legend=categNames,
#                col=brewer.pal(9,"Set1"),
#                lwd=2, cex=0.9)
#       }
#       
#       if(print_to_screen == TRUE) {
#         print("Chi-sq test:")
#         print(survdiff(surv_object~covariable_object,rho = 0))
#         
#         print("Cox PH test:")
#         print(summary(coxph(surv_object~covariable_object)))
#       }
#       coxfit <- coxph(surv_object~covariable_object)
#       
#       # if (plot_graph == TRUE) {
#       #   text(1000,0,labels=paste0("HR=",round(exp(summary(coxfit)$coefficients[1]),2),"; p=",
#       #                             round(summary(coxfit)$logtest[3],3)))
#       #   end of dev...
#       # }
#       
#       ## Save survival stats in table:
#       survival_stats_table[i, 1]<- target_gene
#       survival_stats_table[i, 2]<- colnames(clinical_survival_deletion_category[i+2])
#       survival_stats_table[i, 3]<- round(summary(coxfit)$logtest[3],2)
#       survival_stats_table[i, 4]<- round(summary(coxfit)$waldtest[3],2)
#       survival_stats_table[i, 5]<- round(summary(coxfit)$sctest[3],2)
#       survival_stats_table[i, 6]<- round(summary(coxfit)$coefficients[2],2)
#       survival_stats_table[i, 7]<- sum(grepl(1, covariable_object))
#       survival_stats_table[i, 8]<- sum(grepl(2, covariable_object))
#       survival_stats_table[i, 9]<- sum(grepl(3, covariable_object))
#       survival_stats_table[i, 10]<- sum(grepl(4, covariable_object))
#       ## If all tumours have the same category then a vector is crteated instead of a table:
#       if(is.null(ncol(fittedSurv_mean$matrix))){
#         ## Which category does all the tumours have:
#         category<- which(survival_stats_table[i, 7:10] > 0)
#         ## Get mean survival for that category:
#         survival_stats_table[i, 10+category]<- fittedSurv_mean$matrix[5]
#       }else{
#         list_of_mean_survival<- as.list(fittedSurv_mean$matrix[,5])
#         if(!is.null(list_of_mean_survival$`covariable_object=1`)){
#           survival_stats_table[i, 11]<- list_of_mean_survival$`covariable_object=1`
#         }
#         if(!is.null(list_of_mean_survival$`covariable_object=2`)){
#           survival_stats_table[i, 12]<- list_of_mean_survival$`covariable_object=2`
#         }
#         if(!is.null(list_of_mean_survival$`covariable_object=3`)){
#           survival_stats_table[i, 13]<- list_of_mean_survival$`covariable_object=3`
#         }
#         if(!is.null(list_of_mean_survival$`covariable_object=4`)){
#           survival_stats_table[i, 14]<- list_of_mean_survival$`covariable_object=4`
#         }
#       }
#     }
#     #   ## If less than one category in survival analysis need to adjust values saved:
#     #   if(is.null(ncol(fittedSurv_mean$matrix))){
#     #     #survival_stats_table[i, 7]<- unique(covariable_object) #categories
#     #     survival_stats_table[i, 8]<- paste(fittedSurv_mean$matrix[5], collapse = " ") #mean
#     #   } else {
#     #     #survival_stats_table[i, 7]<- paste(rownames(fittedSurv_mean$matrix), collapse = " ") #categories
#     #     #survival_stats_table[i, 8]<- paste(fittedSurv_mean$matrix[,5], collapse = " ") #mean
#     #     if (survival_stats_table[i, 7] != 0){
#     #       survival_stats_table[i, 11]<- 
#     #   }
#     #   
#     #   
#     
#   }
#   
#   return(survival_stats_table)
#   
# }



##############
### Create a for loop to go through all cancers and perform survival analysis on target genes and co-deletions:
##Cat 1 and 2 only:

##Create an empty list
# survival_stats_cancer_list3<- vector("list", length(threshold_selected_cnv_list_plus_all_loc))
# 
# for (i in 14:length(threshold_selected_cnv_list_plus_all_loc)){
#   
#   ##Get CNV table
#   cnv.table<- threshold_selected_cnv_list_plus_all_loc[[i]]
#   
#   
#   ##Get survival data
#   survival_time_list<- clinical_long_plus_all_survival[[i]]
#   
#   
#   ## Perform survival analysis of co-deletions:
#   print(paste(names(threshold_selected_cnv_list_plus_all_loc[i]), "Cancer type:"))
#   survival_stats_list<- lapply(gene_information_long_list, function(x) survival_analysis_of_gene_list_cat1_and_2(target_gene_list = x, survival_time_list = survival_time_list, cnv.table = cnv.table, distance = 2.5e+06,
#                                                                                                                  threshold = -2, deletion = TRUE, time_of_death_column = 5, 
#                                                                                                                  death_event_column = 6, column_start = 11, start = TRUE, 
#                                                                                                                  remove_NA = TRUE, Cytoband = FALSE, print_to_screen = FALSE, 
#                                                                                                                  plot_graph = FALSE))
#   
#   ##Combine survival stats together and order by log-rank test p-value
#   survival_stats_list_table<- do.call(rbind, survival_stats_list)
#   survival_stats_list_table<- survival_stats_list_table[order(survival_stats_list_table$`p-value_logrank_test`, decreasing = FALSE),]
#   #survival_stats_list_table[1:100, c(1,2,5,6,8,9)]
#   
#   ##Save file as .csv
#   write.csv(survival_stats_list_table, file = paste0(names(threshold_selected_cnv_list_plus_all_loc[i]),"_co-deletion_overall_survival_stats_cat1_and_2.csv"), quote = FALSE)
#   
#   ##Add to list
#   survival_stats_cancer_list3[[i]]<- survival_stats_list_table
#   
#   
#   
# }
# 
# #sapply(survival_stats_cancer_list, function(x) nrow(x))
# #length(survival_stats_cancer_list)
# 
# 
# ## Save object:
# saveRDS(survival_stats_cancer_list3, file = "./R workspaces/survival_stats_overall_survival_cancer_list_cat1_and_2")

######################
### Fixed problem with NA values for tumour types with no disease free survival data:


# t<-sapply(gene_information_long_list, function(x) x[[1]])
# t
# grep("NF1", t)
# 
# 
# target_gene_list<- gene_information_long_list[[48]]
# survival_time_list<- clinical_long_plus_all_survival[[2]]
# cnv.table<- threshold_selected_cnv_list_plus_all_loc[[2]]
# distance = 2.5e+06
# threshold = -2
# deletion = TRUE
# time_of_death_column = 7
# death_event_column = 6
# column_start = 11
# start = TRUE 
# remove_NA = TRUE
# Cytoband = FALSE
# print_to_screen = FALSE 
# plot_graph = FALSE








###########
survival_analysis_of_gene_list_cat1_and_2<- function(target_gene_list,
                                                     survival_time_list,
                                                     cnv.table,
                                                     distance = 2.5e+06,
                                                     threshold = -2,
                                                     deletion = TRUE,
                                                     time_of_death_column, 
                                                     death_event_column,
                                                     column_start = 11,
                                                     start = TRUE, 
                                                     remove_NA = TRUE,
                                                     Cytoband = FALSE,
                                                     print_to_screen = FALSE, 
                                                     plot_graph = FALSE){
  
  #############
  ## Get target gene name, chromosome start and end locations and set interval upstream and downstream 
  #of target gene used to identify genes in close proximety to target gene. 
  target_gene<- target_gene_list[[1]]
  Chromosome<- target_gene_list[[2]]
  selection_criteria<- c(target_gene_list[[4]] - distance, target_gene_list[[5]] + distance)
  print(target_gene)
  
  #############
  ### Get genes surrounding target gene
  
  ##Remove genes with no start site
  if (remove_NA == TRUE){
    
    cnv.table<- cnv.table %>%
      dplyr::filter(!is.na(start))
  } 
  ## Select chromosome
  if (Chromosome[1] > 0 & start == FALSE){
    
    ##Select Chromosome of interest and convert CNV data to matrix:
    matrix<- cnv.table %>% dplyr::filter(CHR %in% Chromosome) 
    cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
    rownames(cnv.matrix)<- matrix$Gene.Symbol
    
  }
  
  ## Select chromosome and region of interest
  if (start == TRUE & Chromosome[1] > 0){
    
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
  
  
  ###########
  ###Create a table with genes in rows and tumours on columns and each entry has a number from 1 to 4 depending on which 
  #deletion/amplification catagory it belongs to: 
  #1 = Deletion of target gene only. 
  #2 = Deletion of both target gene and proximal gene (co-deleted gene or interest)
  #3 = Deletion of proximal gene only.
  #4 = Deletion in neither target gene or proximal gene.
  
  ##Add gene names as new column
  gene_names_cnv.matrix<- cbind.data.frame(gene = row.names(cnv.matrix) ,cnv.matrix)
  gene_names_cnv.matrix$gene<-as.character(gene_names_cnv.matrix$gene)
  gene_names<- row.names(cnv.matrix)
  
  ##########
  ###Function to categorise gene deletions:
  ## If I remove it from within this function it no longer works for some strange reason!
  categorise_deletion_type_function<- function(gene_names, gene_names_cnv.matrix){
    
    ##Determine if target gene is deleted
    target_gene_deleted<- gene_names_cnv.matrix %>%
      dplyr::filter(gene == target_gene) %>%
      dplyr::select(-gene)
    
    target_gene_deleted<- as.vector(target_gene_deleted == 1)
    
    ##Determine if proximal gene is deleted
    proximal_gene_deleted<- gene_names_cnv.matrix %>%
      dplyr::filter(gene == gene_names) %>%
      dplyr::select(-gene)
    
    proximal_gene_deleted<- as.vector(proximal_gene_deleted == 1)
    
    ##Determine if target gene and not proximal gene is deleted(category 1)
    deletion_category1<- proximal_gene_deleted - target_gene_deleted < 0
    
    ##Determine if target gene and proximal gene are deleted together (category 2)
    deletion_category2<- target_gene_deleted + proximal_gene_deleted > 1
    
    ##Determine if proximal gene is deleted without target gene (category 3)
    deletion_category3<- target_gene_deleted - proximal_gene_deleted < 0
    
    ##Determine if neither proximal gene or target gene is deleted (category 4)
    deletion_category4<- target_gene_deleted + proximal_gene_deleted == 0
    
    ## Join deletion_category vectors into table then take the non-NA value for each tumour:
    all_deletion_category_table<- cbind(deletion_category1, deletion_category2, deletion_category3, deletion_category4)
    
    final_deletion_category<- max.col(all_deletion_category_table,ties.method="first")
    
    return(final_deletion_category)
    
  }
  
  ## Use categorise_deletion_type_function to loop over all genes around target gene and 
  #categorise each tumour into group 1,2,3 or 4.
  deletion_category_table<- lapply(gene_names, function(x) categorise_deletion_type_function(x, gene_names_cnv.matrix))
  
  
  deletion_category_table<- do.call(rbind, deletion_category_table)
  colnames(deletion_category_table)<- colnames(cnv.matrix)
  rownames(deletion_category_table)<- rownames(cnv.matrix)
  
  
  #########################################
  ### Join tumour category table above to tumour survival data:
  
  deletion_category<-t(deletion_category_table)
  
  ##Convert patient IDs in deletion_category so they match with the patient IDS in the clinical table
  deletion_category_patient_ID<- rownames(deletion_category) %>%
    substr(0, 12) %>%
    gsub("\\.", "-", .) %>%
    cbind.data.frame(patient_IDs = ., deletion_category)
  
  deletion_category_patient_ID$patient_IDs<- as.character(deletion_category_patient_ID$patient_IDs)
  
  ##Extract patient ID, time of death and if death occured columns from survival table:
  short_survival_time_list<- survival_time_list[,c(1, time_of_death_column, death_event_column)]
  
  ##Join clinical table with tumour deletion category table:
  clinical_survival_deletion_category<- dplyr::full_join(short_survival_time_list, deletion_category_patient_ID, by = "patient_IDs") %>%
    dplyr::select(-patient_IDs)
  
  ## Remove entries with NA deletion categories due to missing clinical data
  vars<- colnames(deletion_category_patient_ID[-1])
  deletion_category_proximal_genes<- clinical_survival_deletion_category %>%
    dplyr::select(one_of(vars))
  
  clinical_survival_deletion_category<- clinical_survival_deletion_category[complete.cases(deletion_category_proximal_genes), ]
  
  
  ################################################
  ### Survival analysis for categories 1 and 2 only:
  
  # ##Create Surv object and remove NA values
  # death_time<- clinical_survival_deletion_category[,1]
  # death_event<- clinical_survival_deletion_category[,2]
  # surv_object<- Surv(death_time, death_event==1)
  # ##Remove entries with NA:
  # NA_object<- !is.na(surv_object)
  # death_time<- death_time[NA_object]
  # death_event<- death_event[NA_object]
  # surv_object<- surv_object[NA_object]
  
  ##Create empty table to store stats:
  survival_stats_table<- data.frame(matrix(NA, ncol = 14, nrow = ncol(clinical_survival_deletion_category) -2))
  names(survival_stats_table)<-c("target_gene", "proximal_gene", "p-value_Likelihood_ratio_test",
                                 "p-value_Wald_test", "p-value_logrank_test", "Hazard_ratio",
                                 "number_of_samples_cat1", "number_of_samples_cat2", "number_of_samples_cat3",
                                 "number_of_samples_cat4","mean_survival_cat_1","mean_survival_cat_2","mean_survival_cat_3"
                                 ,"mean_survival_cat_4")
  
  ## Loop through genes around target gene and calculate survival stats:
  for (i in 1:nrow(survival_stats_table)){
    
    ##Get covariable object for survfit:
    covariable_object<- clinical_survival_deletion_category[,i+2]
    ##Keep only entries with 1 and 2
    covariable_object_cat1_and_2<- which(covariable_object %in% c(1,2))
    covariable_object<- covariable_object[covariable_object_cat1_and_2]
    ##Remove values with death_time = NA
    death_time<- clinical_survival_deletion_category[covariable_object_cat1_and_2,1]
    death_event<- clinical_survival_deletion_category[covariable_object_cat1_and_2,2]
    death_time<- death_time[which(!is.na(death_time))]
    death_event<- death_event[which(!is.na(death_time))]
    covariable_object<- covariable_object[which(!is.na(death_time))]
    ##Remove values with death_event = NA
    death_time<- death_time[which(!is.na(death_event))]
    death_event<- death_event[which(!is.na(death_event))]
    covariable_object<- covariable_object[which(!is.na(death_event))]
    
    ##Create Surv object if there are any tumours in catagory 1 and 2
    if(length(covariable_object) > 1){
      # death_time<- clinical_survival_deletion_category[covariable_object_cat1_and_2,1]
      # death_event<- clinical_survival_deletion_category[covariable_object_cat1_and_2,2]
      surv_object<- Surv(death_time, death_event==1)
      
      ##Fit Kaplain meier graph to one co-variable to compare data with :
      fittedSurv <- survfit(surv_object~covariable_object, na.action = na.exclude)
      fittedSurv_mean<- survival:::survmean(fittedSurv, rmean="individual") 
      #survival:::survmean(fittedSurv, rmean="common")
      
      ##Get names of categories:
      df.categ <- cbind(sapply(names(fittedSurv$strata), function(x) strsplit(x,"=")[[1]][2]),
                        fittedSurv$n)
      df.categ <- data.frame(df.categ)
      
      ##Get category names for one variable
      categNames <- apply(df.categ, 1, function(x) paste0(x[1]," (",x[2],")"))
      
      if (plot_graph == TRUE) {
        #save as tiff in folder....
        plot(fittedSurv, main=plotTitle,
             xlab="Time (days)", ylab=ylabel,
             col=brewer.pal(9,"Set1"), mark.time=T)
        legend("topright", legend=categNames,
               col=brewer.pal(9,"Set1"),
               lwd=2, cex=0.9)
      }
      
      if(print_to_screen == TRUE) {
        print("Chi-sq test:")
        print(survdiff(surv_object~covariable_object,rho = 0))
        
        print("Cox PH test:")
        print(summary(coxph(surv_object~covariable_object)))
      }
      coxfit <- coxph(surv_object~covariable_object)
      
      # if (plot_graph == TRUE) {
      #   text(1000,0,labels=paste0("HR=",round(exp(summary(coxfit)$coefficients[1]),2),"; p=",
      #                             round(summary(coxfit)$logtest[3],3)))
      #   end of dev...
      # }
      
      ## Save survival stats in table:
      survival_stats_table[i, 1]<- target_gene
      survival_stats_table[i, 2]<- colnames(clinical_survival_deletion_category[i+2])
      survival_stats_table[i, 3]<- round(summary(coxfit)$logtest[3],2)
      survival_stats_table[i, 4]<- round(summary(coxfit)$waldtest[3],2)
      survival_stats_table[i, 5]<- round(summary(coxfit)$sctest[3],2)
      survival_stats_table[i, 6]<- round(summary(coxfit)$coefficients[2],2)
      survival_stats_table[i, 7]<- sum(grepl(1, covariable_object))
      survival_stats_table[i, 8]<- sum(grepl(2, covariable_object))
      survival_stats_table[i, 9]<- sum(grepl(3, covariable_object))
      survival_stats_table[i, 10]<- sum(grepl(4, covariable_object))
      ## If all tumours have the same category then a vector is crteated instead of a table:
      if(is.null(ncol(fittedSurv_mean$matrix))){
        ## Which category does all the tumours have:
        category<- which(survival_stats_table[i, 7:10] > 0)
        ## Get mean survival for that category:
        survival_stats_table[i, 10+category]<- fittedSurv_mean$matrix[5]
      }else{
        list_of_mean_survival<- as.list(fittedSurv_mean$matrix[,5])
        if(!is.null(list_of_mean_survival$`covariable_object=1`)){
          survival_stats_table[i, 11]<- list_of_mean_survival$`covariable_object=1`
        }
        if(!is.null(list_of_mean_survival$`covariable_object=2`)){
          survival_stats_table[i, 12]<- list_of_mean_survival$`covariable_object=2`
        }
        if(!is.null(list_of_mean_survival$`covariable_object=3`)){
          survival_stats_table[i, 13]<- list_of_mean_survival$`covariable_object=3`
        }
        if(!is.null(list_of_mean_survival$`covariable_object=4`)){
          survival_stats_table[i, 14]<- list_of_mean_survival$`covariable_object=4`
        }
      }
    }
    #   ## If less than one category in survival analysis need to adjust values saved:
    #   if(is.null(ncol(fittedSurv_mean$matrix))){
    #     #survival_stats_table[i, 7]<- unique(covariable_object) #categories
    #     survival_stats_table[i, 8]<- paste(fittedSurv_mean$matrix[5], collapse = " ") #mean
    #   } else {
    #     #survival_stats_table[i, 7]<- paste(rownames(fittedSurv_mean$matrix), collapse = " ") #categories
    #     #survival_stats_table[i, 8]<- paste(fittedSurv_mean$matrix[,5], collapse = " ") #mean
    #     if (survival_stats_table[i, 7] != 0){
    #       survival_stats_table[i, 11]<- 
    #   }
    #   
    #   
    
  }
  
  return(survival_stats_table)
  
}






#########################
### Repeat loop but perform disease free survival analysis:

##Create an empty list
survival_stats_cancer_list4<- vector("list", length(threshold_selected_cnv_list_plus_all_loc))

for (i in 1:length(threshold_selected_cnv_list_plus_all_loc)){
  
  ##Get CNV table
  cnv.table<- threshold_selected_cnv_list_plus_all_loc[[i]]
  
  
  ##Get survival data
  survival_time_list<- clinical_long_plus_all_survival[[i]]
  
  
  ## Perform survival analysis of co-deletions:
  ##Disease free survival = column 7!!!
  print(paste(names(threshold_selected_cnv_list_plus_all_loc[i]), "Cancer type:"))
  survival_stats_list<- lapply(gene_information_long_list, function(x) survival_analysis_of_gene_list_cat1_and_2(target_gene_list = x, survival_time_list = survival_time_list, cnv.table = cnv.table, distance = 2.5e+06,
                                                                                                                 threshold = -2, deletion = TRUE, time_of_death_column = 7, 
                                                                                                                 death_event_column = 6, column_start = 11, start = TRUE, 
                                                                                                                 remove_NA = TRUE, Cytoband = FALSE, print_to_screen = FALSE, 
                                                                                                                 plot_graph = FALSE))
  
  ##Combine survival stats together and order by log-rank test p-value
  survival_stats_list_table<- do.call(rbind, survival_stats_list)
  survival_stats_list_table<- survival_stats_list_table[order(survival_stats_list_table$`p-value_logrank_test`, decreasing = FALSE),]
  #survival_stats_list_table[1:100, c(1,2,5,6,8,9)]
  
  ##Save file as .csv
  write.csv(survival_stats_list_table, file = paste0(names(threshold_selected_cnv_list_plus_all_loc[i]),"_co-deletion_disease_free_survival_stats_cat1_and_2.csv"), quote = FALSE)
  
  ##Add to list
  survival_stats_cancer_list4[[i]]<- survival_stats_list_table
  
  
  
}

## Save object:
saveRDS(survival_stats_cancer_list4, file = "./R workspaces/survival_stats_disease_free_cancer_list_cat1_and_2")


###########
###Delete cancers from list that dont have correct disease free survival data:
# UCS, PCPG, MESO, ACC










