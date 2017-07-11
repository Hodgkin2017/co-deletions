#################
###Survival statistics function long
#################

##########
###Parameters:
# gene_information_list
# target_gene_list<- gene_information_list[1:2]
# target_gene_list
# cnv.table<- threshold_short_cnv_list_loc[[1]]
# dim(cnv.table)
# gene_information<- gene_information_list[[1]]
# gene_information
# 
# remove_NA = TRUE
# start = TRUE
# Chromosome<- gene_information[[2]]
# distance<- 2.5e+06
# selection_criteria<- c(gene_information[[4]]-distance, gene_information[[5]]+distance)
# column_start = 11
# deletion = TRUE
# threshold = -2
# target_gene = gene_information[[1]]
# survival_time_list<- clinical_survival_list[[1]]
# time_of_death_column<- 5
# death_event_column<- 6
# print_to_screen<- FALSE
# 
# #########
# ##Functions used in main function:
# 
# ##########
# ## Function to categorise the type of deletion(1, 2, 3 of 4):
# 
# # categorise_deletion_type_function<- function(gene_names, gene_names_cnv.matrix){
# #   
# #   ##Determine if target gene is deleted
# #   target_gene_deleted<- gene_names_cnv.matrix %>%
# #     dplyr::filter(gene == target_gene) %>%
# #     dplyr::select(-gene)
# #   
# #   target_gene_deleted<- as.vector(target_gene_deleted == 1)
# #   
# #   ##Determine if proximal gene is deleted
# #   proximal_gene_deleted<- gene_names_cnv.matrix %>%
# #     dplyr::filter(gene == gene_names) %>%
# #     dplyr::select(-gene)
# #   
# #   proximal_gene_deleted<- as.vector(proximal_gene_deleted == 1)
# #   
# #   ##Determine if target gene and not proximal gene is deleted(category 1)
# #   deletion_category1<- proximal_gene_deleted - target_gene_deleted < 0
# #   
# #   ##Determine if target gene and proximal gene are deleted together (category 2)
# #   deletion_category2<- target_gene_deleted + proximal_gene_deleted > 1
# #   
# #   ##Determine if proximal gene is deleted without target gene (category 3)
# #   deletion_category3<- target_gene_deleted - proximal_gene_deleted < 0
# #   
# #   ##Determine if neither proximal gene or target gene is deleted (category 4)
# #   deletion_category4<- target_gene_deleted + proximal_gene_deleted == 0
# #   
# #   ##Join categories together into one table:
# #   all_deletion_category_table<- cbind(deletion_category1, deletion_category2, deletion_category3, deletion_category4)
# #   
# #   ##Find one category per tumour:
# #   final_deletion_category<- max.col(all_deletion_category_table,ties.method="first")
# #   
# #   return(final_deletion_category)
# #   
# # }
# 
# 
# 
# 
# 
# ##########
# ###Function to take a target gene and surrounding genes and characterise each tumour based on if it had a deletion 
# #of the target gene and each suurounbding gene 
# 
# categorise_deletion_type_around_target_gene<- function(cnv.table, target_gene, Chromosome, selection_criteria, threshold = -2, 
#                                                        deletion = TRUE, column_start = 11, start = TRUE, remove_NA = TRUE, Cytoband = FALSE) {
#   
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
#   ###Create a new table with genes in rows and tumours on columns and each entry has a number from 1 to 4 depending on which 
#   #deletion/amplification catagory it belongs to:
#   
#   ##Add gene names as new column
#   gene_names_cnv.matrix<- cbind.data.frame(gene = row.names(cnv.matrix) ,cnv.matrix)
#   gene_names_cnv.matrix$gene<-as.character(gene_names_cnv.matrix$gene)
#   gene_names<- row.names(cnv.matrix)
#   
#   ##########
#   ###Function to categorise gene deletions:
#   ## If I remove it from function the function no longer works!
#   categorise_deletion_type_function<- function(gene_names, gene_names_cnv.matrix){
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
#   deletion_category_table<- lapply(gene_names, function(x) categorise_deletion_type_function(x, gene_names_cnv.matrix))
#   length(deletion_category_table)
#   #deletion_category_table[[22]]
#   
#   deletion_category_table<- do.call(rbind, deletion_category_table)
#   colnames(deletion_category_table)<- colnames(cnv.matrix)
#   rownames(deletion_category_table)<- rownames(cnv.matrix)
#   #deletion_category_table[20:24, 169:173]
#   #cnv.matrix[20:24, 169:173]
#   
#   return(deletion_category_table)
#   
# }
# 
# 
# ################
# ### Function to append deletion category to clinical data and remove any tumours without a deletion category:
# 
# deletion_category_target_gene<-deletion_category_target_gene_list[[2]]
# 
# join_clinical_deletion_category_table<- function(deletion_category_target_gene, survival_time_list, time_of_death_column, death_event_column){
#   # deletion_category_gene_name_table<- cbind.data.frame(gene = row.names(x) ,x)
#   # deletion_category_gene_name_table$gene<-as.character(deletion_category_gene_name_table$gene)
#   # deletion_category_gene_name_table[1:2, 1:5]
#   # 
#   # ##Filter deletion category table by current gene of interest i.e. proximal gene
#   # deletion_category<- deletion_category_gene_name_table %>% 
#   #   dplyr::filter(gene == proximal_gene) %>%
#   #   dplyr::select(-gene)
#   # 
#   # dim(deletion_category)
#   # deletion_category[1, 1:5]
#   deletion_category<-t(deletion_category_target_gene)
#   dim(deletion_category)
#   
#   colnames(deletion_category)
#   rownames(deletion_category)
#   
#   ##Convert patient IDs in deletion_category so they match with the patient IDS in the clinical table
#   deletion_category_patient_ID<- rownames(deletion_category) %>%
#     substr(0, 12) %>%
#     gsub("\\.", "-", .) %>%
#     cbind.data.frame(patient_IDs = ., deletion_category)
#   
#   deletion_category_patient_ID$patient_IDs<- as.character(deletion_category_patient_ID$patient_IDs)
#   
#   class(deletion_category_patient_ID$patient_IDs)
#   deletion_category_patient_ID
#   clinical_survival_list[[1]]$patient_IDs
#   
#   short_survival_time_list<- survival_time_list[,c(1, time_of_death_column, death_event_column)]
#   
#   clinical_survival_deletion_category<- dplyr::full_join(short_survival_time_list, deletion_category_patient_ID, by = "patient_IDs") %>%
#     dplyr::select(-patient_IDs)
#   head(clinical_survival_deletion_category, 40)
#   dim(clinical_survival_deletion_category)
#   dim(clinical_survival_list[[1]])
#   colnames(clinical_survival_deletion_category)
#   rownames(clinical_survival_deletion_category)
#   
#   
#   #######
#   ### Remove entries with NA deletion categories due to missing clinical data
#   vars<- colnames(deletion_category_patient_ID[-1])
#   deletion_category_proximal_genes<- clinical_survival_deletion_category %>%
#     dplyr::select(one_of(vars))
#   dim(clinical_survival_deletion_category)
#   
#   clinical_survival_deletion_category<- clinical_survival_deletion_category[complete.cases(deletion_category_proximal_genes), ]
#   dim(clinical_survival_deletion_category)
#   
#   # apply(test, 2, function(x) !is.na(x))
#   # row.has.na <- apply(test, 1, function(x){any(is.na(x))})
#   # row.has.na
#   # test[c(row.has.na), ]
#   # test[complete.cases(test), ]
#   
#   return(clinical_survival_deletion_category)
# }
# 
# 
# 
# 
# 
# 
# 
# ################
# ### Another function:
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #########
# ### Main Function:
# 
# ##Parameters:
# # gene_information_list
# # target_gene_list<- gene_information_list[1:2]
# # target_gene_list
# # cnv.table<- threshold_short_cnv_list_loc[[1]]
# # survival_time_list<- clinical_survival_list[[1]]
# 
# 
# survival_analysis_of_gene_list<- function(target_gene_list, survival_time_list, cnv.table, distance = 2.5e+06, 
#                                           threshold = -2, deletion = TRUE, time_of_death_column, 
#                                           death_event_column, column_start = 11, start = TRUE, 
#                                           remove_NA = TRUE, Cytoband = FALSE, print_to_screen = FALSE, 
#                                           plot_graph = FALSE){
# 
#   
# ##########
# ###Create a list of tables were each table contains deletion categories(1-4) for each tumour(columns) per 
# #gene(rows) surrounding an interval around a target gene.
# 
# deletion_category_target_gene_list <- lapply(target_gene_list, function(x) categorise_deletion_type_around_target_gene(cnv.table = cnv.table, target_gene = x[[1]], Chromosome = x[[2]], 
#                                                                                 selection_criteria = c(x[[4]] - distance, x[[5]]+ distance), threshold = threshold, deletion = deletion, column_start = column_start, start = start, remove_NA = remove_NA, Cytoband = Cytoband))
# 
# # deletion_category_target_gene_list[[1]][,1:3]
# # deletion_category_target_gene_list[[2]][,1:3]
# 
# 
# 
# 
# 
# 
# ########
# ### Take a deletion category table for a target gene and join to survival table
# ## Add extra column with patient ID:
# # deletion_category_target_gene<-deletion_category_target_gene_list[[2]]
# # survival_time_list<- clinical_survival_list[[1]]
# # 
# # join_clinical_deletion_category_table(deletion_category_target_gene = deletion_category_target_gene, survival_time_list = survival_time_list)
# 
# deletion_category_survival_target_gene_list<- lapply(deletion_category_target_gene_list, function(x) join_clinical_deletion_category_table(x, survival_time_list, time_of_death_column = time_of_death_column, death_event_column = death_event_column))
# # deletion_category_survival_target_gene_list[[1]][1:3, 1:12]
# 
# #########################From here!
# return(deletion_category_survival_target_gene_list)
# 
# }
# 
# test_function<- survival_analysis_of_gene_list(target_gene_list,
#                                                survival_time_list,
#                                                cnv.table,
#                                                distance = 2.5e+06,
#                                                threshold = -2,
#                                                deletion = TRUE,
#                                                column_start = 11,
#                                                time_of_death_column = 5,
#                                                death_event_column = 6,
#                                                start = TRUE,
#                                                remove_NA = TRUE,
#                                                Cytoband = FALSE)
# length(test_function)
# test_function[[1]]
# dim(test_function[[1]])
# dim(test_function[[2]])
# #########################Above here!
# 
# #######
# ### Get Survival statistics
# deletion_category_survival_target_gene<- test_function[[1]]
# 
# survival_statistics<- function(deletion_category_survival_target_gene, print_to_screen = print_to_screen, plot_graph = plot_graph){
#   
# death_time<- deletion_category_survival_target_gene[,1]
# death_event<- deletion_category_survival_target_gene[,2]
# 
# ##Create Surv object
# surv_object<- Surv(death_time, death_event==1)
# 
# ##Create empty table to store stats:
# survival_stats_table<- data.frame(matrix(NA, ncol = 9, nrow = ncol(deletion_category_survival_target_gene) -2))
# names(survival_stats_table)<-c("target_gene", "proximal_gene", "p-value_Likelihood_ratio_test",
#                          "p-value_Wald_test", "p-value_logrank_test", "Hazard_ratio", "Categories",
#                          "mean_survival", "number_of_samples_per_category")
# 
# for (i in 1:nrow(survival_stats_table)){
#   
#   ##Get covariable object for survfit:
#   covariable_object<- deletion_category_survival_target_gene[,i+2]
#   
#   ##Fit Kaplain meier graph to one co-variable to compare data with :
#     fittedSurv <- survfit(surv_object~covariable_object, na.action = na.exclude)
#     fittedSurv_mean<- survival:::survmean(fittedSurv, rmean="individual") 
#     #survival:::survmean(fittedSurv, rmean="common") 
#     df.categ <- cbind(sapply(names(fittedSurv$strata), function(x) strsplit(x,"=")[[1]][2]),
#                       fittedSurv$n)
#     df.categ <- data.frame(df.categ)
# 
#   ##Get category names for one variable
#   categNames <- apply(df.categ, 1, function(x) paste0(x[1]," (",x[2],")"))
#   
#   # if (plot_graph == TRUE) {
#   #   plot(fittedSurv, main=plotTitle,
#   #        xlab="Time (days)", ylab=ylabel, 
#   #        col=brewer.pal(9,"Set1"), mark.time=T)
#   #   legend("topright", legend=categNames, 
#   #          col=brewer.pal(9,"Set1"), 
#   #          lwd=2, cex=0.9)
#   # }
#   if(print_to_screen == TRUE) {
#   print("Chi-sq test:")
#     print(survdiff(surv_object~covariable_object,rho = 0))
#     
#     print("Cox PH test:")
#     print(summary(coxph(surv_object~covariable_object)))
# }
#     coxfit <- coxph(surv_object~covariable_object)
#   
#   # if (plot_graph == TRUE) {
#   #   text(1000,0,labels=paste0("HR=",round(exp(summary(coxfit)$coefficients[1]),2),"; p=",
#   #                             round(summary(coxfit)$logtest[3],3)))
#   # }
#   
#     survival_stats_table[i, 1]<- target_gene
#     survival_stats_table[i, 2]<- colnames(deletion_category_survival_target_gene[i+2])
#     survival_stats_table[i, 3]<- round(summary(coxfit)$logtest[3],2)
#     survival_stats_table[i, 4]<- round(summary(coxfit)$waldtest[3],2)
#     survival_stats_table[i, 5]<- round(summary(coxfit)$sctest[3],2)
#     survival_stats_table[i, 6]<- round(summary(coxfit)$coefficients[2],2)
#   
#     if(is.null(ncol(fittedSurv_mean$matrix))){
#       survival_stats_table[i, 7]<- unique(covariable_object) #categories
#       survival_stats_table[i, 8]<- paste(fittedSurv_mean$matrix[5], collapse = " ") #mean
#       survival_stats_table[i, 9]<-  paste(fittedSurv_mean$matrix[3], collapse = " ") #number of samples
#     } else {
#     survival_stats_table[i, 7]<- paste(rownames(fittedSurv_mean$matrix), collapse = " ") #categories
#     survival_stats_table[i, 8]<- paste(fittedSurv_mean$matrix[,5], collapse = " ") #mean
#     survival_stats_table[i, 9]<-  paste(fittedSurv_mean$matrix[,3], collapse = " ") #number of samples
#     }
#     
#  
# }
#   
#   return(survival_stats_table)
#   
# }
# 
# ########
# ### test function:
# deletion_category_survival_target_gene<- test_function[[2]]
# target_gene<- gene_information_list[[2]][[1]]
# target_gene
# 
# test_stats<- survival_statistics(deletion_category_survival_target_gene, print_to_screen = print_to_screen, plot_graph = plot_graph)
# test_stats
# 
# ## Test function in apply:
# #test_stats<-lapply(test_function, function(x) survival_statistics(x, print_to_screen = print_to_screen, plot_graph = plot_graph))

####################################################################################
########################
### Survival statistics function long
#########################


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

#target_gene, Chromosome, selection_criteria,

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
  ###Create a new table with genes in rows and tumours on columns and each entry has a number from 1 to 4 depending on which 
  #deletion/amplification catagory it belongs to:
  
  ##Add gene names as new column
  gene_names_cnv.matrix<- cbind.data.frame(gene = row.names(cnv.matrix) ,cnv.matrix)
  gene_names_cnv.matrix$gene<-as.character(gene_names_cnv.matrix$gene)
  gene_names<- row.names(cnv.matrix)
  
  ##########
  ###Function to categorise gene deletions:
  ## If I remove it from function the function no longer works!
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
  
  deletion_category_table<- lapply(gene_names, function(x) categorise_deletion_type_function(x, gene_names_cnv.matrix))
  length(deletion_category_table)
  #deletion_category_table[[22]]
  
  deletion_category_table<- do.call(rbind, deletion_category_table)
  colnames(deletion_category_table)<- colnames(cnv.matrix)
  rownames(deletion_category_table)<- rownames(cnv.matrix)
  #deletion_category_table[20:24, 169:173]
  #cnv.matrix[20:24, 169:173]
  
  #return(deletion_category_table)
  
  
  
  
  
  
  #########################################
  
  deletion_category<-t(deletion_category_table)
  dim(deletion_category)
  
  colnames(deletion_category)
  rownames(deletion_category)
  
  ##Convert patient IDs in deletion_category so they match with the patient IDS in the clinical table
  deletion_category_patient_ID<- rownames(deletion_category) %>%
    substr(0, 12) %>%
    gsub("\\.", "-", .) %>%
    cbind.data.frame(patient_IDs = ., deletion_category)
  
  deletion_category_patient_ID$patient_IDs<- as.character(deletion_category_patient_ID$patient_IDs)
  
  # class(deletion_category_patient_ID$patient_IDs)
  # deletion_category_patient_ID
  # clinical_survival_list[[1]]$patient_IDs
  
  short_survival_time_list<- survival_time_list[,c(1, time_of_death_column, death_event_column)]
  
  clinical_survival_deletion_category<- dplyr::full_join(short_survival_time_list, deletion_category_patient_ID, by = "patient_IDs") %>%
    dplyr::select(-patient_IDs)
  head(clinical_survival_deletion_category, 40)
  dim(clinical_survival_deletion_category)
  dim(clinical_survival_list[[1]])
  colnames(clinical_survival_deletion_category)
  rownames(clinical_survival_deletion_category)
  
  
  #######
  ### Remove entries with NA deletion categories due to missing clinical data
  vars<- colnames(deletion_category_patient_ID[-1])
  deletion_category_proximal_genes<- clinical_survival_deletion_category %>%
    dplyr::select(one_of(vars))
  dim(clinical_survival_deletion_category)
  
  clinical_survival_deletion_category<- clinical_survival_deletion_category[complete.cases(deletion_category_proximal_genes), ]
  dim(clinical_survival_deletion_category)
  
  # apply(test, 2, function(x) !is.na(x))
  # row.has.na <- apply(test, 1, function(x){any(is.na(x))})
  # row.has.na
  # test[c(row.has.na), ]
  # test[complete.cases(test), ]
  
  #return(clinical_survival_deletion_category)
  
  
  
  
  
  ################################################
  
  #survival_statistics<- function(deletion_category_survival_target_gene, print_to_screen = print_to_screen, plot_graph = plot_graph)
  
  death_time<- clinical_survival_deletion_category[,1]
  death_event<- clinical_survival_deletion_category[,2]
  
  ##Create Surv object
  surv_object<- Surv(death_time, death_event==1)
  
  ##Create empty table to store stats:
  survival_stats_table<- data.frame(matrix(NA, ncol = 9, nrow = ncol(clinical_survival_deletion_category) -2))
  names(survival_stats_table)<-c("target_gene", "proximal_gene", "p-value_Likelihood_ratio_test",
                                 "p-value_Wald_test", "p-value_logrank_test", "Hazard_ratio", "Categories",
                                 "mean_survival", "number_of_samples_per_category")
  
  for (i in 1:nrow(survival_stats_table)){
    
    ##Get covariable object for survfit:
    covariable_object<- clinical_survival_deletion_category[,i+2]
    
    ##Fit Kaplain meier graph to one co-variable to compare data with :
    fittedSurv <- survfit(surv_object~covariable_object, na.action = na.exclude)
    fittedSurv_mean<- survival:::survmean(fittedSurv, rmean="individual") 
    #survival:::survmean(fittedSurv, rmean="common") 
    df.categ <- cbind(sapply(names(fittedSurv$strata), function(x) strsplit(x,"=")[[1]][2]),
                      fittedSurv$n)
    df.categ <- data.frame(df.categ)
    
    ##Get category names for one variable
    categNames <- apply(df.categ, 1, function(x) paste0(x[1]," (",x[2],")"))
    
    # if (plot_graph == TRUE) {
    #   plot(fittedSurv, main=plotTitle,
    #        xlab="Time (days)", ylab=ylabel, 
    #        col=brewer.pal(9,"Set1"), mark.time=T)
    #   legend("topright", legend=categNames, 
    #          col=brewer.pal(9,"Set1"), 
    #          lwd=2, cex=0.9)
    # }
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
    # }
    
    survival_stats_table[i, 1]<- target_gene
    survival_stats_table[i, 2]<- colnames(clinical_survival_deletion_category[i+2])
    survival_stats_table[i, 3]<- round(summary(coxfit)$logtest[3],2)
    survival_stats_table[i, 4]<- round(summary(coxfit)$waldtest[3],2)
    survival_stats_table[i, 5]<- round(summary(coxfit)$sctest[3],2)
    survival_stats_table[i, 6]<- round(summary(coxfit)$coefficients[2],2)
    
    if(is.null(ncol(fittedSurv_mean$matrix))){
      survival_stats_table[i, 7]<- unique(covariable_object) #categories
      survival_stats_table[i, 8]<- paste(fittedSurv_mean$matrix[5], collapse = " ") #mean
      survival_stats_table[i, 9]<-  paste(fittedSurv_mean$matrix[3], collapse = " ") #number of samples
    } else {
      survival_stats_table[i, 7]<- paste(rownames(fittedSurv_mean$matrix), collapse = " ") #categories
      survival_stats_table[i, 8]<- paste(fittedSurv_mean$matrix[,5], collapse = " ") #mean
      survival_stats_table[i, 9]<-  paste(fittedSurv_mean$matrix[,3], collapse = " ") #number of samples
    }
    
    
  }
  
  return(survival_stats_table)
  
}







###########################################







#   ##########
#   ###Create a list of tables were each table contains deletion categories(1-4) for each tumour(columns) per 
#   #gene(rows) surrounding an interval around a target gene.
#   
#   deletion_category_target_gene_list <- lapply(target_gene_list, function(x) categorise_deletion_type_around_target_gene(cnv.table = cnv.table, target_gene = x[[1]], Chromosome = x[[2]], 
#                                                                                                                          selection_criteria = c(x[[4]] - distance, x[[5]]+ distance), threshold = threshold, deletion = deletion, column_start = column_start, start = start, remove_NA = remove_NA, Cytoband = Cytoband))
#   
#   # deletion_category_target_gene_list[[1]][,1:3]
#   # deletion_category_target_gene_list[[2]][,1:3]
#   
#   
#   
#   
#   
#   
#   ########
#   ### Take a deletion category table for a target gene and join to survival table
#   ## Add extra column with patient ID:
#   # deletion_category_target_gene<-deletion_category_target_gene_list[[2]]
#   # survival_time_list<- clinical_survival_list[[1]]
#   # 
#   # join_clinical_deletion_category_table(deletion_category_target_gene = deletion_category_target_gene, survival_time_list = survival_time_list)
#   
#   deletion_category_survival_target_gene_list<- lapply(deletion_category_target_gene_list, function(x) join_clinical_deletion_category_table(x, survival_time_list, time_of_death_column = time_of_death_column, death_event_column = death_event_column))
#   # deletion_category_survival_target_gene_list[[1]][1:3, 1:12]
#   
#   #########################From here!
#   return(deletion_category_survival_target_gene_list)
#   
# }
# 
# test_function<- survival_analysis_of_gene_list(target_gene_list,
#                                                survival_time_list,
#                                                cnv.table,
#                                                distance = 2.5e+06,
#                                                threshold = -2,
#                                                deletion = TRUE,
#                                                column_start = 11,
#                                                time_of_death_column = 5,
#                                                death_event_column = 6,
#                                                start = TRUE,
#                                                remove_NA = TRUE,
#                                                Cytoband = FALSE)
# length(test_function)
# test_function[[1]]
# dim(test_function[[1]])
# dim(test_function[[2]])



##############
test<- survival_analysis_of_gene_list(target_gene_list = target_gene_list, survival_time_list = survival_time_list, cnv.table = cnv.table, distance = 2.5e+06,
                                      threshold = -2, deletion = TRUE, time_of_death_column = 5, 
                                      death_event_column = 6, column_start = 11, start = TRUE, 
                                      remove_NA = TRUE, Cytoband = FALSE, print_to_screen = FALSE, 
                                      plot_graph = FALSE)

test_apply<- lapply(gene_information_list, function(x) survival_analysis_of_gene_list(target_gene_list = x, survival_time_list = survival_time_list, cnv.table = cnv.table, distance = 2.5e+06,
                                                                                      threshold = -2, deletion = TRUE, time_of_death_column = 5, 
                                                                                      death_event_column = 6, column_start = 11, start = TRUE, 
                                                                                      remove_NA = TRUE, Cytoband = FALSE, print_to_screen = FALSE, 
                                                                                      plot_graph = FALSE))

identical(test, test_apply[[2]])
all.equal(test, test_apply[[2]])
length(test_apply)
test[1:5, 1:5]
test_apply[[2]][1:5, 1:5]
dim(test)
dim( test_apply[[2]])
dim( test_apply[[23]])

test_apply_table<- do.call(rbind, test_apply)
test_apply_table
dim(test_apply_table)
sapply(test_apply, function(x) nrow(x))
sum(sapply(test_apply, function(x) nrow(x)))

test_apply_table_sort1<- test_apply_table[order(test_apply_table$`p-value_logrank_test`, decreasing = FALSE),]
test_apply_table_sort1[1:100, c(1,2,5,6,8,9)]

#################
###Next extract n value information for each catagory: 1,2,3,4



### Plot survival graphs for significant genes with sufficient n numbers.













#####
### Plot survival graph if wanted:








#########
### Create a table for each gene v's target gene containing survival statistics e.g. p-values







