#####################
### Survival analysis stats table function...temp
#####################

#################
##Parameters:

target_gene_list<- gene_information_list[[1]]
target_gene_list
survival_time_list<- clinical_survival_list[[1]]
cnv.table<- threshold_short_cnv_list_loc[[1]]

distance = 2.5e+06
threshold = -2
deletion = TRUE
time_of_death_column = 5 
death_event_column = 6
column_start = 11
start = TRUE 
remove_NA = TRUE
Cytoband = FALSE
print_to_screen = FALSE 
plot_graph = FALSE

#########
### Function to Create a table of survival p-value statistics for co-deletions or co-amplifications

survival_analysis_of_gene_list<- function(target_gene_list,
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
  ### Survival analysis:
  
  ##Create Surv object and remove NA values
  death_time<- clinical_survival_deletion_category[,1]
  death_event<- clinical_survival_deletion_category[,2]
  surv_object<- Surv(death_time, death_event==1)
  ##Remove entries with NA:
  NA_object<- !is.na(surv_object)
  death_time<- death_time[NA_object]
  death_event<- death_event[NA_object]
  surv_object<- surv_object[NA_object]
  
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
    covariable_object<- covariable_object[NA_object]
    
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

###########################################
### Test function:
target_gene_list<- gene_information_list[[1]]
survival_time_list<- clinical_survival_list[[1]]
cnv.table<- threshold_short_cnv_list_loc[[1]]


##O1:
test<- survival_analysis_of_gene_list(target_gene_list = target_gene_list, survival_time_list = survival_time_list, cnv.table = cnv.table, distance = 2.5e+06,
                                      threshold = -2, deletion = TRUE, time_of_death_column = 5, 
                                      death_event_column = 6, column_start = 11, start = TRUE, 
                                      remove_NA = TRUE, Cytoband = FALSE, print_to_screen = FALSE, 
                                      plot_graph = FALSE)

##O2:
test_apply<- lapply(gene_information_list, function(x) survival_analysis_of_gene_list(target_gene_list = x, survival_time_list = survival_time_list, cnv.table = cnv.table, distance = 2.5e+06,
                                                                                      threshold = -2, deletion = TRUE, time_of_death_column = 5, 
                                                                                      death_event_column = 6, column_start = 11, start = TRUE, 
                                                                                      remove_NA = TRUE, Cytoband = FALSE, print_to_screen = FALSE, 
                                                                                      plot_graph = FALSE))








