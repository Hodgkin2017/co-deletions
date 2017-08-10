################
### Survival stratification
################

#######################
#######################
### Stratify survival by cancer stage


###################
### Get pathological grade from all clinical data files

##which cancers contain "pathologic_stage" information

length(clinical_fbget_long_list)
names(clinical_fbget_long_list)
test<- clinical_fbget_long_list[[3]]
grep("pathologic_stage", colnames(test))
any(colnames(test) == "pathologic_stage")

test<- clinical_fbget_long_list[[37]]
colnames(test) == "pathologic_stage"
#any(colnames(test) == "pathologic_stage")
sum(colnames(test) == "pathologic_stage")

sapply(clinical_fbget_long_list, function(x) sum(colnames(x) == "pathologic_stage"))
indicies<- sapply(clinical_fbget_long_list, function(x) any(colnames(x) == "pathologic_stage"))

sum(sapply(clinical_fbget_long_list, function(x) any(colnames(x) == "pathologic_stage")))

##Create list of cancers with pathological stage informations
clinical_fbget_long_list_short<- clinical_fbget_long_list[indicies]
length(clinical_fbget_long_list_short)
clinical_fbget_long_list_short[[1]]

## Create function to extract tcga_participant_barcode and pathologic_stage
#clinical_table<- clinical_fbget_long_list[[3]]

extract_information<- function(clinical_table){
  
  extracted_results<- clinical_table %>%
    dplyr::select(tcga_participant_barcode, pathologic_stage)
  
  return(extracted_results)
}

##Apply function to shorter clinical list:
clinical_fbget_long_list_short_pathologic_stage<- lapply(clinical_fbget_long_list_short, function (x) extract_information(x))
length(clinical_fbget_long_list_short_pathologic_stage)
clinical_fbget_long_list_short_pathologic_stage[[2]]

##Join tables for different cancers together:
clinical_fbget_long_list_short_pathologic_stage<- do.call(rbind, clinical_fbget_long_list_short_pathologic_stage)
class(clinical_fbget_long_list_short_pathologic_stage)
dim(clinical_fbget_long_list_short_pathologic_stage)

##Get unique patients only
sum(duplicated(clinical_fbget_long_list_short_pathologic_stage))
clinical_fbget_long_list_short_pathologic_stage_unique<- clinical_fbget_long_list_short_pathologic_stage[!duplicated(clinical_fbget_long_list_short_pathologic_stage), ]
length(unique(clinical_fbget_long_list_short_pathologic_stage_unique$tcga_participant_barcode))
nrow(clinical_fbget_long_list_short_pathologic_stage_unique)

## What tumour stages exist in dataset?
unique(clinical_fbget_long_list_short_pathologic_stage_unique$pathologic_stage)
table(clinical_fbget_long_list_short_pathologic_stage_unique$pathologic_stage)

##Create function that converts pathological stage to early or late:
pathologic_stage<- clinical_fbget_long_list_short_pathologic_stage_unique$pathologic_stage 
pathologic_stage

pathologic_stage_conversion<- function(pathologic_stage){
  new_pathologic_stage<- gsub(" ", "", pathologic_stage)
  new_pathologic_stage<- gsub("i/iinos", "early", new_pathologic_stage)
  new_pathologic_stage<- gsub("is", "early", new_pathologic_stage)
  new_pathologic_stage<- gsub("stage0", "early", new_pathologic_stage)
  new_pathologic_stage<- gsub("\\<stagei\\>", "early", new_pathologic_stage)
  new_pathologic_stage<- gsub("\\<stageia\\>", "early", new_pathologic_stage)
  new_pathologic_stage<- gsub("\\<stageib\\>", "early", new_pathologic_stage)
  new_pathologic_stage<- gsub("\\<stageii\\>", "early", new_pathologic_stage)
  new_pathologic_stage<- gsub("\\<stageiia\\>", "early", new_pathologic_stage)
  new_pathologic_stage<- gsub("\\<stageiib\\>", "early", new_pathologic_stage)
  new_pathologic_stage<- gsub("\\<stageiic\\>", "early", new_pathologic_stage)
  new_pathologic_stage<- gsub("\\<stageiii\\>", "late", new_pathologic_stage)
  new_pathologic_stage<- gsub("\\<stageiiia\\>", "late", new_pathologic_stage)
  new_pathologic_stage<- gsub("\\<stageiiib\\>", "late", new_pathologic_stage)
  new_pathologic_stage<- gsub("\\<stageiiic\\>", "late", new_pathologic_stage)
  new_pathologic_stage<- gsub("\\<stageiv\\>", "late", new_pathologic_stage)
  new_pathologic_stage<- gsub("\\<stageiva\\>", "late", new_pathologic_stage)
  new_pathologic_stage<- gsub("\\<stageivb\\>", "late", new_pathologic_stage)
  new_pathologic_stage<- gsub("\\<stageivc\\>", "late", new_pathologic_stage)
  new_pathologic_stage<- gsub("\\<stagex\\>", "early", new_pathologic_stage)
  
  
}
table(pathologic_stage)
table(new_pathologic_stage)

clinical_pathologic_stage<- pathologic_stage_conversion(clinical_fbget_long_list_short_pathologic_stage_unique$pathologic_stage )
clinical_pathologic_stage
table(clinical_pathologic_stage)
clinical_fbget_long_list_short_pathologic_stage_unique$binary_pathologic_stage<- clinical_pathologic_stage
dim(clinical_fbget_long_list_short_pathologic_stage_unique)

##Keep only complete cases
clinical_fbget_long_list_short_pathologic_stage_unique<- clinical_fbget_long_list_short_pathologic_stage_unique[complete.cases(clinical_fbget_long_list_short_pathologic_stage_unique),]
dim(clinical_fbget_long_list_short_pathologic_stage_unique)
colnames(clinical_fbget_long_list_short_pathologic_stage_unique)<- c("patient_IDs", "pathologic_stage", "binary_pathologic_stage")
head(clinical_fbget_long_list_short_pathologic_stage_unique)

## Join pathologic_stage table with the clinical table for all tumours
names(clinical_long_plus_all_survival)
length(clinical_long_plus_all_survival)
dim(clinical_long_plus_all_survival$ALL)
head(clinical_long_plus_all_survival$ALL)
colnames(clinical_long_plus_all_survival$ALL)

clinical_all_survival_pathologic<- dplyr::full_join(clinical_long_plus_all_survival$ALL, clinical_fbget_long_list_short_pathologic_stage_unique, by = "patient_IDs")
dim(clinical_all_survival_pathologic)
head(clinical_all_survival_pathologic)

##Remove patients without pathological stage:
clinical_all_survival_pathologic_complete<- clinical_all_survival_pathologic[complete.cases(clinical_all_survival_pathologic[,9]),]
dim(clinical_all_survival_pathologic_complete)
head(clinical_all_survival_pathologic_complete)






###################
### Use pathologic_stage for survival stratification.


cnv.table<- threshold_selected_cnv_list_plus_all_loc$ALL
##Get survival data
survival_time_list<- clinical_all_survival_pathologic_complete
## Get Target gene
target_gene_names<- sapply(gene_information_long_list, function(x) x[[1]])
index<- which(target_gene_names %in%  "RB1" )
target_gene_list<- gene_information_long_list[[index]]

distance = 2.5e+06
threshold = -2
deletion = TRUE
time_of_death_column = 5 
death_event_column = 6
stratification_column = 9
column_start = 11
start = TRUE 
remove_NA = TRUE
Cytoband = FALSE
print_to_screen = TRUE 
plot_graph = TRUE
ylabel = "Overall survival"
path = "./"
#path = ("../../Output/Survival analysis/survival_curves")





#########
###F2: Function to Create a table of survival p-value statistics for co-deletions or co-amplifications

stratified_survival_cat1V2<- function(target_gene_list,
                                                     survival_time_list,
                                                     cnv.table,
                                                     distance = 2.5e+06,
                                                     threshold = -2,
                                                     deletion = TRUE,
                                                     time_of_death_column, 
                                                     death_event_column,
                                                     column_start = 11,
                                                    stratification_column,
                                                     start = TRUE, 
                                                     remove_NA = TRUE,
                                                     Cytoband = FALSE,
                                                     print_to_screen = FALSE, 
                                                     plot_graph = FALSE,
                                                     ylabel = "Overall survival",
                                                     path = "./",
                                                    output_deletion_category = FALSE){
  
  ##########
  ##Get and change working directory
  current_wd<- getwd()
  setwd(path)
  
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
  
  #############################
  ### Instead of performing survival analysis the function outputs the deletion 
  #catagories for each tumour and gene
  
  if(output_deletion_category == TRUE){
    
    return(deletion_category_patient_ID)
    
  }
  
  ##############
  ##Extract patient ID, time of death and if death occured columns from survival table:
  short_survival_time_list<- survival_time_list[,c(1, time_of_death_column, death_event_column, stratification_column)]
  
  ##Join clinical table with tumour deletion category table:
  clinical_survival_deletion_category<- dplyr::full_join(short_survival_time_list, deletion_category_patient_ID, by = "patient_IDs") %>%
    dplyr::select(-patient_IDs)
  
  ## Remove entries with NA deletion categories due to missing clinical data
  # vars<- colnames(deletion_category_patient_ID[-1])
  # deletion_category_proximal_genes<- clinical_survival_deletion_category %>%
  #   dplyr::select(one_of(vars))
  # 
  # clinical_survival_deletion_category<- clinical_survival_deletion_category[complete.cases(deletion_category_proximal_genes), ]
  # 
  ## Remove entries with NA in any column
  clinical_survival_deletion_category<- clinical_survival_deletion_category[complete.cases(clinical_survival_deletion_category), ]
  
  
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
  survival_stats_table<- data.frame(matrix(NA, ncol = 14, nrow = ncol(clinical_survival_deletion_category) -3))
  names(survival_stats_table)<-c("target_gene", "proximal_gene", "p-value_Likelihood_ratio_test",
                                 "p-value_Wald_test", "p-value_logrank_test", "Hazard_ratio",
                                 "number_of_samples_cat1", "number_of_samples_cat2", "number_of_samples_cat3",
                                 "number_of_samples_cat4","mean_survival_cat_1","mean_survival_cat_2","mean_survival_cat_3"
                                 ,"mean_survival_cat_4")
  
  ## Loop through genes around target gene and calculate survival stats:
  for (i in 1:nrow(survival_stats_table)){
    
    ################here!!!!!
    
    ##Get covariable object for survfit:
    covariable_object<- clinical_survival_deletion_category[,i+3]
    ##Keep only entries with 1 and 2
    covariable_object_cat1_and_2<- which(covariable_object %in% c(1,2))
    covariable_object<- covariable_object[covariable_object_cat1_and_2]
    ##Remove values with death_time = NA
    death_time<- clinical_survival_deletion_category[covariable_object_cat1_and_2,1]
    death_event<- clinical_survival_deletion_category[covariable_object_cat1_and_2,2]
    stratification<- clinical_survival_deletion_category[covariable_object_cat1_and_2,3]
    # death_time<- death_time[which(!is.na(death_time))]
    # death_event<- death_event[which(!is.na(death_time))]
    # covariable_object<- covariable_object[which(!is.na(death_time))]
    # ##Remove values with death_event = NA
    # death_time<- death_time[which(!is.na(death_event))]
    # death_event<- death_event[which(!is.na(death_event))]
    # covariable_object<- covariable_object[which(!is.na(death_event))]
    # 
    ##Create Surv object if there are any tumours in catagory 1 and 2
    if(length(covariable_object) > 1){
      # death_time<- clinical_survival_deletion_category[covariable_object_cat1_and_2,1]
      # death_event<- clinical_survival_deletion_category[covariable_object_cat1_and_2,2]
      surv_object<- Surv(death_time, death_event==1)
      
      ##Fit Kaplain meier graph to one co-variable to compare data with :
      #fittedSurv <- survfit(surv_object~covariable_object, na.action = na.exclude)
      fittedSurv <- survfit(surv_object~covariable_object+stratification, na.action = na.exclude)
      #fittedSurv <- survfit(surv~dfCov[,1]+dfCov[,2], na.action = na.exclude)
      fittedSurv_mean<- survival:::survmean(fittedSurv, rmean="individual") 
      #survival:::survmean(fittedSurv, rmean="common")
      
      ##Get names of categories:
      if(length(unique(covariable_object)) > 1){
        #df.categ <- sapply(names(fittedSurv$strata), function(x) strsplit(x,"="))
        df.categ <- cbind(sapply(names(fittedSurv$strata), function(x) strsplit(x,"=")[[1]][3]),
                          fittedSurv$n)
        #df.categ <- data.frame(df.categ)
      }else {
        df.categ <- cbind(unique(covariable_object), fittedSurv$n)
        
      }
      temp_categ<- df.categ[,1]
      codel_cat<- c("deletion", "deletion", "co-deletion", "co-deletion")
      temp_categ<- paste(codel_cat, temp_categ, sep = " ")
      # temp_categ<- sub("\\<1\\>", paste(target_gene, "deletion", sep = " "), temp_categ)
      # temp_categ<- sub("\\<2\\>", paste(target_gene, "and", colnames(clinical_survival_deletion_category)[i+2], "co-deletion", sep = " "), temp_categ)
      # df.categ<- cbind(temp_categ, df.categ[,2]) 
      temp_categ<- paste(target_gene, temp_categ, sep = " ")
      #temp_categ<- sub("\\<2\\>", paste(target_gene, "and", colnames(clinical_survival_deletion_category)[i+2], "co-deletion", sep = " "), temp_categ)
      df.categ<- cbind(temp_categ, df.categ[,2])           
      ## Change name to 1 = Tumour suppressor deletion
      #2 = Co-deletion
      
      
      ##Get category names for one variable
      categNames <- apply(df.categ, 1, function(x) paste0(x[1]," (",x[2],")"))
      
      if (plot_graph == TRUE) {
        tiff(paste(target_gene, "_",  colnames(clinical_survival_deletion_category)[i+3],".tiff", sep =""), width = 7, height = 7, units = 'in', res = 100)
        plot(fittedSurv, main=paste(target_gene, "_",  colnames(clinical_survival_deletion_category)[i+3], " Survival", sep =""),
             xlab="Time (days)", ylab=ylabel,
             col=brewer.pal(9,"Set1"), 
             mark.time=T)
        legend("topright", legend=categNames,
               col=brewer.pal(9,"Set1"),
               lwd=2, cex=0.7)
        
      }
      
      if(print_to_screen == TRUE) {
        print("Chi-sq test:")
        print(survdiff(surv_object~covariable_object,rho = 0))
        
        print("Cox PH test:")
        print(summary(coxph(surv_object~covariable_object)))
      }
      coxfit <- coxph(surv_object~covariable_object)
      
      if (plot_graph == TRUE) {
        text(1000,0,labels=paste0("HR=",round(exp(summary(coxfit)$coefficients[1]),2),"; p=",
                                  round(summary(coxfit)$sctest[3],2)))
        dev.off()
      }
      
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
  
  setwd(current_wd)
  
  return(survival_stats_table)
  
}

####################
### Test function
cnv.table<- threshold_selected_cnv_list_plus_all_loc$ALL
##Get survival data
survival_time_list<- clinical_all_survival_pathologic_complete
## Get Target gene
target_gene_names<- sapply(gene_information_long_list, function(x) x[[1]])
index<- which(target_gene_names %in%  "RB1" )
target_gene_list<- gene_information_long_list[[index]]

##RB1
test<- stratified_survival_cat1V2(cnv.table =  cnv.table,
                                  survival_time_list =  survival_time_list,
                                  target_gene_list =  target_gene_list,
                           distance = 2.5e+06,
                           threshold = -2,
                           deletion = TRUE,
                           time_of_death_column = 5, 
                           death_event_column = 6,
                           stratification_column = 9,
                           column_start = 11,
                           start = TRUE, 
                           remove_NA = TRUE,
                           Cytoband = FALSE,
                           print_to_screen = FALSE, 
                           plot_graph = TRUE,
                           ylabel = "Overall survival",
                           path = "../../Output/Survival analysis/")

##CDKN1B
index<- which(target_gene_names %in%  "CDKN1B" )
target_gene_list<- gene_information_long_list[[index]]

test<- stratified_survival_cat1V2(cnv.table =  cnv.table,
                                  survival_time_list =  survival_time_list,
                                  target_gene_list =  target_gene_list,
                                  distance = 2.5e+06,
                                  threshold = -2,
                                  deletion = TRUE,
                                  time_of_death_column = 5, 
                                  death_event_column = 6,
                                  stratification_column = 9,
                                  column_start = 11,
                                  start = TRUE, 
                                  remove_NA = TRUE,
                                  Cytoband = FALSE,
                                  print_to_screen = FALSE, 
                                  plot_graph = TRUE,
                                  ylabel = "Overall survival",
                                  path = "../../Output/Survival analysis/")

##LRP1B
index<- which(target_gene_names %in%  "LRP1B" )
target_gene_list<- gene_information_long_list[[index]]
target_gene_list

test<- stratified_survival_cat1V2(cnv.table =  cnv.table,
                                  survival_time_list =  survival_time_list,
                                  target_gene_list =  target_gene_list,
                                  distance = 2.5e+06,
                                  threshold = -2,
                                  deletion = TRUE,
                                  time_of_death_column = 5, 
                                  death_event_column = 6,
                                  stratification_column = 9,
                                  column_start = 11,
                                  start = TRUE, 
                                  remove_NA = TRUE,
                                  Cytoband = FALSE,
                                  print_to_screen = FALSE, 
                                  plot_graph = TRUE,
                                  ylabel = "Overall survival",
                                  path = "../../Output/Survival analysis/")

##TGFBR2
index<- which(target_gene_names %in%  "TGFBR2" )
target_gene_list<- gene_information_long_list[[index]]
target_gene_list

test<- stratified_survival_cat1V2(cnv.table =  cnv.table,
                                  survival_time_list =  survival_time_list,
                                  target_gene_list =  target_gene_list,
                                  distance = 2.5e+06,
                                  threshold = -2,
                                  deletion = TRUE,
                                  time_of_death_column = 5, 
                                  death_event_column = 6,
                                  stratification_column = 9,
                                  column_start = 11,
                                  start = TRUE, 
                                  remove_NA = TRUE,
                                  Cytoband = FALSE,
                                  print_to_screen = FALSE, 
                                  plot_graph = TRUE,
                                  ylabel = "Overall survival",
                                  path = "../../Output/Survival analysis/")
##TP53
index<- which(target_gene_names %in%  "TP53" )
target_gene_list<- gene_information_long_list[[index]]
target_gene_list

test<- stratified_survival_cat1V2(cnv.table =  cnv.table,
                                  survival_time_list =  survival_time_list,
                                  target_gene_list =  target_gene_list,
                                  distance = 2.5e+06,
                                  threshold = -2,
                                  deletion = TRUE,
                                  time_of_death_column = 5, 
                                  death_event_column = 6,
                                  stratification_column = 9,
                                  column_start = 11,
                                  start = TRUE, 
                                  remove_NA = TRUE,
                                  Cytoband = FALSE,
                                  print_to_screen = FALSE, 
                                  plot_graph = TRUE,
                                  ylabel = "Overall survival",
                                  path = "../../Output/Survival analysis/")

##ZFHX3
index<- which(target_gene_names %in%  "ZFHX3" )
target_gene_list<- gene_information_long_list[[index]]

test<- stratified_survival_cat1V2(cnv.table =  cnv.table,
                                  survival_time_list =  survival_time_list,
                                  target_gene_list =  target_gene_list,
                                  distance = 2.5e+06,
                                  threshold = -2,
                                  deletion = TRUE,
                                  time_of_death_column = 5, 
                                  death_event_column = 6,
                                  stratification_column = 9,
                                  column_start = 11,
                                  start = TRUE, 
                                  remove_NA = TRUE,
                                  Cytoband = FALSE,
                                  print_to_screen = FALSE, 
                                  plot_graph = TRUE,
                                  ylabel = "Overall survival",
                                  path = "../../Output/Survival analysis/")



