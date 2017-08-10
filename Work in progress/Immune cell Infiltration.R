#####################
### Immune cell infiltration
######################

#####################################
###How many immune genes have I identified using my statistical tests?


#####################################
###Import Immune cell infiltrate table
tcia_immune_infiltrate<- read.delim("../../Input data/TCIA/cellTypeFractionsAll.tsv", header = TRUE, stringsAsFactors = FALSE)
dim(tcia_immune_infiltrate)
head(tcia_immune_infiltrate)

#######
##Order table by cancer type then by immune cell type
tcia_immune_infiltrate<- dplyr::arrange(tcia_immune_infiltrate, disease, cell_type)
head(tcia_immune_infiltrate)

########
##How many unique TCGA barcodes are there?
unique(tcia_immune_infiltrate$patientBarcode)
length(unique(tcia_immune_infiltrate$patientBarcode))

##########
##How many cancer types are there?
unique(tcia_immune_infiltrate$disease)
length(unique(tcia_immune_infiltrate$disease))

##########
##How many immune cell types are there?
unique(tcia_immune_infiltrate$cell_type)
length(unique(tcia_immune_infiltrate$cell_type))

##########
## How many of my CNV data samples do I have Immune infiltration data for?
#threshold_selected_cnv_list_plus_all_loc contains all the CNV data plus the last 
#item table the list is all samples combined
##Get CNV patient barcode
CNV_patient_ids<- colnames(threshold_selected_cnv_list_plus_all_loc$ALL)
head(CNV_patient_ids,11)
CNV_patient_ids<- CNV_patient_ids[11:length(CNV_patient_ids)]
head(CNV_patient_ids,11)

##Convert CNV_patient_ids patient barcode to the same format as tcia_immune_infiltrate
CNV_patient_ids<- CNV_patient_ids %>%
  substr(0, 12) %>%
  gsub("\\.", "-", .)

CNV_patient_ids

length(which(CNV_patient_ids %in% tcia_immune_infiltrate$patientBarcode))
# Total = 7045!
length(which( unique(tcia_immune_infiltrate$patientBarcode) %in% CNV_patient_ids))
# Total = 7045! All tcia TCGA samples have a matching CNV bar code.

###########
###Is cancer type CRC = COADREAD? (coloreactal cancer)
CNV_patient_ids<- colnames(threshold_selected_cnv_list_plus_all_loc$COADREAD)
head(CNV_patient_ids,11)
CNV_patient_ids<- CNV_patient_ids[11:length(CNV_patient_ids)]
head(CNV_patient_ids,11)

##Convert CNV_patient_ids patient barcode to the same format as tcia_immune_infiltrate
CNV_patient_ids<- CNV_patient_ids %>%
  substr(0, 12) %>%
  gsub("\\.", "-", .)

CNV_patient_ids

CNV_patient_ids_CRC<- tcia_immune_infiltrate %>%
  dplyr::filter(disease == "CRC") %>%
  dplyr::select(patientBarcode) %>%
  .$patientBarcode

length(which(CNV_patient_ids %in% CNV_patient_ids_CRC))
# Total = 613
length(which( unique(CNV_patient_ids_CRC) %in% unique(CNV_patient_ids)))
# Total = 613!





##########
##Create a list of dataframes?


######################################
### Function to perform Immune cell infiltration analysis.
## Function takes a CNV table, immune_infiltrate table and a Tumour suppressor gene


#####
##Parameters:
# target_gene_list<- gene_information_list[[2]]
# immune_cell_infiltrate_table<- tcia_immune_infiltrate
# cnv.table<- threshold_selected_cnv_list_plus_all_loc$ALL
# dim(cnv.table)
# distance = 2.5e+06
# threshold = -2
# deletion = TRUE
# infiltrate_column = 4  
# column_start = 11
# start = TRUE 
# remove_NA = TRUE
# Cytoband = FALSE
# print_to_screen = FALSE 
# plot_graph = FALSE
# path = "./"
# 
# 
# target_gene_list = list1
# immune_cell_infiltrate_table = table1
# cnv.table = table2

immune_cell_infiltrate<- function(target_gene_list,
                                          immune_cell_infiltrate_table,
                                          cnv.table,
                                          distance = 2.5e+06,
                                          threshold = -2,
                                          deletion = TRUE,
                                          infiltrate_column, 
                                          column_start = 11,
                                          start = TRUE, 
                                          remove_NA = TRUE,
                                          Cytoband = FALSE,
                                          print_to_screen = FALSE, 
                                          plot_graph = FALSE,
                                  path = "./"){
  
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
  ### Join tumour category table above to immune infiltrate table:
  
  deletion_category<-t(deletion_category_table)
  
  ##Convert patient IDs in deletion_category so they match with the patient IDS in the clinical table
  deletion_category_patient_ID<- rownames(deletion_category) %>%
    substr(0, 12) %>%
    gsub("\\.", "-", .) %>%
    cbind.data.frame(patientBarcode = ., deletion_category)
  
  deletion_category_patient_ID$patientBarcode<- as.character(deletion_category_patient_ID$patientBarcode)
  
  ##Get names of infiltrating cell types:
  tcia_immune_infiltrate_cell_type<- unique(tcia_immune_infiltrate$cell_type)
  tcia_immune_infiltrate_list<- vector("list", length(tcia_immune_infiltrate_cell_type))
  
  for (i in 1: length(tcia_immune_infiltrate_cell_type)){
    
    rows_to_keep<- immune_cell_infiltrate_table$cell_type == tcia_immune_infiltrate_cell_type[i]
    small_immune_cell_infiltrate_table<- immune_cell_infiltrate_table[rows_to_keep,]
    ##Join small_immune_cell_infiltrate_table with CNV data
    small_immune_cell_infiltrate_table<- dplyr::full_join(small_immune_cell_infiltrate_table, deletion_category_patient_ID, by = "patientBarcode")
    ##Remove NA entries
    #column_names<- colnames(small_immune_cell_infiltrate_table)
    small_immune_cell_infiltrate_table<- small_immune_cell_infiltrate_table[!is.na(small_immune_cell_infiltrate_table[,ncol(small_immune_cell_infiltrate_table)]),]
    ##Add to list
    tcia_immune_infiltrate_list[[i]]<- small_immune_cell_infiltrate_table
  }
  
  ## Join immun cell type list together and add column for target gene
  tcia_immune_infiltrate_table<- do.call(rbind, tcia_immune_infiltrate_list)
  target_gene_column<-rep(target_gene, nrow(tcia_immune_infiltrate_table))
  tcia_immune_infiltrate_table<- cbind(target_gene = target_gene_column, tcia_immune_infiltrate_table)

  return(tcia_immune_infiltrate_table)
  
}

#########
### Test function

list1<- gene_information_list[[2]]
table1<- tcia_immune_infiltrate
table2<- threshold_selected_cnv_list_plus_all_loc$ALL

# target_gene_list = gene_information_list[[2]]
# immune_cell_infiltrate_table = tcia_immune_infiltrate
# cnv.table = threshold_selected_cnv_list_plus_all_loc$ALL


test<- immune_cell_infiltrate(target_gene_list = list1,
                              immune_cell_infiltrate_table = table1,
                              cnv.table = table2, 
                              distance = 2.5e+06,
                              threshold = -2,
                              deletion = TRUE,
                              infiltrate_column = 4,  
                              column_start = 11,
                              start = TRUE, 
                              remove_NA = TRUE,
                              Cytoband = FALSE,
                              print_to_screen = FALSE, 
                              plot_graph = FALSE,
                              path = "./")

dim(test)
#identical(test, tcia_immune_infiltrate_table)
test[1:10,]

table2<- threshold_selected_cnv_list_plus_all_loc$BRCA

test2<- immune_cell_infiltrate(target_gene_list = list1,
                              immune_cell_infiltrate_table = table1,
                              cnv.table = table2, 
                              distance = 2.5e+06,
                              threshold = -2,
                              deletion = TRUE,
                              infiltrate_column = 4,  
                              column_start = 11,
                              start = TRUE, 
                              remove_NA = TRUE,
                              Cytoband = FALSE,
                              print_to_screen = FALSE, 
                              plot_graph = FALSE,
                              path = "./")

dim(test2)
test2[1:10,]

list2<- gene_information_list[1:2]
table2<- threshold_selected_cnv_list_plus_all_loc$ALL
table2<- threshold_selected_cnv_list_plus_all_loc$BRCA

test_lapply<- lapply(list2, function(x) immune_cell_infiltrate(target_gene_list = x,
                                                          immune_cell_infiltrate_table = table1,
                                                          cnv.table = table2, 
                                                          distance = 2.5e+06,
                                                          threshold = -2,
                                                          deletion = TRUE,
                                                          infiltrate_column = 4,  
                                                          column_start = 11,
                                                          start = TRUE, 
                                                          remove_NA = TRUE,
                                                          Cytoband = FALSE,
                                                          print_to_screen = FALSE, 
                                                          plot_graph = FALSE,
                                                          path = "./")
                     )

length(test_lapply)
dim(test_lapply[[1]])
test_lapply[[1]][1:5, 1:10]
identical(test2, test_lapply[[2]])


###################
###ANOVA analysis of immune cell infiltrate:
########################




col_start = 7
#ALL_immune_infiltrate_table<- test_lapply[[2]]

immune_infiltrate_table<- test_lapply[[2]]
immune_infiltrate_table<- ALL_immune_infiltrate_table
immune_infiltrate_table<- test_lapply[[1]]

immune_cell_infiltrate_annova<- function(immune_infiltrate_table, col_start, join_genes = TRUE){

## Get immune cell types
cell_types_names<- unique(immune_infiltrate_table$cell_type) 

##Remove NA
cell_types_names<- cell_types_names[complete.cases(cell_types_names)]

##Create an empty list:
annova_immune_infiltration_list<- vector("list", ncol(immune_infiltrate_table) - col_start + 1)

for (j in col_start:ncol(immune_infiltrate_table)){


##Create an empty vector:
annova_immune_infiltration_results<- data.frame(matrix(NA, ncol = 19, nrow = length(cell_types_names)))
colnames(annova_immune_infiltration_results)<- c("cancer", "target_gene", "proximal_gene", "cell_type",
                                                 "number_cat1", "number_cat2",
                                                 "number_cat3", "number_cat4",
                                                 "mean_cibersort_cat1", "mean_cibersort_cat2",
                                                 "mean_cibersort_cat3", "mean_cibersort_cat4",
                                                 "ANOVA_p_value", "p_value_cat2_1", "p_value3_1", 
                                                 "p_value4_1", "p_value3_2", "p_value4_2", "p_value4_3")

for (i in 1: length(cell_types_names)){
  
  ##Filter by immune cell type
  cell_type_name<- cell_types_names[i]
  
  cell_type_immune_table<- immune_infiltrate_table %>%
    dplyr::filter(cell_type == cell_type_name)
  
  ##Perform ANNOVA
  if(length(levels(as.factor(cell_type_immune_table[,j]))) > 1){
  anova_test<- aov(cell_type_immune_table$cibersort_LM22 ~ as.factor(cell_type_immune_table[,j]))
  anova_summary<- summary(anova_test)
  
  ##Perform TUKEY test
  tukey_test<- TukeyHSD(anova_test)
  
  }
  
  ##Get mean cibersort value per deletion category
  means<- round(tapply(cell_type_immune_table$cibersort_LM22, factor(cell_type_immune_table[,7], levels = c(1,2,3,4)), mean), digits=4)
  
  ##number per catagory
  number_per_category<- table(factor(cell_type_immune_table[,j], levels = c(1,2,3,4)))
  
  ##Fill in table
  annova_immune_infiltration_results[i,1]<- immune_infiltrate_table$disease[1]
  annova_immune_infiltration_results[i,2]<- as.character(immune_infiltrate_table$target_gene[1])
  annova_immune_infiltration_results[i,3]<- colnames(immune_infiltrate_table)[j]
  annova_immune_infiltration_results[i,4]<- cell_type_name
  annova_immune_infiltration_results[i,5]<- number_per_category[1]
  annova_immune_infiltration_results[i,6]<- number_per_category[2]
  annova_immune_infiltration_results[i,7]<- number_per_category[3]
  annova_immune_infiltration_results[i,8]<- number_per_category[4]
  annova_immune_infiltration_results[i,9]<- means[1]
  annova_immune_infiltration_results[i,10]<- means[2]
  annova_immune_infiltration_results[i,11]<- means[3]
  annova_immune_infiltration_results[i,12]<- means[4]
  
  if(length(levels(as.factor(cell_type_immune_table[,j]))) > 1){
  annova_immune_infiltration_results[i,13]<- anova_summary[[1]]$`Pr(>F)`[1]

  index<- grep("2-1", rownames(tukey_test$`as.factor(cell_type_immune_table[, j])`))
  if(length(index) == 1){
  annova_immune_infiltration_results[i,14]<- tukey_test$`as.factor(cell_type_immune_table[, j])`["2-1",4]
  }
  index<- grep("3-1", rownames(tukey_test$`as.factor(cell_type_immune_table[, j])`))
  if(length(index) == 1){
  annova_immune_infiltration_results[i,15]<- tukey_test$`as.factor(cell_type_immune_table[, j])`["3-1",4]
  }
  index<- grep("4-1", rownames(tukey_test$`as.factor(cell_type_immune_table[, j])`))
  if(length(index) == 1){
  annova_immune_infiltration_results[i,16]<- tukey_test$`as.factor(cell_type_immune_table[, j])`["4-1",4]
  }
  index<- grep("3-2", rownames(tukey_test$`as.factor(cell_type_immune_table[, j])`))
  if(length(index) == 1){
  annova_immune_infiltration_results[i,17]<- tukey_test$`as.factor(cell_type_immune_table[, j])`["3-2",4]
  }
  index<- grep("4-2", rownames(tukey_test$`as.factor(cell_type_immune_table[, j])`))
  if(length(index) == 1){
  annova_immune_infiltration_results[i,18]<- tukey_test$`as.factor(cell_type_immune_table[, j])`["4-2",4]
  }
  index<- grep("4-3", rownames(tukey_test$`as.factor(cell_type_immune_table[, j])`))
  if(length(index) == 1){
  annova_immune_infiltration_results[i,19]<- tukey_test$`as.factor(cell_type_immune_table[, j])`["4-3",4]
  
  }
  
  }
  
}

print(paste("Cancer:", annova_immune_infiltration_results[i,1],
            "Target gene:", annova_immune_infiltration_results[i,2],
            "Proximal gene:", annova_immune_infiltration_results[i,3]), sep = " ")


annova_immune_infiltration_list[[j]]<- annova_immune_infiltration_results

}

if(join_genes == TRUE){
  
  annova_immune_infiltration_list<- do.call(rbind, annova_immune_infiltration_list)
  
}

return(annova_immune_infiltration_list)

}

###################
### test function

test_annova<- immune_cell_infiltrate_annova(immune_infiltrate_table = test_lapply[[2]], col_start = 7, join_genes = TRUE)

test_annova_apply<- lapply(test_lapply, function(x) immune_cell_infiltrate_annova(x, col_start = 7, join_genes = TRUE))

length(test_annova_apply)
test_annova_apply[[1]]

###########################
### Test all cancers and all Tumour suppressor genes to see if they have 
#significant immune cell infiltration 
#########################


##########
##How many cancer types are there?
unique(tcia_immune_infiltrate$disease)
length(unique(tcia_immune_infiltrate$disease))

names(threshold_selected_cnv_list_plus_all_loc)

##CRC = COADREAD

############
### Create list of CNV tables that I have immune infiltration data for

threshold_immune_selected_cnv_list_plus_all_loc<- threshold_selected_cnv_list_plus_all_loc[c(2,3,4,6,9,10,11,12,13,16,17,18,20,21,23,25,26,28,30,33)]
length(threshold_immune_selected_cnv_list_plus_all_loc)
identical(names(threshold_immune_selected_cnv_list_plus_all_loc[1:19]), unique(tcia_immune_infiltrate$disease))

#############

###create empty list
immune_cell_infiltrate_annova_per_cancer_list<- vector("list", length(threshold_immune_selected_cnv_list_plus_all_loc))

###loop through each cancer
for ( i in 1: length(threshold_immune_selected_cnv_list_plus_all_loc)){
  
  cnv_table<-threshold_immune_selected_cnv_list_plus_all_loc[[i]]
  
  ##Determine immune cell infiltration and type od deletion present per 
  #tumour and per target gene - proximal gene pair
  immune_cell_infiltrate_list<- lapply(gene_information_long_list[1:2], function(x) immune_cell_infiltrate(target_gene_list = x,
                                                                 immune_cell_infiltrate_table = tcia_immune_infiltrate,
                                                                 cnv.table = cnv_table, 
                                                                 distance = 2.5e+06,
                                                                 threshold = -2,
                                                                 deletion = TRUE,
                                                                 infiltrate_column = 4,  
                                                                 column_start = 11,
                                                                 start = TRUE, 
                                                                 remove_NA = TRUE,
                                                                 Cytoband = FALSE,
                                                                 print_to_screen = FALSE, 
                                                                 plot_graph = FALSE,
                                                                 path = "./")
  )
  
  ##Name each item in list with target gene name
  names(immune_cell_infiltrate_list)<- lapply(gene_information_long_list[1:2], function(x) x[[1]])
  
  ##Perform ANNOVA within each immune cell type for each target gene and proximal gene:
  immune_cell_infiltrate_annova_list<- lapply(immune_cell_infiltrate_list, function(x) immune_cell_infiltrate_annova(x, col_start = 7, join_genes = TRUE))
  
  ##Join target gene ANNOVA tables in list together and add to final results list
  immune_cell_infiltrate_annova_long_table<- do.call(rbind, immune_cell_infiltrate_annova_list)
  immune_cell_infiltrate_annova_per_cancer_list[[i]]<- immune_cell_infiltrate_annova_long_table
  
}

names(immune_cell_infiltrate_annova_per_cancer_list)<- names(threshold_immune_selected_cnv_list_plus_all_loc)

## Save object:
saveRDS(immune_cell_infiltrate_annova_per_cancer_list, file = "./R workspaces/immune_cell_infiltrate_annova_per_cancer_list")
















