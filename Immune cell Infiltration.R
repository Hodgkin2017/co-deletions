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

##########
##Create a list of dataframes?


######################################
### Function to perform Immune cell infiltration analysis.

#####
##Parameters:
target_gene_list<- gene_information_list[[2]]
immune_cell_infiltrate_table<- tcia_immune_infiltrate
cnv.table<- threshold_selected_cnv_list_plus_all_loc$ALL
dim(cnv.table)
distance = 2.5e+06
threshold = -2
deletion = TRUE
infiltrate_column = 4  
column_start = 11
start = TRUE 
remove_NA = TRUE
Cytoband = FALSE
print_to_screen = FALSE 
plot_graph = FALSE
path = "./"


target_gene_list = list1
immune_cell_infiltrate_table = table1
cnv.table = table2

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
  
  tcia_immune_infiltrate_table<- do.call(rbind, tcia_immune_infiltrate_list)

  
  ################################################
  ###  analysis:
  

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
identical(test, tcia_immune_infiltrate_table)

###########
###ANOVA analysis?....


