#################
### Survival Analysis of co-deleted genes
#################

#################
### Create a list of tables for each cancer containing overall survival and disease free survival
#################








###############
### Characterise Tumours based on their type of deletion and perform survival analysis:
##############

##Parameters:
gene_information_list
gene_information<- gene_information_list[[2]]
gene_information
cnv.table<- threshold_short_cnv_list_loc[[1]]
dim(cnv.table)

remove_NA = TRUE
start = TRUE
Chromosome<- gene_information[[2]]
distance<- 2.5e+06
selection_criteria<- c(gene_information[[4]]-distance, gene_information[[5]]+distance)
column_start = 11
deletion = TRUE
threshold = -2
target_gene = gene_information[[1]]

##########
###Function to take a target gene and surrounding genes and characterise each tumour based on if it had a deletion 
#of the target gene and each suurounbding gene 

categorise_deletion_type_around_target_gene<- function(cnv.table, target_gene, Chromosome, selection_criteria, threshold = -2, 
                                                       deletion = TRUE, column_start = 11, start = TRUE, remove_NA = TRUE, Cytoband = FALSE) {


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

dim(cnv.matrix)
cnv.matrix[,10:12]

##Create a binary matrix of CNV data such that deletions (deletion = TRUE) or 
#amplifications (deletion = FALSE) below (for deletions) or above (for amplifications) a threshold = 1 
if (deletion == TRUE) {
  
  cnv.matrix<- ifelse(cnv.matrix <= threshold, 1, 0)
  
} else {
  
  cnv.matrix<- ifelse(cnv.matrix >= threshold, 1, 0)
}

dim(cnv.matrix)
which(colSums(cnv.matrix) > 0, arr.ind = T)
cnv.matrix[,170:172]

###########
###Create a new table with genes in rows and tumours on columns and each entry has a number from 1 to 4 depending on which 
#deletion/amplification catagory it belongs to:

##Add gene names as new column
gene_names_cnv.matrix<- cbind.data.frame(gene = row.names(cnv.matrix) ,cnv.matrix)
gene_names_cnv.matrix$gene<-as.character(gene_names_cnv.matrix$gene)
gene_names_cnv.matrix[,c(1, 170:172)]
gene_names<- row.names(cnv.matrix)

##########
categorise_deletion_type_function<- function(gene_names, gene_names_cnv.matrix){
##Determine if target gene is deleted
target_gene_deleted<- gene_names_cnv.matrix %>%
  dplyr::filter(gene == target_gene) %>%
  dplyr::select(-gene)

dim(target_gene_deleted)
target_gene_deleted[,c(1, 170:172)]

target_gene_deleted<- as.vector(target_gene_deleted == 1)
target_gene_deleted[170:173]

##Determine if proximal gene is deleted
proximal_gene_deleted<- gene_names_cnv.matrix %>%
  dplyr::filter(gene == gene_names) %>%
  dplyr::select(-gene)

  dim(proximal_gene_deleted)
  proximal_gene_deleted[,c(1, 170:172)]
  
  proximal_gene_deleted<- as.vector(proximal_gene_deleted == 1)
  proximal_gene_deleted[170:173]

##Determine if target gene and not proximal gene is deleted(category 1)
deletion_category1<- proximal_gene_deleted - target_gene_deleted < 0
deletion_category1

##Determine if target gene and proximal gene are deleted together (category 2)
deletion_category2<- target_gene_deleted + proximal_gene_deleted > 1
deletion_category2

##Determine if proximal gene is deleted without target gene (category 3)
deletion_category3<- target_gene_deleted - proximal_gene_deleted < 0
deletion_category3

##Determine if neither proximal gene or target gene is deleted (category 4)
deletion_category4<- target_gene_deleted + proximal_gene_deleted == 0
deletion_category4

all_deletion_category_table<- cbind(deletion_category1, deletion_category2, deletion_category3, deletion_category4)
all_deletion_category_table
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

return(deletion_category_table)

}
#################
### Test function

# cnv.table, target_gene, Chromosome, selection_criteria, threshold = -2, 
# column_start = 11, start = TRUE, remove_NA = TRUE, Cytoband = FALSE)

x<- gene_information_list[[2]]
distance<- 2.5e+06
##BRCA CNV table:
cnv.table<- threshold_short_cnv_list_loc[[1]]

test3<- categorise_deletion_type_around_target_gene(cnv.table = cnv.table, target_gene = x[[1]], Chromosome = x[[2]], 
                                            selection_criteria = c(x[[4]] - distance, x[[5]]+ distance), threshold = -2, deletion = TRUE)

#identical(deletion_category_table, test3)

##Test function in lapply with a list of target genes:
my_list<- gene_information_list[1:2]
test4<- lapply(my_list, function(x) categorise_deletion_type_around_target_gene(cnv.table = cnv.table, target_gene = x[[1]], Chromosome = x[[2]], 
                                                                             selection_criteria = c(x[[4]] - distance, x[[5]]+ distance), threshold = -2))


identical(test4[[2]], test3)

test4[[1]][,1:3]
test4[[2]][,1:3]










####################
###Perform survival analysis for CDKN2A and MTAP
####################

#Parameters:
test4[[2]][,1:3]
proximal_gene<- "MTAP"
#survival data
clinical_survival_list[[1]]

####################
###Get gene deletion data for MTAP
## Add extra column with patient ID:
deletion_category_gene_name_table<- cbind.data.frame(gene = row.names(test4[[2]]) ,test4[[2]])
deletion_category_gene_name_table$gene<-as.character(deletion_category_gene_name_table$gene)
deletion_category_gene_name_table[1:2, 1:5]


deletion_category<- deletion_category_gene_name_table %>% 
  dplyr::filter(gene == proximal_gene) %>%
  dplyr::select(-gene)

dim(deletion_category)
deletion_category[1, 1:5]
deletion_category<-t(deletion_category)
dim(deletion_category)
table(deletion_category[,1])
deletion_category

rownames(deletion_category)

deletion_category_patient_ID<- rownames(deletion_category) %>%
  substr(0, 12) %>%
  gsub("\\.", "-", .) %>%
  cbind.data.frame(patient_IDs = ., deletion_category)

deletion_category_patient_ID$patient_IDs<- as.character(deletion_category_patient_ID$patient_IDs)

class(deletion_category_patient_ID$patient_IDs)
deletion_category_patient_ID
clinical_survival_list[[1]]$patient_IDs



clinical_survival_deletion_category<- dplyr::full_join(clinical_survival_list[[1]], deletion_category_patient_ID, by = "patient_IDs")
head(clinical_survival_deletion_category, 40)
dim(clinical_survival_deletion_category)
dim(clinical_survival_list[[1]])
##Comment 18 tumours do not match between clinical and cnv data
sum(is.na(clinical_survival_deletion_category$deletion_category))
sum(is.na(test4[[2]][2,]))




########
### Survival analysis

##Overall survival:
s<- survfit(Surv(new_death, event = death_event == 1)~1, data = clinical_survival_deletion_category)
s
summary(s)
plot(s)

##Disease free survival
s1<- survfit(Surv(disease_free_survival, event = death_event == 1)~1, data = clinical_survival_deletion_category)
s1
summary(s1)
plot(s1)
plot(s1, mark.time = T)

##Overall survival by deletion_category

s<- survfit(Surv(new_death, event = death_event == 1)~deletion_category, data = clinical_survival_deletion_category)
s
summary(s)
plot(s)
print(survfit(Surv(new_death, event = death_event == 1)~deletion_category, data = clinical_survival_deletion_category), print.rmean=TRUE)

##Disease free survival by deletion_category
s1<- survfit(Surv(disease_free_survival, event = death_event == 1)~deletion_category, data = clinical_survival_deletion_category)
s1
summary(s1)
plot(s1)
plot(s1, mark.time = T)

###############
## log-rank test:

#surv_object<- Surv(disease_free_survival, event = death_event == 1, data = clinical_survival_deletion_category)
survdiff(Surv(disease_free_survival, event = death_event == 1)~deletion_category, data = clinical_survival_deletion_category)

###############
##Cox-ph model:

coxph(Surv(disease_free_survival, event = death_event == 1)~deletion_category, data = clinical_survival_deletion_category)
summary(coxph(Surv(disease_free_survival, event = death_event == 1)~deletion_category, data = clinical_survival_deletion_category))


#############
##Cox-ph model using just deletion_categories 1 and 2 only!

clinical_survival_deletion_category_1_2<- clinical_survival_deletion_category %>%
  dplyr::filter(deletion_category == 1 | deletion_category == 2)

table(clinical_survival_deletion_category_1_2$deletion_category)

##Overall survival:
summary(coxph(Surv(new_death, event = death_event == 1)~deletion_category, data = clinical_survival_deletion_category_1_2))
##Disease free survival
summary(coxph(Surv(disease_free_survival, event = death_event == 1)~deletion_category, data = clinical_survival_deletion_category_1_2))





############
## Select first gene and catagorise tumours 1,2,3 or 4 depending on what type of deletion each tumour has
# Make a vector for each tumour in BRCA if it contains a deletion in the target gene only (1), 
#target gene and proximal gene (2), proximal gene only (3) or deletion in neither gene (4)

###########
## Perform survival analysis and record p-value and mean/median survival, number of samples tested
#print(survfit(my.surv ~ 1), print.rmean=TRUE)




##########







#################
###
################





















