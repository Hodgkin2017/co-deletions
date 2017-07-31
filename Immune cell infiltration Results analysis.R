########################
### Immune cell infiltration Results analysis
#######################

##########################
### All target genes

###############
###Import data

immune_cell_infiltrate_annova_per_cancer_list<- readRDS(file = "./R workspaces/immune_cell_infiltrate_annova_per_cancer_list")
head(immune_cell_infiltrate_annova_per_cancer_list[[1]])

length(immune_cell_infiltrate_annova_per_cancer_list)
names(immune_cell_infiltrate_annova_per_cancer_list)
colnames(immune_cell_infiltrate_annova_per_cancer_list[[1]])


###############
###Perform BH multiple testing correction

##Add a column to each table in list with the cancer name

dim(immune_cell_infiltrate_annova_per_cancer_list[[1]])
immune_cell_infiltrate_annova_per_cancer_list[[20]][1:10,]
names(immune_cell_infiltrate_annova_per_cancer_list[1])

##It appears the name of the cancer for ALL was incorrectly imputted in table.
sapply(immune_cell_infiltrate_annova_per_cancer_list, function(x) x$cancer[1])

immune_cell_infiltrate_annova_per_cancer_list[[20]] %>% 
  dplyr::filter(target_gene == "CDKN2A") %>%
  dplyr::select(number_cat1, number_cat2, number_cat3, number_cat4)

##Correct cancer name in immune_cell_infiltrate_annova_per_cancer_list$ALL table
immune_cell_infiltrate_annova_per_cancer_list[[20]]$cancer<- rep("ALL", nrow(immune_cell_infiltrate_annova_per_cancer_list[[20]]))
sapply(immune_cell_infiltrate_annova_per_cancer_list, function(x) x$cancer[1])
immune_cell_infiltrate_annova_per_cancer_list[[20]][1:10,]

## Join tables together except for 'all'
immune_cell_infiltrate_annova_per_cancer_table<- do.call(rbind, immune_cell_infiltrate_annova_per_cancer_list[1:19])
dim(immune_cell_infiltrate_annova_per_cancer_list)
dim(immune_cell_infiltrate_annova_per_cancer_table)
unique(immune_cell_infiltrate_annova_per_cancer_table$cancer)

########
###Split each cell type into a different table and add to a list

##Get the number of cell types:
immune_cell_types<- unique(immune_cell_infiltrate_annova_per_cancer_table$cell_type)
length(immune_cell_types)
immune_cell_types

## Create an empty list
immune_cell_infiltrate_annova_per_cell_type_list<- vector("list", length(immune_cell_types))
immune_cell_infiltrate_annova_per_cell_type_list

## For loop to create list:
for (i in 1: length(immune_cell_types)){
  
  immune_cell_infiltrate_annova_per_cell_type_list[[i]]<- immune_cell_infiltrate_annova_per_cancer_table %>%
    dplyr::filter(cell_type == immune_cell_types[i])
}

#############
###BH multiple testing for cancers 1:19 in list per cell type:
for (i in 1: length(immune_cell_infiltrate_annova_per_cell_type_list)){
  
  immune_table<- immune_cell_infiltrate_annova_per_cell_type_list[[i]]
  immune_table$BH_adjust_ANOVA<- p.adjust(immune_table$ANOVA_p_value, method = "BH")
  immune_table$BH_adjust_cat2_1<- p.adjust(immune_table$p_value_cat2_1, method = "BH")
  immune_table$BH_adjust_cat3_1<- p.adjust(immune_table$p_value3_1, method = "BH")
  immune_table$BH_adjust_cat4_1<- p.adjust(immune_table$p_value4_1, method = "BH")
  immune_table$BH_adjust_cat3_2<- p.adjust(immune_table$p_value3_2, method = "BH")
  immune_table$BH_adjust_cat4_2<- p.adjust(immune_table$p_value4_2, method = "BH")
  immune_table$BH_adjust_Cat4_3<- p.adjust(immune_table$p_value4_3, method = "BH")
  immune_table<- dplyr::arrange(immune_table, BH_adjust_cat2_1)
  immune_cell_infiltrate_annova_per_cell_type_list[[i]]<- immune_table
}

############
###BH multiple testing for all cancers (item 20 in list)
immune_cell_infiltrate_annova_per_cancer_all_table<- immune_cell_infiltrate_annova_per_cancer_list[[20]]
dim(immune_cell_infiltrate_annova_per_cancer_all_table)

########
###Split each cell type into a different table and add to a list

##Get the number of cell types:
immune_cell_types<- unique(immune_cell_infiltrate_annova_per_cancer_all_table$cell_type)
length(immune_cell_types)
immune_cell_types

## Create an empty list
immune_cell_infiltrate_annova_per_cell_type_all_list<- vector("list", length(immune_cell_types))
immune_cell_infiltrate_annova_per_cell_type_all_list

## For loop to create list:
for (i in 1: length(immune_cell_types)){
  
  immune_cell_infiltrate_annova_per_cell_type_all_list[[i]]<- immune_cell_infiltrate_annova_per_cancer_all_table %>%
    dplyr::filter(cell_type == immune_cell_types[i])
}

length(immune_cell_infiltrate_annova_per_cell_type_all_list)
dim(immune_cell_infiltrate_annova_per_cell_type_all_list[[1]])

#############
###BH multiple testing for all cancers table in list per cell type:
for (i in 1: length(immune_cell_infiltrate_annova_per_cell_type_list)){
  
  immune_table<- immune_cell_infiltrate_annova_per_cell_type_all_list[[i]]
  immune_table$BH_adjust_ANOVA<- p.adjust(immune_table$ANOVA_p_value, method = "BH")
  immune_table$BH_adjust_cat2_1<- p.adjust(immune_table$p_value_cat2_1, method = "BH")
  immune_table$BH_adjust_cat3_1<- p.adjust(immune_table$p_value3_1, method = "BH")
  immune_table$BH_adjust_cat4_1<- p.adjust(immune_table$p_value4_1, method = "BH")
  immune_table$BH_adjust_cat3_2<- p.adjust(immune_table$p_value3_2, method = "BH")
  immune_table$BH_adjust_cat4_2<- p.adjust(immune_table$p_value4_2, method = "BH")
  immune_table$BH_adjust_Cat4_3<- p.adjust(immune_table$p_value4_3, method = "BH")
  immune_table<- dplyr::arrange(immune_table, BH_adjust_cat2_1)
  immune_cell_infiltrate_annova_per_cell_type_all_list[[i]]<- immune_table
}

##############
###Join all cancer list with individual cancer list:
unique(immune_cell_infiltrate_annova_per_cell_type_list[[1]]$cancer)
dim(immune_cell_infiltrate_annova_per_cell_type_list[[1]])
dim(immune_cell_infiltrate_annova_per_cell_type_all_list[[1]])

for (i in 1: length(immune_cell_infiltrate_annova_per_cell_type_list)) {
  
  immune_cell_infiltrate_annova_per_cell_type_list[[i]]<- rbind(immune_cell_infiltrate_annova_per_cell_type_list[[i]],
                                                                immune_cell_infiltrate_annova_per_cell_type_all_list[[i]])
  
  
  
  
}

dim(immune_cell_infiltrate_annova_per_cell_type_list[[1]])
unique(immune_cell_infiltrate_annova_per_cell_type_list[[1]]$cancer)
sapply(immune_cell_infiltrate_annova_per_cell_type_list, function(x) unique(x$cell_type))
lapply(immune_cell_infiltrate_annova_per_cell_type_list, function(x) unique(x$cancer))

#####################
###Identify significant genes per immune cell type
sapply(immune_cell_infiltrate_annova_per_cell_type_list, function(x) unique(x$cell_type))
##  "Resting Natural killer cell" = 12, "Activated Natural killer cell" = 13, "CD8 T cell" = 19

###########
### Find the number of significant genes per cell type:
colnames(immune_cell_infiltrate_annova_per_cell_type_list[[1]])
##p<= 0.05
sapply(immune_cell_infiltrate_annova_per_cell_type_list, function(x) sum(x$BH_adjust_cat2_1 <= 0.05, na.rm = TRUE))
##p<= 0.1
sapply(immune_cell_infiltrate_annova_per_cell_type_list, function(x) sum(x$BH_adjust_cat2_1 <= 0.1, na.rm = TRUE))
##p<= 0.1 n>20

##Make a function to select significant gene pairs using sapply
input_table = immune_cell_infiltrate_annova_per_cell_type_list[[1]]
p_value = 0.05
# filter1 = "BH_adjust_cat2_1"
# select1 = "number_cat1"
# select2 = "number_cat2"
#significant_selection<- function(input_table, p_value, select1, select2)
significant_selection<- function(input_table, p_value){
  
  output_table<- input_table %>%
    #dplyr::filter(as.name(filter1) <= p_value) %>%
    dplyr::filter(BH_adjust_cat2_1 <= p_value) %>%
    dplyr::filter(number_cat1 >= 20 & number_cat2 >= 20) %>%
    dplyr::arrange(BH_adjust_cat2_1)
}

immune_cell_infiltrate_annova_per_cell_type_significant_list<- lapply(immune_cell_infiltrate_annova_per_cell_type_list, function(x) significant_selection(x, 0.05))
sapply(immune_cell_infiltrate_annova_per_cell_type_significant_list, function(x) nrow(x))
sapply(immune_cell_infiltrate_annova_per_cell_type_list, function(x) unique(x$cell_type))







#########################
### Genes identified in previous analysis
##"CD8 T cell" = 19 CDKN2A
colnames(immune_cell_infiltrate_annova_per_cell_type_list[[19]])
immune_cell_infiltrate_annova_per_cell_type_list[[19]] %>%
  dplyr::filter(target_gene == "CDKN2A") %>%
  dplyr::filter(p_value_cat2_1 <= 0.05) %>%
  dplyr::filter(number_cat1 >= 20 & number_cat2 >= 20)

## Create a function to use with lapply:
significant_gene_selection<- function(x, gene){
x %>%
  dplyr::filter(target_gene == gene) %>%
  dplyr::filter(p_value_cat2_1 <= 0.05) %>%
  dplyr::filter(number_cat1 >= 20 & number_cat2 >= 20)
}

CDKN2A_significant_immune_infiltration<- lapply(immune_cell_infiltrate_annova_per_cell_type_list, function(x)significant_gene_selection(x, "CDKN2A"))
names(CDKN2A_significant_immune_infiltration)<- immune_cell_types
sapply(CDKN2A_significant_immune_infiltration, function(x) nrow(x))

CDKN2A_significant_immune_infiltration$`Activated Natural killer cell`
CDKN2A_significant_immune_infiltration$`CD8 T cell`







