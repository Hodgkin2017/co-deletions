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


#########here!!


## Identify the number of significantly co-deleted genes NOT BH adjusted with p-value < 0.05
sum(survival_stats_overall_survival_cancer_cat1_and_2_table$`p-value_logrank_test` <= 0.05, na.rm = TRUE)
sum(survival_stats_overall_survival_cancer_cat1_and_2_all_table$`p-value_logrank_test` <= 0.05, na.rm = TRUE)

## Identify the number of significantly co-deleted genes WITH BH adjusted with p-value < 0.05
sum(survival_stats_overall_survival_cancer_cat1_and_2_table$BH_adjust_logrank <= 0.05, na.rm = TRUE)
sum(survival_stats_overall_survival_cancer_cat1_and_2_all_table$BH_adjust_logrank <= 0.05, na.rm = TRUE)

survival_stats_ovsurv_cat1_2_significant_table_p0_05<- survival_stats_overall_survival_cancer_cat1_and_2_table %>%
  dplyr::filter(BH_adjust_logrank <= 0.05)
# survival_stats_ovsurv_cat1_2_significant_table_p0_05<- survival_stats_overall_survival_cancer_cat1_and_2_table$BH_adjust_logrank <= 0.05
# survival_stats_ovsurv_cat1_2_significant_table_p0_05<- survival_stats_overall_survival_cancer_cat1_and_2_table[survival_stats_ovsurv_cat1_2_significant_table_p0_05,]
dim(survival_stats_overall_survival_cancer_cat1_and_2_table)
dim(survival_stats_ovsurv_cat1_2_significant_table_p0_05)
head(survival_stats_ovsurv_cat1_2_significant_table_p0_05)

##number of significance genes per cancer type:
table(survival_stats_ovsurv_cat1_2_significant_table_p0_05$cancer)

## Identify the number of significantly co-deleted genes with p-value = 0.1
sum(survival_stats_overall_survival_cancer_cat1_and_2_table$BH_adjust_logrank <= 0.1, na.rm = TRUE)
sum(survival_stats_overall_survival_cancer_cat1_and_2_all_table$BH_adjust_logrank <= 0.1, na.rm = TRUE)

survival_stats_ovsurv_cat1_2_significant_table_p0_1<- survival_stats_overall_survival_cancer_cat1_and_2_table %>%
  dplyr::filter(BH_adjust_logrank <= 0.1)
# survival_stats_ovsurv_cat1_2_significant_table_p0_05<- survival_stats_overall_survival_cancer_cat1_and_2_table$BH_adjust_logrank <= 0.05
# survival_stats_ovsurv_cat1_2_significant_table_p0_05<- survival_stats_overall_survival_cancer_cat1_and_2_table[survival_stats_ovsurv_cat1_2_significant_table_p0_05,]
dim(survival_stats_overall_survival_cancer_cat1_and_2_table)
dim(survival_stats_ovsurv_cat1_2_significant_table_p0_1)
head(survival_stats_ovsurv_cat1_2_significant_table_p0_1)

##number of significance genes per cancer type:
table(survival_stats_ovsurv_cat1_2_significant_table_p0_1$cancer)

####################
### Identify significant genes with greater than 20 values in cat 1 and 2.

########
##pvalue <= 0.05
survival_stats_ovsurv_cat1_2_significant_table_p0_05[1:10,]
dim(survival_stats_ovsurv_cat1_2_significant_table_p0_05)

survival_stats_ovsurv_cat1_2_significant_table_p0_05_more_than_20<- 
  survival_stats_ovsurv_cat1_2_significant_table_p0_05 %>% 
  dplyr::filter(number_of_samples_cat1 >= 20 & number_of_samples_cat2 >= 20)

dim(survival_stats_ovsurv_cat1_2_significant_table_p0_05_more_than_20)

##number of significance genes per cancer type:
table(survival_stats_ovsurv_cat1_2_significant_table_p0_05_more_than_20$cancer)

survival_stats_ovsurv_cat1_2_significant_table_p0_05_more_than_20 %>%
  dplyr::group_by(cancer) %>%
  dplyr::summarise(total = n())


survival_stats_ovsurv_cat1_2_significant_table_p0_05_more_than_20 %>%
  dplyr::filter(cancer == "ALL") %>%
  dplyr::select(target_gene) %>%
  unique()

survival_stats_ovsurv_cat1_2_significant_table_p0_05_more_than_20 %>%
  dplyr::filter(cancer == "LUAD") %>%
  dplyr::select(target_gene) %>%
  unique() 

survival_stats_ovsurv_cat1_2_significant_table_p0_05_more_than_20 %>%
  dplyr::filter(cancer == "LUAD")

survival_stats_ovsurv_cat1_2_significant_table_p0_05_more_than_20 %>%
  dplyr::filter(cancer == "ALL")

##Sort by BH p-value and cancer type?

write.csv(survival_stats_ovsurv_cat1_2_significant_table_p0_05_more_than_20, file = "survival_stats_ovsurv_cat1_2_significant_table_p0_05_more_than_20.csv", quote = FALSE)

##########
### Export data tables:

## Create csv for LUAD
# survival_stats_ovsurv_cat1_2_LUAD<- survival_stats_ovsurv_cat1_2_significant_table_p0_05_more_than_20 %>%
#   dplyr::filter(cancer == "LUAD") %>%
#   dplyr::arrange(BH_adjust_logrank)
# 
# survival_stats_ovsurv_cat1_2_LUAD
# 
# write.csv(survival_stats_ovsurv_cat1_2_LUAD, file = "survival_stats_ovsurv_cat1_2_LUAD.csv", quote = FALSE)

## Create csv for ALL
survival_stats_ovsurv_cat1_2_ALL<- survival_stats_ovsurv_cat1_2_significant_table_p0_05_more_than_20 %>%
  dplyr::filter(cancer == "ALL") %>%
  dplyr::arrange(BH_adjust_logrank)

survival_stats_ovsurv_cat1_2_ALL

write.csv(survival_stats_ovsurv_cat1_2_ALL, file = "survival_stats_ovsurv_cat1_2_ALL.csv", quote = FALSE)

##########################

########
##pvalue <= 0.1
survival_stats_ovsurv_cat1_2_significant_table_p0_1[1:10,]
dim(survival_stats_ovsurv_cat1_2_significant_table_p0_1)

survival_stats_ovsurv_cat1_2_significant_table_p0_1_more_than_20<- 
  survival_stats_ovsurv_cat1_2_significant_table_p0_1 %>% 
  dplyr::filter(number_of_samples_cat1 >= 20 & number_of_samples_cat2 >= 20)

dim(survival_stats_ovsurv_cat1_2_significant_table_p0_1_more_than_20)

##number of significance genes per cancer type:
table(survival_stats_ovsurv_cat1_2_significant_table_p0_1_more_than_20$cancer)

survival_stats_ovsurv_cat1_2_significant_table_p0_1_more_than_20 %>%
  dplyr::group_by(cancer) %>%
  dplyr::summarise(total = n())


survival_stats_ovsurv_cat1_2_significant_table_p0_1_more_than_20 %>%
  dplyr::filter(cancer == "ALL") %>%
  dplyr::select(target_gene) %>%
  unique()

survival_stats_ovsurv_cat1_2_significant_table_p0_1_more_than_20 %>%
  dplyr::filter(cancer == "LUAD") %>%
  dplyr::select(target_gene) %>%
  unique() 

survival_stats_ovsurv_cat1_2_significant_table_p0_1_more_than_20 %>%
  dplyr::filter(cancer == "LUAD")

survival_stats_ovsurv_cat1_2_significant_table_p0_1_more_than_20 %>%
  dplyr::filter(cancer == "ALL")

##Sort by BH p-value and cancer type?

write.csv(survival_stats_ovsurv_cat1_2_significant_table_p0_1_more_than_20, file = "survival_stats_ovsurv_cat1_2_significant_table_p0_1_more_than_20.csv", quote = FALSE)

##########
### Export data tables:

## Create csv for LUAD
# survival_stats_ovsurv_cat1_2_LUAD<- survival_stats_ovsurv_cat1_2_significant_table_p0_05_more_than_20 %>%
#   dplyr::filter(cancer == "LUAD") %>%
#   dplyr::arrange(BH_adjust_logrank)
# 
# survival_stats_ovsurv_cat1_2_LUAD
# 
# write.csv(survival_stats_ovsurv_cat1_2_LUAD, file = "survival_stats_ovsurv_cat1_2_LUAD.csv", quote = FALSE)

## Create csv for ALL
survival_stats_ovsurv_cat1_2_ALL<- survival_stats_ovsurv_cat1_2_significant_table_p0_1_more_than_20 %>%
  dplyr::filter(cancer == "ALL") %>%
  dplyr::arrange(BH_adjust_logrank)

survival_stats_ovsurv_cat1_2_ALL

write.csv(survival_stats_ovsurv_cat1_2_ALL, file = "survival_stats_ovsurv_cat1_2_ALL_p=0_1.csv", quote = FALSE)





















#########################
### Genes identified in previous analysis
















