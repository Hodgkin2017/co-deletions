####################
### Survival results analysis
####################

###################
### Maria's and Christophe's gene list: Total 169 genes
#################

#########################
### BRCA survival results analysis

##Import gene survival analysis results:

# BRCA_co_deletion_overall_survival_stats_long<- read.csv("../../Output/Survival analysis/tables/Maria's-Christophe gene list/BRCA_co-deletion_overall_survival_stats_long.csv", header = TRUE, stringsAsFactors = FALSE)
# dim(BRCA_co_deletion_overall_survival_stats_long)
# 
# BRCA_co_deletion_overall_survival_stats_long_significant<- BRCA_co_deletion_overall_survival_stats_long %>% dplyr::filter(p.value_logrank_test <= 0.05)
# dim(BRCA_co_deletion_overall_survival_stats_long_significant)
# 
# 
# BRCA_co_deletion_overall_survival_stats_long_significant_n10<- BRCA_co_deletion_overall_survival_stats_long_significant %>% dplyr::filter(number_of_samples_cat1 >= 10) %>%
#   dplyr::filter(number_of_samples_cat2 >= 10)
# 
# dim(BRCA_co_deletion_overall_survival_stats_long_significant_n10)
# 
# write.csv(BRCA_co_deletion_overall_survival_stats_long_significant_n10, file = "BRCA_co_deletion_overall_survival_stats_long_significant_n10.csv")


###################
### Overall survival analysis between deletion catagories 1 and 2
####################

###############
###Import data
##old data rounded to 2.d.p
# survival_stats_overall_survival_cancer_list_cat1_and_2<- readRDS(file = "./R workspaces/survival_stats_overall_survival_cancer_list_cat1_and_2_169genes")
# survival_stats_overall_survival_cancer_list_cat1_and_2[[1]]


survival_stats_overall_survival_cancer_list_cat1_and_2<- readRDS(file = "./R workspaces/survival_stats_overall_survival_cancer_list_cat1_and_2_169genes_2")
survival_stats_overall_survival_cancer_list_cat1_and_2[[1]]

length(survival_stats_overall_survival_cancer_list_cat1_and_2)
names(survival_stats_overall_survival_cancer_list_cat1_and_2)<- names(threshold_selected_cnv_list_plus_all_loc)
names(survival_stats_overall_survival_cancer_list_cat1_and_2)
colnames(survival_stats_overall_survival_cancer_list_cat1_and_2[[1]])


###############
###Perform BH multiple testing correction

##Add a column to each table in list with the cancer name

dim(survival_stats_overall_survival_cancer_list_cat1_and_2[[1]])
survival_stats_overall_survival_cancer_list_cat1_and_2[[1]][1:10,]
names(survival_stats_overall_survival_cancer_list_cat1_and_2[1])

for (i in 1: length(survival_stats_overall_survival_cancer_list_cat1_and_2)){
  
  column_name<- names(survival_stats_overall_survival_cancer_list_cat1_and_2[i])
  column_name<- rep(column_name, nrow(survival_stats_overall_survival_cancer_list_cat1_and_2[[i]]))
  survival_stats_overall_survival_cancer_list_cat1_and_2[[i]]<-cbind(cancer = column_name, survival_stats_overall_survival_cancer_list_cat1_and_2[[i]])
  
}
head(survival_stats_overall_survival_cancer_list_cat1_and_2[[1]])
head(survival_stats_overall_survival_cancer_list_cat1_and_2[[33]])

## Join tables together except for 'all'
survival_stats_overall_survival_cancer_cat1_and_2_table<- do.call(rbind, survival_stats_overall_survival_cancer_list_cat1_and_2[1:32])
dim(survival_stats_overall_survival_cancer_list_cat1_and_2)
dim(survival_stats_overall_survival_cancer_cat1_and_2_table)

##BH multiple testing for cancers 1:32 in list
survival_stats_overall_survival_cancer_cat1_and_2_table$BH_adjust_logrank<- p.adjust(survival_stats_overall_survival_cancer_cat1_and_2_table$`p-value_logrank_test`, method = "BH")
head(survival_stats_overall_survival_cancer_cat1_and_2_table)

##BH multiple testing for all cancers (item 33 in list)
survival_stats_overall_survival_cancer_cat1_and_2_all_table<- survival_stats_overall_survival_cancer_list_cat1_and_2[[33]]
dim(survival_stats_overall_survival_cancer_cat1_and_2_all_table)
survival_stats_overall_survival_cancer_cat1_and_2_all_table$BH_adjust_logrank<- p.adjust(survival_stats_overall_survival_cancer_cat1_and_2_all_table$`p-value_logrank_test`, method = "BH")
head(survival_stats_overall_survival_cancer_cat1_and_2_all_table)

##Join BH tables together:
survival_stats_overall_survival_cancer_cat1_and_2_table<- rbind(survival_stats_overall_survival_cancer_cat1_and_2_table, survival_stats_overall_survival_cancer_cat1_and_2_all_table)
dim(survival_stats_overall_survival_cancer_cat1_and_2_table)


#####################
###Identify significant genes

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



###################################################################

###################
### Disease free survival analysis between deletion catagories 1 and 2
####################

###############
###Import data
##old data rounded to 2.d.p
# survival_stats_DiseFreeSurv_cancer_list_cat1_and_2<- readRDS(file = "./R workspaces/survival_stats_disease_free_cancer_list_cat1_and_2_169genes")
# survival_stats_DiseFreeSurv_cancer_list_cat1_and_2[[1]]

survival_stats_DiseFreeSurv_cancer_list_cat1_and_2<- readRDS(file = "./R workspaces/survival_stats_disease_free_cancer_list_cat1_and_2_169genes_2")
survival_stats_DiseFreeSurv_cancer_list_cat1_and_2[[1]]

length(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2)
names(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2)<- names(threshold_selected_cnv_list_plus_all_loc)
names(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2)
colnames(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2[[1]])

#############
### Remove tables for cancer types without disease free survival data:
#UCS = 31, PCPG = 22, MESO = 19, ACC = 1
survival_stats_DiseFreeSurv_cancer_list_cat1_and_2<- survival_stats_DiseFreeSurv_cancer_list_cat1_and_2[c(2:18,20, 21, 23:30, 32, 33)]
names(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2)

###############
###Perform BH multiple testing correction

##Add a column to each table in list with the cancer name

dim(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2[[1]])
survival_stats_DiseFreeSurv_cancer_list_cat1_and_2[[1]][1:10,]
names(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2[1])

for (i in 1: length(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2)){
  
  column_name<- names(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2[i])
  column_name<- rep(column_name, nrow(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2[[i]]))
  survival_stats_DiseFreeSurv_cancer_list_cat1_and_2[[i]]<-cbind(cancer = column_name, survival_stats_DiseFreeSurv_cancer_list_cat1_and_2[[i]])
  
}
head(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2[[1]])
head(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2[[29]])

## Join tables together except for 'all'
survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table<- do.call(rbind, survival_stats_DiseFreeSurv_cancer_list_cat1_and_2[1:32])
dim(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2)
dim(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table)

##BH multiple testing for cancers 1:32 in list
survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table$BH_adjust_logrank<- p.adjust(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table$`p-value_logrank_test`, method = "BH")
head(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table)

##BH multiple testing for all cancers table (item 29 in list)
survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_all_table<- survival_stats_DiseFreeSurv_cancer_list_cat1_and_2[[29]]
dim(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_all_table)
survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_all_table$BH_adjust_logrank<- p.adjust(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_all_table$`p-value_logrank_test`, method = "BH")
head(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_all_table)

##Join BH tables together:
survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table<- rbind(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table, survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_all_table)
dim(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table)


#####################
###Identify significant genes

## Identify the number of significantly co-deleted genes NOT BH adjusted with p-value < 0.05
sum(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table$`p-value_logrank_test` <= 0.05, na.rm = TRUE)
sum(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_all_table$`p-value_logrank_test` <= 0.05, na.rm = TRUE)

## Identify the number of significantly co-deleted genes WITH BH adjusted with p-value < 0.05
sum(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table$BH_adjust_logrank <= 0.05, na.rm = TRUE)
sum(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_all_table$BH_adjust_logrank <= 0.05, na.rm = TRUE)

survival_stats_DiseFreeSurv_cancer_list_cat1_2_table_p0_05<- survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table %>%
  dplyr::filter(BH_adjust_logrank <= 0.05)
# survival_stats_ovsurv_cat1_2_significant_table_p0_05<- survival_stats_overall_survival_cancer_cat1_and_2_table$BH_adjust_logrank <= 0.05
# survival_stats_ovsurv_cat1_2_significant_table_p0_05<- survival_stats_overall_survival_cancer_cat1_and_2_table[survival_stats_ovsurv_cat1_2_significant_table_p0_05,]
dim(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table)
dim(survival_stats_DiseFreeSurv_cancer_list_cat1_2_table_p0_05)
head(survival_stats_DiseFreeSurv_cancer_list_cat1_2_table_p0_05)

##number of significance genes per cancer type:
table(survival_stats_DiseFreeSurv_cancer_list_cat1_2_table_p0_05$cancer)

## Identify the number of significantly co-deleted genes with p-value = 0.1
sum(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table$BH_adjust_logrank <= 0.1, na.rm = TRUE)
sum(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_all_table$BH_adjust_logrank <= 0.1, na.rm = TRUE)

survival_stats_DiseFreeSurv_cancer_list_cat1_2_table_p0_1<- survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table %>%
  dplyr::filter(BH_adjust_logrank <= 0.1)
# survival_stats_ovsurv_cat1_2_significant_table_p0_05<- survival_stats_overall_survival_cancer_cat1_and_2_table$BH_adjust_logrank <= 0.05
# survival_stats_ovsurv_cat1_2_significant_table_p0_05<- survival_stats_overall_survival_cancer_cat1_and_2_table[survival_stats_ovsurv_cat1_2_significant_table_p0_05,]
dim(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table)
dim(survival_stats_DiseFreeSurv_cancer_list_cat1_2_table_p0_1)
head(survival_stats_DiseFreeSurv_cancer_list_cat1_2_table_p0_1)


####################
### Identify significant genes with greater than 20 values in cat 1 and 2.

########
##pvalue <= 0.05
survival_stats_DiseFreeSurv_cancer_list_cat1_2_table_p0_05[1:10,]
dim(survival_stats_DiseFreeSurv_cancer_list_cat1_2_table_p0_05)

survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20<- 
  survival_stats_DiseFreeSurv_cancer_list_cat1_2_table_p0_05 %>% 
  dplyr::filter(number_of_samples_cat1 >= 20 & number_of_samples_cat2 >= 20)

dim(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20)

##number of significance genes per cancer type:
table(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20$cancer)

survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20 %>%
  dplyr::group_by(cancer) %>%
  dplyr::summarise(total = n())


survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20 %>%
  dplyr::filter(cancer == "ALL") %>%
  dplyr::select(target_gene) %>%
  unique()

# survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20 %>%
#   dplyr::filter(cancer == "GBM") %>%
#   dplyr::select(target_gene) %>%
#   unique() 
# 
# survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20 %>%
#   dplyr::filter(cancer == "GBM")

survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20 %>%
  dplyr::filter(cancer == "ALL")

##Sort by BH p-value and cancer type?

write.csv(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20, file = "survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20.csv", quote = FALSE)

##########
### Export data tables:

## Create csv for GBM
# survival_stats_disefreeSurv_cat1_2_GBM<- survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20 %>%
#   dplyr::filter(cancer == "GBM") %>%
#   dplyr::arrange(BH_adjust_logrank)
# 
# survival_stats_disefreeSurv_cat1_2_GBM
# 
# write.csv(survival_stats_disefreeSurv_cat1_2_GBM, file = "survival_stats_disefreeSurv_cat1_2_GBM.csv", quote = FALSE)

## Create csv for ALL
survival_stats_disefreeSurv_cat1_2_ALL<- survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20 %>%
  dplyr::filter(cancer == "ALL") %>%
  dplyr::arrange(BH_adjust_logrank)

survival_stats_disefreeSurv_cat1_2_ALL

write.csv(survival_stats_disefreeSurv_cat1_2_ALL, file = "survival_stats_disefreeSurv_cat1_2_ALL.csv", quote = FALSE)


###################
########
##pvalue <= 0.1
survival_stats_DiseFreeSurv_cancer_list_cat1_2_table_p0_1[1:10,]
dim(survival_stats_DiseFreeSurv_cancer_list_cat1_2_table_p0_1)

survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20<- 
  survival_stats_DiseFreeSurv_cancer_list_cat1_2_table_p0_1 %>% 
  dplyr::filter(number_of_samples_cat1 >= 20 & number_of_samples_cat2 >= 20)

dim(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20)

##number of significance genes per cancer type:
table(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20$cancer)

survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20 %>%
  dplyr::group_by(cancer) %>%
  dplyr::summarise(total = n())


survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20 %>%
  dplyr::filter(cancer == "ALL") %>%
  dplyr::select(target_gene) %>%
  unique()

# survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20 %>%
#   dplyr::filter(cancer == "GBM") %>%
#   dplyr::select(target_gene) %>%
#   unique() 
# 
# survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20 %>%
#   dplyr::filter(cancer == "GBM")

survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20 %>%
  dplyr::filter(cancer == "ALL")

##Sort by BH p-value and cancer type?

write.csv(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20, file = "survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20_p=0_1.csv", quote = FALSE)

##########
### Export data tables:

## Create csv for GBM
# survival_stats_disefreeSurv_cat1_2_GBM<- survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20 %>%
#   dplyr::filter(cancer == "GBM") %>%
#   dplyr::arrange(BH_adjust_logrank)
# 
# survival_stats_disefreeSurv_cat1_2_GBM
# 
# write.csv(survival_stats_disefreeSurv_cat1_2_GBM, file = "survival_stats_disefreeSurv_cat1_2_GBM.csv", quote = FALSE)

## Create csv for ALL
survival_stats_disefreeSurv_cat1_2_ALL_0_1<- survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20 %>%
  dplyr::filter(cancer == "ALL") %>%
  dplyr::arrange(BH_adjust_logrank)

survival_stats_disefreeSurv_cat1_2_ALL

write.csv(survival_stats_disefreeSurv_cat1_2_ALL_0_1, file = "survival_stats_disefreeSurv_cat1_2_ALL_0_1.csv", quote = FALSE)





##Draw venn diagrams to compare gene pairs identified on overall 
#survival and disease free survival?


########################################################

#####################
###Stratify survival curves for significance hits
####################

##################
### Overall survival
###################

#################
### What genes do I need to stratify?
unique(survival_stats_ovsurv_cat1_2_LUAD$target_gene)
unique(survival_stats_ovsurv_cat1_2_LUAD$proximal_gene)

unique(survival_stats_ovsurv_cat1_2_ALL$target_gene)
unique(survival_stats_ovsurv_cat1_2_ALL$proximal_gene)

################
### ADD Tumour stage to survival_time_list

names(survival_time_list)






#################
### Stratify LUAD: CDKN2A

##############
### Create a for loop to go through all cancers and perform survival analysis on target genes and co-deletions:


# ##Create an empty list
# survival_stats_cancer_list_LUAD<- vector("list", 1)
#   
#   ##Get CNV table
#   cnv.table<- threshold_selected_cnv_list_plus_all_loc$LUAD
#   
#   
#   ##Get survival data
#   survival_time_list<- clinical_long_plus_all_survival$LUAD_clinical.tsv
#   
#   # target_gene_names<- sapply(gene_information_long_list, function(x) x[[1]])
#   # which(target_gene_names %in%  "CDKN2A" )
#   
#   ## Perform survival analysis of co-deletions:
#   survival_stats<- survival_analysis_of_gene_list_cat1_and_2(target_gene_list = gene_information_long_list[[16]], survival_time_list = survival_time_list, cnv.table = cnv.table, distance = 2.5e+06,
#                                                                                                                  threshold = -2, deletion = TRUE, time_of_death_column = 5, 
#                                                                                                                  death_event_column = 6, column_start = 11, start = TRUE, 
#                                                                                                                  remove_NA = TRUE, Cytoband = FALSE, print_to_screen = FALSE, 
#                                                                                                                  plot_graph = FALSE)
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   ##Combine survival stats together and order by log-rank test p-value
#   survival_stats_list_table<- do.call(rbind, survival_stats_list)
#   survival_stats_list_table<- survival_stats_list_table[order(survival_stats_list_table$`p-value_logrank_test`, decreasing = FALSE),]
#   #survival_stats_list_table[1:100, c(1,2,5,6,8,9)]
#   
#   
#   
#   
#   
#   
#   
#   ##Save file as .csv
#   write.csv(survival_stats_list_table, file = paste0(names(threshold_selected_cnv_list_plus_all_loc[i]),"_co-deletion_overall_survival_stats_cat1_and_2.csv"), quote = FALSE)
#   
# 
#   
# 
# 
# 
# ## Save object:
# saveRDS(survival_stats_cancer_list3, file = "./R workspaces/survival_stats_overall_survival_cancer_list_cat1_and_2")
# 


##################
### Draw survival curves
###################

##############
### Overall survival
#############

###Parameters:
##Get CNV table
cnv.table<- threshold_selected_cnv_list_plus_all_loc$LUAD
##Get survival data
survival_time_list<- clinical_long_plus_all_survival$LUAD_clinical.tsv
## Get Target gene
target_gene_names<- sapply(gene_information_long_list, function(x) x[[1]])
which(target_gene_names %in%  "CDKN2A" )
target_gene_list<- gene_information_long_list[[16]]

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
ylabel = "Overall survival"
path = "./"
path = ("../../Output/Survival analysis/survival_curves")








#########
###F2: Function to Create a table of survival p-value statistics for co-deletions or co-amplifications

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
                                                     plot_graph = FALSE,
                                                     ylabel = "Overall survival",
                                                     path = "./"){
  
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
    
    ################here!!!!!
    
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
      if(length(unique(covariable_object)) > 1){
        df.categ <- cbind(sapply(names(fittedSurv$strata), function(x) strsplit(x,"=")[[1]][2]),
                          fittedSurv$n)
        #df.categ <- data.frame(df.categ)
      }else {
        df.categ <- cbind(unique(covariable_object), fittedSurv$n)
        
      }
      temp_categ<- df.categ[,1]
      temp_categ<- sub("\\<1\\>", paste(target_gene, "deletion", sep = " "), temp_categ)
      temp_categ<- sub("\\<2\\>", paste(target_gene, "and", colnames(clinical_survival_deletion_category)[i+2], "co-deletion", sep = " "), temp_categ)
      df.categ<- cbind(temp_categ, df.categ[,2])           
      
      ## Change name to 1 = Tumour suppressor deletion
      #2 = Co-deletion
      
      
      ##Get category names for one variable
      categNames <- apply(df.categ, 1, function(x) paste0(x[1]," (",x[2],")"))
      
      if (plot_graph == TRUE) {
        tiff(paste(target_gene, "_",  colnames(clinical_survival_deletion_category)[i+2],".tiff", sep =""), width = 7, height = 7, units = 'in', res = 100)
        plot(fittedSurv, main=paste(target_gene, "_",  colnames(clinical_survival_deletion_category)[i+2], " Survival", sep =""),
             xlab="Time (days)", ylab=ylabel,
             col=brewer.pal(9,"Set1"), 
             mark.time=T)
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



###########################################
### Test function:

setwd("../../../Code/co-deletions/")

##Get CNV table
cnv.table<- threshold_selected_cnv_list_plus_all_loc$LUAD
##Get survival data
survival_time_list<- clinical_long_plus_all_survival$LUAD_clinical.tsv
## Get Target gene
target_gene_names<- sapply(gene_information_long_list, function(x) x[[1]])
which(target_gene_names %in%  "CDKN2A" )
target_gene_list<- gene_information_long_list[[16]]

##O3:
test_plot1<- survival_analysis_of_gene_list_cat1_and_2(target_gene_list = target_gene_list, survival_time_list = survival_time_list, cnv.table = cnv.table, distance = 2.5e+06,
                                                       threshold = -2, deletion = TRUE, time_of_death_column = 5, 
                                                       death_event_column = 6, column_start = 11, start = TRUE, 
                                                       remove_NA = TRUE, Cytoband = FALSE, print_to_screen = FALSE, 
                                                       plot_graph = TRUE, ylabel = "Overall survival",
                                                       path = ("../../Output/Survival analysis/survival_curves"))


####################


##Get CNV table
cnv.table<- threshold_selected_cnv_list_plus_all_loc$ALL
##Get survival data
survival_time_list<- clinical_long_plus_all_survival$ALL
## Get Target gene
target_gene_names<- sapply(gene_information_long_list, function(x) x[[1]])
selected_target_genes<- which(target_gene_names %in%  c("NF1", "RB1", "TP53", "CAMTA1", "CDKN1B", "LRP1B", "TGFBR2", "ZFHX3") )
target_gene_list<- gene_information_long_list[selected_target_genes]
length(target_gene_list)

##O3:
test_plot2<- lapply(target_gene_list, function (x) survival_analysis_of_gene_list_cat1_and_2(target_gene_list = x, survival_time_list = survival_time_list, cnv.table = cnv.table, distance = 2.5e+06,
                                                                                             threshold = -2, deletion = TRUE, time_of_death_column = 5, 
                                                                                             death_event_column = 6, column_start = 11, start = TRUE, 
                                                                                             remove_NA = TRUE, Cytoband = FALSE, print_to_screen = FALSE, 
                                                                                             plot_graph = TRUE, ylabel = "Overall survival",
                                                                                             path = ("../../Output/Survival analysis/survival_curves/ALL")
)
)



setwd("../../../../Code/co-deletions/")


##############
### Disease Free survival
#############

##################
### GBM significant gene


##Get CNV table
cnv.table<- threshold_selected_cnv_list_plus_all_loc$GBM
##Get survival data
survival_time_list<- clinical_long_plus_all_survival$GBM_clinical.tsv
## Get Target gene
target_gene_names<- sapply(gene_information_long_list, function(x) x[[1]])
which(target_gene_names %in%  "CDKN2A" )
target_gene_list<- gene_information_long_list[[16]]

##O3:
test_plot1<- survival_analysis_of_gene_list_cat1_and_2(target_gene_list = target_gene_list, survival_time_list = survival_time_list, cnv.table = cnv.table, distance = 2.5e+06,
                                                       threshold = -2, deletion = TRUE, time_of_death_column = 5, 
                                                       death_event_column = 6, column_start = 11, start = TRUE, 
                                                       remove_NA = TRUE, Cytoband = FALSE, print_to_screen = FALSE, 
                                                       plot_graph = TRUE, ylabel = "Disease Free Survival",
                                                       path = ("../../Output/Survival analysis/survival_curves/Disease_free/GBM"))


##################
### pan-cancer significant genes


##Get CNV table
cnv.table<- threshold_selected_cnv_list_plus_all_loc$ALL
##Get survival data
survival_time_list<- clinical_long_plus_all_survival$ALL
## Get Target gene
target_gene_names<- sapply(gene_information_long_list, function(x) x[[1]])
selected_target_genes<- which(target_gene_names %in%  c("RB1", "TP53", "NCOR1","CAMTA1", "CDKN1B", "LRP1B", "TGFBR2", "ZFHX3") )
target_gene_list<- gene_information_long_list[selected_target_genes]
length(target_gene_list)
target_gene_list


##O3:
test_plot2<- lapply(target_gene_list, function (x) survival_analysis_of_gene_list_cat1_and_2(target_gene_list = x, survival_time_list = survival_time_list, cnv.table = cnv.table, distance = 2.5e+06,
                                                                                             threshold = -2, deletion = TRUE, time_of_death_column = 5, 
                                                                                             death_event_column = 6, column_start = 11, start = TRUE, 
                                                                                             remove_NA = TRUE, Cytoband = FALSE, print_to_screen = FALSE, 
                                                                                             plot_graph = TRUE, ylabel = "Disease Free Survival",
                                                                                             path = ("../../Output/Survival analysis/survival_curves/Disease_free/ALL")
)
)



setwd("../../../../Code/co-deletions/")

#######################
##Venn diagram to show intersection of genes between overall survival and disease free survival






