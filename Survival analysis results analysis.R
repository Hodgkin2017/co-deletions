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
survival_stats_overall_survival_cancer_list_cat1_and_2<- readRDS(file = "./R workspaces/survival_stats_overall_survival_cancer_list_cat1_and_2_169genes")
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



###################################################################

###################
### Disease free survival analysis between deletion catagories 1 and 2
####################

###############
###Import data
survival_stats_DiseFreeSurv_cancer_list_cat1_and_2<- readRDS(file = "./R workspaces/survival_stats_disease_free_cancer_list_cat1_and_2_169genes")
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

survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20 %>%
  dplyr::filter(cancer == "GBM") %>%
  dplyr::select(target_gene) %>%
  unique() 

survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20 %>%
  dplyr::filter(cancer == "GBM")

survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20 %>%
  dplyr::filter(cancer == "ALL")

##Sort by BH p-value and cancer type?

write.csv(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20, file = "survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20.csv", quote = FALSE)

##Draw venn diagrams to compare gene pairs identified on overall 
#survival and disease free survival


