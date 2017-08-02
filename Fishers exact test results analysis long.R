##################
###Fishers exact test results analysis
###################

####################
### Fishers exact test per gene between cat 1,2,3 & 4 for short 23 gene list
#####################

# fishers_co_deletion_test_results_list<- readRDS(file = "./R workspaces/fishers_co_deletion_test_results_list")
# fishers_co_deletion_test_results_list
# 
# length(fishers_co_deletion_test_results_list)
# names(fishers_co_deletion_test_results_list)<- names(threshold_selected_cnv_list_plus_all_loc)
# names(fishers_co_deletion_test_results_list)
# 
# fishers_co_deletion_test_results_list[[1]]$`co-deletion_p-value`<= 0.05
# 
# significant_hits<- lapply(fishers_co_deletion_test_results_list, function(x) x$`co-deletion_p-value`<= 0.05)
# 
# lapply(fishers_co_deletion_test_results_list, function(x) x$`co-deletion_p-value`<= 0.05)
# 
# ###Loop to create new list with co-deletions with p-value <= 0.5
# fishers_co_deletion_test_results_list_significant<- vector("list", length(fishers_co_deletion_test_results_list))
# 
# names(fishers_co_deletion_test_results_list)<- names(threshold_selected_cnv_list_plus_all_loc)
# 
# for (i in 1: length(fishers_co_deletion_test_results_list)) {
#   
#   df<- fishers_co_deletion_test_results_list[[i]]
#   
#   fishers_co_deletion_test_results_list_significant[[i]]<- df %>% dplyr::filter(`co-deletion_p-value`<= 0.05)
#   
# }
# 
# names(fishers_co_deletion_test_results_list_significant)<- names(threshold_selected_cnv_list_plus_all_loc)
# 
# sapply(fishers_co_deletion_test_results_list, function(x) nrow(x))
# sapply(fishers_co_deletion_test_results_list_significant, function(x) nrow(x))
# 
# BRCA_fishers_exact_test<- fishers_co_deletion_test_results_list_significant[[3]]
# BRCA_fishers_exact_test<- BRCA_fishers_exact_test[order(BRCA_fishers_exact_test$`co-deletion_p-value`),]
# head(BRCA_fishers_exact_test, 40)
# BRCA_fishers_exact_test[41:80,]

####################
### Fishers exact test per gene between cat 1,2,3 & 4 for long 169 gene list
#####################

###############
###Import data
fishers_co_deletion_per_gene_list<- readRDS(file = "./R workspaces/fishers_co_deletion_test_results_list_169genes")
fishers_co_deletion_per_gene_list[[1]]

length(fishers_co_deletion_per_gene_list)
names(fishers_co_deletion_per_gene_list)<- names(threshold_selected_cnv_list_plus_all_loc)
names(fishers_co_deletion_per_gene_list)
colnames(fishers_co_deletion_per_gene_list[[1]])

#############
### Add number of cat 1, 2, 3, 4 values to tables.

survival_stats_overall_survival_cancer_list_169genes<- readRDS(file = "./R workspaces/survival_stats_overall_survival_cancer_list_169genes")
length(survival_stats_overall_survival_cancer_list_169genes)
survival_stats_overall_survival_cancer_list_169genes[[1]][1:10,1:3]
fishers_co_deletion_per_gene_list[[1]][1:10,]
survival_stats_overall_survival_cancer_list_169genes[[33]][1:10,1:3]
fishers_co_deletion_per_gene_list[[33]][1:10,]

## for loop to add cat 1, 2, 3, 4 values to fishers exact test tables in list
for (i in 1: length(survival_stats_overall_survival_cancer_list_169genes)){
  
  number_of_samples_table<- survival_stats_overall_survival_cancer_list_169genes[[i]][,7:10]
  fishers_co_deletion_per_gene_list[[i]]<- cbind(fishers_co_deletion_per_gene_list[[i]], number_of_samples_table)
  
  
}
fishers_co_deletion_per_gene_list[[1]][1:10,]

###############
###Perform BH multiple testing correction

##Add a column to each table in list with the cancer name

dim(fishers_co_deletion_per_gene_list[[1]])
fishers_co_deletion_per_gene_list[[1]][1:10,]
names(fishers_co_deletion_per_gene_list[1])

for (i in 1: length(fishers_co_deletion_per_gene_list)){
  
  column_name<- names(fishers_co_deletion_per_gene_list[i])
  column_name<- rep(column_name, nrow(fishers_co_deletion_per_gene_list[[i]]))
  fishers_co_deletion_per_gene_list[[i]]<-cbind(cancer = column_name, fishers_co_deletion_per_gene_list[[i]])
  
}
head(fishers_co_deletion_per_gene_list[[1]])
head(fishers_co_deletion_per_gene_list[[33]])

## Join tables together except for 'all'
fishers_co_deletion_per_gene_long_table<- do.call(rbind, fishers_co_deletion_per_gene_list[1:32])
dim(fishers_co_deletion_per_gene_list)
dim(fishers_co_deletion_per_gene_long_table)

##BH multiple testing for cancers 1:32 in list
fishers_co_deletion_per_gene_long_table$BH_adjust<- p.adjust(fishers_co_deletion_per_gene_long_table$`co-deletion_p-value`, method = "BH")
head(fishers_co_deletion_per_gene_long_table)

##BH multiple testing for all cancers (item 33 in list)
fishers_co_deletion_per_gene_all_table<- fishers_co_deletion_per_gene_list[[33]]
dim(fishers_co_deletion_per_gene_all_table)
fishers_co_deletion_per_gene_all_table$BH_adjust<- p.adjust(fishers_co_deletion_per_gene_all_table$`co-deletion_p-value`, method = "BH")
head(fishers_co_deletion_per_gene_all_table)

##Join BH tables together:
fishers_co_deletion_per_gene_long_table<- rbind(fishers_co_deletion_per_gene_long_table, fishers_co_deletion_per_gene_all_table)
dim(fishers_co_deletion_per_gene_long_table)


#####################
###Identify significant genes

## Identify the number of significantly co-deleted genes with p-value = 0.05
sum(fishers_co_deletion_per_gene_long_table$BH_adjust <= 0.05)

fishers_co_deletion_per_gene_long_table_significanct_p0.05<- fishers_co_deletion_per_gene_long_table$BH_adjust <= 0.05
fishers_co_deletion_per_gene_long_table_significanct_p0.05<- fishers_co_deletion_per_gene_long_table[fishers_co_deletion_per_gene_long_table_significanct_p0.05,]
dim(fishers_co_deletion_per_gene_long_table)
dim(fishers_co_deletion_per_gene_long_table_significanct_p0.05)
head(fishers_co_deletion_per_gene_long_table_significanct_p0.05)

##number of significance genes per cancer type:
table(fishers_co_deletion_per_gene_long_table_significanct_p0.05$cancer)

## Identify the number of significantly co-deleted genes with p-value = 0.1
sum(fishers_co_deletion_per_gene_long_table$BH_adjust <= 0.1)

fishers_co_deletion_per_gene_long_table_significanct_p0.1<- fishers_co_deletion_per_gene_long_table$BH_adjust <= 0.1
fishers_co_deletion_per_gene_long_table_significanct_p0.1<- fishers_co_deletion_per_gene_long_table[fishers_co_deletion_per_gene_long_table_significanct_p0.1,]
dim(fishers_co_deletion_per_gene_long_table)
dim(fishers_co_deletion_per_gene_long_table_significanct_p0.1)
head(fishers_co_deletion_per_gene_long_table_significanct_p0.1)

##number of significance genes per cancer type:
table(fishers_co_deletion_per_gene_long_table_significanct_p0.1$cancer)

####################
### Identify significant genes with greater than 20 values in cat 1 and 2.

########
##pvalue <= 0.05
fishers_co_deletion_per_gene_long_table_significanct_p0.05[1:10,]
dim(fishers_co_deletion_per_gene_long_table_significanct_p0.05)

fishers_co_deletion_per_gene_long_table_significanct_p0.05_more_than_20<- 
  fishers_co_deletion_per_gene_long_table_significanct_p0.05 %>% 
  dplyr::filter(number_of_samples_cat1 >= 20 & number_of_samples_cat2 >= 20)

dim(fishers_co_deletion_per_gene_long_table_significanct_p0.05_more_than_20)

##number of significance genes per cancer type:
table(fishers_co_deletion_per_gene_long_table_significanct_p0.05_more_than_20$cancer)

fishers_co_deletion_per_gene_long_table_significanct_p0.05_more_than_20 %>%
  dplyr::group_by(cancer) %>%
  dplyr::summarise(total = n())

########
##pvalue <= 0.1
fishers_co_deletion_per_gene_long_table_significanct_p0.1[1:10,]
dim(fishers_co_deletion_per_gene_long_table_significanct_p0.1)

fishers_co_deletion_per_gene_long_table_significanct_p0.1_more_than_20<- 
  fishers_co_deletion_per_gene_long_table_significanct_p0.1 %>% 
  dplyr::filter(number_of_samples_cat1 >= 20 & number_of_samples_cat2 >= 20)

dim(fishers_co_deletion_per_gene_long_table_significanct_p0.1_more_than_20)

##number of significance genes per cancer type:
table(fishers_co_deletion_per_gene_long_table_significanct_p0.1_more_than_20$cancer)

fishers_co_deletion_per_gene_long_table_significanct_p0.1_more_than_20 %>%
  dplyr::group_by(cancer) %>%
  dplyr::summarise(total = n())

fishers_co_deletion_per_gene_long_table_significanct_p0.1_more_than_20 %>%
  dplyr::filter(cancer == "ALL") %>%
  dplyr::select(target_gene) %>%
  unique() %>%
  dim()

fishers_co_deletion_per_gene_long_table_significanct_p0.1_more_than_20 %>%
  dplyr::filter(cancer == "ALL") %>%
  dplyr::select(target_gene) %>%
  unique()

fishers_co_deletion_per_gene_ALL<- fishers_co_deletion_per_gene_long_table_significanct_p0.1_more_than_20 %>%
  dplyr::filter(cancer == "ALL") %>%
  dplyr::arrange(BH_adjust)

fishers_co_deletion_per_gene_ALL

#write.csv(fishers_co_deletion_per_gene_ALL, file = "fishers_co_deletion_per_gene_ALL.csv", quote = FALSE)


fishers_co_deletion_per_gene_long_table_significanct_p0.1_more_than_20 %>%
  dplyr::filter(cancer == "BRCA") %>%
  dplyr::select(target_gene) %>%
  unique() %>%
  dim()

fishers_co_deletion_per_gene_BRCA<- fishers_co_deletion_per_gene_long_table_significanct_p0.1_more_than_20 %>%
  dplyr::filter(cancer == "BRCA") %>%
  dplyr::arrange(BH_adjust)

fishers_co_deletion_per_gene_BRCA

#write.csv(fishers_co_deletion_per_gene_BRCA, file = "fishers_co_deletion_per_gene_BRCA.csv", quote = FALSE)


#write.csv(fishers_co_deletion_per_gene_long_table_significanct_p0.1_more_than_20, file = "fishers_co_deletion_per_gene_long_table_significanct_p0.1_more_than_20.csv", quote = FALSE)


########################
### Bar plot of number of significant genes:
##https://www.r-bloggers.com/make-a-bar-plot-with-ggplot/

## p<= 0.05
dim(fishers_co_deletion_per_gene_long_table_significanct_p0.05_more_than_20)
head(fishers_co_deletion_per_gene_long_table_significanct_p0.05_more_than_20)

table(fishers_co_deletion_per_gene_long_table_significanct_p0.05_more_than_20$cancer)
bar_plot_3<- fishers_co_deletion_per_gene_long_table_significanct_p0.05_more_than_20 %>%
  dplyr::group_by(cancer) %>%
  dplyr::summarise(total = n())

bar_plot_1<- fishers_co_deletion_per_gene_long_table_significanct_p0.05_more_than_20 %>%
  dplyr::group_by(cancer, target_gene) %>%
  dplyr::summarise(total = n())

bar_plot_2<- bar_plot_1 %>%
  dplyr::group_by(cancer) %>%
  dplyr::summarise(total = n())

bar_plot<- data.frame(cancer = rep(bar_plot_2$cancer,2), genes =  rbind(bar_plot_2[,2], bar_plot_3[,2]))
bar_plot$cancer<- factor(bar_plot$cancer, levels = bar_plot$cancer)
bar_plot
dim(bar_plot)

Key<- c(rep("Significant Tumour Suppressors",14), rep("Significant co-deletions",14))
Key

ggplot(bar_plot, aes(cancer, c(target_genes, total))) + 
  geom_bar(stat = "identity", aes(fill = type), position = "dodge") +
  xlab("Months") + ylab("Count") +
  ggtitle("Chickens & Eggs") +
  theme_bw()

############
months <-rep(c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec"), 2)
chickens <-c(1, 2, 3, 3, 3, 4, 5, 4, 3, 4, 2, 2)
eggs <-c(0, 8, 10, 13, 16, 20, 25, 20, 18, 16, 10, 8)
values <-c(chickens, eggs)
type <-c(rep("chickens", 12), rep("eggs", 12))
mydata <-data.frame(months, values)
mydata
type <-c(rep("chickens", 14), rep("eggs", 14))

mydata<- bar_plot
colnames(mydata)<- c("months", "values")
mydata$months<- factor(mydata$months, levels = mydata$months)

p <-ggplot(mydata, aes(months, values))
p +geom_bar(stat = "identity", aes(fill = type), position = "dodge") +
  xlab("Months") + ylab("Count") +
  ggtitle("Chickens & Eggs") +
  theme_bw()

bar_plot<- data.frame(cancer = bar_plot_2$cancer, as.numeric(bar_plot_2[,2]), bar_plot_3[,2])
bar_plot$cancer<- factor(bar_plot$cancer, levels = bar_plot$cancer)
colnames(bar_plot)<- c("cancer", "Significant_Tumour_Suppressor", "Significant_codeletions")
bar_plot
dim(bar_plot)
class(bar_plot$Significant_Tumour_Suppressor)

p <-ggplot(bar_plot, aes(x = cancer, y = Significant_Tumour_Suppressor))
p +geom_bar(stat = "identity") +
  xlab("Cancer") + ylab("Number of Significant Tumour suppressors") +
  ggtitle("Number of Significant Tumour suppressors per Cancer") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)
        ) +
##Save plot
ggsave("bar_fishers_SignifCodel_tumour_suppressors.tiff")


p <-ggplot(bar_plot, aes(x = cancer, y = Significant_codeletions))
p +geom_bar(stat = "identity") +
  xlab("Cancer") + ylab("Number of Significant Co-deletions") +
  ggtitle("Number of Significant Co-deletions per Cancer") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)
  ) +
  
  ##Save plot
  ggsave("bar_fishers_SignifCodel_codeletion.tiff")





###################
## Table of significant genes




####################
### Fishers exact test per cancer between cat 2 and 1, 2 and 3 and 2 and 4 for long 169 gene list
#####################
## Is a co-deletion found in a particular cancer more than any other cancer?

###############
###Import data
fishers_test_per_cancer_list<- readRDS(file = "./R workspaces/fishers_test_per_cancer_results_list_169genes")
fishers_test_per_cancer_list[[1]]
fishers_test_per_cancer_list[[33]]

## Remove empty value in list item 33
fishers_test_per_cancer_list<- fishers_test_per_cancer_list[1:32]

length(fishers_test_per_cancer_list)
names(fishers_test_per_cancer_list)<- names(threshold_selected_cnv_list_plus_all_loc[1:32])
names(fishers_test_per_cancer_list)
colnames(fishers_test_per_cancer_list[[1]])

###############
###Perform BH multiple testing correction

##Add a column to each table in list with the cancer name

dim(fishers_test_per_cancer_list[[1]])
fishers_test_per_cancer_list[[1]][1:10,]
names(fishers_test_per_cancer_list[1])

for (i in 1: length(fishers_test_per_cancer_list)){
  
  column_name<- names(fishers_test_per_cancer_list[i])
  column_name<- rep(column_name, nrow(fishers_test_per_cancer_list[[i]]))
  fishers_test_per_cancer_list[[i]]<-cbind(cancer = column_name, fishers_test_per_cancer_list[[i]])
  
}
head(fishers_test_per_cancer_list[[1]])
head(fishers_test_per_cancer_list[[32]])

## Join tables together
fishers_test_per_cancer_long_table<- do.call(rbind, fishers_test_per_cancer_list)
dim(fishers_test_per_cancer_list)
dim(fishers_test_per_cancer_long_table)

##BH multiple testing for each p-value in list:
colnames(fishers_test_per_cancer_long_table)
fishers_test_per_cancer_long_table$BH_adjust_cat2vs1<- p.adjust(fishers_test_per_cancer_long_table$`co-deletion_p-value_cat2vs1`, method = "BH")
fishers_test_per_cancer_long_table$BH_adjust_cat2vs3<- p.adjust(fishers_test_per_cancer_long_table$`co-deletion_p-value_cat2vs3`, method = "BH")
fishers_test_per_cancer_long_table$BH_adjust_cat2vs4<- p.adjust(fishers_test_per_cancer_long_table$`co-deletion_p-value_cat2vs4`, method = "BH")
head(fishers_test_per_cancer_long_table)

#####################
###Identify significant genes
#cat2vs1 = co-deletions v's deletion of tumour suppressor alone
#cat2vs3 = co-deletions v's deletion of proximal gene alone
#cat2vs4 = co-deletions v's deletion of no genes

## Cat2vs1: Identify the number of significantly co-deleted genes with p-value = 0.05
sum(fishers_test_per_cancer_long_table$BH_adjust_cat2vs1 <= 0.05)

fishers_test_per_cancer_long_table_significant_p0_05<- fishers_test_per_cancer_long_table$BH_adjust_cat2vs1 <= 0.05
fishers_test_per_cancer_long_table_significant_p0_05<- fishers_test_per_cancer_long_table[fishers_test_per_cancer_long_table_significant_p0_05,]
dim(fishers_test_per_cancer_long_table)
dim(fishers_test_per_cancer_long_table_significant_p0_05)
head(fishers_test_per_cancer_long_table_significant_p0_05)

##number of significance genes per cancer type:
table(fishers_test_per_cancer_long_table_significant_p0_05$cancer)
#e.g. In BRCA, 196 genes are more significantly co-deleted than deleted compared to the other cancers.
## Identify BRCA genes:
test<- "BRCA"
fishers_test_per_cancer_long_table_significant_p0_05 %>%
  dplyr::filter(cancer == test) %>%
  dplyr::select(target_gene) %>%
  unique

fishers_test_per_cancer_long_table_significant_p0_05 %>%
  dplyr::filter(cancer == test) %>%
  dplyr::group_by(target_gene) %>%
  dplyr::select(proximal_gene) %>%
  do.call(rbind, .)

## Identify ACC genes:
fishers_test_per_cancer_long_table_significant_p0_05 %>%
  dplyr::filter(cancer == "ACC") %>%
  dplyr::select(target_gene) %>%
  unique

## Identify BLCA genes:
fishers_test_per_cancer_long_table_significant_p0_05 %>%
  dplyr::filter(cancer == "BLCA") %>%
  dplyr::select(target_gene) %>%
  unique

## Cat2vs3: Identify the number of significantly co-deleted genes with p-value = 0.05
sum(fishers_test_per_cancer_long_table$BH_adjust_cat2vs3 <= 0.05)

fishers_test_per_cancer_long_table_significant_p0_05<- fishers_test_per_cancer_long_table$BH_adjust_cat2vs3 <= 0.05
fishers_test_per_cancer_long_table_significant_p0_05<- fishers_test_per_cancer_long_table[fishers_test_per_cancer_long_table_significant_p0_05,]
dim(fishers_test_per_cancer_long_table)
dim(fishers_test_per_cancer_long_table_significant_p0_05)
head(fishers_test_per_cancer_long_table_significant_p0_05)

##number of significance genes per cancer type:
table(fishers_test_per_cancer_long_table_significant_p0_05$cancer)
#e.g. In BRCA, 196 genes are more significantly co-deleted than deleted compared to the other cancers.
## Identify BRCA genes:
test<- "BRCA"
fishers_test_per_cancer_long_table_significant_p0_05 %>%
  dplyr::filter(cancer == test) %>%
  dplyr::select(target_gene) %>%
  unique

fishers_test_per_cancer_long_table_significant_p0_05 %>%
  dplyr::filter(cancer == test) %>%
  dplyr::group_by(target_gene) %>%
  dplyr::select(proximal_gene) %>%
  do.call(rbind, .)

## Identify ACC genes:
fishers_test_per_cancer_long_table_significant_p0_05 %>%
  dplyr::filter(cancer == "ACC") %>%
  dplyr::select(target_gene) %>%
  unique

## Identify BLCA genes:
fishers_test_per_cancer_long_table_significant_p0_05 %>%
  dplyr::filter(cancer == "BLCA") %>%
  dplyr::select(target_gene) %>%
  unique

## Cat2vs4: Identify the number of significantly co-deleted genes with p-value = 0.05
sum(fishers_test_per_cancer_long_table$BH_adjust_cat2vs4 <= 0.05)

fishers_test_per_cancer_long_table_significant_p0_05<- fishers_test_per_cancer_long_table$BH_adjust_cat2vs4 <= 0.05
fishers_test_per_cancer_long_table_significant_p0_05<- fishers_test_per_cancer_long_table[fishers_test_per_cancer_long_table_significant_p0_05,]
dim(fishers_test_per_cancer_long_table)
dim(fishers_test_per_cancer_long_table_significant_p0_05)
head(fishers_test_per_cancer_long_table_significant_p0_05)

##number of significance genes per cancer type:
table(fishers_test_per_cancer_long_table_significant_p0_05$cancer)
#e.g. In BRCA, 196 genes are more significantly co-deleted than deleted compared to the other cancers.
## Identify BRCA genes:
test<- "BRCA"
fishers_test_per_cancer_long_table_significant_p0_05 %>%
  dplyr::filter(cancer == test) %>%
  dplyr::select(target_gene) %>%
  unique

fishers_test_per_cancer_long_table_significant_p0_05 %>%
  dplyr::filter(cancer == test) %>%
  dplyr::select(target_gene) %>%
  unique %>%
  dim()

fishers_test_per_cancer_long_table_significant_p0_05 %>%
  dplyr::filter(cancer == test) %>%
  dplyr::group_by(target_gene) %>%
  dplyr::select(proximal_gene) %>%
  do.call(rbind, .)

## Identify ACC genes:
fishers_test_per_cancer_long_table_significant_p0_05 %>%
  dplyr::filter(cancer == "ACC") %>%
  dplyr::select(target_gene) %>%
  unique

## Identify BLCA genes:
fishers_test_per_cancer_long_table_significant_p0_05 %>%
  dplyr::filter(cancer == "BLCA") %>%
  dplyr::select(target_gene) %>%
  unique


####################
### Identify significant genes with greater than 20 values in cat 1 and 2 for the cancer being investigated.

########
##pvalue <= 0.05

fishers_test_per_cancer_long_table_significant_p0_05<- fishers_test_per_cancer_long_table$BH_adjust_cat2vs1 <= 0.05
fishers_test_per_cancer_long_table_significant_p0_05<- fishers_test_per_cancer_long_table[fishers_test_per_cancer_long_table_significant_p0_05,]

fishers_test_per_cancer_long_table_significant_p0_05[1:10,]
dim(fishers_test_per_cancer_long_table_significant_p0_05)

fishers_test_per_cancer_long_table_significant_p0_05_more_than_20<- 
  fishers_test_per_cancer_long_table_significant_p0_05 %>% 
  dplyr::filter(n_cat1_per_cancer >= 20 & n_cat2_per_cancer >= 20)

dim(fishers_test_per_cancer_long_table_significant_p0_05_more_than_20)

##number of significance genes per cancer type:
table(fishers_test_per_cancer_long_table_significant_p0_05_more_than_20$cancer)

fishers_test_per_cancer_long_table_significant_p0_05_more_than_20 %>%
  dplyr::group_by(cancer) %>%
  dplyr::summarise(total = n())

## Identify BRCA genes: target gene
test<- "BRCA"
fishers_test_per_cancer_long_table_significant_p0_05_more_than_20 %>%
  dplyr::filter(cancer == test) %>%
  dplyr::select(target_gene) %>%
  unique

## Identify BRCA genes: proximal genes
fishers_test_per_cancer_long_table_significant_p0_05_more_than_20 %>%
  dplyr::filter(cancer == test) %>%
  dplyr::group_by(target_gene) %>%
  dplyr::select(proximal_gene) %>%
  do.call(rbind, .)

## Identify BRCA genes: target gene
test<- "SKCM"
fishers_test_per_cancer_long_table_significant_p0_05_more_than_20 %>%
  dplyr::filter(cancer == test) %>%
  dplyr::select(target_gene) %>%
  unique

## Identify BRCA genes: proximal genes
fishers_test_per_cancer_long_table_significant_p0_05_more_than_20 %>%
  dplyr::filter(cancer == test) %>%
  dplyr::group_by(target_gene) %>%
  dplyr::select(proximal_gene) %>%
  do.call(rbind, .)

## Create csv for BRCA
fishers_co_deletion_per_cancer_BRCA<- fishers_test_per_cancer_long_table_significant_p0_05_more_than_20 %>%
  dplyr::filter(cancer == "BRCA") %>%
  dplyr::arrange(BH_adjust_cat2vs1)

fishers_co_deletion_per_cancer_BRCA

#write.csv(fishers_co_deletion_per_cancer_BRCA, file = "fishers_co_deletion_per_cancer_BRCA.csv", quote = FALSE)

## Create csv for SKCM
fishers_co_deletion_per_cancer_SKCM<- fishers_test_per_cancer_long_table_significant_p0_05_more_than_20 %>%
  dplyr::filter(cancer == "SKCM") %>%
  dplyr::arrange(BH_adjust_cat2vs1)

fishers_co_deletion_per_cancer_SKCM

#write.csv(fishers_co_deletion_per_cancer_SKCM, file = "fishers_co_deletion_per_cancer_SKCM.csv", quote = FALSE)


#write.csv(fishers_test_per_cancer_long_table_significant_p0_05_more_than_20, file = "fishers_test_per_cancer_long_table_significant_p0_05_more_than_20.csv", quote = FALSE)

###################
### Bar plot of Fishers exact test results for co-deletions between cancers
##bar_fishers_SignifCodel_perCancer

## p<= 0.05
dim(fishers_test_per_cancer_long_table_significant_p0_05_more_than_20)
head(fishers_test_per_cancer_long_table_significant_p0_05_more_than_20)

table(fishers_test_per_cancer_long_table_significant_p0_05_more_than_20$cancer)
bar_plot_3<- fishers_test_per_cancer_long_table_significant_p0_05_more_than_20 %>%
  dplyr::group_by(cancer) %>%
  dplyr::summarise(total = n())

bar_plot_3

bar_plot_1<- fishers_test_per_cancer_long_table_significant_p0_05_more_than_20 %>%
  dplyr::group_by(cancer, target_gene) %>%
  dplyr::summarise(total = n())

bar_plot_1

bar_plot_2<- bar_plot_1 %>%
  dplyr::group_by(cancer) %>%
  dplyr::summarise(total = n())

bar_plot_2

bar_plot<- data.frame(cancer = rep(bar_plot_2$cancer,2), genes =  rbind(bar_plot_2[,2], bar_plot_3[,2]))
bar_plot$cancer<- factor(bar_plot$cancer, levels = bar_plot$cancer)
bar_plot
dim(bar_plot)

# Key<- c(rep("Significant Tumour Suppressors",14), rep("Significant co-deletions",14))
# Key

bar_plot<- data.frame(cancer = bar_plot_2$cancer, bar_plot_2[,2], bar_plot_3[,2])
bar_plot$cancer<- factor(bar_plot$cancer, levels = bar_plot$cancer)
colnames(bar_plot)<- c("cancer", "Significant_Tumour_Suppressor", "Significant_codeletions")
bar_plot
dim(bar_plot)
class(bar_plot$Significant_Tumour_Suppressor)

p <-ggplot(bar_plot, aes(x = cancer, y = Significant_Tumour_Suppressor))
p +geom_bar(stat = "identity") +
  xlab("Cancer") + ylab("Number of Significant Tumour suppressors") +
  ggtitle("Number of Significant Tumour suppressors per Cancer") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)
  ) +
  ##Save plot
  ggsave("bar_fishers_SignifCodel_perCancer_tumour_suppressors.tiff")


p <-ggplot(bar_plot, aes(x = cancer, y = Significant_codeletions))
p +geom_bar(stat = "identity") +
  xlab("Cancer") + ylab("Number of Significant Co-deletions") +
  ggtitle("Number of Significant Co-deletions per Cancer") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)
  ) +
  
  ##Save plot
  ggsave("bar_fishers_SignifCodel_perCancer_codeletion.tiff")









##########################
### Venn diagrams to compare results between Fishers exact tests

## Intersection of Tumour suppressors

venn.plot <- venn.diagram(
  x = list(
    fisher_codeletions = fishers_co_deletion_per_gene_long_table_significanct_p0.05_more_than_20$target_gene,
    fisher_cancer = fishers_test_per_cancer_long_table_significant_p0_05_more_than_20$target_gene
  ),
  euler.d = TRUE,
  scaled = TRUE,
  filename = "Euler_3set_scaled.tiff",
  cex = 2.5,
  cat.cex = 1.5,
  fill=c("white", "white"),
  col= c("blue", "green"), 
  alpha=c(0.5,0.5),
  cat.col = c("blue", "green"),
  cat.pos = c(0,0),
  cat.fontface=4,
  ext.text = TRUE,
  category.names=c("Significant co-deletions", "Significant co-deletions between cancers"),
  main="Significant Tumour suppressors")

###########
##Intersection of co-deletions

venn.plot <- venn.diagram(
  x = list(
    fisher_codeletions = fishers_co_deletion_per_gene_long_table_significanct_p0.05_more_than_20$proximal_gene,
    fisher_cancer = fishers_test_per_cancer_long_table_significant_p0_05_more_than_20$proximal_gene
  ),
  euler.d = TRUE,
  scaled = TRUE,
  filename = "Euler_3set_scaled.tiff",
  cex = 2.5,
  cat.cex = 1.5,
  fill=c("white", "white"),
  col= c("blue", "green"), 
  alpha=c(0.5,0.5),
  cat.col = c("blue", "green"),
  cat.pos = c(0,0),
  cat.fontface=4,
  ext.text = TRUE,
  category.names=c("Significant co-deletions", "Significant co-deletions between cancers"),
  main="Significant Co-deletions")






