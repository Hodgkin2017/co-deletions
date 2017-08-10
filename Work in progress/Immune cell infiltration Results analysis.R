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
sapply(immune_cell_infiltrate_annova_per_cell_type_list, function(x) sum(x$BH_adjust_cat2_1 <= 0.1 & x$number_cat1 >= 20 & x$number_cat2 >= 20, na.rm = TRUE))

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

immune_cell_infiltrate_annova_per_cell_type_significant_list<- lapply(immune_cell_infiltrate_annova_per_cell_type_list, function(x) significant_selection(x, 0.1))
sapply(immune_cell_infiltrate_annova_per_cell_type_significant_list, function(x) nrow(x))
sapply(immune_cell_infiltrate_annova_per_cell_type_list, function(x) unique(x$cell_type))
sapply(immune_cell_infiltrate_annova_per_cell_type_significant_list, function(x) unique(x$cancer))

immune_cell_infiltrate_annova_per_cell_type_significant_table<- do.call(rbind, immune_cell_infiltrate_annova_per_cell_type_significant_list)
immune_cell_infiltrate_annova_per_cell_type_significant_table<- immune_cell_infiltrate_annova_per_cell_type_significant_table %>%
  dplyr::arrange(cancer, target_gene, cell_type, BH_adjust_cat2_1)

head(immune_cell_infiltrate_annova_per_cell_type_significant_table)
write.csv(immune_cell_infiltrate_annova_per_cell_type_significant_table, file = "immune_cell_infiltrate_annova_per_cell_type_significant_table_p0_1_greather_than_20.csv")

immune_cell_infiltrate_annova_per_cell_type_significant_table<- immune_cell_infiltrate_annova_per_cell_type_significant_table %>%
  dplyr::arrange(cell_type, cancer, target_gene, BH_adjust_cat2_1)

head(immune_cell_infiltrate_annova_per_cell_type_significant_table)
write.csv(immune_cell_infiltrate_annova_per_cell_type_significant_table, file = "immune_cell_infiltrate_annova_per_cell_type_significant_table_p0_1_greather_than_20_by_celltype.csv")


#########################
### Bar plot of number of significant co-deletions per cell type:


bar_plot1<- sapply(immune_cell_infiltrate_annova_per_cell_type_significant_list, function(x) nrow(x))
bar_plot2<- sapply(immune_cell_infiltrate_annova_per_cell_type_list, function(x) unique(x$cell_type))

barplot<-data.frame(cell_type = bar_plot2, value = bar_plot1)
barplot<- barplot %>%
  dplyr::filter(value > 0)
barplot

p <-ggplot(barplot, aes(x = cell_type, y = value))
p +geom_bar(stat = "identity") +
  xlab("Immune Cell Type") + ylab("Number of Significant Co-deletions") +
  ggtitle("Number of Significant Co-deletions per Immune Cell Type") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)
  ) +
  ##Save plot
  ggsave("bar_all_immune_codeletions.tiff")

#######
#Bar plot of intersections with significant co-deleted genes?

immune_cell_infiltrate_intersection<- read.csv("../../Output/Immune_Infiltration/immune_cell_infiltration_intersection.csv", header = TRUE, stringsAsFactors = FALSE)
dim(immune_cell_infiltrate_intersection)
head(immune_cell_infiltrate_intersection)

bar_plot1<- unique(immune_cell_infiltrate_intersection$cell_type)
bar_plot1
bar_plot2<- immune_cell_infiltrate_intersection %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(n())
bar_plot2
bar_plot2[,2]

barplot<-data.frame(cell_type = bar_plot1, value = bar_plot2[,2])
colnames(barplot)<- c("cell_type", "value")
barplot

p <-ggplot(barplot, aes(x = cell_type, y = value))
p +geom_bar(stat = "identity") +
  xlab("Immune Cell Type") + ylab("Number of Significant Co-deletions") +
  ggtitle("Number of Significant Co-deletions per Immune Cell Type") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)
  ) +
  ##Save plot
  ggsave("bar_significant_selected_immune_codeletions.tiff")





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


###############
##CDKN2A

CDKN2A_significant_immune_infiltration<- lapply(immune_cell_infiltrate_annova_per_cell_type_list, function(x)significant_gene_selection(x, "CDKN2A"))
names(CDKN2A_significant_immune_infiltration)<- immune_cell_types
CDKN2A_significant_immune_infiltration_table<- do.call(rbind, CDKN2A_significant_immune_infiltration) %>%
  dplyr::arrange(cell_type, cancer, target_gene, p_value_cat2_1)

head(CDKN2A_significant_immune_infiltration_table)
write.csv(CDKN2A_significant_immune_infiltration_table, file = "CDKN2A_significant_immune_infiltration_table.csv")

sapply(CDKN2A_significant_immune_infiltration, function(x) nrow(x))

CDKN2A_significant_immune_infiltration$`Activated Natural killer cell`
CDKN2A_significant_immune_infiltration$`CD8 T cell`



##########
##CDKN1B
CDKN1B_significant_immune_infiltration<- lapply(immune_cell_infiltrate_annova_per_cell_type_list, function(x)significant_gene_selection(x, "CDKN1B"))
names(CDKN1B_significant_immune_infiltration)<- immune_cell_types
CDKN1B_significant_immune_infiltration_table<- do.call(rbind, CDKN1B_significant_immune_infiltration) %>%
  dplyr::arrange(cell_type, cancer, target_gene, p_value_cat2_1)

sapply(CDKN1B_significant_immune_infiltration, function(x) nrow(x))
head(CDKN1B_significant_immune_infiltration_table)
write.csv(CDKN1B_significant_immune_infiltration_table, file = "CDKN1B_significant_immune_infiltration_table.csv")

##########
##LRP1B
LRP1B_significant_immune_infiltration<- lapply(immune_cell_infiltrate_annova_per_cell_type_list, function(x)significant_gene_selection(x, "LRP1B"))
names(LRP1B_significant_immune_infiltration)<- immune_cell_types
LRP1B_significant_immune_infiltration_table<- do.call(rbind, LRP1B_significant_immune_infiltration) %>%
  dplyr::arrange(cell_type, cancer, target_gene, p_value_cat2_1)

sapply(LRP1B_significant_immune_infiltration, function(x) nrow(x))
head(LRP1B_significant_immune_infiltration_table)
write.csv(LRP1B_significant_immune_infiltration_table, file = "LRP1B_significant_immune_infiltration_table.csv")


##########
##RB1
RB1_significant_immune_infiltration<- lapply(immune_cell_infiltrate_annova_per_cell_type_list, function(x)significant_gene_selection(x, "RB1"))
names(RB1_significant_immune_infiltration)<- immune_cell_types
RB1_significant_immune_infiltration_table<- do.call(rbind, RB1_significant_immune_infiltration) %>%
  dplyr::arrange(cell_type, cancer, target_gene, p_value_cat2_1)

sapply(RB1_significant_immune_infiltration, function(x) nrow(x))
head(RB1_significant_immune_infiltration_table)
write.csv(RB1_significant_immune_infiltration_table, file = "RB1_significant_immune_infiltration_table.csv")

##########
##TGFBR2
TGFBR2_significant_immune_infiltration<- lapply(immune_cell_infiltrate_annova_per_cell_type_list, function(x)significant_gene_selection(x, "TGFBR2"))
names(TGFBR2_significant_immune_infiltration)<- immune_cell_types
TGFBR2_significant_immune_infiltration_table<- do.call(rbind, TGFBR2_significant_immune_infiltration) %>%
  dplyr::arrange(cell_type, cancer, target_gene, p_value_cat2_1)

sapply(TGFBR2_significant_immune_infiltration, function(x) nrow(x))
head(TGFBR2_significant_immune_infiltration_table)
write.csv(TGFBR2_significant_immune_infiltration_table, file = "TGFBR2_significant_immune_infiltration_table.csv")

##########
##TP53
TP53_significant_immune_infiltration<- lapply(immune_cell_infiltrate_annova_per_cell_type_list, function(x)significant_gene_selection(x, "TP53"))
names(TP53_significant_immune_infiltration)<- immune_cell_types
TP53_significant_immune_infiltration_table<- do.call(rbind, TP53_significant_immune_infiltration) %>%
  dplyr::arrange(cell_type, cancer, target_gene, p_value_cat2_1)

sapply(TP53_significant_immune_infiltration, function(x) nrow(x))
head(TP53_significant_immune_infiltration_table)
write.csv(TP53_significant_immune_infiltration_table, file = "TP53_significant_immune_infiltration_table.csv")


##########
##ZFHX3
ZFHX3_significant_immune_infiltration<- lapply(immune_cell_infiltrate_annova_per_cell_type_list, function(x)significant_gene_selection(x, "ZFHX3"))
names(ZFHX3_significant_immune_infiltration)<- immune_cell_types
ZFHX3_significant_immune_infiltration_table<- do.call(rbind, ZFHX3_significant_immune_infiltration) %>%
  dplyr::arrange(cell_type, cancer, target_gene, p_value_cat2_1)

sapply(ZFHX3_significant_immune_infiltration, function(x) nrow(x))
head(ZFHX3_significant_immune_infiltration_table)
write.csv(ZFHX3_significant_immune_infiltration_table, file = "ZFHX3_significant_immune_infiltration_table.csv")





#################
###Boxplots
##http://ggplot2.tidyverse.org/reference/geom_boxplot.html
## Will I need to create tables of co-deletion catagories again?
mpg
p <- ggplot(mpg, aes(class, hwy))
p + geom_boxplot()


immune_cell_infiltrate_signif<- read.csv("../../Output/Immune_Infiltration/immune_cell_infiltrate_annova_per_cell_type_significant_table_p0_1_greather_than_20.csv", header = TRUE, stringsAsFactors = FALSE)
immune_cell_infiltrate_signif

df<- immune_cell_infiltrate_signif %>%
  dplyr::select(cell_type, mean_cibersort_cat1, mean_cibersort_cat2)
df3<- cbind(df[,1:2], rep("cat_1", nrow(df)))
df4<- cbind(df[,c(1,3)], rep("cat_2", nrow(df)))
colnames(df3)<- c("cell_type", "mean_CIBERSORT", "category")
colnames(df4)<- c("cell_type", "mean_CIBERSORT", "category")
df2<-rbind(df3, df4)
  
##############
p <- ggplot(df, aes(cell_type, mean_cibersort_cat1))
p + geom_boxplot()
p + geom_boxplot() + geom_jitter(width = 0.2)
################
ggplot(df2, aes(cell_type, mean_CIBERSORT)) +
  geom_boxplot(aes(group = category))
#############
ggplot(aes(y = mean_CIBERSORT, x = cell_type, fill = category), data = df2) + 
  geom_boxplot() + 
  geom_jitter(width = 0.5) +
  theme(legend.position="bottom") +
  xlab("Immune cell type") +
  ylab("CIBERSORT value") +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5, size = 10))

  
ggsave("box_immune_infiltration.tiff",width = 14, height = 10, dpi = 300)

#############
ggplot(aes(y = mean_CIBERSORT, x = cell_type, fill=factor(category,labels=c("Single deletion","Co-deletion"))), data = df2) + 
  geom_boxplot() + 
  #geom_jitter(width = 0.5) +
  theme(legend.position="bottom") +
  xlab("Immune cell type") +
  ylab("CIBERSORT value") +
  labs(fill="Deletion category") #+
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5, size = 10))
  
  ggsave("box_immune_infiltration.tiff",width = 12, height = 8, dpi = 300)
  

####################
# immune_cell_infiltrate_signif_intersection<- read.csv("../../Output/Immune_Infiltration/", header = TRUE, stringsAsFactors = FALSE)
#   immune_cell_infiltrate_signif_intersection
#   
#   df<- immune_cell_infiltrate_signif_all %>%
#     dplyr::select(cell_type, mean_cibersort_cat1, mean_cibersort_cat2)
#   df3<- cbind(df[,1:2], rep("cat_1", nrow(df)))
#   df4<- cbind(df[,c(1,3)], rep("cat_2", nrow(df)))
#   colnames(df3)<- c("cell_type", "mean_CIBERSORT", "category")
#   colnames(df4)<- c("cell_type", "mean_CIBERSORT", "category")
#   df2<-rbind(df3, df4)
#   
#   
#   
#   #############
#   ggplot(aes(y = mean_CIBERSORT, x = cell_type, fill=factor(category,labels=c("Single deletion","Co-deletion"))), data = df2) + 
#     geom_boxplot() + 
#     geom_jitter(width = 0.5) +
#     theme(legend.position="bottom") +
#     xlab("Immune cell type") +
#     ylab("CIBERSORT value") +
#     labs(fill="Deletion category") #+
#   theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5, size = 10))
#   
#   ggsave("box_immune_infiltration.tiff",width = 12, height = 8, dpi = 300)
#   
##########
  #######################
  ##Venn diagram to show intersection of genes between overall survival and disease free survival
  ## Intersection of Tumour suppressors
  
  venn.plot <- venn.diagram(
    x = list(
      overall_survival = survival_stats_ovsurv_cat1_2_significant_table_p0_1_more_than_20$target_gene,
      immune = immune_cell_infiltrate_intersection$target_gene
    ),
    euler.d = TRUE,
    scaled = TRUE,
    filename = "venn_immune_intersectio_OvSurv_TS.tiff",
    cex = 2.5,
    cat.cex = 1.5,
    fill=c("white", "white"),
    col= c("blue", "green"), 
    alpha=c(0.5,0.5),
    cat.col = c("blue", "green"),
    cat.pos = c(0,180),
    cat.fontface=4,
    ext.text = TRUE,
    category.names=c("Overall Survival", "Immune cell Infiltration"),
    main="Significant Tumour suppressors")
  
  
  venn.plot <- venn.diagram(
    x = list(
      overall_survival = survival_stats_ovsurv_cat1_2_significant_table_p0_1_more_than_20$proximal_gene,
      immune = immune_cell_infiltrate_intersection$proximal_gene
    ),
    euler.d = TRUE,
    scaled = TRUE,
    filename = "venn_immune_intersectio_OvSurv_codel.tiff",
    cex = 2.5,
    cat.cex = 1.5,
    fill=c("white", "white"),
    col= c("blue", "green"), 
    alpha=c(0.5,0.5),
    cat.col = c("blue", "green"),
    cat.pos = c(0,180),
    cat.fontface=4,
    ext.text = TRUE,
    category.names=c("Overall Survival", "Immune cell Infiltration"),
    main="Significant Co-deletions")
  