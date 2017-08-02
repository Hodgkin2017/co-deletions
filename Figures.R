####################
### Figures
####################

######################################
###From: Number of deletions or amplifications per cytoband long.R

### heatmap_deletion_perCytoband

CNV.data<- threshold_selected_cnv_list_plus_all_loc
names(threshold_selected_cnv_list_plus_all_loc)
acc.cnv.chr.location<- threshold_selected_cnv_list_plus_all_loc[[1]]
colnames(acc.cnv.chr.location)[1:11]

cytoband.del.matrix<- function(x,y){
  
  z<- x
  cytoband.list<- events.per.cytoband(z, threshold = -2, cytoband_column = 10, column_data_start = 11, chromosome_interval = 0,  deletion = TRUE)
  return(cytoband.list[[2]]$proportion.of.deletions)
  
}

system.time(test3<- lapply(CNV.data, function(x) {cytoband.del.matrix(x, acc.cnv.chr.location)}))
test4<- do.call(cbind, test3) %>%as.matrix

##Add row and column names
rownames(test4)<- test1[[2]]$cytoband
#colnames(test4)<- c(names(cnv.list), "ALL")
colnames(test4)<- names(CNV.data)

##Make annotation row dataframe
annotation_row<- data.frame(chromosome = test1[[2]]$chromosome)
rownames(annotation_row)<- test1[[2]]$cytoband
dim(annotation_row)
head(annotation_row)

##Make heat map

col.pal<- colorRampPalette(c( "white","navy", "firebrick3"))(1000)

pheatmap(test4,
         cluster_row = F,
         cluster_cols = F,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = col.pal,
         fontsize_row=1,
         #cellwidth = 10,
         annotation_row = annotation_row,
         annotation_legend = FALSE
)


######
### heatmap_amplification_perCytoband

CNV.data<- threshold_selected_cnv_list_plus_all_loc
names(threshold_selected_cnv_list_plus_all_loc)
acc.cnv.chr.location<- threshold_selected_cnv_list_plus_all_loc[[1]]
colnames(acc.cnv.chr.location)[1:11]

cytoband.del.matrix<- function(x,y){
  
  z<- x
  cytoband.list<- events.per.cytoband(z, threshold = 2, cytoband_column = 10, column_data_start = 11, chromosome_interval = 0,  deletion = FALSE)
  return(cytoband.list[[2]]$proportion.of.deletions)
  
}

system.time(test3<- lapply(CNV.data, function(x) {cytoband.del.matrix(x, acc.cnv.chr.location)}))
test4<- do.call(cbind, test3) %>%as.matrix

##Add row and column names
rownames(test4)<- test1[[2]]$cytoband
#colnames(test4)<- c(names(cnv.list), "ALL")
colnames(test4)<- names(CNV.data)

##Make annotation row dataframe
annotation_row<- data.frame(chromosome = test1[[2]]$chromosome)
rownames(annotation_row)<- test1[[2]]$cytoband
dim(annotation_row)
head(annotation_row)

##Make heat map

col.pal<- colorRampPalette(c( "white","navy", "firebrick3"))(1000)

pheatmap(test4,
         cluster_row = F,
         cluster_cols = F,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = col.pal,
         fontsize_row=1,
         #cellwidth = 10,
         annotation_row = annotation_row,
         annotation_legend = FALSE
)


###########################
### heatmap_deletion_perMB

##################
### Sum of deletions for all tumour types together:

acc.cnv.chr.location<- threshold_selected_cnv_list_plus_all_loc[[1]]

test2<- events.per.cytoband(acc.cnv.chr.location, threshold = -1, cytoband_column = 10, column_data_start = 11 , select_chromosome = 1, chromosome_interval = 10,  deletion = TRUE)

chromosome_interval_1MB<- ceiling(test2[[3]]$intervals_for_Mb)
chromosome_interval_1MB
##Create one large dataframe with all CNV data in it:
#x<-join.cnv.datasets(cnv.list, 4)
#x<- chromosomal_location(CNV.all.table)
names(threshold_selected_cnv_list_plus_all_loc)
x<- threshold_selected_cnv_list_plus_all_loc[[33]]
dim(x)

##Matrix to store results per chromosome:
heatmap.matrix.chr.interval250<- data.frame(matrix(NA, ncol = 1, nrow = 2500))
dim(heatmap.matrix.chr.interval250)
heatmap.matrix.chr.interval250[,1]<-seq(1:2500)
colnames(heatmap.matrix.chr.interval250)<-c("intervals")
colnames(heatmap.matrix.chr.interval250)
chr<-c(seq(1:22), "X")
chr

##Calculate proportion of deletions per cytoband and add to matrix
for (i in 1:23) {
  j<-chr[i]
  chromosome_interval<-chromosome_interval_1MB[i]
  cytoband.list<- events.per.cytoband(x, threshold = -2, cytoband_column = 10, column_data_start = 11, select_chromosome = j , chromosome_interval = chromosome_interval,  deletion = TRUE)
  y<-data.frame(intervals = seq(1:length(cytoband.list[[4]]$proportion.of.deletions)), chr = cytoband.list[[4]]$proportion.of.deletions)
  # glimpse(cytoband.list[[1]])
  # glimpse(cytoband.list[[2]])
  # glimpse(cytoband.list[[3]])
  # glimpse(cytoband.list[[4]])
  heatmap.matrix.chr.interval250<- full_join(heatmap.matrix.chr.interval250, y, by = "intervals")
  # heatmap.matrix.chr.interval250[,i]<- cytoband.list[[4]][,5]
  print(i)
  print(length(cytoband.list[[4]]$proportion.of.deletions))
}


# a.function<- function(object_name, ){
#   
# ###.........make above loop into function...
# }

colnames(heatmap.matrix.chr.interval250)<- c("Interval", paste0("Chromosome ", seq(1:22)), "Chromosome X")

head(heatmap.matrix.chr.interval250)
ncol(heatmap.matrix.chr.interval250)

##############
### Plot heatmaps


#my.matrix<- as.matrix(heatmap.matrix.chr.interval250[,2:24])
#is.na(heatmap.matrix.chr.interval250$chr.x)


col.pal<- colorRampPalette(c( "white","navy", "firebrick3"))(1000)

pheatmap(heatmap.matrix.chr.interval250[1:205,2:24],
         cluster_row = F,
         cluster_cols = F,
         show_rownames = FALSE,
         show_colnames = TRUE,
         color = col.pal,
         fontsize_row=1
         #cellwidth = 10,
         #annotation_row = annotation_row,
         #annotation_legend = FALSE
)






###########################
### heatmap_amplification_perMB


##################
### Sum of deletions for all tumour types together:


acc.cnv.chr.location<- threshold_selected_cnv_list_plus_all_loc[[1]]

test2<- events.per.cytoband(acc.cnv.chr.location, threshold = -1, cytoband_column = 10, column_data_start = 11 , select_chromosome = 1, chromosome_interval = 10,  deletion = TRUE)

chromosome_interval_1MB<- ceiling(test2[[3]]$intervals_for_Mb)
chromosome_interval_1MB
##Create one large dataframe with all CNV data in it:
#x<-join.cnv.datasets(cnv.list, 4)
#x<- chromosomal_location(CNV.all.table)
names(threshold_selected_cnv_list_plus_all_loc)
x<- threshold_selected_cnv_list_plus_all_loc[[33]]
dim(x)

##Matrix to store results per chromosome:
heatmap.matrix.chr.interval250<- data.frame(matrix(NA, ncol = 1, nrow = 2500))
dim(heatmap.matrix.chr.interval250)
heatmap.matrix.chr.interval250[,1]<-seq(1:2500)
colnames(heatmap.matrix.chr.interval250)<-c("intervals")
colnames(heatmap.matrix.chr.interval250)
chr<-c(seq(1:22), "X")
chr

##Calculate proportion of deletions per cytoband and add to matrix
for (i in 1:23) {
  j<-chr[i]
  chromosome_interval<-chromosome_interval_1MB[i]
  cytoband.list<- events.per.cytoband(x, threshold = 2, cytoband_column = 10, column_data_start = 11, select_chromosome = j , chromosome_interval = chromosome_interval,  deletion = FALSE)
  y<-data.frame(intervals = seq(1:length(cytoband.list[[4]]$proportion.of.deletions)), chr = cytoband.list[[4]]$proportion.of.deletions)
  # glimpse(cytoband.list[[1]])
  # glimpse(cytoband.list[[2]])
  # glimpse(cytoband.list[[3]])
  # glimpse(cytoband.list[[4]])
  heatmap.matrix.chr.interval250<- full_join(heatmap.matrix.chr.interval250, y, by = "intervals")
  # heatmap.matrix.chr.interval250[,i]<- cytoband.list[[4]][,5]
  print(i)
  print(length(cytoband.list[[4]]$proportion.of.deletions))
}


# a.function<- function(object_name, ){
#   
# ###.........make above loop into function...
# }

colnames(heatmap.matrix.chr.interval250)<- c("Interval", paste0("Chromosome ", seq(1:22)), "Chromosome X")

head(heatmap.matrix.chr.interval250)
ncol(heatmap.matrix.chr.interval250)

##############
### Plot heatmaps


#my.matrix<- as.matrix(heatmap.matrix.chr.interval250[,2:24])
#is.na(heatmap.matrix.chr.interval250$chr.x)


col.pal<- colorRampPalette(c( "white","navy", "firebrick3"))(1000)

pheatmap(heatmap.matrix.chr.interval250[1:205,2:24],
         cluster_row = F,
         cluster_cols = F,
         show_rownames = FALSE,
         show_colnames = TRUE,
         color = col.pal,
         fontsize_row=1
         #cellwidth = 10,
         #annotation_row = annotation_row,
         #annotation_legend = FALSE
)





#########################
###heatmap_deletion_chromosome9
##From: Co-deletion and co-amplification heatmaps with thresholded data.R

heatmap_deletion_chromosome9<- co.deletion_co.amplification_matrix(threshold_selected_cnv_list_plus_all_loc[[33]], column_start = 11, threshold = -2,Chromosome = 9, 
                                                                        deletion = TRUE, normalisation = "total.tumour.number")


chromosome<- "9"
plot_data<- heatmap_deletion_chromosome9

tiff(paste("Chromosome_", chromosome, "_deletion = true.tiff", sep =""), width = 13, height = 11, units = 'in', res = 100)
pheatmap(plot_data,
         cluster_row = F,
         cluster_cols = F,
         show_rownames = FALSE,
         show_colnames = FALSE
)
dev.off()



#########################
###heatmap_amplification_chromosome9
##From: Co-deletion and co-amplification heatmaps with thresholded data.R

heatmap_amplification_chromosome9<- co.deletion_co.amplification_matrix(threshold_selected_cnv_list_plus_all_loc[[33]], column_start = 11, threshold = 2,Chromosome = 9, 
                                                                   deletion = FALSE, normalisation = "total.tumour.number")


chromosome<- "9"
plot_data<- heatmap_amplification_chromosome9

tiff(paste("Chromosome_", chromosome, "_deletion = true.tiff", sep =""), width = 13, height = 11, units = 'in', res = 100)
pheatmap(plot_data,
         cluster_row = F,
         cluster_cols = F,
         show_rownames = FALSE,
         show_colnames = FALSE
)
dev.off()


########################
######################
###table_top_codeleted
#Co-deletion and co-amplification heatmaps with thresholded data




#######################
######################
###Circle_perCancer_perCytoband

cytoband.table<- unique(threshold_selected_cnv_list_plus_all_loc[[1]]$Cytoband)
length(cytoband.table)
cytoband.table

co.deletions.per.cytoband.circle.plots.table<- data.frame(cytoband = NA, cancer = NA, proportion_deletions = NA)
co.deletions.per.cytoband.circle.plots.table

for (i in 1:length(threshold_selected_cnv_list_plus_all_loc)){
  
  cancer.type<- names(threshold_selected_cnv_list_plus_all_loc)[i]
  cancer.table<- threshold_selected_cnv_list_plus_all_loc[[i]]
  cancer.list<- sapply(cytoband.table, function(x) co.deletion_co.amplification_matrix(cancer.table, column_start = 11, threshold = -2, Cytoband = TRUE, selection_criteria = x, deletion = TRUE, normalisation = "frequency.for.whole.sample", remove_NA = FALSE))
  #cancer.list2<-unlist(cancer.list)
  cancer.deletions.per.cytoband<- cbind(cytoband = cytoband.table, cancer = rep(cancer.type, 806), proportion_deletions = cancer.list)
  co.deletions.per.cytoband.circle.plots.table<-rbind(co.deletions.per.cytoband.circle.plots.table, cancer.deletions.per.cytoband) 
  print(cancer.type)
}

## delete first row of deletions.per.cytoband.circle.plots.table
dim(co.deletions.per.cytoband.circle.plots.table)
co.deletions.per.cytoband.circle.plots.table<-co.deletions.per.cytoband.circle.plots.table[-1,]
head(co.deletions.per.cytoband.circle.plots.table)
tail(co.deletions.per.cytoband.circle.plots.table)
nrow(co.deletions.per.cytoband.circle.plots.table)

## Plot circle plot

#df2<- co.deletions.per.cytoband.circle.plots.table[c(1:10, 807:817),]

ggplot(co.deletions.per.cytoband.circle.plots.table, aes(y = factor(cancer),
                                                         x = factor(cytoband))) +
  xlab("Cytoband") +
  ylab("Cancer") +
  geom_point(aes(size = as.numeric(proportion_deletions))) + 
  scale_size_area("Proportion of\nco-deletion\nevents") +
  
  ##Save plot
  ggsave("co_deletions_per_cytoband_circle_plots.tiff")









#####################
######################
###Circle_perCancer_perCytoband_target_genes

#############
### Repeat co-deltion plot for target gene cytobands only:

##Get cytobands of target genes
# target.genes.cytoband <- short.cnv.list[[1]] %>% 
#   dplyr::filter(Gene.Symbol %in% target.genes) %>%
#   dplyr::select(Cytoband) %>%
#   t() %>% #required to convert output from data.frame to vector
#   as.character()

length(gene_information_long_list)
gene_information_long_list[[1]]
target_genes_names<- sapply(gene_information_long_list, function(x) x[[1]])
target_genes_names
target_genes_names_order<- order(target_genes_names)
target_genes_names_order
target_genes_names<- target_genes_names[target_genes_names_order]
target_genes_names
target_genes_cytobands<- sapply(gene_information_long_list, function(x) x[[3]])
target_genes_cytobands
length(target_genes_cytobands)
target_genes_cytobands<- target_genes_cytobands[target_genes_names_order]
target_genes_cytobands


head(co.deletions.per.cytoband.circle.plots.table)

##create table of proportion of tumours with deletions per cytoband for target genes
target.genes.co.deletions.per.cytoband.circle.plots.table<- co.deletions.per.cytoband.circle.plots.table %>% 
  dplyr::filter(cytoband %in% target_genes_cytobands)

dim(target.genes.co.deletions.per.cytoband.circle.plots.table)
# (target_genes_cytobands %in% co.deletions.per.cytoband.circle.plots.table$cytoband)


# target.genes.co.deletions.per.cytoband.circle.plots.table<- co.deletions.per.cytoband.circle.plots.table[co.deletions.per.cytoband.circle.plots.table$cytoband %in% target_genes_cytobands,]
# dim(target.genes.co.deletions.per.cytoband.circle.plots.table)

##Plot:
##Comment: in order to get the cytobands to be in numerical order need to use levels =
ggplot(target.genes.co.deletions.per.cytoband.circle.plots.table, aes(y = factor(cancer, levels = names(threshold_selected_cnv_list_plus_all_loc)),
                                                                      x = factor(cytoband, levels = mixedsort(target_genes_cytobands)))) +
  xlab("Cytoband") +
  ylab("Cancer") +
  geom_point(aes(size = as.numeric(proportion_deletions))) + 
  scale_size_area("Proportion of\nco-deletion\nevents")+
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
  ) + #theme_bw()
  
  ##Save plot
  ggsave("target_genes_co_deletions_per_cytoband_circle_plots.tiff")

########

##Plot:
##Comment: in order to get the cytobands to be in numerical order need to use levels =
ggplot(target.genes.co.deletions.per.cytoband.circle.plots.table, aes(y = factor(cytoband, levels = mixedsort(target_genes_cytobands)) ,
                                                                      x = factor(cancer, levels = names(threshold_selected_cnv_list_plus_all_loc)))) +
  xlab("Cancer") +
  ylab("Cytoband") +
  geom_point(aes(size = as.numeric(proportion_deletions))) + 
  scale_size_area("Proportion of\nco-deletion\nevents")+
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
  ) + #theme_bw()
  
  ##Save plot
  ggsave("target_genes_co_deletions_per_cytoband_circle_plots.tiff", width = 8, height = 16)


##comment: can not replace cytoband with Tumour suppressor name for some reason.





####################
####################
### Table of significant co-deletions by fishers exact test








####################
####################
### Bar plot of Fishers exact test results for co-deletions
##bar_fishers_SignifCodel_codeletion

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


####################
####################
###Odds ratio plot for significant co-deletions
### Example code for odds ration forest map:
##http://www.surefoss.org/dataanalysis/plotting-odds-ratios-aka-a-forrestplot-with-ggplot2/







####################
####################
### Circle plots between cancers
###Plot: Circle plots of significant tumour suppressors on x axis and significant co-deleted genes 
#on y-axis. Circle size = number of cancers co-deletion is found in. Colour of circle = all (red), 
#unique (blue), not unique (black).
##Need long table with 4 columns, 1st column = target gene, second column = proximal gene, 
#third column = number of cancers co-deletion present, firth column = colour of circle catagory






####################
####################
### Table of significant co-deletions between Cancers by fishers exact test








####################
####################
### Bar plot of Fishers exact test results for co-deletions between cancers
##bar_fishers_SignifCodel_perCancer_tumour_suppressors

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








####################
####################
### Odds ratio figure for fishers exact test between cancers:
##http://www.surefoss.org/dataanalysis/plotting-odds-ratios-aka-a-forrestplot-with-ggplot2/





####################
####################
### Circle plots between cancers
###Plot: Circle plots of significant tumour suppressors on x axis and significant co-deleted genes 
#on y-axis. Circle size = number of cancers co-deletion is found in. Colour of circle = all (red), 
#unique (blue), not unique (black).










###########################
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









##############################################################################
####################
####################
### Bar plot of significant genes in overall survival analysis
## Survival analysis results analysis.R

## p<= 0.1

head(survival_stats_ovsurv_cat1_2_significant_table_p0_1_more_than_20)
dim(survival_stats_ovsurv_cat1_2_significant_table_p0_1_more_than_20)
table(survival_stats_ovsurv_cat1_2_significant_table_p0_1_more_than_20$target_gene)

bar_plot<- survival_stats_ovsurv_cat1_2_ALL %>%
  dplyr::group_by(target_gene) %>%
  dplyr::summarise(total = n())


bar_plot<- data.frame(cbind(bar_plot[,1], cbind(bar_plot[,2])))
bar_plot
colnames(bar_plot)<- c("target_gene", "values")

sum(as.numeric(bar_plot$values))

p <-ggplot(bar_plot, aes(x = factor(target_gene), y = as.numeric(values)))
p +geom_bar(stat = "identity") +
  xlab("Tummour Suppressor") + ylab("Number of Significant co-deletions") +
  ggtitle("Number of Significant Co-deletions per Tummour Suppressor") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)
  ) +
  ##Save plot
  ggsave("bar_Overall_survival_codeletions.tiff")








####################
####################
### Bar plot of significant genes in disease free survival analysis
## p<= 0.1

head(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20)
dim(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20)
table(survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20$target_gene)

bar_plot<- survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20 %>%
  dplyr::group_by(target_gene) %>%
  dplyr::summarise(total = n())


bar_plot<- data.frame(cbind(bar_plot[,1], cbind(bar_plot[,2])))
bar_plot
colnames(bar_plot)<- c("target_gene", "values")

sum(as.numeric(bar_plot$values))

p <-ggplot(bar_plot, aes(x = factor(target_gene), y = as.numeric(values)))
p +geom_bar(stat = "identity") +
  xlab("Tummour Suppressor") + ylab("Number of Significant co-deletions") +
  ggtitle("Number of Significant Co-deletions per Tummour Suppressor") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)
  ) +
  ##Save plot
  ggsave("bar_Disease_Free_survival_codeletions.tiff")




########################
###Survival curves














############################
###Venn diagrans to compare Fishers exact test, overall survival and disease free survival.
## Intersection of Tumour suppressors

venn.plot <- venn.diagram(
  x = list(
    fisher_codeletions = fishers_co_deletion_per_gene_long_table_significanct_p0.05_more_than_20$target_gene,
    overall_survival = survival_stats_ovsurv_cat1_2_significant_table_p0_1_more_than_20$target_gene,
    diseasefree_survival = survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20$target_gene
  ),
  euler.d = TRUE,
  scaled = TRUE,
  filename = "Euler_3set_scaled.tiff",
  cex = 2.5,
  cat.cex = 1.5,
  fill=c("white", "white", "white"),
  col= c("blue", "green", "purple"), 
  alpha=c(0.5,0.5, 0.5),
  cat.col = c("blue", "green", "purple"),
  cat.pos = c(0,0,180),
  cat.fontface=4,
  ext.text = TRUE,
  category.names=c("Fisher's Exact Test", "Overall Survival", "Disease Free Survival"),
  main="Significant Tumour suppressors")

###########
##Intersection of co-deletions

venn.plot <- venn.diagram(
  x = list(
    fisher_codeletions = fishers_co_deletion_per_gene_long_table_significanct_p0.05_more_than_20$proximal_gene,
    overall_survival = survival_stats_ovsurv_cat1_2_significant_table_p0_1_more_than_20$proximal_gene,
    diseasefree_survival = survival_stats_DiseFreeSurv_cancer_list_cat1_and_2_table_more_than_20$proximal_gene
  ),
  euler.d = TRUE,
  scaled = TRUE,
  filename = "Euler_3set_scaled.tiff",
  cex = 2.5,
  cat.cex = 1.5,
  fill=c("white", "white", "white"),
  col= c("blue", "green", "purple"), 
  alpha=c(0.5,0.5, 0.5),
  cat.col = c("blue", "green", "purple"),
  cat.pos = c(0,0,180),
  cat.fontface=4,
  ext.text = TRUE,
  category.names=c("Fisher's Exact Test", "Overall Survival", "Disease Free Survival"),
  main="Significant Co-deletions")


