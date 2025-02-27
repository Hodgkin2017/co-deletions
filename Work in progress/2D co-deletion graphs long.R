##################
### 2D co-deletion graphs
#################

##Make a line plot of distance from target gene and proportion of deletions using ggplot2.
##Draw a second line of average co-deletion score of itself and surrounding genes and plot this too.
##Facet wrap in such a way that genes will be on the horizontal axis and cancer along the vertical axis.

################
### Get distance from target gene and proportion of deletions from 
#Scatter plot- co-deletions co-amplifications against distance between genes.R

## Comment: Make scatter plot method into  function so it is easier to run and change attributes e.g. distance and normalisation
#Make it have two outputs a list of just target gene co-deletions and a list of all v's all target gene co-deletions

##Get Tumours suppressors that have significant co-deletions in their loci.
# length(unique(fishers_co_deletion_per_gene_long_table_significanct_p0.05_more_than_20$target_gene))
# unique(fishers_co_deletion_per_gene_long_table_significanct_p0.05_more_than_20$target_gene)
# target_gene_names<- sapply(gene_information_long_list, function(x) x[[1]])
# target_gene_names
# target_gene_index<- which(target_gene_names %in% unique(fishers_co_deletion_per_gene_long_table_significanct_p0.05_more_than_20$target_gene))
# length(target_gene_index)
# target_gene_index<- target_gene_index[c(1:30,32:34, 36:49)]
# target_gene_index
# length(target_gene_index)

grep("ZFHX3" ,sapply(gene_information_long_list, function(x) x[[1]]))
gene_information_list_significant2<-gene_information_long_list[c(142, 157, 60, 166, 76, 168)]
cnv.table<- threshold_selected_cnv_list_plus_all_loc$ALL
cnv.table<- threshold_selected_cnv_list_plus_all_loc$BRCA


# test<- sapply(gene_information_long_list[target_gene_index], function(x) x[[1]])
# test2<- unique(fishers_co_deletion_per_gene_long_table_significanct_p0.05_more_than_20$target_gene)
# test<- test[order(test)]
# test2<- test2[order(test2)]
# test
# test2
# sapply(gene_information_long_list, function(x) x[[1]])
# identical(gene_information_long_list[[101]], gene_information_long_list[[103]])
# identical(gene_information_long_list[[93]], gene_information_long_list[[125]])





# gene_information_long_significant_list<- gene_information_long_list[target_gene_index]
# length(gene_information_long_significant_list)
# 
# #cnv.table<-threshold_short_cnv_list_loc[[1]]
# cnv.table<-threshold_selected_cnv_list_plus_all_loc$ALL

#co_deletions_distance_from_target_gene_plot_table3<- distance_from_target_gene_co_deletion_co_amplification_function(cnv.table = cnv.table, gene_information_list = gene_information_list, distance = 2.5e+06, deletion = TRUE, threshold = -2, compare_all_genes = FALSE, normalisation = "tumours.with.event")
co_deletions_distance_from_target_gene_plot_table<- distance_from_target_gene_co_deletion_co_amplification_function(cnv.table = cnv.table, gene_information_list = gene_information_list_significant2, distance = 2.5e+06, deletion = TRUE, threshold = -2, compare_all_genes = FALSE, normalisation = "tumours.with.event")

# cnv.table<-threshold_selected_cnv_list_plus_all_loc$ACC
# gene_information_list<- gene_information_long_significant_list
# co_deletions_distance_from_target_gene_plot_table2<- distance_from_target_gene_co_deletion_co_amplification_function(cnv.table = cnv.table, gene_information_list = gene_information_list, distance = 2.5e+06, deletion = TRUE, threshold = -2, compare_all_genes = FALSE, normalisation = "tumours.with.event")
# dim(co_deletions_distance_from_target_gene_plot_table2)

dim(co_deletions_distance_from_target_gene_plot_table)
head(co_deletions_distance_from_target_gene_plot_table)
tail(co_deletions_distance_from_target_gene_plot_table)
class(co_deletions_distance_from_target_gene_plot_table$proportion_co_del_amp)
unique(co_deletions_distance_from_target_gene_plot_table$Target_gene)

#saveRDS(co_deletions_distance_from_target_gene_plot_table, file = "temp.data.rds")

###########
###Plot data:
ggplot(co_deletions_distance_from_target_gene_plot_table[1:100,], aes(Comparison_gene, as.numeric(proportion_co_del_amp))) +
  geom_bar(stat = "identity", aes(fill = Target_gene)) +
  facet_wrap(~Target_gene, nrow = 1, scales = "free_x") + 
  theme(legend.position="none") +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))


ggplot(co_deletions_distance_from_target_gene_plot_table[1:100,], 
       aes(Comparison_gene, as.numeric(proportion_co_del_amp), colour = Target_gene)) +
  geom_line(stat = "identity", aes(group = Target_gene)) +
  facet_wrap(~Target_gene, nrow = 1, scales = "free_x") + 
  theme(legend.position="none") +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))


# ggplot(co_deletions_distance_from_target_gene_plot_table[1:100,], 
#        aes(Comparison_gene, as.numeric(proportion_co_del_amp), colour = Target_gene)) +
#   geom_point(size = 0.5, shape = 1) +
#   #geom_area(aes(fill=Target_gene)) +
#   geom_line(stat = "identity", aes(group = Target_gene)) +
#   facet_wrap(~Target_gene, nrow = 1, scales = "free_x") + 
#   theme(legend.position="none") +
#   theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))

# ggplot(co_deletions_distance_from_target_gene_plot_table[1:100,], 
#        aes(Comparison_gene, as.numeric(proportion_co_del_amp), colour = Target_gene)) +
#   geom_area(mapping = aes(x = ), fill=Target_gene)

########### From: https://stackoverflow.com/questions/30343494/ggplot2-geom-area-with-factorial-x-axis
# library(directlabels)
# ggplot(foo, aes(x=as.numeric(factor(age)), y = n, fill=diag)) +
#   geom_area() +
#   scale_x_discrete(labels = levels(foo$age)) + 
#   geom_dl(aes(label = diag), list("top.points", cex = .6)) + 
#   guides(fill = FALSE)
# 
# ownload.file("https://dl.dropboxusercontent.com/u/12226044/admissions.Rdata", 
#              destfile = fn <- file.path(tempdir(), "admissions.Rdata"), 
#              mode = "wb")
# load(fn)
# library(ggplot2)
# library(directlabels)
# library(dplyr)
# shaped %>%
#   group_by(diag, age) %>%
#   summarise(n = mean(n)) %>%
#   ggplot(aes(x=as.numeric(factor(age)), y = n, fill=diag)) +
#   geom_area(position = "stack") +
#   scale_x_discrete(labels = levels(shaped$age), expand = c(.1, .1)) + 
#   geom_dl(aes(label = diag, color = diag), position = "stack", list("last.bumpup", rot=-30, cex = .5)) + 
#   guides(fill = FALSE, colour = FALSE)

#########
# dummydata2<-data.frame(Comparison_gene = as.factor(letters[1:8]), Target_gene = rep(c("b", "g"), c(4,4)), proportion_co_del_amp = c(1.5, 2.2, 3.7, 5.2, 3.2, 1.2, 4.3, 3.9))
# dummydata2
# 
# ggplot(dummydata2, aes(Comparison_gene, as.numeric(proportion_co_del_amp))) +
#   geom_bar(stat = "identity", aes(fill = Target_gene)) +
#   facet_wrap(~Target_gene, scales = "free_x")

###################
### Calculate average co-deletion in neighbouring genes

##Practice:
# dummydata<- matrix(base::sample(seq(1:20), 64, replace = TRUE), ncol= 8, nrow = 8)
# colnames(dummydata)<- letters[1:8]
# rownames(dummydata)<- letters[1:8]
# dummydata
# 
# i=3
# n=2
# for(i in 1: nrow(dummydata)){
#   
#   if(i < n){
#   start<- i
#   end<- i+n
#   print(dummydata[i, start:end])
#   
#   } else if((i >= n) & (i+n-1 < ncol(dummydata))) {
#     start<- i-n
#     end<- i+n
#     print(dummydata[i, start:end])
#     
#   } else if(i+n > ncol(dummydata)) {
#     
#     start<- i-n
#     end<- ncol(dummydata)
#     print(dummydata[i, start:end])
#     
#   }
# }

##Function to get mean co-deletion of genes +/- n genes away from current gene:

mean_co_deletion_co_amplification_values_around_gene<- function(co_deletion_table, 
                                                                distance_from_gene_to_calculate_mean){
  
  n<- distance_from_gene_to_calculate_mean
  
  result<-rep(NA, nrow(co_deletion_table))
  
  for(i in 1: nrow(co_deletion_table)){
    
    if(i < n){
      start<- i
      end<- i+n
      result[i]<- mean(co_deletion_table[i, start:end])
      
    } else if((i >= n) & (i+n-1 < ncol(co_deletion_table))) {
      start<- i-n
      end<- i+n
      result[i]<- mean(co_deletion_table[i, start:end])
      
    } else if(i+n > ncol(co_deletion_table)) {
      
      start<- i-n
      end<- ncol(co_deletion_table)
      result[i]<- mean(co_deletion_table[i, start:end])
      
    }
  }
  
  return(result)
  
}

##Test function:
# mean_co_deletion_co_amplification_values_around_gene(dummydata, 1)
# mean_co_deletion_co_amplification_values_around_gene(dummydata, 2)

###############
### Get average co-deletion of neighbouring genes 2.5MB around target genes then join to table 
#with co-deletions and distance per target gene

##Create list with one dataframe per target gene of co-deletions
distance<- 2.5e+06
co.deletion.per.target.gene<- lapply(gene_information_list_significant2, function(x) co.deletion_co.amplification_matrix(cnv.table, column_start = 11, threshold = -2, start = TRUE, Chromosome = x[[2]], selection_criteria = c(x[[4]] - distance, x[[5]] + distance), deletion = TRUE, normalisation = "tumours.with.event"))
length(co.deletion.per.target.gene)
co.deletion.per.target.gene[[2]]

##Calculate average co-deletion for target gene +/- one gene.
test<- mean_co_deletion_co_amplification_values_around_gene(co.deletion.per.target.gene[[2]], 1)
test
length(test)
dim(co.deletion.per.target.gene[[2]])

co_deletion_around_target_gene<- lapply(co.deletion.per.target.gene, function(x) mean_co_deletion_co_amplification_values_around_gene(x,1))
length(co_deletion_around_target_gene)
co_deletion_around_target_gene[[2]]

identical(test, co_deletion_around_target_gene[[2]])

##Calculate average co-deletion for target gene +/- two genes.
co_deletion_around_target_gene2<- lapply(co.deletion.per.target.gene, function(x) mean_co_deletion_co_amplification_values_around_gene(x,2))
length(co_deletion_around_target_gene2)
co_deletion_around_target_gene2[[2]]

##Join each item of co.deletion.around.target.gene list into a one column dataframe:
co_deletion_around_target_gene<- unlist(co_deletion_around_target_gene)
co_deletion_around_target_gene

co_deletion_around_target_gene2<- unlist(co_deletion_around_target_gene2)
co_deletion_around_target_gene2

##Join co.deletion.around.target.gene to original dataframe with per gene proportion of co-deletions
#with target gene
length(co_deletion_around_target_gene)
dim(co_deletions_distance_from_target_gene_plot_table)
co_deletions_distance_from_target_gene_plot_table<- cbind(co_deletions_distance_from_target_gene_plot_table, co_deletion_around_target_gene)
dim(co_deletions_distance_from_target_gene_plot_table)
head(co_deletions_distance_from_target_gene_plot_table)
tail(co_deletions_distance_from_target_gene_plot_table)

co_deletions_distance_from_target_gene_plot_table<- cbind(co_deletions_distance_from_target_gene_plot_table, co_deletion_around_target_gene2)


##Plot proportion of co-deletions with target gene and with surounding genes

# ggplot(co_deletions_distance_from_target_gene_plot_table[1:100,], 
#        aes(Comparison_gene, as.numeric(proportion_co_del_amp), colour = Target_gene)) +
#   geom_point(size = 0.5, shape = 1) +
#   geom_line(stat = "identity", aes(group = Target_gene)) +
#   geom_point(aes(y= co_deletion_around_target_gene), size = 0.5, shape = 1) +
#   geom_line(aes(y= co_deletion_around_target_gene), stat = "identity", aes(group = Target_gene)) +
#   facet_wrap(~Target_gene, nrow = 1, scales = "free_x") + 
#   theme(legend.position="none") +
#   theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))

##Plot neighbouring/local co-deletions one gene away
# ggplot(co_deletions_distance_from_target_gene_plot_table, aes(x=Comparison_gene, group = Target_gene)) + 
#   geom_point(aes(y = as.numeric(proportion_co_del_amp), colour = "Co-deletion with target gene"), size = 0.5, shape = 1) +
#   geom_line(aes(y = as.numeric(proportion_co_del_amp), colour = "Co-deletion with target gene")) + 
#   geom_point(aes(y = as.numeric(co_deletion_around_target_gene), colour = "Mean co-deletion for genes 1 gene away"), size = 0.5, shape = 1) +
#   geom_line(aes(y = as.numeric(co_deletion_around_target_gene), colour = "Mean co-deletion for genes 1 gene away")) +
#   #facet_wrap(~Target_gene, nrow = 1, scales = "free_x") + 
#   facet_wrap(~Target_gene, scales = "free_x") + 
#   #theme(legend.position="none") +
#   theme(legend.position="bottom") +
#   theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
# 
# ##Save plot
# ggsave("BRCA_co-deletion_distance_1geneaway_lineplot.tiff",width = 16, height = 12, dpi = 100)

##Plot neighbouring/local co-deletions two genes away
ggplot(co_deletions_distance_from_target_gene_plot_table, aes(x=factor(Comparison_gene, levels = co_deletions_distance_from_target_gene_plot_table$Comparison_gene), group = Target_gene)) + 
  geom_point(aes(y = as.numeric(proportion_co_del_amp), colour = "Co-deletion with target gene"), size = 0.5, shape = 1) +
  geom_line(aes(y = as.numeric(proportion_co_del_amp), colour = "Co-deletion with target gene")) + 
  geom_point(aes(y = as.numeric(co_deletion_around_target_gene2), colour = "Mean co-deletion for genes 2 gene away"), size = 0.5, shape = 1) +
  geom_line(aes(y = as.numeric(co_deletion_around_target_gene2), colour = "Mean co-deletion for genes 2 gene away")) +
  #facet_wrap(~Target_gene, nrow = 1, scales = "free_x") + 
  facet_wrap(~Target_gene, scales = "free_x",  ncol = 2) + 
  #theme(legend.position="none") +
  theme(legend.position="bottom") +
  xlab("Gene") +
  ylab("Proportion of co-deletions in tumours with any deletion within loci") +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5, size = 5)) 
#+ geom_vline(xintercept = Target_gene)


##Save plot
ggsave("ALL_co-deletion_distance_2genesaway_lineplot2.tiff",width = 14, height = 10, dpi = 300)


##Plot neighbouring/local co-deletions one and two genes away
# ggplot(co_deletions_distance_from_target_gene_plot_table, aes(x=Comparison_gene, group = Target_gene)) + 
#   geom_point(aes(y = as.numeric(proportion_co_del_amp), colour = "Co-deletion with target gene"), size = 0.5, shape = 1) +
#   geom_line(aes(y = as.numeric(proportion_co_del_amp), colour = "Co-deletion with target gene")) +
#   geom_point(aes(y = as.numeric(co_deletion_around_target_gene), colour = "Mean co-deletion for genes 1 gene away"), size = 0.5, shape = 1) +
#   geom_line(aes(y = as.numeric(co_deletion_around_target_gene), colour = "Mean co-deletion for genes 1 gene away")) +
#   geom_point(aes(y = as.numeric(co_deletion_around_target_gene2), colour = "Mean co-deletion for genes 2 gene away"), size = 0.5, shape = 1) +
#   geom_line(aes(y = as.numeric(co_deletion_around_target_gene2), colour = "Mean co-deletion for genes 2 gene away")) +
#   #facet_wrap(~Target_gene, nrow = 1, scales = "free_x") + 
#   facet_wrap(~Target_gene, scales = "free_x") + 
#   #theme(legend.position="none") +
#   theme(legend.position="bottom") +
#   theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
# 
# ##Save plot
# ggsave("BRCA_co-deletion_distance_1and2genesaway_lineplot.tiff",width = 16, height = 12, dpi = 100)
