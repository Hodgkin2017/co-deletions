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
co_deletions_distance_from_target_gene_plot_table<- readRDS(file = "temp.data.rds")
dim(co_deletions_distance_from_target_gene_plot_table)
head(co_deletions_distance_from_target_gene_plot_table)
tail(co_deletions_distance_from_target_gene_plot_table)
class(co_deletions_distance_from_target_gene_plot_table$proportion_co_del_amp)

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


ggplot(co_deletions_distance_from_target_gene_plot_table[1:100,], 
       aes(Comparison_gene, as.numeric(proportion_co_del_amp), colour = Target_gene)) +
  geom_point() +
  #geom_area(aes(fill=Target_gene)) +
  geom_line(stat = "identity", aes(group = Target_gene, fill = Target_gene)) +
  facet_wrap(~Target_gene, nrow = 1, scales = "free_x") + 
  theme(legend.position="none") +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))

ggplot(co_deletions_distance_from_target_gene_plot_table[1:100,], 
       aes(Comparison_gene, as.numeric(proportion_co_del_amp), colour = Target_gene)) +
  geom_area(mapping = aes(x = ), fill=Target_gene)

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
dummydata2<-data.frame(Comparison_gene = as.factor(letters[1:8]), Target_gene = rep(c("b", "g"), c(4,4)), proportion_co_del_amp = c(1.5, 2.2, 3.7, 5.2, 3.2, 1.2, 4.3, 3.9))
dummydata2

ggplot(dummydata2, aes(Comparison_gene, as.numeric(proportion_co_del_amp))) +
  geom_bar(stat = "identity", aes(fill = Target_gene)) +
  facet_wrap(~Target_gene, scales = "free_x")

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
mean_co_deletion_co_amplification_values_around_gene(dummydata, 1)
mean_co_deletion_co_amplification_values_around_gene(dummydata, 2)


