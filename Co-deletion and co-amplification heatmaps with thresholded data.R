#######################
### Co-deletion and co-amplification heatmaps
#######################

#########################
### Heatmaps per chromosome
########################

##############
##Heatmaps showing co-deletions and co-amplifications for each chromosome for all tumour types combined.

threshold_CNV_all_table_loc[1:2, 1:12]

##Use a deletion threshold of -1
list.CNV.all.table.co.del<-lapply(c(seq(1:22), "X"), function(x) co.deletion_co.amplification_matrix(threshold_CNV_all_table_loc, column_start = 11, threshold = -1,Chromosome = x, 
                                                                                                     deletion = TRUE, normalisation = "tumours.with.event"))
length(list.CNV.all.table.co.del)
names(list.CNV.all.table.co.del)<- c(seq(1:22), "X")
names(list.CNV.all.table.co.del)
list.CNV.all.table.co.del[[1]][1:6, 1:12]
names(list.CNV.all.table.co.del[1])

##For loop to plot heatmaps for each chromosome (apply functions not working for me)
for (i in 1:length(list.CNV.all.table.co.del)){
  
  chromosome<- names(list.CNV.all.table.co.del[i])
  plot_data<- list.CNV.all.table.co.del[[i]]
  
  tiff(paste("Chromosome_", chromosome, "_deletion = true.tiff", sep =""), width = 25, height = 22, units = 'in', res = 100)
  pheatmap(plot_data,
           cluster_row = F,
           cluster_cols = F,
           show_rownames = FALSE,
           show_colnames = FALSE
  )
  dev.off()
  
  print(i)
}


# plot_heatmap<-function(x){
#   
#   tiff(paste("Chromosome_", x, "_deletion = true.tiff", sep =""), width = 25, height = 22, units = 'in', res = 100)
#   pheatmap(x,
#            cluster_row = F,
#            cluster_cols = F,
#            show_rownames = FALSE,
#            show_colnames = FALSE
#   )
#   dev.off()
#   
#   print(x)
#   
# }
# 
# lapply(names(list.CNV.all.table.co.del[9]), function(x) plot_heatmap(list.CNV.all.table.co.del[[x]]))
# 
# lapply(names(list.CNV.all.table.co.del), function(x) plot_heatmap(list.CNV.all.table.co.del[[x]]))
# 
# lapply(list.CNV.all.table.co.del[9], function(x) plot_heatmap(x))


##Use a deletion threshold of -2
list.CNV.all.table.co.del<-lapply(c(seq(1:22), "X"), function(x) co.deletion_co.amplification_matrix(threshold_CNV_all_table_loc, column_start = 11, threshold = -2,Chromosome = x, 
                                                                                                     deletion = TRUE, normalisation = "tumours.with.event"))
length(list.CNV.all.table.co.del)
names(list.CNV.all.table.co.del)<- c(seq(1:22), "X")
names(list.CNV.all.table.co.del)


##For loop to plot heatmaps for each chromosome (apply functions not working for me)
for (i in 1:length(list.CNV.all.table.co.del)){
  
  chromosome<- names(list.CNV.all.table.co.del[i])
  plot_data<- list.CNV.all.table.co.del[[i]]
  
  tiff(paste("Chromosome_", chromosome, "_deletion = true.tiff", sep =""), width = 25, height = 22, units = 'in', res = 100)
  pheatmap(plot_data,
           cluster_row = F,
           cluster_cols = F,
           show_rownames = FALSE,
           show_colnames = FALSE
  )
  dev.off()
  
  print(i)
}



##Use a amplification threshold of +2
list.CNV.all.table.co.amp<-lapply(c(seq(1:22), "X"), function(x) co.deletion_co.amplification_matrix(threshold_CNV_all_table_loc, 
                                                                                                     column_start = 11, threshold = 2,Chromosome = x, deletion = FALSE, normalisation = "tumours.with.event"))
names(list.CNV.all.table.co.amp)<- c(seq(1:22), "X")
names(list.CNV.all.table.co.amp)


##For loop to plot heatmaps for each chromosome (apply functions not working for me)
for (i in 1:length(list.CNV.all.table.co.amp)){
  
  chromosome<- names(list.CNV.all.table.co.amp[i])
  plot_data<- list.CNV.all.table.co.amp[[i]]
  
  tiff(paste("Chromosome_", chromosome, "_deletion = FALSE.tiff", sep =""), width = 25, height = 22, units = 'in', res = 100)
  pheatmap(plot_data,
           cluster_row = F,
           cluster_cols = F,
           show_rownames = FALSE,
           show_colnames = FALSE
  )
  dev.off()
  
  print(i)
}

##Comment: Maybe also do co-amplifications using a threshold of +1?

####################
###Identification of top co-deletion and co-amplification genes.
####################




#########################
### Co-deletions and co-amplifications per cytoband of target genes
########################



















#########################
### Co-deletions and co-amplifications for genes 2.5MB away from target gene.
########################

##Comment: Should I normalise against number of individuals with deletion in target gene?


