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








###########################
### heatmap_amplification_perMB







#########################
###heatmap example of chromosome 9







####################
### Example code for odds ration forest map:
##http://www.surefoss.org/dataanalysis/plotting-odds-ratios-aka-a-forrestplot-with-ggplot2/






