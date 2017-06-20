#######################
### Co-deletion and co-amplification heatmaps with GISTIC thresholded data
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

### Parameters for function

## Specify Heatmap colours
ann_colors = list(
  strand = c("-" ="palegreen1", "+" ="skyblue1"),
  Distance = colorRampPalette(c("white", "firebrick3"))(100)
)

##Create input file of cytoband.cordinates
cytoband.cordinates<- read.delim("../../Input data/Annotations/Cytoband_start_end_hgTables", header = TRUE, stringsAsFactors = FALSE )

cytoband.cordinates.X.chrom<-cytoband.cordinates$X.chrom
cytoband.cordinates.X.chrom<- gsub(".*chr","",cytoband.cordinates.X.chrom)
cytoband.cordinates.X.chrom<- paste(cytoband.cordinates.X.chrom, cytoband.cordinates$name, sep = "")
cytoband.cordinates<-cbind(cytoband_name = cytoband.cordinates.X.chrom, cytoband.cordinates)
head(cytoband.cordinates)



Plot.target.genes.cytoband.heatmap<- function(cnv.table, target.gene,cytoband.cordinates, deletion = TRUE, threshold = -1){
  
  #############
  ###Get cytoband for cnv.table
  
  cytoband<- cnv.table %>% 
    dplyr::filter(Gene.Symbol %in% target.gene) %>%
    dplyr::select(Cytoband)
  
  cytoband<- cytoband[1,1]
  
  #############
  ###Get co-deletions for cytoband
  matrix<- co.deletion_co.amplification_matrix(cnv.table, column_start =  11, threshold = threshold, selection_criteria = cytoband, Cytoband = TRUE, deletion = deletion, normalisation = "tumours.with.event")
  
  #############
  ##Check matrix contains more than one value otherwise pheatmap wont plot heatmap
  
  #if (sum(matrix) == 0) {
  
  if (ncol(unique(matrix, MARGIN = 2)) ==1){ 
    
    print (paste("No heatmap was plotted for ", target.gene, "as all values in matrix were the same:", matrix[1,1]), sep = " ")
    
  }else{
    
    ################
    ###Make annotation bar for cytoband heatmap
    
    ## get cytoband coordinates for cytoband
    cytoband.start.end<- cytoband.cordinates %>% 
      dplyr::filter(cytoband_name %in% cytoband) %>%
      dplyr::select(chromStart, chromEnd)
    
    ## Get start and end coordinates for genes in cytoband
    ## Remove samples with NA
    matrix.gene.names<- colnames(matrix)
    
    matrix.gene.names.loc<- cnv.table %>%
      dplyr::filter(Gene.Symbol %in% matrix.gene.names) %>%
      dplyr::select(Gene.Symbol, strand, start, end) 
    
    ## calculate distance of start and end of genes to start and end of Cytoband
    ## Calculate which value is smallest
    matrix.gene.names.loc<- matrix.gene.names.loc %>% mutate(distance_from_start = start - cytoband.start.end[1,1],
                                                             distance_from_end = cytoband.start.end[1,2] - end,
                                                             minimum_distance = pmin(distance_from_start, distance_from_end) 
    )
    
    ##comment: Some genes cross cytoband boundary
    ##Give genes that cross cytoband boundary a value of 0
    #ifelse(matrix.gene.names.loc$minimum_distance < 0, 0, matrix.gene.names.loc$minimum_distance)
    matrix.gene.names.loc$minimum_distance[matrix.gene.names.loc$minimum_distance < 0]<- 0
    
    
    ##################
    ## Create annotation table for strand
    
    annotation.table<- as.data.frame(matrix.gene.names.loc$strand)
    rownames(annotation.table)<-  matrix.gene.names.loc$Gene.Symbol
    
    ## Create annotation table for Distance to cytoband edge
    annotation.table<- cbind(annotation.table, matrix.gene.names.loc$minimum_distance)
    
    colnames(annotation.table)<- c("strand", "Distance")
    
    
    ###############
    ###Plot heat maps
    
    ##colour scheme for heat map and annotation table?
    
    #create object to plot correct fontsize?
    
    ##Plot heatmap:
    tiff(paste(target.gene,"_deletion = ", deletion, "_threshold =", threshold, ".tiff", sep =""), width = 25, height = 22, units = 'in', res = 100)
    pheatmap(matrix,
             cluster_row = TRUE,
             cluster_cols = FALSE,
             breaks = NA,
             show_rownames = TRUE,
             show_colnames = TRUE,
             #fontsize_row = 2,
             #fontsize_col = 2,
             annotation_col = annotation.table,
             annotation_colors = ann_colors
    )
    dev.off()
    
    print(target.gene)
    
  }
}


#########
###Test function

cnv.table<- threshold_short_cnv_list_loc[[1]]

Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, target.gene = "TP53", cytoband.cordinates = cytoband.cordinates,threshold = 1, deletion = FALSE)
Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, target.gene = "CDKN2A", cytoband.cordinates = cytoband.cordinates, threshold = -1, deletion = TRUE)
Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, target.gene = "MET", cytoband.cordinates = cytoband.cordinates,threshold = -1, deletion = TRUE)
Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, target.gene = "STK11", cytoband.cordinates = cytoband.cordinates,threshold = -1, deletion = TRUE)


##########################
###Run apply function on co-deletion plotting function created above to plot heat maps automatically
########################

x<- c("MET", "CDKN2A", "RB1", "WWOX", 
      "LRP1B", "PDE4D", "CCNE1", "TP53",
      "FGFR1", "MYC", "EGFR","WHSC1L1",
      "ERBB2", "MCL1", "MDM2", "CCND1", "ATM",
      "NOTCH1", "PPP2R2A", "BRD4", "ARID1A",
      "STK11", "PARK2")

# #x<- c("CDKN2A", "RB1", "WWOX", 
#       "LRP1B", "PDE4D", "CCNE1", "TP53",
#       "FGFR1", "MYC", "EGFR","WHSC1L1",
#       "ERBB2", "MCL1", "MDM2", "CCND1", "ATM",
#       "NOTCH1", "PPP2R2A", "BRD4", "ARID1A",
#       "STK11", "PARK2")
# 
# #x<- c("MET")

# x<- c("MET", "CDKN2A", "RB1")
# 
# lapply(x, function(x) Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, 
#                                                          target.gene = x, 
#                                                          cytoband.cordinates = cytoband.cordinates,
#                                                          threshold = -1, 
#                                                          deletion = TRUE))
# 
# lapply(x, function(x) Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, 
#                                                          target.gene = x, 
#                                                          cytoband.cordinates = cytoband.cordinates,
#                                                          threshold = 1, 
#                                                          deletion = FALSE))
# 
# #################
# ### Perform co-amplification and co-deletion analysis for specific cancer types
# #Also see for loop furthur down
# ###############
# names(cnv.list)
# 
# 
# ##Breast:
# cnv.table<- chromosomal_location(cnv.list$BRCA)
# 
# setwd("../../Output/plots/170607 co-amp co-del (BRCA)/")
# 
# 
# lapply(x, function(x) Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, 
#                                                          target.gene = x, 
#                                                          cytoband.cordinates = cytoband.cordinates,
#                                                          threshold = -1, 
#                                                          deletion = TRUE))
# 
# lapply(x, function(x) Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, 
#                                                          target.gene = x, 
#                                                          cytoband.cordinates = cytoband.cordinates,
#                                                          threshold = 1, 
#                                                          deletion = FALSE))
# 
# 
# 

# cancer.type<- "BRCA"
# target.gene<- "CDKN2A"
# 
# Plot.target.genes.cytoband.heatmap.cancers.list<- function(cancer.type, target.gene){}
# 
# dir.create(paste("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Output/plots/", x, sep = ""))

################
### For loop to create directory and add plots for each cancer type I am most interested in. 
# Normalised to number of tumours with deletion or amplification:
################


##Tables to store heatmap plotting success 
results_del_table_1<- data.frame(matrix(NA, ncol = length(short.cnv.list), nrow = length(x)))
results_del_table_2<- data.frame(matrix(NA, ncol = length(short.cnv.list), nrow = length(x)))
results_amp_table_1<- data.frame(matrix(NA, ncol = length(short.cnv.list), nrow = length(x)))
results_amp_table_2<- data.frame(matrix(NA, ncol = length(short.cnv.list), nrow = length(x)))

##For loop to make new directory and save co-amplification and co-deletion plots in it
for (i in 1: length(threshold_short_cnv_list_loc)){
  
  tumour.type<- names(threshold_short_cnv_list_loc[i])
  dir.create(paste("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Output/plots/170620 co-amp co-del (", tumour.type, ")", sep = ""))
  setwd(paste("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Output/plots/170620 co-amp co-del (", tumour.type, ")", sep = ""))
  
  cnv.table<- threshold_short_cnv_list_loc[[i]]
  
  output1<- lapply(x, function(x) Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, 
                                                                     target.gene = x, 
                                                                     cytoband.cordinates = cytoband.cordinates,
                                                                     threshold = -1, 
                                                                     deletion = TRUE))
  
  output2<- lapply(x, function(x) Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, 
                                                                     target.gene = x, 
                                                                     cytoband.cordinates = cytoband.cordinates,
                                                                     threshold = -2, 
                                                                     deletion = TRUE))
  
  output3<- lapply(x, function(x) Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, 
                                                                     target.gene = x, 
                                                                     cytoband.cordinates = cytoband.cordinates,
                                                                     threshold = 2, 
                                                                     deletion = FALSE))
  
  output4<- lapply(x, function(x) Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, 
                                                                     target.gene = x, 
                                                                     cytoband.cordinates = cytoband.cordinates,
                                                                     threshold = 2, 
                                                                     deletion = FALSE))
  
  ## save results of wether heatmaps were plotted and if not what the value was:
  #Co-deletion with threshold -1
  output1<- gsub("No heatmap was plotted for", "All", output1)
  output1<- gsub("as all values in matrix were the same:", "values = ", output1)
  results_del_table_1[,i]<- output1
  #Co-deletion with threshold -2
  output2<- gsub("No heatmap was plotted for", "All", output2)
  output2<- gsub("as all values in matrix were the same:", "values = ", output2)
  results_del_table_2[,i]<- output2
  #Co-amplification with threshold 1
  output3<- gsub("No heatmap was plotted for", "All", output3)
  output3<- gsub("as all values in matrix were the same:", "values = ", output3)
  results_amp_table_1[,i]<- output3
  #Co-amplification with threshold 1
  output4<- gsub("No heatmap was plotted for", "All", output4)
  output4<- gsub("as all values in matrix were the same:", "values = ", output4)
  results_amp_table_2[,i]<- output4
  
}

## Tables of heatmap plotting success: Add row names and column names:
setwd("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Output/plots/")
rownames(results_del_table_1)<- x
colnames(results_del_table_1)<- names(threshold_short_cnv_list_loc)
rownames(results_del_table_2)<- x
colnames(results_del_table_2)<- names(threshold_short_cnv_list_loc)
rownames(results_amp_table_1)<- x
colnames(results_amp_table_1)<- names(threshold_short_cnv_list_loc)
rownames(results_amp_table_2)<- x
colnames(results_amp_table_2)<- names(threshold_short_cnv_list_loc)

results_del_table_1
results_del_table_2
results_amp_table_1
results_amp_table_2


## Save output files:
write.csv(results_del_table_1, file = "Deletion heatmap output threshold -1.csv")
write.csv(results_del_table_2, file = "Deletion heatmap output threshold -2.csv")

write.csv(results_amp_table_1, file = "Amplification heatmap output threshold +1.csv")
write.csv(results_amp_table_2, file = "Amplification heatmap output threshold +2.csv")







#########################
### Co-deletions and co-amplifications for genes 2.5MB away from target gene.
########################

##Comment: Should I normalise against number of individuals with deletion in target gene?

gene_information_list

##Colour scheme:
ann_colors = list(
  strand = c("-" ="palegreen1", "+" ="skyblue1"),
  Distance = colorRampPalette(c("firebrick3", "white"))(100)
)

##parameters:
distance<- 2.5e+06
cnv.table
x<- gene_information_list[[2]]
target.gene<- x[[1]]
deletion = TRUE
threshold = -1
start<- TRUE
column_start = 11
Chromosome = x[[2]]
selection_criteria = c(x[[4]] - distance, x[[5]] + distance)
# selection_criteria = c(x[[1]][[4]], x[[1]][[5]])
normalisation = "total.tumour.number"

x[[1]]
x[1]

plot_genes_surrounding_target_genes_heatmap<- function(cnv.table, target.gene, column_start =  11, distance, deletion = TRUE, threshold = -1, start = TRUE, Chromosome, selection_criteria, normalisation = "tumours.with.event"){
  
  #############
  ###Get cytoband for cnv.table
  
  # cytoband<- cnv.table %>% 
  #   dplyr::filter(Gene.Symbol %in% target.gene) %>%
  #   dplyr::select(Cytoband)
  # 
  # cytoband<- cytoband[1,1]
  
  ###########
  ###Get genomic interval
  new_selection_criteria = c(x[[4]] - distance, x[[5]] + distance)
  
  #############
  ###Get co-deletions/co-amplifications for geneomic interval
  matrix<- co.deletion_co.amplification_matrix(cnv.table, column_start =  column_start, threshold = threshold, start = TRUE, Chromosome = Chromosome, selection_criteria = new_selection_criteria,  deletion = deletion, normalisation = normalisation)
  
  #############
  ##Check matrix contains more than one value otherwise pheatmap wont plot heatmap
  
  if (ncol(unique(matrix, MARGIN = 2)) ==1){ 
    
    print (paste("No heatmap was plotted for ", target.gene, "as all values in matrix were the same:", matrix[1,1]), sep = " ")
    
  }else{
    
    ################
    ###Make annotation bar for cytoband heatmap
    
    ## get cytoband coordinates for cytoband
    # cytoband.start.end<- cytoband.cordinates %>% 
    #   dplyr::filter(cytoband_name %in% cytoband) %>%
    #   dplyr::select(chromStart, chromEnd)
    
    ##Calcule distance of each gene from target gene
    distance_from_target_gene_dataframe<- distance_from_target_gene_function(cnv.table = cnv.table, x = x, distance = distance)
    
    
    ## Get Gene.Symbol, strand, start and end coordinates for genes in cytoband
    matrix.gene.names<- colnames(matrix)
    
    matrix.gene.names.loc<- cnv.table %>%
      dplyr::filter(Gene.Symbol %in% matrix.gene.names) %>%
      dplyr::select(Gene.Symbol, strand, start, end) 
    
    ## calculate distance of start and end of genes to start and end of Cytoband
    ## Calculate which value is smallest
    # matrix.gene.names.loc<- matrix.gene.names.loc %>% mutate(distance_from_start = start - cytoband.start.end[1,1],
    #                                                          distance_from_end = cytoband.start.end[1,2] - end,
    #                                                          minimum_distance = pmin(distance_from_start, distance_from_end) 
    # )
    
    ##comment: Some genes cross cytoband boundary
    ##Give genes that cross cytoband boundary a value of 0
    #ifelse(matrix.gene.names.loc$minimum_distance < 0, 0, matrix.gene.names.loc$minimum_distance)
    # matrix.gene.names.loc$minimum_distance[matrix.gene.names.loc$minimum_distance < 0]<- 0
    
    
    ##################
    ## Create annotation table for strand
    
    annotation_table<- as.data.frame(matrix.gene.names.loc$strand)
    rownames(annotation_table)<-  matrix.gene.names.loc$Gene.Symbol
    
    ## Create annotation table for Distance to cytoband edge
    annotation_table<- cbind(annotation_table, distance_from_target_gene_dataframe)
    
    colnames(annotation_table)<- c("strand", "Distance")
    
    
    ###############
    ###Plot heat maps
    
    ##colour scheme for heat map and annotation table?
    
    #create object to plot correct fontsize?
    
    ##Plot heatmap:
    tiff(paste(target.gene,"_deletion = ", deletion, "_threshold =", threshold,"_distance = ",distance, ".tiff", sep =""), width = 25, height = 22, units = 'in', res = 100)
    pheatmap(matrix,
             cluster_row = TRUE,
             cluster_cols = FALSE,
             breaks = NA,
             show_rownames = TRUE,
             show_colnames = TRUE,
             #fontsize_row = 2,
             #fontsize_col = 2,
             annotation_col = annotation_table,
             annotation_colors = ann_colors
    )
    dev.off()
    
    print(target.gene)
    
  }
}


#########
###Test function

cnv.table<- threshold_short_cnv_list_loc[[1]]
x<-gene_information_list[[10]]
x
# distance<-2.5e+06

plot_genes_surrounding_target_genes_heatmap(cnv.table, target.gene = x[[1]], distance = 2.5e+06, deletion = TRUE, threshold = -1, start = TRUE, Chromosome = x[[2]], selection_criteria = c(x[[4]], x[[5]]) , normalisation = "tumours.with.event")

#########
### Test ploting multiple heatmaps for target genes using apply function:

#test<- lapply(gene_information_list ,function(x) plot_genes_surrounding_target_genes_heatmap(cnv.table, target.gene = x[[1]], distance = 2.5e+06, deletion = TRUE, threshold = -1, start = TRUE, Chromosome = x[[2]], selection_criteria = c(x[[4]], x[[5]]) , normalisation = "tumours.with.event"))
##Comment: can not get lapply to work on plot_genes_surrounding_target_genes_heatmap so use for loop:


###For loop to plot heatmaps for different target genes:

for (i in 1: length(gene_information_list)){
  x<- gene_information_list[[i]]
  plot_genes_surrounding_target_genes_heatmap(cnv.table, target.gene = x[[1]], distance = 2.5e+06, deletion = TRUE, threshold = -1, start = TRUE, Chromosome = x[[2]], selection_criteria = c(x[[4]]-distance, x[[5]]+distance) , normalisation = "tumours.with.event")
  
}


###############
###For loop to plot heat maps for different cancers for genes 2.5MB 5' and 3' of target gene.

x<- c("MET", "CDKN2A", "RB1", "WWOX", 
      "LRP1B", "PDE4D", "CCNE1", "TP53",
      "FGFR1", "MYC", "EGFR","WHSC1L1",
      "ERBB2", "MCL1", "MDM2", "CCND1", "ATM",
      "NOTCH1", "PPP2R2A", "BRD4", "ARID1A",
      "STK11", "PARK2")



##Tables to store heatmap plotting success 
results_del_table_1<- data.frame(matrix(NA, ncol = length(short.cnv.list), nrow = length(x)))
results_del_table_2<- data.frame(matrix(NA, ncol = length(short.cnv.list), nrow = length(x)))
results_amp_table_1<- data.frame(matrix(NA, ncol = length(short.cnv.list), nrow = length(x)))
results_amp_table_2<- data.frame(matrix(NA, ncol = length(short.cnv.list), nrow = length(x)))

##Vector to tempoarily store output of heatmap plotting success  
output1<- rep(NA, length(gene_information_list))
output2<- rep(NA, length(gene_information_list))
output3<- rep(NA, length(gene_information_list))
output4<- rep(NA, length(gene_information_list))

##For loop to make new directory and save co-amplification and co-deletion plots in it
for (i in 1: length(threshold_short_cnv_list_loc)){
  
  tumour.type<- names(threshold_short_cnv_list_loc[i])
  dir.create(paste("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Output/plots/170620 co-amp co-del (", tumour.type, ")", sep = ""))
  setwd(paste("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Output/plots/170620 co-amp co-del (", tumour.type, ")", sep = ""))
  
  cnv.table<- threshold_short_cnv_list_loc[[i]]
  
  for (j in 1: length(gene_information_list)){
    x<- gene_information_list[[j]]
    output1[j]<- plot_genes_surrounding_target_genes_heatmap(cnv.table, target.gene = x[[1]], distance = 2.5e+06, deletion = TRUE, threshold = -1, start = TRUE, Chromosome = x[[2]], selection_criteria = c(x[[4]]-distance, x[[5]]+distance) , normalisation = "tumours.with.event")
  }
  for (j in 1: length(gene_information_list)){
    x<- gene_information_list[[j]]
    output2[j]<- plot_genes_surrounding_target_genes_heatmap(cnv.table, target.gene = x[[1]], distance = 2.5e+06, deletion = TRUE, threshold = -2, start = TRUE, Chromosome = x[[2]], selection_criteria = c(x[[4]]-distance, x[[5]]+distance) , normalisation = "tumours.with.event")
  }
  for (j in 1: length(gene_information_list)){
    x<- gene_information_list[[j]]
    output3[j]<- plot_genes_surrounding_target_genes_heatmap(cnv.table, target.gene = x[[1]], distance = 2.5e+06, deletion = FALSE, threshold = 1, start = TRUE, Chromosome = x[[2]], selection_criteria = c(x[[4]]-distance, x[[5]]+distance) , normalisation = "tumours.with.event")
  }
  for (j in 1: length(gene_information_list)){
    x<- gene_information_list[[j]]
    output4[j]<- plot_genes_surrounding_target_genes_heatmap(cnv.table, target.gene = x[[1]], distance = 2.5e+06, deletion = FALSE, threshold = 2, start = TRUE, Chromosome = x[[2]], selection_criteria = c(x[[4]]-distance, x[[5]]+distance) , normalisation = "tumours.with.event")
  }
  
 
  ## save results of wether heatmaps were plotted and if not what the value was:
  #Co-deletion with threshold -1
  output1<- gsub("No heatmap was plotted for", "All", output1)
  output1<- gsub("as all values in matrix were the same:", "values = ", output1)
  results_del_table_1[,i]<- output1
  #Co-deletion with threshold -2
  output2<- gsub("No heatmap was plotted for", "All", output2)
  output2<- gsub("as all values in matrix were the same:", "values = ", output2)
  results_del_table_2[,i]<- output2
  #Co-amplification with threshold 1
  output3<- gsub("No heatmap was plotted for", "All", output3)
  output3<- gsub("as all values in matrix were the same:", "values = ", output3)
  results_amp_table_1[,i]<- output3
  #Co-amplification with threshold 1
  output4<- gsub("No heatmap was plotted for", "All", output4)
  output4<- gsub("as all values in matrix were the same:", "values = ", output4)
  results_amp_table_2[,i]<- output4
  
}

## Tables of heatmap plotting success: Add row names and column names:
setwd("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Output/plots/")
x<- c("MET", "CDKN2A", "RB1", "WWOX", 
      "LRP1B", "PDE4D", "CCNE1", "TP53",
      "FGFR1", "MYC", "EGFR","WHSC1L1",
      "ERBB2", "MCL1", "MDM2", "CCND1", "ATM",
      "NOTCH1", "PPP2R2A", "BRD4", "ARID1A",
      "STK11", "PARK2")
rownames(results_del_table_1)<- x
colnames(results_del_table_1)<- names(threshold_short_cnv_list_loc)
rownames(results_del_table_2)<- x
colnames(results_del_table_2)<- names(threshold_short_cnv_list_loc)
rownames(results_amp_table_1)<- x
colnames(results_amp_table_1)<- names(threshold_short_cnv_list_loc)
rownames(results_amp_table_2)<- x
colnames(results_amp_table_2)<- names(threshold_short_cnv_list_loc)

results_del_table_1
results_del_table_2
results_amp_table_1
results_amp_table_2


## Save output files:
write.csv(results_del_table_1, file = "Deletion heatmap output threshold -1.csv")
write.csv(results_del_table_2, file = "Deletion heatmap output threshold -2.csv")

write.csv(results_amp_table_1, file = "Amplification heatmap output threshold +1.csv")
write.csv(results_amp_table_2, file = "Amplification heatmap output threshold +2.csv")


###############
###For loop to plot heat maps for different cancers for genes 5MB 5' and 3' of target gene.

x<- c("MET", "CDKN2A", "RB1", "WWOX", 
      "LRP1B", "PDE4D", "CCNE1", "TP53",
      "FGFR1", "MYC", "EGFR","WHSC1L1",
      "ERBB2", "MCL1", "MDM2", "CCND1", "ATM",
      "NOTCH1", "PPP2R2A", "BRD4", "ARID1A",
      "STK11", "PARK2")



##Tables to store heatmap plotting success 
results_del_table_1<- data.frame(matrix(NA, ncol = length(short.cnv.list), nrow = length(x)))
results_del_table_2<- data.frame(matrix(NA, ncol = length(short.cnv.list), nrow = length(x)))
results_amp_table_1<- data.frame(matrix(NA, ncol = length(short.cnv.list), nrow = length(x)))
results_amp_table_2<- data.frame(matrix(NA, ncol = length(short.cnv.list), nrow = length(x)))

##Vector to tempoarily store output of heatmap plotting success  
output1<- rep(NA, length(gene_information_list))
output2<- rep(NA, length(gene_information_list))
output3<- rep(NA, length(gene_information_list))
output4<- rep(NA, length(gene_information_list))

##For loop to make new directory and save co-amplification and co-deletion plots in it
for (i in 1: length(threshold_short_cnv_list_loc)){
  
  tumour.type<- names(threshold_short_cnv_list_loc[i])
  dir.create(paste("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Output/plots/170620 co-amp co-del_5MB (", tumour.type, ")", sep = ""))
  setwd(paste("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Output/plots/170620 co-amp co-del_5MB (", tumour.type, ")", sep = ""))
  
  cnv.table<- threshold_short_cnv_list_loc[[i]]
  
  for (j in 1: length(gene_information_list)){
    x<- gene_information_list[[j]]
    output1[j]<- plot_genes_surrounding_target_genes_heatmap(cnv.table, target.gene = x[[1]], distance = 5e+06, deletion = TRUE, threshold = -1, start = TRUE, Chromosome = x[[2]], selection_criteria = c(x[[4]]-distance, x[[5]]+distance) , normalisation = "tumours.with.event")
  }
  for (j in 1: length(gene_information_list)){
    x<- gene_information_list[[j]]
    output2[j]<- plot_genes_surrounding_target_genes_heatmap(cnv.table, target.gene = x[[1]], distance = 5e+06, deletion = TRUE, threshold = -2, start = TRUE, Chromosome = x[[2]], selection_criteria = c(x[[4]]-distance, x[[5]]+distance) , normalisation = "tumours.with.event")
  }
  for (j in 1: length(gene_information_list)){
    x<- gene_information_list[[j]]
    output3[j]<- plot_genes_surrounding_target_genes_heatmap(cnv.table, target.gene = x[[1]], distance = 5e+06, deletion = FALSE, threshold = 1, start = TRUE, Chromosome = x[[2]], selection_criteria = c(x[[4]]-distance, x[[5]]+distance) , normalisation = "tumours.with.event")
  }
  for (j in 1: length(gene_information_list)){
    x<- gene_information_list[[j]]
    output4[j]<- plot_genes_surrounding_target_genes_heatmap(cnv.table, target.gene = x[[1]], distance = 5e+06, deletion = FALSE, threshold = 2, start = TRUE, Chromosome = x[[2]], selection_criteria = c(x[[4]]-distance, x[[5]]+distance) , normalisation = "tumours.with.event")
  }
  
  
  ## save results of wether heatmaps were plotted and if not what the value was:
  #Co-deletion with threshold -1
  output1<- gsub("No heatmap was plotted for", "All", output1)
  output1<- gsub("as all values in matrix were the same:", "values = ", output1)
  results_del_table_1[,i]<- output1
  #Co-deletion with threshold -2
  output2<- gsub("No heatmap was plotted for", "All", output2)
  output2<- gsub("as all values in matrix were the same:", "values = ", output2)
  results_del_table_2[,i]<- output2
  #Co-amplification with threshold 1
  output3<- gsub("No heatmap was plotted for", "All", output3)
  output3<- gsub("as all values in matrix were the same:", "values = ", output3)
  results_amp_table_1[,i]<- output3
  #Co-amplification with threshold 1
  output4<- gsub("No heatmap was plotted for", "All", output4)
  output4<- gsub("as all values in matrix were the same:", "values = ", output4)
  results_amp_table_2[,i]<- output4
  
}

## Tables of heatmap plotting success: Add row names and column names:
setwd("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Output/plots/")
rownames(results_del_table_1)<- x
colnames(results_del_table_1)<- names(threshold_short_cnv_list_loc)
rownames(results_del_table_2)<- x
colnames(results_del_table_2)<- names(threshold_short_cnv_list_loc)
rownames(results_amp_table_1)<- x
colnames(results_amp_table_1)<- names(threshold_short_cnv_list_loc)
rownames(results_amp_table_2)<- x
colnames(results_amp_table_2)<- names(threshold_short_cnv_list_loc)

results_del_table_1
results_del_table_2
results_amp_table_1
results_amp_table_2


## Save output files:
write.csv(results_del_table_1, file = "Deletion heatmap output threshold -1.csv")
write.csv(results_del_table_2, file = "Deletion heatmap output threshold -2.csv")

write.csv(results_amp_table_1, file = "Amplification heatmap output threshold +1.csv")
write.csv(results_amp_table_2, file = "Amplification heatmap output threshold +2.csv")





