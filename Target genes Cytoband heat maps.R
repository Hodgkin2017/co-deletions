###############
### Target genes cytoband heat maps
###############

#############
###Get cytoband for CDKN2A (Use ACC data)
acc.cnv.chr.location<- chromosomal_location(cnv.list[[1]])

acc.cnv.chr.location[1:2, 1:12]
acc.cnv.chr.location %>% 
  dplyr::filter(Gene.Symbol == "CDKN2A") %>%
  dplyr::select(Cytoband)

#################
##Obtain genes of interest (in Cytoband 9p21.3) with known start and end sites
##Comment: Maybe modify function to remove NAs
##Comment: Now default in function to remove genes with no start:

# acc.cnv.chr.location.filt<- acc.cnv.chr.location %>%
#   dplyr::filter(!is.na(start))
# 
# dim(acc.cnv.chr.location.filt)

#############
###Get co-deletions for CDKN2A locus (9p21.3)
acc.cnv.9p21.3.co.del<- co.deletion_co.amplification_matrix(acc.cnv.chr.location, column_start =  11, threshold = -1, selection_criteria = "9p21.3", Cytoband = TRUE, deletion = TRUE)
dim(acc.cnv.9p21.3.co.del)
head(acc.cnv.9p21.3.co.del)


#############
###Normalise for number of patients with CDKN2A deletions/number of patients with deletions in locus


##Talk to Christophe.....adjust function



################
###Make annotation bar for cytoband heatmap

##Read in Cytoband cordinates 
cytoband.cordinates<- read.delim("../../Input data/Annotations/Cytoband_start_end_hgTables", header = TRUE, stringsAsFactors = FALSE )
dim(cytoband.cordinates)
head(cytoband.cordinates)

## get cytoband coordinates for 9p21.3
cytoband.cordinates %>% 
  dplyr::filter(X.chrom == "chr9" & name == "p21.3") %>%
  dplyr::select(chromStart, chromEnd)

## Get start and end coordinates for genes in 9p21.3
## Remove samples with NA
gene.names.9p21.3<- colnames(acc.cnv.9p21.3.co.del)
length(gene.names.9p21.3)
gene.names.9p21.3
gene.names.9p21.3.loc<- acc.cnv.chr.location %>%
  dplyr::filter(Gene.Symbol %in% gene.names.9p21.3) %>%
  dplyr::select(Gene.Symbol, strand, start, end) 

#%>% filter(!is.na(start))

##Comment: several genes with NA's!!!....removed now.
gene.names.9p21.3.loc

## calculate distance of start and end of genes to start and end of Cytoband
## Calculate which value is smallest
gene.names.9p21.3.loc.dist<- gene.names.9p21.3.loc %>% mutate(distance_from_start = start - 19900000,
                                 distance_from_end = 25600000 - end,
                                 minimum_distance = pmin(distance_from_start, distance_from_end) 
                                 )

gene.names.9p21.3.loc.dist

##Create table with this value in for each gene (this is the annotation table for distance)

annotation.table.distance<- gene.names.9p21.3.loc.dist %>% 
  dplyr::select(minimum_distance)
class(annotation.table.distance)
dim(annotation.table.distance)

rownames(annotation.table.distance)<-  gene.names.9p21.3.loc.dist$Gene.Symbol
colnames(annotation.table.distance)<-"Distance"

annotation.table.distance
class(annotation.table.distance$Distance)

##################
### Create annotation table for strand

annotation.table.strand<- as.data.frame(gene.names.9p21.3.loc.dist$strand)
rownames(annotation.table.strand)<-  gene.names.9p21.3.loc.dist$Gene.Symbol
colnames(annotation.table.strand)<- "strand"
annotation.table.strand
class(annotation.table.strand$strand)

#################
##Combine annotations
annotation.table<- cbind(annotation.table.strand, annotation.table.distance)
annotation.table


################
### Remove genes from heatmap without start and end coordinates...done before analysis

# genes.to.keep<-as.data.frame(gene.names.9p21.3.loc$Gene.Symbol)
# #rownames(genes.to.keep)<- gene.names.9p21.3.loc$Gene.Symbol
# #
# 
# acc.cnv.9p21.3.co.del
# 
# # test<-dplyr::intersect(acc.cnv.9p21.3.co.del, genes.to.keep)
# # test
# 
# test<- dplyr::right_join(acc.cnv.9p21.3.co.del, genes.to.keep, by = "x1")



###############
###Plot heat map...cluster on rows. Columns = genomic co-ordinates

##colour scheme for heat map and annotation table?

## Specify colors
ann_colors = list(
  strand = c("-" ="palegreen1", "+" ="skyblue1"),
  Distance = colorRampPalette(c("white", "firebrick3"))(100)
)




##Plot heatmap:
pheatmap(acc.cnv.9p21.3.co.del,
         cluster_row = T,
         cluster_cols = F,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 5,
         fontsize_col = 5,
         annotation_col = annotation.table,
         annotation_colors = ann_colors
)
















#####################
### Make fuction which creates heatmap of cytobands for target genes: Co-deletions:
#####################

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




### temp parameters:
cnv.table<- chromosomal_location(cnv.list[[1]])
deletion = TRUE
target.gene = "TP53"
target.gene = "CDKN2A"
target.gene = "ARID1A"
target.gene = "MET"
threshold = -1



Plot.target.genes.cytoband.heatmap<- function(cnv.table, target.gene,cytoband.cordinates, deletion = TRUE, threshold = -1){

#############
###Get cytoband for cnv.table

cytoband<- cnv.table %>% 
  dplyr::filter(Gene.Symbol %in% target.gene) %>%
  dplyr::select(Cytoband)

cytoband<- cytoband[1,1]

#############
###Get co-deletions for cytoband
matrix<- co.deletion_co.amplification_matrix(cnv.table, column_start =  11, threshold = threshold, selection_criteria = cytoband, Cytoband = TRUE, deletion = TRUE, normalisation = "tumours.with.event")

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
tiff(paste(target.gene,"_deletion = ", deletion, ".tiff", sep =""), width = 25, height = 22, units = 'in', res = 100)
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

Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, target.gene = "TP53", cytoband.cordinates = cytoband.cordinates,threshold = -1, deletion = TRUE)
Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, target.gene = "CDKN2A", cytoband.cordinates = cytoband.cordinates, threshold = -1, deletion = TRUE)
Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, target.gene = "MET", cytoband.cordinates = cytoband.cordinates,threshold = -1, deletion = TRUE)
Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, target.gene = "STK11", cytoband.cordinates = cytoband.cordinates,threshold = -1, deletion = TRUE)

##Comment: Normalise by number of individulas with deletions in cytoband?

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

lapply(x, function(x) Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, 
                                                         target.gene = x, 
                                                         cytoband.cordinates = cytoband.cordinates,
                                                         threshold = -1, 
                                                         deletion = TRUE))

lapply(x, function(x) Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, 
                                                         target.gene = x, 
                                                         cytoband.cordinates = cytoband.cordinates,
                                                         threshold = 1, 
                                                         deletion = FALSE))

#################
### Perform co-amplification and co-deletion analysis for specific cancer types
#Also see for loop furthur down
###############
names(cnv.list)


##Breast:
cnv.table<- chromosomal_location(cnv.list$BRCA)

setwd("../../Output/plots/170607 co-amp co-del (BRCA)/")


lapply(x, function(x) Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, 
                                                         target.gene = x, 
                                                         cytoband.cordinates = cytoband.cordinates,
                                                         threshold = -1, 
                                                         deletion = TRUE))

lapply(x, function(x) Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, 
                                                         target.gene = x, 
                                                         cytoband.cordinates = cytoband.cordinates,
                                                         threshold = 1, 
                                                         deletion = FALSE))




# cancer.type<- "BRCA"
# target.gene<- "CDKN2A"
# 
# Plot.target.genes.cytoband.heatmap.cancers.list<- function(cancer.type, target.gene){}
# 
# dir.create(paste("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Output/plots/", x, sep = ""))

################
### For loop to create directory and add plots for each cancer type I am most interested in
################

##Target genes:
x<- c("MET", "CDKN2A", "RB1", "WWOX", 
      "LRP1B", "PDE4D", "CCNE1", "TP53",
      "FGFR1", "MYC", "EGFR","WHSC1L1",
      "ERBB2", "MCL1", "MDM2", "CCND1", "ATM",
      "NOTCH1", "PPP2R2A", "BRD4", "ARID1A",
      "STK11", "PARK2")

##New cnv list of cancer types we are most interested in:
short.cnv.list<- cnv.list[c(3, 7, 9, 12, 20, 21, 23, 24, 26, 29, 30)]
length(short.cnv.list)
names(short.cnv.list)

# -        Breast invasive carcinoma (BRCA)
# -        Esophageal cancer (ESCA)
# -        Head and neck squamous cell carcinoma (HNSC)
# -        Lung adenocarcinoma (LUAD)
# -        Lung squamous cell carcinoma (LUSC)
# -        Ovarian serous cystadenocarcinoma (OV)
# -        Pancreatic ductal adenocarcinoma (PAAD)
# -        Stomach adenocarcinoma (STAD)
# -        Skin cutaneous melanoma (SKCM)
# -        Prostate adenocarcinoma (PRAD)
# -        Colorectal adenocarcinoma (COADREAD)


##For loop to make new directory and save co-amplification and co-deletion plots in it
results.table<- data.frame(matrix(NA, ncol = length(short.cnv.list), nrow = 2*length(x)))

for (i in 1: length(short.cnv.list)){
  
  tumour.type<- names(short.cnv.list[i])
  dir.create(paste("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Output/plots/170615 co-amp co-del (", tumour.type, ") (Thresh = 2)", sep = ""))
  setwd(paste("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Output/plots/170615 co-amp co-del (", tumour.type, ") (Thresh = 2)", sep = ""))
  
  cnv.table<- chromosomal_location(short.cnv.list[[i]])
  
  lapply(x, function(x) Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, 
                                                           target.gene = x, 
                                                           cytoband.cordinates = cytoband.cordinates,
                                                           threshold = -2, 
                                                           deletion = TRUE))
  
  lapply(x, function(x) Plot.target.genes.cytoband.heatmap(cnv.table = cnv.table, 
                                                           target.gene = x, 
                                                           cytoband.cordinates = cytoband.cordinates,
                                                           threshold = 2, 
                                                           deletion = FALSE))
  
  ## save results of wether heatmaps were plotted and if not what the value was....how?
  #results.table[,i]<- 
}






##Comment: TP53 co-deletion heatmap has few TP53 mutations...not correct?! (17p13.1)
# Nope it is ok!!!!! was accidently counting genomic coordinates!
# test1<- cnv.table %>% dplyr::filter(Gene.Symbol == "TP53") 
# dim(test1)
# sum(test1 < -1)
# which(test1 < -1, arr.ind = T)
# test2<- cnv.table %>% dplyr::filter(Gene.Symbol == "TXNDC17") 
# sum(test2 < -1)
# 2/90
# ncol(cnv.table)
# cnv.table[1,1:12]
# 5/90
# matrix<- co.deletion_co.amplification_matrix(cnv.table, column_start =  11, threshold = -1, selection_criteria = "17p13.1", Cytoband = TRUE, deletion = TRUE)
# matrix[1,1:12]
# df<- as.data.frame(matrix)
# df$TP53
# 
# co.deletion_co.amplification_matrix(df, column_start =  1, threshold = -1, deletion = TRUE)

##########################
###Test saving plots
########################


##Try:
# tiff("Plot3.tiff", width = 4, height = 4, units = 'in', res = 300)
# plot(x, y) # Make plot
# dev.off()
# 
# 
# ##Plot heatmap:
# x<- "CDKN2A"
# 
# tiff("CDKN2A.tiff", width = 5, height = 5, units = 'in', res = 300)
# 
# 
# tiff(paste(x, ".tiff", sep =""), width = 5, height = 5, units = 'in', res = 300)
# pheatmap(matrix,
#          cluster_row = T,
#          cluster_cols = F,
#          show_rownames = TRUE,
#          show_colnames = TRUE,
#          fontsize_row = 2,
#          fontsize_col = 2,
#          annotation_col = annotation.table,
#          annotation_colors = ann_colors,
#          
# )
# dev.off()



































