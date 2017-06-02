###############
### Target genes cytoband heat maps
###############

#############
###Get cytoband for CDKN2A (Use ACC data)
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
matrix<- co.deletion_co.amplification_matrix(cnv.table, column_start =  11, threshold = threshold, selection_criteria = cytoband, Cytoband = TRUE, deletion = TRUE)

################
###Make annotation bar for cytoband heatmap

## get cytoband coordinates for 9p21.3
cytoband.start.end<- cytoband.cordinates %>% 
  dplyr::filter(cytoband_name %in% cytoband) %>%
  dplyr::select(chromStart, chromEnd)

## Get start and end coordinates for genes in 9p21.3
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






##Plot heatmap:
pheatmap(matrix,
         cluster_row = T,
         cluster_cols = F,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 2,
         fontsize_col = 2,
         annotation_col = annotation.table,
         annotation_colors = ann_colors
)
}

#########
###Test function

Plot.target.genes.cytoband.heatmap(cnv.table = acc.cnv.chr.location, target.gene = "TP53", cytoband.cordinates = cytoband.cordinates,threshold = -1, deletion = TRUE)

##Comment: Normalise by number of individulas with deletions in cytoband?



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
###Run apply on fuction to plot heat maps automatically
########################


##Try:
tiff("Plot3.tiff", width = 4, height = 4, units = 'in', res = 300)
plot(x, y) # Make plot
dev.off()






































