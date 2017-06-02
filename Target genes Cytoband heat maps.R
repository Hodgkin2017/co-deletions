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
##Maybe modify function to remove NAs
acc.cnv.chr.location.filt<- acc.cnv.chr.location %>%
  dplyr::filter(!is.na(start))

dim(acc.cnv.chr.location.filt)

#############
###Get co-deletions for CDKN2A locus (9p21.3)
acc.cnv.9p21.3.co.del<- co.deletion_co.amplification_matrix(acc.cnv.chr.location.filt, column_start =  11, threshold = -1, selection_criteria = "9p21.3", Cytoband = TRUE, deletion = TRUE)
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

##Comment: several genes with NA's!!!....remove
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
class(annotation.table.distance$minimum_distance)

##################
### Create annotation table for strand

annotation.table.strand<- as.data.frame(gene.names.9p21.3.loc.dist$strand)
rownames(annotation.table.strand)<-  gene.names.9p21.3.loc.dist$Gene.Symbol
colnames(annotation.table.strand)<- "strand"
annotation.table.strand
class(annotation.table.strand$`gene.names.9p21.3.loc.dist$strand`)

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

pheatmap(acc.cnv.9p21.3.co.del,
         cluster_row = T,
         cluster_cols = F,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 5,
         fontsize_col = 5,
         annotation_col = annotation.table.strand)
)

colour<- colorRampPalette(c("firebrick3", "white"))(100)
colour

# Specify colors
ann_colors = list(
  strand = c("-" ="palegreen1", "+" ="skyblue1"),
  Distance = colorRampPalette(c("white", "firebrick3"))(100)
)





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
































