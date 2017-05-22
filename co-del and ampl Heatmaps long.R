#############
## Co-deletion and amplification heatmaps
##############

##Install pacakges
source("https://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")

##Load packages
library(tidyr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

library(org.Hs.eg.db)


##Set working directory
getwd()
setwd("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/")
dir()

## Read 1st CNV file
## Find way to autamatically load each file
acc.cnv<- read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Data/unzipped data/ACC/all_data_by_genes.txt", stringsAsFactors = FALSE, header = TRUE)
dim(acc.cnv)
class(acc.cnv)
#acc.cnv<- tbl_df(acc.cnv)
#acc.cnv
View(acc.cnv)

## function to import all tables automatically and name object after folder name

#....




##Test all cancer CNV files have same gene names 

#....





##Combine all cancer CNV tables into one large dataframe

#.....



## Test to see how to transpose data
df<- data.frame(gene= c("w", "x", "y", "z"), ncbi = 11:14, patient1 = 1:4, patient2 = 2:5, stringsAsFactors = FALSE)
df
data.frame(t(df))

##Test to see how to transform CNV data to proportion of genes with deletion A that have deletion B i.e co-deletions:

##Create dummy CNV data:
A<- c(1, 1, 1, 0, -1)
B<- c(-1, 1, 0, 0, -1)
C<- c(-1, 1, -1, -1, -1)
D<- c(-1, 1, -1, 0, 1)
E<- c(0, 0, -1, -1, 1)
G<- c(1, 1, 0, -1, 0)

df<- data.frame(rbind(A,B,C,D,E,G))
class(df)
df

##Threshold data such that anything with a value of the threshold or lower = 0

df<=-1

df2<- df<=-1
df2

## Convert dataframe to 1 and 0:
df3<- df2*1
df3

rowSums(df3)
ncol(df3)

## Calculate the proportion of individuals with a deletion:
prop.individ.with.del<- rowSums(df3)/ncol(df3)

##loop to create a table with proportion of deletions shared between genes

#How to multiply each row?
#for loop
#apply
#vector based multiplication, i.e one dataframe with another?


df4<-data.frame(matrix(NA, ncol = nrow(df3), nrow = nrow(df3)))
df4


for(i in 1: nrow(df3)){
  
  for (j in 1: nrow(df3)) {
    
    value<- df3[i,]*df3[j,]
    df4[i,j]<- sum(value)/sum(df3[i,])
        
  }
}
df4

for(i in 1:nrow(df4)) {
  
  df4[i,i]<- prop.individ.with.del[i]
}
df4
rownames(df4)<- c("A", "B", "C", "D", "E", "F")
colnames(df4)<- c("A", "B", "C", "D", "E", "F")
df4

class(df4)
df5<- as.matrix(df4)
df5

##############
## Create heatmap with dummy data using pHeatmap

?pheatmap

RColorBrewer::display.brewer.all()

col.pal <- RColorBrewer::brewer.pal(9, "YlGnBu")
col.pal

pheatmap(df5, 
         cluster_row = T,
         cluster_cols = F,
         color = col.pal,
         fontsize = 6.5,
         fontsize_row=10, 
         fontsize_col = 10,
         clustering_distance_rows = "correlation",
         show_rownames = TRUE,
         show_colnames = TRUE
         )

pheatmap(df5, 
         cluster_row = F,
         cluster_cols = F,
         color = col.pal,
         fontsize = 6.5,
         fontsize_row=10, 
         fontsize_col = 10,
         clustering_distance_rows = "correlation",
         show_rownames = TRUE,
         show_colnames = TRUE
)


##########
## Obtain chromosomal locations of genes


keys<- keys(org.Hs.eg.db, keytype = "ENTREZID")
keys
columns<- c("CHR", "CHRLOC", "CHRLOCEND")
sel<- select(org.Hs.eg.db, keys, columns, keytype = "ENTREZID")
View(sel)
dim(sel)

acc.locus.id<- acc.cnv$Locus.ID
length(acc.locus.id)

sel2<- sel[sel$ENTREZID %in% acc.locus.id,]
sel2
dim(sel2)

sel3<- na.omit(sel2[!duplicated(sel2$ENTREZID), ])
sel3
dim(sel3)

sel3$strand<- ifelse(sel3$CHRLOC <0, "-", "+")
sel3
sel3$start<- abs(sel3$CHRLOC)
sel3$end<- abs(sel3$CHRLOCEND)
head(sel3)
dim(sel3)

#################
## Attach CHRLOCCHR and start to acc cnv file

head(acc.cnv)

?dplyr::rename
sel3<- dplyr::rename(sel3, Locus.ID = ENTREZID)

?dplyr::full_join
class(acc.cnv$Locus.ID)
class(sel3$Locus.ID)
sel3$Locus.ID<- as.integer(sel3$Locus.ID)
class(sel3$Locus.ID)

acc.cnv.loc<- full_join(sel3,acc.cnv, by="Locus.ID")
head(acc.cnv.loc)
class(acc.cnv.loc$CHR)
acc.cnv.loc$CHR<- as.integer(acc.cnv.loc$CHR)

is.na(acc.cnv.loc$CHR)
which(is.na(acc.cnv.loc), arr.ind = TRUE)
dim(acc.cnv.loc)
acc.cnv.loc<- dplyr::arrange(acc.cnv.loc, CHR, start)
View(acc.cnv.loc)

## comment: No Y chromosome...check intersection of genes on both X and Y and a lot of NA's!!


#################
##Function to convert CNV file to matrix for heatmap
#Very slow:
#See third function "create.heatmap.matrix.ampl.del.optimised" which is faster
create.heatmap.matrix<- function(x, column_start, threshold){

##Convert CNV data to matrix:

cnv.matrix<- as.matrix(x[,column_start:ncol(x)])

##Threshold data such that anything with a value of the threshold or lower = TRUE

cnv.matrix<- cnv.matrix <= threshold

## Convert dataframe to 1 and 0:
cnv.matrix<- cnv.matrix*1

## Calculate the proportion of individuals with a deletion:
prop.individ.with.del<- rowSums(cnv.matrix)/ncol(cnv.matrix)

##loop to create a table with proportion of deletions shared between genes


heatmap.matrix<-data.frame(matrix(NA, ncol = nrow(cnv.matrix), nrow = nrow(cnv.matrix)))

for(i in 1: nrow(cnv.matrix)){
  
  for (j in 1: nrow(cnv.matrix)) {
    
    value<- cnv.matrix[i,]*cnv.matrix[j,]
    heatmap.matrix[i,j]<- sum(value)/sum(cnv.matrix[i,])
    
  }
}

for(i in 1:nrow(heatmap.matrix)) {
  
  heatmap.matrix[i,i]<- prop.individ.with.del[i]
}

rownames(heatmap.matrix)<- rownames(x)
colnames(heatmap.matrix)<- rownames(x)

heatmap.matrix
}


##function that can create a matrix of proportion of co-amplifications aswell as deletions.
#Very slow:
create.heatmap.matrix.ampl.del<- function(x, column_start, threshold, deletion = TRUE){
  
  ##Convert CNV data to matrix:
  
  cnv.matrix<- as.matrix(x[,column_start:ncol(x)])
  
  if (deletion == TRUE) {
  ##Threshold data such that anything with a value of the threshold or lower = TRUE
  
  cnv.matrix<- cnv.matrix<= threshold
  } else {
    
    cnv.matrix<- cnv.matrix >= threshold
    
  }
  
  ## Convert dataframe to 1 and 0:
  cnv.matrix<- cnv.matrix*1
  
  ## Calculate the proportion of individuals with a deletion:
  prop.individ.with.del<- rowSums(cnv.matrix)/ncol(cnv.matrix)
  
  ##loop to create a table with proportion of deletions shared between genes
  
  
  heatmap.matrix<-data.frame(matrix(NA, ncol = nrow(cnv.matrix), nrow = nrow(cnv.matrix)))
  
  for(i in 1: nrow(cnv.matrix)){
    
    for (j in 1: nrow(cnv.matrix)) {
      
      value<- cnv.matrix[i,]*cnv.matrix[j,]
      heatmap.matrix[i,j]<- sum(value)/sum(cnv.matrix[i,])
      
    }
  }
  
  for(i in 1:nrow(heatmap.matrix)) {
    
    heatmap.matrix[i,i]<- prop.individ.with.del[i]
  }
  
  rownames(heatmap.matrix)<- rownames(x)
  colnames(heatmap.matrix)<- rownames(x)
  
  heatmap.matrix
}


##Faster function to create matrix of proportions of co-deletions
# and co-amplifications

create.heatmap.matrix.ampl.del.optimised<- function(x, column_start, threshold, deletion = TRUE){
  
  ##Convert CNV data to matrix:
  
  cnv.matrix<- as.matrix(x[,column_start:ncol(x)])
  
  ##Create a binary matrix of CNV data such that deletions (deletion = TRUE) or 
  #amplifications (deletion = FALSE) below (for deletions) or above (for amplifications) a threshold = 1 
  if (deletion == TRUE) {
    
    #cnv.matrix[cnv.matrix > threshold] = 0
    #cnv.matrix = -cnv.matrix
    cnv.matrix<- ifelse(cnv.matrix <= threshold, 1, 0)
    
  } else {
    
    #cnv.matrix<- cnv.matrix >= threshold
    #cnv.matrix[cnv.matrix < threshold] = 0
    cnv.matrix<- ifelse(cnv.matrix >= threshold, 1, 0)
  }
  
  ## Calculate the proportion of individuals with a deletion:
  prop.individ.with.del<- rowSums(cnv.matrix)/ncol(cnv.matrix)
  
  ## Calculate total number of co-deletions of amplifications
  
  heatmap.matrix<- cnv.matrix%*%t(cnv.matrix)
  heatmap.matrix<- heatmap.matrix/rowSums(cnv.matrix)
  
  ##NaN caused when rowSums = 0. Replace NaN with 0
  #cnv.matrix[cnv.matrix == NaN] = 0
  heatmap.matrix[is.nan(heatmap.matrix)] = 0
  
  ## Complete diagonals with proportion of individuals with deletions.
  for(i in 1:nrow(heatmap.matrix)) {
    
    heatmap.matrix[i,i]<- prop.individ.with.del[i]
  }
  
  ## Add row and column names
  rownames(heatmap.matrix)<- rownames(x)
  colnames(heatmap.matrix)<- rownames(x)
  
  heatmap.matrix
}







#####################
##Use function to create matrix for pHeatmap

x<- df
column_start<- 1
threshold<- -1
deletion = TRUE

heatmap.matrix<- create.heatmap.matrix(x, column_start, threshold)
heatmap.matrix

rownames(acc.cnv.loc)<- acc.cnv.loc$Gene.Symbol
x<- acc.cnv.loc
column_start<- 11
threshold<- -1

heatmap.matrix<- create.heatmap.matrix(x, column_start, threshold)
View(heatmap.matrix)

##New optimised function
rownames(acc.cnv.loc)<- acc.cnv.loc$Gene.Symbol
x<- acc.cnv.loc

system.time(heatmap.matrix.del<- create.heatmap.matrix.ampl.del.optimised(x, column_start = 11, threshold = -1, deletion = TRUE))
head(heatmap.matrix.del)
dim(heatmap.matrix.del)
View(heatmap.matrix.del)

system.time(heatmap.matrix.amp<- create.heatmap.matrix.ampl.del.optimised(x, column_start = 11, threshold = 1, deletion = FALSE))
head(heatmap.matrix.amp)
dim(heatmap.matrix.amp)
View(heatmap.matrix.amp)


#############
## Create Heatmap of ACC co-deletions data using pHeatmap
##Need to look up different clustering algorithms
##Problem with NA's?

small.heatmap.matrix.del<- heatmap.matrix.del[1:1000,1:1000]
class(small.heatmap.matrix.del)
dim(small.heatmap.matrix.del)
which(is.na(small.heatmap.matrix.del))

which(acc.cnv.loc$CHR == 9, arr.ind = TRUE)
chr9.heatmap.matrix.del<- heatmap.matrix.del[7809:8538,7809:8538]
dim(chr9.heatmap.matrix.del)

col.pal <- RColorBrewer::brewer.pal(9, "YlGnBu")
col.pal

pdf( "mygraph2.pdf", width = 12, height = 12 )
pheatmap(chr9.heatmap.matrix.del,
         cluster_row = F,
         cluster_cols = F,
         fontsize = 6.5,
         fontsize_row=0.5, 
         fontsize_col = 0.5,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cellwidth = 0.8,
         cellheight = 0.8 
         )
dev.off()

pheatmap(chr9.heatmap.matrix.del,
         cluster_row = T,
         cluster_cols = F,
         fontsize = 6.5,
         fontsize_row=0.3, 
         fontsize_col = 0.3,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cellwidth = 0.3,
         cellheight = 0.3 
)

chr9.heatmap.matrix.amp<- heatmap.matrix.amp[7809:8538,7809:8538]
dim(chr9.heatmap.matrix.del)

pheatmap(chr9.heatmap.matrix.amp,
         cluster_row = F,
         cluster_cols = F,
         fontsize = 6.5,
         fontsize_row=0.3, 
         fontsize_col = 0.3,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cellwidth = 0.3,
         cellheight = 0.3 
)

pdf( "mygraph.pdf", width = 12, height = 12 )
pheatmap(chr9.heatmap.matrix.amp,
         cluster_row = T,
         cluster_cols = F,
         fontsize = 6.5,
         fontsize_row=0.5, 
         fontsize_col = 0.5,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cellwidth = 0.8,
         cellheight = 0.8 
)
dev.off()

pheatmap(small.heatmap.matrix.del, 
         cluster_row = T,
         cluster_cols = F,
         color = col.pal,
         fontsize = 6.5,
         #fontsize_row=10, 
         #fontsize_col = 10,
         #clustering_distance_rows = "correlation",
         show_rownames = TRUE,
         show_colnames = TRUE
)

pheatmap(heatmap.matrix.del, 
         cluster_row = F,
         cluster_cols = F,
         color = col.pal,
         fontsize = 6.5,
         fontsize_row=10, 
         fontsize_col = 10,
         #clustering_distance_rows = "correlation",
         show_rownames = TRUE,
         show_colnames = TRUE
)
