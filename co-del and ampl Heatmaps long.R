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
setwd("/Users/Matt/Documents/co-deletions/")
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
## Create heatmap with pHeatmap

?pheatmap



#########
## My heatmap so far.....
########
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

acc.cnv.loc<- full_join(acc.cnv, sel3, by="Locus.ID")
head(acc.cnv.loc)
class(acc.cnv.loc$CHR)
acc.cnv.loc$CHR<- as.integer(acc.cnv.loc$CHR)

is.na(acc.cnv.loc$CHR)
which(is.na(acc.cnv.loc), arr.ind = TRUE)
dim(acc.cnv.loc)
acc.cnv.loc<- dplyr::arrange(acc.cnv.loc, CHR, start)
View(acc.cnv.loc)

## comment: No Y chromosomes and a lot of NA's!!


#################
##Function to convert CNV file to matrix for heatmap

