#############
## Measure the number of deletions or amplifications per cytoband
##############

library(tidyr)
library(dplyr)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)

##dummy data

# dummy.data<- acc.cnv[1:100,]
# head(dummy.data)

##Parameters
# object_name<- dummy.data
# object_name<- acc.cnv
threshold<- -1
#column_gene_ID<- 2
column_data_start<- 11
cytoband_column<-10
#chromosome_column<- 4
#gene_location_column<-5
deletion<- TRUE
chromosome_interval<- 0
chromosome_interval<- 10
select_chromosome<- 1


acc.cnv.chr.location<- chromosomal_location(acc.cnv)
object_name<- acc.cnv.chr.location

events.per.cytoband<- function(object_name, threshold, cytoband_column, column_data_start, chromosome_interval = 0, select_chromosome, deletion = TRUE){
##########
## Obtain chromosomal locations of genes

# keys<- keys(org.Hs.eg.db, keytype = "ENTREZID")
# columns<- c("CHR", "CHRLOC", "CHRLOCEND")
# all.gene.locations<- AnnotationDbi::select(org.Hs.eg.db, keys, columns, keytype = "ENTREZID")
# 
# locus.id<- object_name$Locus.ID
# 
# genes.of.interest<- all.gene.locations[all.gene.locations$ENTREZID %in% locus.id,]
# 
# genes.of.interest<- na.omit(genes.of.interest[!duplicated(genes.of.interest$ENTREZID), ])
# 
# ## Remove genes with poorly defined chromosome
# chr<-c(seq(1:22), "X", "Y")
# 
# rows.of.interest<- which(genes.of.interest$CHR %in% chr)
# 
# genes.of.interest<- genes.of.interest[rows.of.interest,]
# 
# genes.of.interest$strand<- ifelse(genes.of.interest$CHRLOC <0, "-", "+")
# genes.of.interest$start<- abs(genes.of.interest$CHRLOC)
# genes.of.interest$end<- abs(genes.of.interest$CHRLOCEND)
# 
# 
# 
# ##Join gene location dataframe to original CNV data table
# 
# genes.of.interest<- dplyr::rename(genes.of.interest, Locus.ID = ENTREZID)
# 
# genes.of.interest$Locus.ID<- as.integer(genes.of.interest$Locus.ID)
# 
# genes.of.interest<- full_join(genes.of.interest,object_name, by="Locus.ID")
# 
# ##Order genes by chromosome and location:
# genes.of.interest<- dplyr::arrange(genes.of.interest, CHR, start)
# 
# column_data_start<- column_data_start + 7
# cytoband_column<- cytoband_column + 7

## Code to get proportion of deletions per cytoband
cnv.matrix<- as.matrix(object_name[,column_data_start:ncol(object_name)])

if (deletion == TRUE) {
  
  cnv.matrix<- ifelse(cnv.matrix <= threshold, 1, 0)
  
} else {

  cnv.matrix<- ifelse(cnv.matrix >= threshold, 1, 0)
}

##Create list to store data:
if (chromosome_interval == 0){
results <- vector("list", 2)
} else {
  results <- vector("list", 4)
}

##Create matrix with raw deletions/amplification values and number of 
#deletions/amplifications per gene
results[[1]]<- cnv.matrix %>% 
  as.data.frame() %>% 
  dplyr::mutate(sum.of.deletions = rowSums(.)) %>%
  dplyr::mutate(number.of.tumours = ncol(.)-1) %>%
  cbind(cytoband=object_name[,cytoband_column], 
        chromosome = object_name$CHR, 
        start = object_name$start, .) 

##Create dataframe with number of deletions per gene, number of genes and number of 
#potenatial deletion events that could occur
results[[2]]<- results[[1]]%>%
  group_by(cytoband) %>%
  dplyr::summarise(sum.of.genes.deleted = sum(sum.of.deletions),
                   total.number.of.genes=n(),
                   total.number.of.events = sum(number.of.tumours)) %>%
  dplyr::mutate(proportion.of.deletions = sum.of.genes.deleted/total.number.of.events) %>%
  tidyr::separate(cytoband, c("chromosome", "band"), sep = "[p:q]", remove = FALSE, convert = TRUE)
#dplyr::select(-band) %>%

results[[2]]$chromosome<-sub("X", "23", results[[2]]$chromosome)
results[[2]]$chromosome<-sub("Y", "24", results[[2]]$chromosome)
results[[2]]$chromosome<- as.integer(results[[2]]$chromosome)
results[[2]]<- dplyr::arrange(results[[2]], chromosome, cytoband)
results[[2]]$chromosome<-sub("23", "X", results[[2]]$chromosome)
results[[2]]$chromosome<-sub("24", "Y", results[[2]]$chromosome)

## Find row with highest value:
#test3 %>%
  #arrange(desc(proportion.of.deletions)) %>%
  #head()

## Finding chromosome names:
#grep("^1", test$cytoband)
#grepl("^1", test$cytoband)
#dplyr::filter(test, cytoband "^1"==TRUE)

#dplyr::filter(test, grepl("^1", test$cytoband))tid
#rowSums(cnv.matrix)

#############
#### Code to get proportion of deletions per unit size of chromosome:

#seq.test<- seq(1:100)

#table(cut(seq.test, 10, labels = seq(1:10)))
#table(cut(seq.test, chromosome_interval, labels = seq(1:chromosome_interval)))
#dplyr::ntile(seq.test, 10)


## Remove rows with no known start site:?

##Remove genes without chromosomal locations:
#genes.of.interest<- genes.of.interest[!is.na(genes.of.interest$CHRLOC)]
#which(is.na(genes.of.interest), arr.ind = TRUE)

if (chromosome_interval > 0){
## Very rough estimate of lengths of chromosomes:
results[[3]]<- results[[1]] %>%
  group_by(chromosome) %>%
  summarise(chromosome_start = min(start),
            chromosome_end = max(start)) %>%
  mutate(estimated_chromosome_length = chromosome_end - chromosome_start,
         intervals_for_kb = estimated_chromosome_length/1000,
         intervals_for_10kb = estimated_chromosome_length/10000,
         intervals_for_100kb = estimated_chromosome_length/100000,
         intervals_for_Mb = estimated_chromosome_length/1000000,
         intervals_for_10Mb = estimated_chromosome_length/10000000)

  filter.table<- results[[1]] %>%
    dplyr::filter(chromosome == select_chromosome)

  results[[4]]<- filter.table$start %>%
    cut(chromosome_interval, labels = seq(1:chromosome_interval)) %>%
    cbind(Intervals = ., filter.table) %>%
    dplyr::group_by(Intervals) %>%
    dplyr::summarise(sum.of.genes.deleted = sum(sum.of.deletions),
                     total.number.of.genes=n(),
                     total.number.of.potential.events = sum(number.of.tumours)) %>%
    dplyr::mutate(proportion.of.deletions = sum.of.genes.deleted/total.number.of.potential.events)
}

return(results)
}
    
## Test function:

# test<- events.per.cytoband(object_name, threshold, cytoband_column, column_data_start, chromosome_interval = 0,  deletion = TRUE)
# glimpse(test[[1]])
# glimpse(test[[2]])

test1<- events.per.cytoband(acc.cnv.chr.location, threshold = -1, cytoband_column = 10, column_data_start = 11, chromosome_interval = 0,  deletion = TRUE)
glimpse(test1[[1]])
glimpse(test1[[2]])
glimpse(test1[[3]])
glimpse(test1[[4]])

test2<- events.per.cytoband(acc.cnv.chr.location, threshold = -1, cytoband_column = 10, column_data_start = 11 , select_chromosome = 1, chromosome_interval = 10,  deletion = TRUE)
glimpse(test2[[1]])
glimpse(test2[[2]])
glimpse(test2[[3]])
glimpse(test2[[4]])

####################
### Heatmap of the number of deletions per cytoband per chromosome

test1[[2]]$proportion.of.deletions

my.matrix<- as.matrix(test1[[2]]$proportion.of.deletions)
rownames(my.matrix)<-test1[[2]]$cytoband
colnames(my.matrix)<- "ACC"

annotation_row<- data.frame(chromosome = test1[[2]]$chromosome)
rownames(annotation_row)<- test1[[2]]$cytoband
dim(annotation_row)
annotation_row[,1]<- as.character(annotation_row[,1])
class(annotation_row$chromosome)

pheatmap(my.matrix,
         cluster_row = F,
         cluster_cols = F,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row=1,
         cellwidth = 10,
         annotation_row = annotation_row,
         annotation_legend = FALSE
         )

##Have side heatmap that contains chromosome

################
#### Create matrix of proportion of deletions/amplifications per cytoband per cancer type and plot heatmap.
##################

################
##Import all cancer data

## Choose directory and name of files to import from each directory:

# x<-"/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/unzipped original broad TCGA CNV data"
# file.to.import<-"all_data_by_genes.txt"

##O1: Run function to obtain a list object containing all CNV dataframes

# cnv.list<- import.files.from.directories(x, file.to.import)

#################
###Deletion heatmap:

##apply that creates matrix

CNV.data<-cnv.list

CNV.data[[38]]<- CNV.all.table

lapply(1:3, function(x) c(x, x^2, x^3))

new.function<- function(x){
  x %>% +1 %>% *2
}

lapply(1:3, function(x) new.function(x))

###########
##Loop that creates matrix

cancer.type<- names(cnv.list)
cancer.type

heatmap.matrix.cytoband.del<- matrix(NA, ncol = length(cancer.type)+1, nrow = 806)
dim(heatmap.matrix.cytoband.del)

for (i in 1:length(cancer.type)){
  
  x<-cnv.list[[i]]
  x<- dplyr::full_join(acc.cnv.chr.location[,1:8],x, by = "Locus.ID")
  cytoband.list<- events.per.cytoband(x, threshold = -1, cytoband_column = 10, column_data_start = 11, chromosome_interval = 0,  deletion = TRUE)
  heatmap.matrix.cytoband.del[,i]<- cytoband.list[[2]]$proportion.of.deletions
  print(cancer.type[i])
}

head(heatmap.matrix.cytoband.del)

##Final column contains all cancer types:

##Create one large dataframe with all CNV data in it:
#x<-join.cnv.datasets(cnv.list, 4)

##Calculate proportion of deletions per cytoband and add to matrix
x<- dplyr::full_join(acc.cnv.chr.location[,1:8],CNV.all.table, by = "Locus.ID")
cytoband.list<- events.per.cytoband(x, threshold = -1, cytoband_column = 10, column_data_start = 11, chromosome_interval = 0,  deletion = TRUE)
heatmap.matrix.cytoband.del[,38]<- cytoband.list[[2]]$proportion.of.deletions

head(heatmap.matrix.cytoband.del)
class(heatmap.matrix.cytoband.del)
dim(x)

heatmap.matrix.cytoband.del[200:250,4:5]

##Add row and column names

rownames(heatmap.matrix.cytoband.del)<- cytoband.list[[2]]$cytoband
colnames(heatmap.matrix.cytoband.del)<- c(cancer.type, "ALL")
  
  
##Make annotation row dataframe

annotation_row<- data.frame(chromosome = cytoband.list[[2]]$chromosome)
rownames(annotation_row)<- cytoband.list[[2]]$cytoband
dim(annotation_row)
head(annotation_row)



##Make heat map

# pheatmap(heatmap.matrix.cytoband.del[,1:38],
#          cluster_row = F,
#          cluster_cols = F,
#          show_rownames = TRUE,
#          show_colnames = TRUE,
#          fontsize_row=1,
#          #cellwidth = 10,
#          annotation_row = annotation_row,
#          annotation_legend = FALSE
# )

# ?rainbow
# col.pal<-rainbow(9)
# col.pal<-heat.colors(9)
# col.pal<-terrain.colors(9)
# col.pal<-topo.colors(9)
# col.pal<-cm.colors(9)

col.pal<- colorRampPalette(c("white", "navy", "firebrick3"))(1000)


# col.pal<- rev(col.pal)
# 
# col.pal <- RColorBrewer::brewer.pal(9, "YlGnBu")
# col.pal

pheatmap(heatmap.matrix.cytoband.del[,1:38],
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

pheatmap(heatmap.matrix.cytoband.del[,1:37],
         cluster_row = F,
         cluster_cols = T,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = col.pal,
         fontsize_row=1,
         #cellwidth = 10,
         annotation_row = annotation_row,
         annotation_legend = FALSE
)


#################
###Amplification heatmap:

###########
##Loop that creates matrix

cancer.type<- names(cnv.list)
cancer.type

heatmap.matrix.cytoband.ampl<- matrix(NA, ncol = length(cancer.type)+1, nrow = 806)
dim(heatmap.matrix.cytoband.ampl)

for (i in 1:length(cancer.type)){
  
  x<-cnv.list[[i]]
  cytoband.list<- events.per.cytoband(x, threshold = 1, cytoband_column = 3, column_data_start = 4, chromosome_interval = 0,  deletion = FALSE)
  heatmap.matrix.cytoband.ampl[,i]<- cytoband.list[[2]]$proportion.of.deletions
  print(cancer.type[i])
}

head(heatmap.matrix.cytoband.ampl)

##Final column contains all cancer types:

##Create one large dataframe with all CNV data in it:
x<-join.cnv.datasets(cnv.list, 4)

##Calculate proportion of deletions per cytoband and add to matrix
cytoband.list<- events.per.cytoband(x, threshold = 1, cytoband_column = 3, column_data_start = 4, chromosome_interval = 0,  deletion = FALSE)
heatmap.matrix.cytoband.ampl[,38]<- cytoband.list[[2]]$proportion.of.deletions

head(heatmap.matrix.cytoband.ampl)
class(heatmap.matrix.cytoband.ampl)
dim(x)

heatmap.matrix.cytoband.ampl[250:350,37]

##Add row and column names

rownames(heatmap.matrix.cytoband.ampl)<- cytoband.list[[2]]$cytoband
colnames(heatmap.matrix.cytoband.ampl)<- c(cancer.type, "ALL")


##Make annotation row dataframe

annotation_row<- data.frame(chromosome = cytoband.list[[2]]$chromosome)
rownames(annotation_row)<- cytoband.list[[2]]$cytoband
dim(annotation_row)
head(annotation_row)



##Make heat map

col.pal<- colorRampPalette(c( "white","navy", "firebrick3"))(1000)

pheatmap(heatmap.matrix.cytoband.ampl[,1:38],
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

pheatmap(heatmap.matrix.cytoband.ampl[,1:37],
         cluster_row = F,
         cluster_cols = T,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = col.pal,
         fontsize_row=1,
         #cellwidth = 10,
         annotation_row = annotation_row,
         annotation_legend = FALSE
)




################
#### Create matrix of proportion of deletions/amplifications per MB per chromosome and plot heatmap.
##################

##################
### Sum of deletions for all tumour types together:
##Comment: Could make an array to store all chromosomes for each tumour.

chromosome_interval_1MB<- ceiling(cytoband.list[[3]]$intervals_for_Mb)

##Create one large dataframe with all CNV data in it:
#x<-join.cnv.datasets(cnv.list, 4)
dim(x)

##Matrix to store results per chromosome:
heatmap.matrix.chr.interval250<- data.frame(matrix(NA, ncol = 1, nrow = 2500))
dim(heatmap.matrix.chr.interval250)
heatmap.matrix.chr.interval250[,1]<-seq(1:2500)
colnames(heatmap.matrix.chr.interval250)<-c("intervals")
colnames(heatmap.matrix.chr.interval250)
chr<-c(seq(1:22), "X")

##Calculate proportion of deletions per cytoband and add to matrix
for (i in 1:23) {
  j<-chr[i]
  chromosome_interval<-chromosome_interval_1MB[i]
cytoband.list<- events.per.cytoband(x, threshold = -1, cytoband_column = 3, column_data_start = 4, select_chromosome = j , chromosome_interval = chromosome_interval,  deletion = TRUE)
y<-data.frame(intervals = seq(1:length(cytoband.list[[4]]$proportion.of.deletions)), chr = cytoband.list[[4]]$proportion.of.deletions)
# glimpse(cytoband.list[[1]])
# glimpse(cytoband.list[[2]])
# glimpse(cytoband.list[[3]])
# glimpse(cytoband.list[[4]])
heatmap.matrix.chr.interval250<- full_join(heatmap.matrix.chr.interval250, y, by = "intervals")
# heatmap.matrix.chr.interval250[,i]<- cytoband.list[[4]][,5]
print(i)
print(length(cytoband.list[[4]]$proportion.of.deletions))
}

colnames(heatmap.matrix.chr.interval250)<- c("Interval", paste0("Chromosome ", seq(1:22)), "Chromosome X")

head(heatmap.matrix.chr.interval250)
ncol(heatmap.matrix.chr.interval250)

##############
### Plot heatmaps


#my.matrix<- as.matrix(heatmap.matrix.chr.interval250[,2:24])
#is.na(heatmap.matrix.chr.interval250$chr.x)


col.pal<- colorRampPalette(c( "white","navy", "firebrick3"))(1000)

pheatmap(heatmap.matrix.chr.interval250[1:205,2:24],
         cluster_row = F,
         cluster_cols = F,
         show_rownames = FALSE,
         show_colnames = TRUE,
         color = col.pal,
         fontsize_row=1
         #cellwidth = 10,
         #annotation_row = annotation_row,
         #annotation_legend = FALSE
)

pheatmap(heatmap.matrix.chr.interval250[1:205,2:24],
         cluster_row = F,
         cluster_cols = T,
         show_rownames = FALSE,
         show_colnames = TRUE,
         color = col.pal,
         fontsize_row=1
         #cellwidth = 10,
         #annotation_row = annotation_row,
         #annotation_legend = FALSE
)




################
#### Create matrix of ratio of deletions to amplifications per cytoband and plot heatmap.
##################