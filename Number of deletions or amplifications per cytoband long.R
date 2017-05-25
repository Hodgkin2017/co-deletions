#############
## Measure the number of deletions or amplifications per cytoband
##############

library(tidyr)
library(dplyr)
library(org.Hs.eg.db)


##dummy data

dummy.data<- acc.cnv[1:100,]
head(dummy.data)

##Parameters
x<- dummy.data
x<- acc.cnv
threshold<- -1
#column_gene_ID<- 2
column_data_start<- 4
cytoband_column<-3
#chromosome_column<- 4
#gene_location_column<-5
deletion<- TRUE
chromosome_interval<- 0
chromosome_interval<- 10
select_chromosome<- 1

events.per.cytoband<- function(x, threshold, cytoband_column, column_data_start, chromosome_interval = 0, select_chromosome, deletion = TRUE){
##########
## Obtain chromosomal locations of genes

keys<- keys(org.Hs.eg.db, keytype = "ENTREZID")
columns<- c("CHR", "CHRLOC", "CHRLOCEND")
all.gene.locations<- AnnotationDbi::select(org.Hs.eg.db, keys, columns, keytype = "ENTREZID")

locus.id<- x$Locus.ID

genes.of.interest<- all.gene.locations[all.gene.locations$ENTREZID %in% locus.id,]

genes.of.interest<- na.omit(genes.of.interest[!duplicated(genes.of.interest$ENTREZID), ])

## Remove genes with poorly defined chromosome
chr<-c(seq(1:22), "X", "Y")

rows.of.interest<- which(genes.of.interest$CHR %in% chr)

genes.of.interest<- genes.of.interest[rows.of.interest,]

genes.of.interest$strand<- ifelse(genes.of.interest$CHRLOC <0, "-", "+")
genes.of.interest$start<- abs(genes.of.interest$CHRLOC)
genes.of.interest$end<- abs(genes.of.interest$CHRLOCEND)



##Join gene location dataframe to original CNV data table

genes.of.interest<- dplyr::rename(genes.of.interest, Locus.ID = ENTREZID)

genes.of.interest$Locus.ID<- as.integer(genes.of.interest$Locus.ID)

genes.of.interest<- full_join(genes.of.interest,x, by="Locus.ID")

##Order genes by chromosome and location:
genes.of.interest<- dplyr::arrange(genes.of.interest, CHR, start)

column_data_start<- column_data_start + 7
cytoband_column<- cytoband_column + 7

## Code to get proportion of deletions per cytoband
cnv.matrix<- as.matrix(genes.of.interest[,column_data_start:ncol(genes.of.interest)])

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
  cbind(cytoband=genes.of.interest[,cytoband_column], 
        chromosome = genes.of.interest$CHR, 
        start = genes.of.interest$start, .) 

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

results
}
    
## Test function:

test<- events.per.cytoband(x, threshold, cytoband_column, column_data_start, chromosome_interval = 0,  deletion = TRUE)
glimpse(test[[1]])
glimpse(test[[2]])

test1<- events.per.cytoband(acc.cnv, threshold = -1, cytoband_column = 3, column_data_start = 4, chromosome_interval = 0,  deletion = TRUE)
glimpse(test1[[1]])
glimpse(test1[[2]])
glimpse(test1[[3]])
glimpse(test1[[4]])

test2<- events.per.cytoband(acc.cnv, threshold = -1, cytoband_column = 3, column_data_start = 4 , select_chromosome = 1, chromosome_interval = 10,  deletion = TRUE)
glimpse(test2[[1]])
glimpse(test2[[2]])
glimpse(test2[[3]])
glimpse(test2[[4]])

