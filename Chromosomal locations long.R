#############
## Measure the number of deletions or amplifications per cytoband or chromosome interval
##############

###Contents:


##Packages
library(tidyr)
library(dplyr)
library(org.Hs.eg.db)
library(biomaRt)

#######
###Load example data:
acc.cnv<- read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/unzipped original broad TCGA CNV data/ACC/all_data_by_genes.txt", stringsAsFactors = FALSE, header = TRUE)
object_name<- acc.cnv
object_name<- CNV.all.table

##########
### Obtain chromosomal locations of genes

# keys<- keys(org.Hs.eg.db, keytype = "ENTREZID")
# keys
# columns<- c("CHR", "CHRLOC", "CHRLOCEND")
# all.gene.locations<- AnnotationDbi::select(org.Hs.eg.db, keys, columns, keytype = "ENTREZID")
# 
# locus.id<- x$Locus.ID
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
# head(genes.of.interest)
# dim(genes.of.interest)

###################
## Quicker method?

chromosomal_location<- function(object_name){
keys<- as.character(object_name$Locus.ID)
columns<- c("CHR", "CHRLOC", "CHRLOCEND")

genes.of.interest<- AnnotationDbi::select(org.Hs.eg.db, keys, columns, keytype = "ENTREZID")

genes.of.interest<- na.omit(genes.of.interest[!duplicated(genes.of.interest$ENTREZID), ])

genes.of.interest$strand<- ifelse(genes.of.interest$CHRLOC <0, "-", "+")
genes.of.interest$start<- abs(genes.of.interest$CHRLOC)
genes.of.interest$end<- abs(genes.of.interest$CHRLOCEND)
dim(genes.of.interest)

##########
###Join gene location dataframe to original CNV data table

genes.of.interest<- dplyr::rename(genes.of.interest, Locus.ID = ENTREZID)

genes.of.interest$Locus.ID<- as.integer(genes.of.interest$Locus.ID)

genes.of.interest<- full_join(genes.of.interest,object_name, by="Locus.ID")
dim(genes.of.interest)

##########
## Order genes by chromosome and location:
#genes.of.interest<- dplyr::arrange(genes.of.interest, CHR, start)

genes.of.interest$CHR<-sub("X", "23", genes.of.interest$CHR)
genes.of.interest$CHR<-sub("Y", "24", genes.of.interest$CHR)
genes.of.interest$CHR<- as.integer(genes.of.interest$CHR)
genes.of.interest<- dplyr::arrange(genes.of.interest, CHR, start)
genes.of.interest$CHR<-sub("23","X",  genes.of.interest$CHR)
genes.of.interest$CHR<-sub("24", "Y", genes.of.interest$CHR)
dim(genes.of.interest)

return(genes.of.interest)
}


dim(CNV.all.table)
test<-chromosomal_location(CNV.all.table)
head(test)
dim(test)
##Comment: Add if statement that allows the function to obtain chromosomal locations via 
#org.Hs.eg.db or Biomart using either gene name or entrez ID????










############################################
### Old:
############################################

##########
## Obtain chromosomal locations of genes


keys<- keys(org.Hs.eg.db, keytype = "ENTREZID")
keys
columns<- c("CHR", "CHRLOC", "CHRLOCEND")
sel<- AnnotationDbi::select(org.Hs.eg.db, keys, columns, keytype = "ENTREZID")
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

## Remove genes with poorly defined chromosome
chr<-c(seq(1:22), "X", "Y")
chr

rows.of.interest<- which(sel3$CHR %in% chr)
head(rows.of.interest)
length(rows.of.interest)

sel3.final<- sel3[rows.of.interest,]
dim(sel3.final)


sel3.final$strand<- ifelse(sel3.final$CHRLOC <0, "-", "+")
sel3.final
sel3.final$start<- abs(sel3.final$CHRLOC)
sel3.final$end<- abs(sel3.final$CHRLOCEND)
head(sel3.final)
dim(sel3.final)

##which genes did not have defined chromosome
rows.not.of.interest<-which(!sel3$CHR %in% chr)
length(which(!sel3$CHR %in% chr))
genes.removed.from.sel3<- sel3[rows.not.of.interest,]
head(genes.removed.from.sel3)

##Biomart method:

myattributes <- c("ensembl_gene_id",
                  "entrezgene",
                  "external_gene_name",
                  "chromosome_name",
                  "start_position",
                  "end_position",
                  "strand",
                  "gene_biotype",
                  "description")

# Human
grch38 <- useMart("ensembl") %>% 
  useDataset(mart=., dataset="hsapiens_gene_ensembl") %>% 
  getBM(mart=., attributes=myattributes) #%>% 
#fix_genes
View(grch38)


## find rows in downloaded human database whose entrez IDs match those from the CNV files
rows.of.interest<- which(grch38$entrezgene %in% acc.locus.id)
head(rows.of.interest)

genes.in.cnv.file<- grch38[rows.of.interest,]
dim(genes.in.cnv.file)
rownames(genes.in.cnv.file)<- seq(1:nrow(genes.in.cnv.file))

##Remove duplicate genes
sel4<- na.omit(genes.in.cnv.file[!duplicated(genes.in.cnv.file$entrezgene), ])
head(sel4)
dim(sel4)
## Comment: not as good as org.Hs.eg.db method

length(which(duplicated(genes.in.cnv.file$entrezgene)))

##Use gene names:

rows.of.interest<- which(grch38$external_gene_name %in% acc.cnv$Gene.Symbol)
head(rows.of.interest)

genes.in.cnv.file<- grch38[rows.of.interest,]
dim(genes.in.cnv.file)
#rownames(genes.in.cnv.file)<- seq(1:nrow(genes.in.cnv.file))
#genes.in.cnv.file<- arrange(genes.in.cnv.file, chromosome_name)


##Remove duplicate genes
sel4<- na.omit(genes.in.cnv.file[!duplicated(genes.in.cnv.file$entrezgene), ])
head(sel4)
dim(sel4)
sel5<- na.omit(genes.in.cnv.file[!duplicated(genes.in.cnv.file$external_gene_name), ])
head(sel5)
dim(sel5)
## Comment: Better than org.Hs.eg.db method but is it accurate?Number of genes that dont have correct chromosome annotation


##Remove genes without correct chr name:
chr<-c(seq(1:22), "X", "Y")
chr

rows.of.interest<- which(sel4$chromosome_name %in% chr)
head(rows.of.interest)

sel6<- sel4[rows.of.interest,]
dim(sel6)

rows.of.interest<- which(sel5$chromosome_name %in% chr)
head(rows.of.interest)

sel7<- sel5[rows.of.interest,]
dim(sel7)

length(which(duplicated(sel7$entrezgene)))
##Comment remove replicates of both gene ID and gene name?

sel8<- na.omit(sel7[!duplicated(sel7$entrezgene), ])
head(sel8)
dim(sel8)


##Check for duplicates in my data
length(which(duplicated(acc.cnv$Gene.Symbol)))
length(which(duplicated(acc.cnv$Locus.ID)))


################
## What are the genes that do not have unique entrex IDs (using org.Hs.eg.db)?

## Can not do it automatically if they are not in the database!!!!
#Look manually...

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























