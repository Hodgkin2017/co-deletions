#############
## Importing tables from multiple folders
#############

setwd("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/unzipped original broad TCGA CNV data/")
setwd("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/unzipped original broad TCGA CNV data/test")
directory.names<-dir()
directory.names

temp = list.files(pattern="*")
temp

myfiles = lapply(directory.names, read.delim, stringsAsFactors = FALSE, header = TRUE)
names(myfiles)<- directory.names
names(myfiles)
head(myfiles$acc_all_data_by_genes.txt)

x<- directory.names
y<- "all_data_by_genes.txt"



function(x,y){
  list.files(path= ~/x, pattern=y)
  read.delim(x, stringsAsFactors = FALSE, header = TRUE)
  all_data_by_genes.txt
  
  
}

x
y

d<-"/"
c<-"ACC"

#########
x<-"/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/unzipped original broad TCGA CNV data"
file.to.import<-"all_data_by_genes.txt" 

import.files.from.directories<-function(x,file.to.import){
currentwd<- getwd()
setwd(x)
directory.names<-dir()
my.list <- vector("list", length(directory.names))
for (i in 1: length(directory.names)){
  move.to.directory<-paste0(x,"/",directory.names[i])
  setwd(move.to.directory)
  my.list[[i]]<-read.delim(file.to.import, stringsAsFactors = FALSE, header = TRUE)
  print(directory.names[i])
}
names(my.list)<- directory.names
setwd(currentwd)
my.list
}

cnv.list<- import.files.from.directories(x, file.to.import)
cnv.list[[1]][[1]]
cnv.list[[1]][,1]
cnv.list[[1]][1,]
lapply(cnv.list, '[[', 1)

##Make table of gene names
gene.names<- sapply(cnv.list, '[[',1)
lapply(cnv.list, dim)


head(gene.names)
class(gene.names)
dim(gene.names)

##Are the gene names columns identical between CNV datasets
gene.names<- as.data.frame(gene.names)
apply(gene.names, 2, identical, gene.names$ACC) 
apply(gene.names, 2, function(x) identical(x, gene.names$ACC)) 
apply(gene.names, 2, function(x) (x==gene.names$ACC)) 

identical(gene.names$ACC, gene.names$BLCA)

result<- rep(NA, ncol(gene.names))
for (i in 1: ncol(gene.names)){
  
  result[i]<- identical(gene.names[,1], gene.names[,i])
}
result
 
##Test which genes and how many are not identical between datasets

test<- lapply(cnv.list, '[[',1)
dim(test)
length(test)
class(test)

cnv.list$ACC[1:3,1:3]
test2<-cnv.list$ACC[,1:2]
test3<-cnv.list$BLCA[,1:2]
test4<-cnv.list$CESC[,1:2]
dim(test3)
head(test2)

test.compare<- full_join(test2, test3, by="Gene.Symbol")
dim(test.compare)
test.compare2<- full_join(test.compare, test4, by="Gene.Symbol")
dim(test.compare2)
head(test.compare2)
test.compare2[1,4]<- NA
which(is.na(test.compare2), arr.ind = TRUE)


## Maybe gene names are not ordered properly, sort and then check gene names are identical:
gene.names.sorted <- apply(gene.names,2,sort,decreasing=F)
head(gene.names)
head(gene.names.sorted)

result2<- rep(NA, ncol(gene.names.sorted))
for (i in 1: ncol(gene.names.sorted)){
  
  result2[i]<- identical(gene.names.sorted[,1], gene.names.sorted[,i])
}
result2

all.equal(gene.names.sorted[,1], gene.names.sorted[,4])

###Check Locus.ID is the same between CNV data sets:

##Make table of entrez IDs
Locus.IDs<- sapply(cnv.list, '[[',2)
head(Locus.IDs)
class(Locus.IDs)
dim(Locus.IDs)


Locus.IDs.sorted <- apply(Locus.IDs,2,sort,decreasing=F)

result3<- rep(NA, ncol(Locus.IDs.sorted))
for (i in 1: ncol(Locus.IDs.sorted)){
  
  result3[i]<- identical(Locus.IDs.sorted[,1], Locus.IDs.sorted[,i])
}
result3

###############
### Combine dataframes within my list

new.list<- cnv.list[1:2]
names(new.list)
length(new.list)
dim(new.list$ACC)

lapply(new.list,function(x) select(x, c(-Cytoband, -Gene.Symbol)))

cnv.list$ACC %>% dplyr::select(c(-Cytoband, -Gene.Symbol)) %>% colnames
new.list[[1]] %>% dplyr::select(c(-Cytoband, -Gene.Symbol))

dummylist<- list(ACC=cnv.list$ACC[1:10,1:10], BRCA=cnv.list$BRCA[1:10,1:10], STAD=cnv.list$STAD[1:10,1:10])
dummylist

#Paramaters:
x<- cnv.list
x<- dummylist
column<- 4
data.sets<- "all data sets"
data.sets<-c("ACC", "STAD")

##Function to identify which tables in list to join

join.cnv.datasets<- function(x, column, data.sets = "all data sets"){

if(data.sets == "all data sets"){
  
  index<- seq(1:length(x))
  
} else {
names.of.tables<-names(x)
names.of.tables
index<-which(names.of.tables %in% data.sets)
print("DO NOT WORRY ABOUT THE FOLLOWING WARNING MESSAGE:")

}

df<- x[[1]] %>% dplyr::select(c(Gene.Symbol,Locus.ID,Cytoband))

for (i in index){
  
  df2<-x[[i]][,c(1,column:ncol(x[[i]]))]
  
  df<- full_join(df, df2, by = "Gene.Symbol")

}
df
}



df<-join.cnv.datasets(x, column, data.sets = c("ACC", "BRCA", "STAD"))
dim(df)
which(is.na(df), arr.ind = T)

df<-join.cnv.datasets(x, column)
dim(df)
which(is.na(df), arr.ind = T)


