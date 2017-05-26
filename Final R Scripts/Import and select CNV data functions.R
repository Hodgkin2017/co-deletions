#############
## Import and select CNV data functions
#############

###Contents:

##Functions:
#F1. Function to import multiple files from multiple folders
#F2. Function to select some or all CNV datasets from CNV list object 
    #and combining data into one large dataframe.

##Actions and or Loops:
#A1. For loop to check all CNV datasets contain the same gene names.
#A2. For loop to check all CNV datasets contain the same entrez gene IDs.

##Selected Objects:
#O1. cnv.list                  #List object containing all cancer CNV dataframes
#O2. identical.gene.names      #Vector confirming that all CNV dataframes contain the same gene names
#O3. identical.Locus.IDs       #Vector confirming that all CNV dataframes contain the same Entrez gene 
                              #IDs (Locus.ID)
#O4. CNV.ACC.BRCA.STAD.table   #Dataframe containing ACC BRCA and STAD CNV data sets 
#O5. CNV.all.table             #Dataframe containing all CNV data from all cancer types 
                              #dim= 24776 genes x 10234 tumours


#############
## Importing tables from multiple folders
#############

##F1: Function to import multiple files from multiple folders

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



## Choose directory and name of files to import from each directory:

x<-"/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/unzipped original broad TCGA CNV data"
file.to.import<-"all_data_by_genes.txt"

##O1: Run function to obtain a list object containing all CNV dataframes

cnv.list<- import.files.from.directories(x, file.to.import)

#################
## Check that each CNV dataset contains the same genes and Entrez IDs (Locus.ID column):

##Make table of gene names
gene.names<- sapply(cnv.list, '[[',1)
lapply(cnv.list, dim)

##Sort gene names before checking each column is identical
gene.names.sorted <- apply(gene.names,2,sort,decreasing=F)

## A1: Check the gene names for each CNV dataset is identical
identical.gene.names<- rep(NA, ncol(gene.names.sorted))
for (i in 1: ncol(gene.names.sorted)){
  
  identical.gene.names[i]<- identical(gene.names.sorted[,1], gene.names.sorted[,i])
}
##O2:
identical.gene.names


##Make table of entrez IDs
Locus.IDs<- sapply(cnv.list, '[[',2)

##Sort entrez IDs before checking each column is identical
Locus.IDs.sorted <- apply(Locus.IDs,2,sort,decreasing=F)

## A2: Check the entrez IDs for each CNV dataset are identical
identical.Locus.IDs<- rep(NA, ncol(Locus.IDs.sorted))
for (i in 1: ncol(Locus.IDs.sorted)){
  
  identical.Locus.IDs[i]<- identical(Locus.IDs.sorted[,1], Locus.IDs.sorted[,i])
}
##O3:
identical.Locus.IDs


#############
## Selecting some or all CNV datasets from CNV list object and combining data into one large dataframe.
#############

##F2: Function to select some or all CNV datasets from CNV list object and combining data into one large dataframe.

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


##O4: Create dataframe containing ACC, BRCA and STAD CNV data:
CNV.ACC.BRCA.STAD.table<-join.cnv.datasets(cnv.list, column, data.sets = c("ACC", "BRCA", "STAD"))
dim(CNV.ACC.BRCA.STAD.table)
which(is.na(CNV.ACC.BRCA.STAD.table), arr.ind = T)

##O5: Create one large dataframe with all CNV data in it:
CNV.all.table<-join.cnv.datasets(cnv.list, column = 4)
dim(CNV.all.table)
which(is.na(CNV.all.table), arr.ind = T)












