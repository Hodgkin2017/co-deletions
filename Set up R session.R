####################
## Set up R session for co-deletions and amplifications project
####################

## This scripts installs packages, data and functions that I might need for the R session.

#################
###Load packages
library(tidyr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(org.Hs.eg.db)
library(biomaRt)

###################
###Set working directory
setwd("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/")


#################
### Load data
load("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/R workspaces/co-deletions workspace.RData")

################
### Load Functions.

###F1: Function to import multiple files from multiple folders

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


##F2: Function to select some or all CNV datasets from CNV list object 
#and combining data into one large dataframe.
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

