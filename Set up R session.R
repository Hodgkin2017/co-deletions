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

##########
### F1: Obtain chromosomal locations of genes

chromosomal_location<- function(object_name){
  
  ##Variables I want:
  keys<- as.character(object_name$Locus.ID)
  columns<- c("CHR", "CHRLOC", "CHRLOCEND")
  
  ## Search for gene chromosomal locations:
  genes.of.interest<- AnnotationDbi::select(org.Hs.eg.db, keys, columns, keytype = "ENTREZID")
  
  ## Remove duplicated entries:
  genes.of.interest<- na.omit(genes.of.interest[!duplicated(genes.of.interest$ENTREZID), ])
  
  ## Add additional columns to dataframe including strand, gene start and end:
  genes.of.interest$strand<- ifelse(genes.of.interest$CHRLOC <0, "-", "+")
  genes.of.interest$start<- abs(genes.of.interest$CHRLOC)
  genes.of.interest$end<- abs(genes.of.interest$CHRLOCEND)
  
  ##########
  ###Join gene location dataframe to original CNV data table
  
  ##Rename ENTREZID column so it can be joined to oringinal dataframe:
  genes.of.interest<- dplyr::rename(genes.of.interest, Locus.ID = ENTREZID)
  
  genes.of.interest$Locus.ID<- as.integer(genes.of.interest$Locus.ID)
  
  genes.of.interest<- full_join(genes.of.interest,object_name, by="Locus.ID")
  
  ##########
  ### Order genes by chromosome and location:
  
  ## Rename X and Y chromosomes to integer so they can be properly sorted.
  genes.of.interest$CHR<-sub("X", "23", genes.of.interest$CHR)
  genes.of.interest$CHR<-sub("Y", "24", genes.of.interest$CHR)
  genes.of.interest$CHR<- as.integer(genes.of.interest$CHR)
  genes.of.interest<- dplyr::arrange(genes.of.interest, CHR, start)
  genes.of.interest$CHR<-sub("23","X",  genes.of.interest$CHR)
  genes.of.interest$CHR<-sub("24", "Y", genes.of.interest$CHR)
  
  return(genes.of.interest)
}
