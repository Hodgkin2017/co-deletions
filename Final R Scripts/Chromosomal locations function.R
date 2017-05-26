#############
## Measure the number of deletions or amplifications per cytoband or chromosome interval
##############

###Contents:

##Functions:
#F1. chromosomal_location   #Function to obtain the chromosomal location of genes in a table using  
                            #their ENTREZ gene IDs  
                            #Returns a dataframe.
                            #Dataframe 1: Contains original input data plus gene chromosomal locations
                            #Function arguments:
                            #object_name                  # CNV file to attch gene chromosomal  
                                                          # locations to.

##Actions and or Loops outside of functions:
#None

##Selected Objects:
#O1. test                  #Test function using all tumour CNV data (CNV.all.data)

##Packages
library(tidyr)
library(dplyr)
library(org.Hs.eg.db)
library(biomaRt)

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

#############
###O1: Test function


test<-chromosomal_location(CNV.all.table)
head(test)
dim(test)
dim(CNV.all.table)



