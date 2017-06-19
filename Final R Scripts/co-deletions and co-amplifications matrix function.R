#############
## Measure the proportion of co-deletions or co-amplifications
##############

###Contents:

##Functions:
#F1. co.deletion_co.amplification_matrix    #Function to obtain the proportion of co-deletions or co-amplifications
                                            #Returns a matrix containing proportion of co-deletions or co-amplifications.
                                            #Function arguments:
                                            #cnv.table          # CNV file
                                            #column_start         # First column in file containing tumour CNV data
                                            #threshold            # Value above (deletion = FALSE) or below 
                                                                  #(deletion = TRUE) which CNVs will be included in analysis
                                            #selection_criteria   #vector of terms to select appropriate CNV table rows
                                                                  #e.g. Gene.Symbol, start and end chromosomal intervals or Cytobands
                                            #Gene.Symbol          #If Gene.Symbol = TRUE: Select for genes matching selection_criteria
                                                                  #If Gene.Symbol = FALSE: Do not select for genes
                                            #start                #If start = TRUE: Select for chromosomal interval matching 
                                                                  #selection_criteria = c(start, end). Can be used with Chromosome > 0.
                                                                  #If start = FALSE: Do not select for chromosomal interval
                                            #Chromosome           # Select chromosomes of interset e.g. Chromosome = c(3,9)
                                                                  #If Chromosome = 0 then no chromosome is selected.
                                            #Cytoband             # If Cytoband = TRUE: Select for Cytobands matching selection_criteria
                                                                  # If deletion = FALSE: Do not select for Cytobands
                                            #remove_NA = TRUE     # If remove_NA = TRUE: Remove all genes without a start position.
                                            #deletion = TRUE      # If deletion = TRUE: Count number of deletion events
                                                                  # If deletion = FALSE: Count number of amplification events
                                            #normalisation        #If normalisation = "total.tumour.number" then normalise by total number of tumours
                                                                  #If normalisation = "tumours.with.event" then normalise by 
                                                                  # number of tumours with deletion or amplification event. e.g. if only 7 out of 90 
                                                                  # individuals had a deletion in a cytoband then divide number of co-deletions by 7
                                                                  #If normalisation = "none" then raw number of co-deletions or co-amplifications returned.
                                                                  #If normalisation = "frequency.for.whole.sample" then normalise number of events by total number of possible events
                                                                  # reurns a single value and not a matrix
##Actions and or Loops outside of functions:
#None

##Selected Objects:
#O1. test1                  #Output all co-deletions
#O2. test2                  #Output co-deletions for some genes only
#O3. test2a                 #Output co-amplifications for some genes only
#O4. test3                  #Output co-deletions for chromosome 9 only
#O5. test4                  #Output co-deletions for Cytobands "9p21.2" and "9p21.3" only
#O6. test5                  #Output co-deletions between residues 100,000 and 250,000 only
#O7. test6                  #Output co-deletions on chromosome 9 between 100,000 and 250,000 only

################
###Packages
library(dplyr)

####################
###F1: Function to create matrix containing proportions of co-deletions or amplifications 
#normalised by the total number of tumours compared

co.deletion_co.amplification_matrix<- function(cnv.table, column_start = 11, threshold = -1, 
                                               selection_criteria, Gene.Symbol = FALSE, start = FALSE, 
                                               Chromosome = 0, Cytoband = FALSE, remove_NA = TRUE, deletion = TRUE, normalisation = "total.tumour.number"){
  if (remove_NA == TRUE){
    
    cnv.table<- cnv.table %>%
      dplyr::filter(!is.na(start))
  }
  
  if (Gene.Symbol == TRUE){
    
    ##Select Genes of interest by Gene Symbol and convert CNV data to matrix:
    matrix<- cnv.table %>% dplyr::filter(Gene.Symbol %in% selection_criteria) 
    cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
    rownames(cnv.matrix)<- matrix$Gene.Symbol
    
  } else if (Chromosome[1] > 0 & start == FALSE){
    
    ##Select Chromosome of interest and convert CNV data to matrix:
    matrix<- cnv.table %>% dplyr::filter(CHR %in% Chromosome) 
    cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
    rownames(cnv.matrix)<- matrix$Gene.Symbol
    
  } else if (start == TRUE & Chromosome[1] == 0){
    
    ##Select Chromosomal region of interest and convert CNV data to matrix:
    matrix<- cnv.table$start %>%
      dplyr::between(selection_criteria[1], selection_criteria[2]) %>%
      cnv.table[.,]
    
    ## Remove NAs which are maintained by dplyr::between function
    matrix <- matrix %>% filter(!is.na(start))
    
    cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
    rownames(cnv.matrix)<- matrix$Gene.Symbol
    
  }  else if (start == TRUE & Chromosome[1] > 0){
    
    ##Select Chromosome of interest and convert CNV data to matrix:
    matrix<- cnv.table %>% dplyr::filter(CHR %in% Chromosome) 
    
    ##Select Chromosomal region of interest and convert CNV data to matrix:
    # matrix<- matrix$start %>%
    # dplyr::between(selection_criteria[1], selection_criteria[2]) %>%
    # matrix[.,]
    matrix<- matrix %>%
      dplyr::filter(start >= selection_criteria[1], end <= selection_criteria[2])
    
    ##Convert Chromosome and region of interest into matrix:
    cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
    rownames(cnv.matrix)<- matrix$Gene.Symbol
    
  } else if (Cytoband == TRUE){
    ##Select Cytobands of interest and convert CNV data to matrix:
    matrix<- cnv.table %>% dplyr::filter(Cytoband %in% selection_criteria) 
    cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
    rownames(cnv.matrix)<- matrix$Gene.Symbol
    
  } else {
    ##Convert ALL CNV data to matrix:
    cnv.matrix<- as.matrix(cnv.table[,column_start:ncol(cnv.table)])
    rownames(cnv.matrix)<- cnv.table$Gene.Symbol
  }
  
  ##Create a binary matrix of CNV data such that deletions (deletion = TRUE) or 
  #amplifications (deletion = FALSE) below (for deletions) or above (for amplifications) a threshold = 1 
  if (deletion == TRUE) {
    
    cnv.matrix<- ifelse(cnv.matrix <= threshold, 1, 0)
    
  } else {
    
    cnv.matrix<- ifelse(cnv.matrix >= threshold, 1, 0)
  }
  
  
  ## Calculate total number of co-deletions or amplifications
  heatmap.matrix<- cnv.matrix%*%t(cnv.matrix)
  
  if (normalisation == "total.tumour.number") {
    
    heatmap.matrix<- heatmap.matrix/ncol(cnv.matrix)
    
  } else if (normalisation == "tumours.with.event") {
    tumours.with.del.or.amp<- colSums(cnv.matrix)
    tumours.with.del.or.amp<-sum(tumours.with.del.or.amp >0)
    heatmap.matrix<- heatmap.matrix/tumours.with.del.or.amp
    
  } else if (normalisation == "none") {
    
  } else if (normalisation == "frequency.for.whole.sample"){
    heatmap.matrix<- sum(heatmap.matrix)/(ncol(cnv.matrix)*(nrow(cnv.matrix)^2))
    
  }
  
  ## Remove NAs caused by divideing by 0.
  heatmap.matrix[is.nan(heatmap.matrix)] = 0
  
  return(heatmap.matrix)
}

##############
### Examples of function in use:

## Use CNV data for cancer type ACC. Table included chromosomal location information:
acc.cnv.chr.location<- chromosomal_location(cnv.list[[1]])

##O1: Output all co-deletions
test1<- co.deletion_co.amplification_matrix(acc.cnv.chr.location, column_start = 11, threshold = -1, 
                                            deletion = TRUE)

##O2: Output co-deletions for some genes only
test2<- co.deletion_co.amplification_matrix(acc.cnv.chr.location, column_start = 11, threshold = -1, 
                                            Gene.Symbol = TRUE, selection_criteria = c("MET", "CDKN2A", 
                                                                                       "RB1", "WWOX", "LRP1B", "PDE4D", "CCNE1", "TP53","FGFR1", 
                                                                                       "MYC", "EGFR","WHSC1L1", "ERBB2", "MCL1", "MDM2", "CCND1",
                                                                                       "ATM","NOTCH1","PPP2R2A", "BRD4", "ARID1A", "STK11", 
                                                                                       "PARK2"), deletion = TRUE)

##O3: Output co-amplifications for some genes only
test2a<- co.deletion_co.amplification_matrix(acc.cnv.chr.location, column_start = 11, threshold = 1, 
                                             Gene.Symbol = TRUE, selection_criteria = c("MET", "CDKN2A",
                                                                                        "RB1", "WWOX", "LRP1B", "PDE4D", "CCNE1", "TP53","FGFR1", 
                                                                                        "MYC", "EGFR","WHSC1L1", "ERBB2", "MCL1", "MDM2", "CCND1", 
                                                                                        "ATM","NOTCH1", "PPP2R2A", "BRD4", "ARID1A","STK11", 
                                                                                        "PARK2"), deletion = FALSE)                                                                                                                                                     

##O4: Output co-deletions for chromosome 9 only
test3<- co.deletion_co.amplification_matrix(acc.cnv.chr.location, column_start = 11, threshold = -1, Chromosome = 9, deletion = TRUE)

##O5: Output co-deletions for Cytobands "9p21.2" and "9p21.3" only
test4<- co.deletion_co.amplification_matrix(acc.cnv.chr.location, column_start = 11, threshold = -1, Cytoband = TRUE, selection_criteria = c("9p21.2", "9p21.3"), deletion = TRUE)

##O6: Output co-deletions between residues 100,000 and 250,000 only
test5<- co.deletion_co.amplification_matrix(acc.cnv.chr.location, column_start = 11, threshold = -1, start = TRUE, selection_criteria = c(100000, 250000), deletion = TRUE)

##O7: Output co-deletions on chromosome 9 between 100,000 and 250,000 only
test6<-co.deletion_co.amplification_matrix(acc.cnv.chr.location, column_start = 11, threshold = -1, Chromosome = 9, start = TRUE, selection_criteria = c(100000, 250000), deletion = TRUE)
