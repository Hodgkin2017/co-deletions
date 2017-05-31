A<- c(1, 1, 1, 0, -1)
B<- c(-1, 1, 0, 0, -1)
C<- c(-1, 1, -1, -1, -1)
D<- c(-1, 1, -1, 0, 1)
E<- c(0, 0, -1, -1, 1)
G<- c(1, 1, 0, -1, 0)

df<- data.frame(rbind(A,B,C,D,E,G))
class(df)
df
object_name<- df
column_start = 1
threshold = -1
deletion = TRUE


Gene.Symbol<- c("A", "B", "C", "D", "E", "G")
start<- c(100, 200, 250, 300, 400, 500)
Cytoband<- c("1p", "2p", "3p", "3q", "4q", "5p")
CHR<- c(1, 2, 3, 3, 4, 5)
df2<- cbind(Gene.Symbol, start, CHR, Cytoband, df)
df2
object_name<- df2

column_start = 5
threshold = -1
deletion = TRUE


#########
object_name<- acc.cnv.chr.location
head(object_name)
column_start = 11
threshold = -1
deletion = TRUE
Chromosome = 9

##Faster function to create matrix of proportions of co-deletions
# and co-amplifications

create.heatmap.matrix.ampl.del.optimised<- function(object_name, column_start = 11, threshold = -1, 
                                                    selection_criteria, Gene.Symbol = FALSE, start = FALSE, 
                                                    Chromosome = 0, Cytoband = FALSE, deletion = TRUE){
  if (Gene.Symbol == TRUE){

    matrix<- object_name %>% dplyr::filter(Gene.Symbol %in% selection_criteria) 
    cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
    rownames(cnv.matrix)<- matrix$Gene.Symbol
    
  } else if (Chromosome[1] > 0){
    
    matrix<- object_name %>% dplyr::filter(CHR %in% Chromosome) 
    cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
    rownames(cnv.matrix)<- matrix$Gene.Symbol
    
  } else if (start == TRUE & Chromosome[1] == 0){
    ##filter by chromosome 1st?
    
    # matrix<- matrix$start %>%
    #   dplyr::between(selection_criteria[1], selection_criteria[2]) %>%
    #   object_name[.,]
    
    matrix<- object_name$start %>%
      dplyr::between(selection_criteria[1], selection_criteria[2]) %>%
      object_name[.,]
    
    cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
    rownames(cnv.matrix)<- matrix$Gene.Symbol
    
  }  else if (start == TRUE & Chromosome[1] > 0){
  
    matrix<- matrix$start %>%
      dplyr::between(selection_criteria[1], selection_criteria[2]) %>%
      object_name[.,]
    
    cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
    rownames(cnv.matrix)<- matrix$Gene.Symbol
    
  } else if (Cytoband == TRUE){
    
    matrix<- object_name %>% dplyr::filter(Cytoband %in% selection_criteria) 
    cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
    rownames(cnv.matrix)<- matrix$Gene.Symbol
      
  } else {
  ##Convert CNV data to matrix:
  
  cnv.matrix<- as.matrix(object_name[,column_start:ncol(object_name)])
  rownames(cnv.matrix)<- object_name$Gene.Symbol
  }
  
  ##Create a binary matrix of CNV data such that deletions (deletion = TRUE) or 
  #amplifications (deletion = FALSE) below (for deletions) or above (for amplifications) a threshold = 1 
  if (deletion == TRUE) {
    
    #cnv.matrix[cnv.matrix > threshold] = 0
    #cnv.matrix = -cnv.matrix
    cnv.matrix<- ifelse(cnv.matrix <= threshold, 1, 0)
    
  } else {
    
    #cnv.matrix<- cnv.matrix >= threshold
    #cnv.matrix[cnv.matrix < threshold] = 0
    cnv.matrix<- ifelse(cnv.matrix >= threshold, 1, 0)
  }
  
  ## Calculate total number of co-deletions or amplifications
  
  heatmap.matrix<- cnv.matrix%*%t(cnv.matrix)
  heatmap.matrix<- heatmap.matrix/ncol(cnv.matrix)
  
  ##NaN caused when rowSums = 0. Replace NaN with 0
  #cnv.matrix[cnv.matrix == NaN] = 0
  heatmap.matrix[is.nan(heatmap.matrix)] = 0
  
  ## Add row and column names
  #rownames(heatmap.matrix)<- rownames(object_name)
  #colnames(heatmap.matrix)<- rownames(object_name)
  
  return(heatmap.matrix)
}

create.heatmap.matrix.ampl.del.optimised(df, column_start = 1, threshold = -1, deletion = TRUE)

test<-create.heatmap.matrix.ampl.del.optimised(acc.cnv.chr.location, column_start = 11, threshold = -1, deletion = TRUE)
rownames(test)

create.heatmap.matrix.ampl.del.optimised(df2, column_start = 4, threshold = -1, deletion = TRUE)

create.heatmap.matrix.ampl.del.optimised(df2, column_start = 5, threshold = -1, selection_criteria = c("A", "B"), deletion = TRUE)
create.heatmap.matrix.ampl.del.optimised(df2, column_start = 5, threshold = -1, Gene.Symbol = TRUE, selection_criteria = c("A", "B"), deletion = TRUE)
create.heatmap.matrix.ampl.del.optimised(df2, column_start = 5, threshold = -1, Chromosome = c(3,4), selection_criteria = c(1,2), deletion = TRUE)
create.heatmap.matrix.ampl.del.optimised(df2, column_start = 5, threshold = -1, Cytoband = TRUE, selection_criteria = c("1p", "3p"), deletion = TRUE)
create.heatmap.matrix.ampl.del.optimised(df2, column_start = 5, threshold = -1, start = TRUE, selection_criteria = c(0,400), deletion = TRUE)
create.heatmap.matrix.ampl.del.optimised(df2, column_start = 5, threshold = -1, start = TRUE, selection_criteria = c(101,300), deletion = TRUE)
create.heatmap.matrix.ampl.del.optimised(df2, column_start = 5, threshold = -1, Chromosome = c(2,3), start = TRUE, selection_criteria = c(0,400), deletion = TRUE)

##co-deletions ... no filtering
test1<- create.heatmap.matrix.ampl.del.optimised(acc.cnv.chr.location, column_start = 11, threshold = -1, deletion = TRUE)

test2<- create.heatmap.matrix.ampl.del.optimised(acc.cnv.chr.location, column_start = 11, threshold = -1, Gene.Symbol = TRUE, selection_criteria = c("A", "B"), deletion = TRUE)
test3<- create.heatmap.matrix.ampl.del.optimised(acc.cnv.chr.location, column_start = 11, threshold = -1, Chromosome = 9, deletion = TRUE)
test4<- create.heatmap.matrix.ampl.del.optimised(acc.cnv.chr.location, column_start = 11, threshold = -1, Cytoband = TRUE, selection_criteria = c("9p21.2", "9p21.3"), deletion = TRUE)
test5<- create.heatmap.matrix.ampl.del.optimised(acc.cnv.chr.location, column_start = 11, threshold = -1, start = TRUE, selection_criteria = c(4000,100000), deletion = TRUE)
test6<-create.heatmap.matrix.ampl.del.optimised(acc.cnv.chr.location, column_start = 11, threshold = -1, Chromosome = 9, start = TRUE, selection_criteria = c(4000,100000), deletion = TRUE)
dim(test6)

pheatmap(test1[1:100,1:100],
         cluster_row = F,
         cluster_cols = F,
         show_rownames = FALSE,
         show_colnames = FALSE
)
pheatmap(test6,
         cluster_row = F,
         cluster_cols = F,
         show_rownames = FALSE,
         show_colnames = FALSE
)
