##Create dummy CNV data:
A<- c(1, 1, 1, 0, -1)
B<- c(-1, 1, 0, 0, -1)
C<- c(-1, 1, -1, -1, -1)
D<- c(-1, 1, -1, 0, 1)
E<- c(0, 0, -1, -1, 1)
G<- c(1, 1, 0, -1, 0)

df<- data.frame(rbind(A,B,C,D,E,G))
class(df)
df

###############
##Function

create.heatmap.matrix.ampl.del<- function(x, column_start, threshold, deletion = TRUE){
  
  ##Convert CNV data to matrix:
  
  cnv.matrix<- as.matrix(x[,column_start:ncol(x)])
  
  if (deletion == TRUE) {
    ##Threshold data such that anything with a value of the threshold or lower = TRUE
    
    cnv.matrix<- cnv.matrix<= threshold
    
  } else {
    
    cnv.matrix<- cnv.matrix >= threshold
    
  }
  
  ## Convert dataframe to 1 and 0:
  cnv.matrix<- cnv.matrix*1
  
  ## Calculate the proportion of individuals with a deletion:
  prop.individ.with.del<- rowSums(cnv.matrix)/ncol(cnv.matrix)
  
  ##loop to create a table with proportion of deletions shared between genes
  
  heatmap.matrix<-data.frame(matrix(NA, ncol = nrow(cnv.matrix), nrow = nrow(cnv.matrix)))
  
  for(i in 1: nrow(cnv.matrix)){
    
    for (j in 1: nrow(cnv.matrix)) {
      
      value<- cnv.matrix[i,]*cnv.matrix[j,]
      heatmap.matrix[i,j]<- sum(value)/sum(cnv.matrix[i,])
      
    }
  }
  
  for(i in 1:nrow(heatmap.matrix)) {
    
    heatmap.matrix[i,i]<- prop.individ.with.del[i]
  }
  
  rownames(heatmap.matrix)<- rownames(x)
  colnames(heatmap.matrix)<- rownames(x)
  
  heatmap.matrix
}

####################
##Use function


system.time(test.matrix<- create.heatmap.matrix.ampl.del(df, column_start = 1, threshold = -1, deletion = TRUE))
test.matrix

system.time(test.matrix2<- create.heatmap.matrix.ampl.del(df, column_start = 1, threshold = 1, deletion = FALSE))
test.matrix2

##################
## Increase speed of function:
x<- df
column_start<- 1
threshold<- -1
deletion = TRUE

 
create.heatmap.matrix.ampl.del.optimised<- function(x, column_start, threshold, deletion = TRUE){
  
  ##Convert CNV data to matrix:
  
  cnv.matrix<- as.matrix(x[,column_start:ncol(x)])
  
  ##Create a binary matrix of CNV data such that deletions (deletion = TRUE) or 
  #amplifications (deletion = FALSE) below (for deletions) or above (for amplifications) a threshold = 1 
  if (deletion == TRUE) {
    
    cnv.matrix[cnv.matrix > threshold] = 0
    cnv.matrix = -cnv.matrix
    
  } else {
    
    #cnv.matrix<- cnv.matrix >= threshold
    cnv.matrix[cnv.matrix < threshold] = 0
    
  }
  
  ## Calculate the proportion of individuals with a deletion:
  prop.individ.with.del<- rowSums(cnv.matrix)/ncol(cnv.matrix)
  
  ## Calculate total number of co-deletions of amplifications
  
  heatmap.matrix<- cnv.matrix%*%t(cnv.matrix)
  heatmap.matrix<- heatmap.matrix/rowSums(cnv.matrix)
  
  ## Complete diagonals with proportion of individuals with deletions.
  for(i in 1:nrow(heatmap.matrix)) {
    
    heatmap.matrix[i,i]<- prop.individ.with.del[i]
  }
  
  ## Add row and column names
  rownames(heatmap.matrix)<- rownames(x)
  colnames(heatmap.matrix)<- rownames(x)
  
  heatmap.matrix
}


############
## Run function

system.time(test.matrix3<- create.heatmap.matrix.ampl.del.optimised(df, column_start = 1, threshold = -1, deletion = TRUE))
test.matrix3
test.matrix
identical(test.matrix3, test.matrix)
table(test.matrix3 == test.matrix, useNA = 'ifany')
all.equal(test.matrix3,test.matrix)

system.time(test.matrix4<- create.heatmap.matrix.ampl.del.optimised(df, column_start = 1, threshold = 1, deletion = FALSE))
test.matrix4
test.matrix2


