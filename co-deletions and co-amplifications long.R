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
acc.cnv.chr.location<- chromosomal_location(cnv.list[[1]])
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
    
  } else if (Chromosome[1] > 0 & start == FALSE){
    
    matrix<- object_name %>% dplyr::filter(CHR %in% Chromosome) 
    cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
    rownames(cnv.matrix)<- matrix$Gene.Symbol
    
  } else if (start == TRUE & Chromosome[1] == 0){

    matrix<- object_name$start %>%
      dplyr::between(selection_criteria[1], selection_criteria[2]) %>%
      object_name[.,]
    
    matrix <- matrix %>% filter(!is.na(start))
    
    cnv.matrix<- as.matrix(matrix[,column_start:ncol(matrix)])
    rownames(cnv.matrix)<- matrix$Gene.Symbol
    
  }  else if (start == TRUE & Chromosome[1] > 0){
  
    #object_name %>% filter(CHR == 9) %>% dplyr::select(start)
    
    
    matrix<- object_name %>% dplyr::filter(CHR %in% Chromosome) 
    
    matrix<- matrix$start %>%
      dplyr::between(selection_criteria[1], selection_criteria[2]) %>%
      matrix[.,]
    
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

test2<- create.heatmap.matrix.ampl.del.optimised(acc.cnv.chr.location, column_start = 11, threshold = -1, Gene.Symbol = TRUE, selection_criteria = c("MET", "CDKN2A", "RB1", "WWOX", 
                                                                                                                                                     "LRP1B", "PDE4D", "CCNE1", "TP53",
                                                                                                                                                     "FGFR1", "MYC", "EGFR","WHSC1L1",
                                                                                                                                                     "ERBB2", "MCL1", "MDM2", "CCND1", "ATM",
                                                                                                                                                     "NOTCH1", "PPP2R2A", "BRD4", "ARID1A",
                                                                                                                                                     "STK11", "PARK2"), deletion = TRUE)
test2a<- create.heatmap.matrix.ampl.del.optimised(acc.cnv.chr.location, column_start = 11, threshold = 1, Gene.Symbol = TRUE, selection_criteria = c("MET", "CDKN2A", "RB1", "WWOX", 
                                                                                                                                                     "LRP1B", "PDE4D", "CCNE1", "TP53",
                                                                                                                                                     "FGFR1", "MYC", "EGFR","WHSC1L1",
                                                                                                                                                     "ERBB2", "MCL1", "MDM2", "CCND1", "ATM",
                                                                                                                                                     "NOTCH1", "PPP2R2A", "BRD4", "ARID1A",
                                                                                                                                                     "STK11", "PARK2"), deletion = FALSE)                                                                                                                                                     
test3<- create.heatmap.matrix.ampl.del.optimised(acc.cnv.chr.location, column_start = 11, threshold = -1, Chromosome = 9, deletion = TRUE)
test4<- create.heatmap.matrix.ampl.del.optimised(acc.cnv.chr.location, column_start = 11, threshold = -1, Cytoband = TRUE, selection_criteria = c("9p21.2", "9p21.3"), deletion = TRUE)
test5<- create.heatmap.matrix.ampl.del.optimised(acc.cnv.chr.location, column_start = 11, threshold = -1, start = TRUE, selection_criteria = c(100000, 250000), deletion = TRUE)

test6<-create.heatmap.matrix.ampl.del.optimised(acc.cnv.chr.location, column_start = 11, threshold = -1, Chromosome = 9, start = TRUE, selection_criteria = c(100000, 250000), deletion = TRUE)
test7<-create.heatmap.matrix.ampl.del.optimised(acc.cnv.chr.location, column_start = 11, threshold = -1, Chromosome = 9, start = TRUE, selection_criteria = c(100000, 50000000), deletion = TRUE)
object_name %>% filter(CHR == 9) %>% dplyr::select(start) %>% range()
dim(test2)

pheatmap(test1[16000:16100,16000:16100],
         cluster_row = F,
         cluster_cols = F,
         show_rownames = FALSE,
         show_colnames = FALSE
)
pheatmap(test4,
         cluster_row = F,
         cluster_cols = F,
         show_rownames = FALSE,
         show_colnames = FALSE
)
dim(test2)
pheatmap(test2a,
         cluster_row = F,
         cluster_cols = F,
         show_rownames = TRUE,
         show_colnames = TRUE
)


##############
##Heatmaps showing co-deletions and co-amplifications for each chromosome for all tumour types combined.

CNV.all.table.chr.location<- chromosomal_location(CNV.all.table)

# CNV.all.table.co.del<- co.deletion_co.amplification_matrix(CNV.all.table.chr.location, column_start = 11, threshold = -1, 
#                                                            deletion = TRUE)
# 
# CNV.all.table.co.amp<- co.deletion_co.amplification_matrix(CNV.all.table.chr.location, column_start = 11, threshold = 1, 
#                                                            deletion = FALSE)


list.CNV.all.table.co.del<-lapply(c(seq(1:23), "X"), function(x) co.deletion_co.amplification_matrix(CNV.all.table.chr.location, column_start = 11, threshold = -1,Chromosome = x, 
                                                                                                     deletion = TRUE))

# test10<-lapply(c(seq(1:23), "X"), function(x) co.deletion_co.amplification_matrix(acc.cnv.chr.location, column_start = 11, threshold = 1, Chromosome = x, 
#                                                                                                      deletion = FALSE))
# names(test10)<- paste("chromosome", c(seq(1:23), "X"),  sep = " ")
# test11<-test10[9:10]
# names(test11)
# #lapply(test11, function(x) head(x))
# 
# lapply(test11, function(x) plot_heatmap(x))
# lapply(test11, function(x) plot_heatmap.dev(x))

plot_heatmap<-function(x){

  pheatmap(x,
           cluster_row = F,
           cluster_cols = F,
           show_rownames = FALSE,
           show_colnames = FALSE
  )


}
lapply(list.CNV.all.table.co.del, function(x) plot_heatmap(x))

list.CNV.all.table.co.amp<-lapply(c(seq(1:23), "X"), function(x) co.deletion_co.amplification_matrix(CNV.all.table.chr.location, column_start = 11, threshold = 1,Chromosome = x, 
                                                                                                     deletion = FALSE))
lapply(list.CNV.all.table.co.amp, function(x) plot_heatmap(x))




##Does not work:
plot_heatmap.dev<-function(x){
  
  pdf(filename=paste(names(x),".pdf",sep=""), width=6, height=6)
  pheatmap(x,
           cluster_row = F,
           cluster_cols = F,
           show_rownames = FALSE,
           show_colnames = FALSE
  )
  dev.off()
  
}

##Try?:
plot_heatmap.dev<-function(x){
  
  tiff(file=paste(names(x),".pdf",sep=""), width=6, height=6)
  pheatmap(x,
           cluster_row = F,
           cluster_cols = F,
           show_rownames = FALSE,
           show_colnames = FALSE
  )
  dev.off()
  
}



# plot_heatmap.dev<-function(x){
#   
#   pdf(filename=paste(names(x),".pdf",sep=""), width=6, height=6)
#   pheatmap(x,
#            cluster_row = F,
#            cluster_cols = F,
#            show_rownames = FALSE,
#            show_colnames = FALSE
#   )
#   dev.off()
#   
# }

#list.CNV.all.table.co.del<-lapply(c(seq(1:23), "X"), function(x) dplyr::filter(CNV.all.table.co.del, CHR == x))
# 
# test11<-test10[1:2]
# lapply(test11, plot_heatmap(x))
# 
# plot_heatmap<-function(x){
#   
#   pheatmap(x,
#            cluster_row = F,
#            cluster_cols = F,
#            show_rownames = FALSE,
#            show_colnames = FALSE
#   )
#   
#   
# }
###############
##Co-deletion and co-amplification of target genes

head(colnames(CNV.all.table.chr.location),15)

christophe.genes.del<- co.deletion_co.amplification_matrix(CNV.all.table.chr.location, column_start = 11, threshold = -1, Gene.Symbol = TRUE, selection_criteria = c("MET", "CDKN2A", "RB1", "WWOX", 
                                                                                                                                                                     "LRP1B", "PDE4D", "CCNE1", "TP53",
                                                                                                                                                                     "FGFR1", "MYC", "EGFR","WHSC1L1",
                                                                                                                                                                     "ERBB2", "MCL1", "MDM2", "CCND1", "ATM",
                                                                                                                                                                     "NOTCH1", "PPP2R2A", "BRD4", "ARID1A",
                                                                                                                                                                     "STK11", "PARK2"), deletion = TRUE)
christophe.genes.amp<- co.deletion_co.amplification_matrix(CNV.all.table.chr.location, column_start = 11, threshold = 1, Gene.Symbol = TRUE, selection_criteria = c("MET", "CDKN2A", "RB1", "WWOX", 
                                                                                                                                                                    "LRP1B", "PDE4D", "CCNE1", "TP53",
                                                                                                                                                                    "FGFR1", "MYC", "EGFR","WHSC1L1",
                                                                                                                                                                    "ERBB2", "MCL1", "MDM2", "CCND1", "ATM",
                                                                                                                                                                    "NOTCH1", "PPP2R2A", "BRD4", "ARID1A",
                                                                                                                                                                    "STK11", "PARK2"), deletion = FALSE)

pdf(file="christophe.genes.del.pdf")
pheatmap(christophe.genes.del,
         cluster_row = F,
         cluster_cols = F,
         show_rownames = TRUE,
         show_colnames = TRUE
)
dev.off()

pdf(file="christophe.genes.amp.pdf")
pheatmap(christophe.genes.amp,
         cluster_row = F,
         cluster_cols = F,
         show_rownames = TRUE,
         show_colnames = TRUE
)
dev.off()







############
##Identification of top co-deletion and co-amplification genes.


CNV.all.table.co.del<- co.deletion_co.amplification_matrix(CNV.all.table.chr.location, column_start = 11, threshold = -1, deletion = TRUE)
CNV.all.table.co.amp<- co.deletion_co.amplification_matrix(CNV.all.table.chr.location, column_start = 11, threshold = 1, deletion = FALSE)

dim(CNV.all.table.co.del)
CNV.all.table.co.del2<-cbind(Gene.Symbol = rownames(CNV.all.table.co.del), CNV.all.table.co.del)
dim(CNV.all.table.co.del2)
head(colnames(CNV.all.table.co.del2),12)
# CNV.all.table.co.del.names<- rownames(CNV.all.table.co.del)
# CNV.all.table.co.del3<-dplyr::bind_cols(CNV.all.table.co.del.names, CNV.all.table.co.del)
# dim(CNV.all.table.co.del3)
CNV.all.table.co.del2<- as.data.frame(CNV.all.table.co.del2)
dim(CNV.all.table.co.del2)
CNV.all.table.co.del.long<- tidyr::gather(CNV.all.table.co.del2, key = Gene.Symbol )
dim(CNV.all.table.co.del.long)
head(CNV.all.table.co.del.long)
colnames(CNV.all.table.co.del.long)<- c("Gene.Symbol1", "Gene.Symbol2", "value")
head(CNV.all.table.co.del.long, 20)
CNV.all.table.co.del.long.sort<- dplyr::arrange(CNV.all.table.co.del.long, desc(value))
CNV.all.table.co.del.long.sort2<- dplyr::arrange(CNV.all.table.co.del.long, value)
head(CNV.all.table.co.del.long.sort, 20)
head(CNV.all.table.co.del.long.sort2, 20)

CNV.all.table.co.del.long.ungroup<- dplyr::ungroup(CNV.all.table.co.del.long)
CNV.all.table.co.del.long.sort3<- dplyr::arrange(CNV.all.table.co.del.long.ungroup, desc(value))
head(CNV.all.table.co.del.long.sort3, 20)
max(CNV.all.table.co.del.long.sort3$value)
tail(CNV.all.table.co.del.long.sort3, 20)
head(CNV.all.table.co.del.long.sort3, 100)
tail(CNV.all.table.co.del.long.sort3, 20)

top1000.co.deletions<- head(CNV.all.table.co.del.long.sort3, 1000)
write.csv(top1000.co.deletions, file = "top1000.co.deletions.csv", quote = FALSE)

CNV.all.table.co.del.long.sort3.CDKN2A<- dplyr::filter(CNV.all.table.co.del.long.sort3, Gene.Symbol2 == "CDKN2A")

pheatmap(CNV.all.table.co.del[1:100,1:100],
         cluster_row = F,
         cluster_cols = F,
         show_rownames = TRUE,
         show_colnames = TRUE
)


# d1<- c(2,3,4)
# d2<- c(1,2,3)
# d3<- c(10,11,12)
# 
# dummy<- data.frame(A=d1, B=d2, C=d3)
# rownames(dummy)<- c("A", "B", "C")
# 
# tidyr::gather(dummy)
# 
# dummy2<-cbind(name = c("A", "B", "C"), dummy)
# dummy2
# 
# df2<- tidyr::gather(dummy2, key = name )
# df2
# colnames(df2)<- c("Gene.Symbol1", "Gene.Symbol2", "value")
# df2
# dplyr::arrange(df2, desc(value))













