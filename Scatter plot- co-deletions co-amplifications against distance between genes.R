###################
### Plot propotion of tumours with co-deletion or co-amplification as a function of distance.
###################

#######
##Test of making long dataframe with gather function:
i=1
cancer.table<- chromosomal_location(short.cnv.list[[i]])
cancer.table[1:2,1:11]
ncol(cancer.table) -10
co.deletion_co.amplification_matrix(cancer.table, column_start = 11, threshold = -1, Cytoband = TRUE, selection_criteria = "9p21.3", deletion = TRUE, normalisation = "none")
test<-co.deletion_co.amplification_matrix(cancer.table, column_start = 11, threshold = -1, Cytoband = TRUE, selection_criteria = "9p21.3", deletion = TRUE, normalisation = "total.tumour.number")
dim(test)
# test2<-data.frame(gene1 = 1:4, gene2 = 1:4, gene3 = 1:4, gene4 = 1:4)
# rownames(test2)<-c("gene1", "gene2", "gene3", "gene4")
# test2<- cbind(rownames(test2), test2)
# test2

tidyr::gather(test2, gene.id1, gene.id1, 2:5)

test<- as.data.frame(cbind(Gene.Symbol.row = rownames(test), test))
test
test3<- tidyr::gather(test, Gene.Symbol.col,proportion, 2:ncol(test))

nrow(test3)

##Maybe plot per chromosome arm rather than per chromosome
##Try plotting for each cytoband.
test4<- co.deletion_co.amplification_matrix(cancer.table, column_start = 11, threshold = -1, Chromosome = 9, deletion = TRUE, normalisation = "total.tumour.number")

###############
##Create three column dataframe of proportion of co-deletions normalised by number of tumours. One row = one gene pair


##Test:
# cancer.table[1:2, 1:12]
# 
# dummydata<-data.frame(Gene.Symbol = c("MET", "CDKN2A", "CDKN2B", "RB", "MET2", "CDKN2A2", "CDKN2B2", "RB2"), start = c(10, 100, 120, 300, 20, 200, 220, 400), CHR = c(9, 9, 9, 9, 10, 10, 10, 10), sample1 = c(-1,-1,-1,-1, -1,-1,-1,-1), sample2 = c(-1,-1,-1,-1, 0,0,0,0), sample3 = c(0,0,0,0, 0,0,0,0), sample4 = c(0, 0, 0, 0, 0,0,0,0))
# dummydata
# co.deletion_co.amplification_matrix(dummydata, column_start = 4, threshold = -1, Chromosome = 9, deletion = TRUE, normalisation = "total.tumour.number")
# co.deletion_co.amplification_matrix(dummydata, column_start = 4, threshold = -1, Chromosome = 10, deletion = TRUE, normalisation = "total.tumour.number")
# 
# co.deletion.per.chromosome<- lapply(c(9,10), function(x) co.deletion_co.amplification_matrix(dummydata, column_start = 4, threshold = -1, Chromosome = x, deletion = TRUE, normalisation = "total.tumour.number"))
# co.deletion.per.chromosome
# dim(co.deletion.per.chromosome[[2]])
# co.deletion.per.chromosome<- lapply(co.deletion.per.chromosome, function(x) as.data.frame(cbind(Gene.Symbol.row = rownames(x), x)))
# co.deletion.per.chromosome
# gathered<- lapply(co.deletion.per.chromosome, function(x) tidyr::gather(x, Gene.Symbol.col,proportion, 2:ncol(x)))
# gathered
# 
# # ?do.call
# do.call(rbind, gathered)








co.deletion.per.chromosome.table<- data.frame(Gene.Symbol.row = NA, Gene.Symbol.col = NA, proportion_deletions = NA)
co.deletion.per.chromosome.table

c((1:22), "X")

i=1
cancer.table<- chromosomal_location(short.cnv.list[[i]])

co.deletion.per.chromosome<- lapply(c((1:22), "X"), function(x) co.deletion_co.amplification_matrix(cancer.table, column_start = 11, threshold = -1, Chromosome = x, deletion = TRUE, normalisation = "total.tumour.number"))
dim(co.deletion.per.chromosome[[2]])
co.deletion.per.chromosome<- lapply(co.deletion.per.chromosome, function(x) as.data.frame(cbind(Gene.Symbol.row = rownames(x), x)))

gathered<- lapply(co.deletion.per.chromosome, function(x) tidyr::gather(x, Gene.Symbol.col,proportion, 2:ncol(x)))
glimpse(gathered)
gathered.co.deletion.per.chromosome<- do.call(rbind, gathered)
dim(gathered.co.deletion.per.chromosome)

#test3<- tidyr::gather(test, Gene.Symbol.col,proportion, 2:ncol(test))








############
## Create a datafame of inter gene distances using the start of the gene.


# cancer.table[1:2,1:11]
# chromosome9<- cancer.table %>%
#   dplyr::filter(CHR == 9)
# 
# chromosome9.start<- cancer.table %>%
#   dplyr::filter(CHR == 9) %>%
#   dplyr::select(Gene.Symbol, start)
# 
# intergene.distance<- chromosome9.start$start[2:nrow(chromosome9.start)] - chromosome9.start$start[1:(nrow(chromosome9.start)-1)]
# length(intergene.distance)
# 
# 
# chromosome9[1:2,1:11]


## Create for loop function to loop through each each chromosome and generate intergene distance

# apply()
# 
# dim(cancer.table)
# 
# interval.values<- c((1:22), "X")
# interval.values
# 
# for (i in 1:23) {}
# i=1
# 
# selected.interval<- interval.values[i]
# selected.interval

intergene_distance_function<- function(cancer.table, selected.interval){

selected.genes.start<- cancer.table %>%
  dplyr::filter(CHR == selected.interval) %>%
  dplyr::select(start)

selected.genes.start

intergene.distance.table <- data.frame(matrix(NA, ncol = nrow(selected.genes.start), nrow = nrow(selected.genes.start)))
dim(intergene.distance.table)
intergene.distance.table

for (j in 1:nrow(selected.genes.start)) {

intergene.distance.table[,j]<- abs(selected.genes.start[,1] - selected.genes.start[j,1])

}

Gene_Symbol<- cancer.table %>%
  dplyr::filter(CHR == selected.interval) %>%
  dplyr::select(Gene.Symbol) %>%
  t() %>% #required to convert output from data.frame to vector
  as.character()

class(Gene_Symbol)
length(Gene_Symbol)

colnames(intergene.distance.table)<- Gene_Symbol
rownames(intergene.distance.table)<- Gene_Symbol

return(intergene.distance.table)

}


##Test function:
##Make a list of intergene distances:
# x<- selected.interval
# 
# test<-intergene_distance_function(cancer.table, selected.interval = x)
# 
# identical(test, intergene.distance.table)

intergene_distance<-lapply(c((1:22), "X"), function(x) intergene_distance_function(cancer.table, selected.interval = x))

dim(intergene_distance[[2]])


#write.csv(intergene.distance.table, file = "intergene.distance.table.csv")

## Convert dataframes in list to dataframes with one pairwise distance per row and join dataframes together

intergene_distance<- lapply(intergene_distance, function(x) as.data.frame(cbind(Gene.Symbol.row = rownames(x), x)))
dim(intergene_distance[[2]])
gathered_intergene_distance<- lapply(intergene_distance, function(x) tidyr::gather(x, Gene.Symbol.col,proportion, 2:ncol(x)))
dim(gathered_intergene_distance[[2]])
glimpse(gathered_intergene_distance)
##add chromosome name as extra column to each dataframe in list


selected_intergene_distance<- do.call(rbind, gathered_intergene_distance)
dim(selected_intergene_distance)


##join intergene distances to co-deletions list

co_deletions_plot_table<- cbind(inter_gene_distance = selected_intergene_distance, co_deletions = gathered.co.deletion.per.chromosome$proportion)
head(co_deletions_plot_table)
dim(co_deletions_plot_table)
## repeat for other chromosomes and have chromosome as a column? 
#actually I have done this..repeat for other cancers

##Remove all entries with 0 distance or 0 co-deletions

# dummy<-data.frame(gene1 = co_deletions_plot_table[1:10,1], 
#                   gene2 = co_deletions_plot_table[1:10,2],
#                   x = c(1,2,0,4,5,6,7,8,0, 10),
#                   y = c(1,2,3,4,0,6,0,8,9, 10)
#                   )
# dummy %>%
#   dplyr::filter(x != 0) %>%
#   dplyr::filter(y != 0)


co_deletions_removed_zeros_plot_table<-co_deletions_plot_table %>%
  dplyr::filter(inter_gene_distance.proportion != 0) %>%
  dplyr::filter(co_deletions != 0)
dim(co_deletions_removed_zeros_plot_table)

##Make scatter plot of propotion of tumours with co-deletion or co-amplification as a function of distance.
small.co_deletions_plot_table<- co_deletions_removed_zeros_plot_table[1:10000,]
ggplot(small.co_deletions_plot_table, aes(x = as.numeric(inter_gene_distance.proportion), y = as.numeric(co_deletions))) +
  geom_point(size = 1, shape = 1) +
  geom_smooth()

ggplot(co_deletions_removed_zeros_plot_table, aes(x = inter_gene_distance.proportion, y = as.integer(co_deletions))) +
  geom_point(size = 1, shape = 1) +
  geom_smooth()

# geom_point(colour = "red", size = 3)
# geom_smooth()...data processing intensive

##Comment: do not plot any with a score of 0 for x or y?
#colour and shape by chromosome or cancer, ?


