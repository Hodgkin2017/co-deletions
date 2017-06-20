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
##Create long three column dataframe of proportion of co-deletions normalised by number of tumours. One row = one gene pair


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
glimpse(gathered_intergene_distance[[2]])
names(gathered_intergene_distance)<-c((1:22), "X")

################
##add chromosome name as extra column to each dataframe in list
#gathered_intergene_distance2<- lapply(gathered_intergene_distance, function(x) as.data.frame(cbind(x, chromosome = rep(names(x),nrow(x)))))

#head(gathered_intergene_distance2[[2]])

gathered_intergene_distance2<-Map(cbind, gathered_intergene_distance, chromosome = names(gathered_intergene_distance))

# head(test[[1]])
# head(test[[2]])
# head(test[[23]])
# tail(test[[23]])
# list2env(my_list, .GlobalEnv) 
# my_list <- Map(cbind, my_list, new_clumn = names(my_list))

#############

selected_intergene_distance2<- do.call(rbind, gathered_intergene_distance2)
dim(selected_intergene_distance2)

##join intergene distances to co-deletions list

co_deletions_plot_table2<- cbind(inter_gene_distance = selected_intergene_distance2, co_deletions = gathered.co.deletion.per.chromosome$proportion)
head(co_deletions_plot_table2)
tail(co_deletions_plot_table2)
dim(co_deletions_plot_table2)
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


co_deletions_removed_zeros_plot_table2<-co_deletions_plot_table2 %>%
  dplyr::filter(inter_gene_distance.proportion != 0) %>%
  dplyr::filter(co_deletions != 0)
dim(co_deletions_removed_zeros_plot_table)

##Convert co-deletion from factor to numeric
class(co_deletions_removed_zeros_plot_table2[,5])
co_deletions_removed_zeros_plot_table2$co_deletions<- as.numeric(levels(co_deletions_removed_zeros_plot_table2$co_deletions))[co_deletions_removed_zeros_plot_table2$co_deletions]
class(co_deletions_removed_zeros_plot_table2[,5])

## Save object
saveRDS(co_deletions_removed_zeros_plot_table2, "/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/BRCA_co_deletion_distance_plot_table.rds")


##Make scatter plot of propotion of tumours with co-deletion or co-amplification as a function of distance.
range(co_deletions_removed_zeros_plot_table2$co_deletions)

## first 10,000 co-deletions
small.co_deletions_plot_table<- co_deletions_removed_zeros_plot_table2[1:10000,]

ggplot(small.co_deletions_plot_table, aes(x = as.numeric(inter_gene_distance.proportion), y = as.numeric(co_deletions))) +
  geom_point(size = 1, shape = 1) +
  geom_smooth()

## Plot all data:
ggplot(co_deletions_removed_zeros_plot_table2, aes(x = inter_gene_distance.proportion, 
                                                   y = co_deletions)) +
  geom_point(size = 1, shape = 1) +
  geom_smooth() +
  xlab("Inter gene distance") +
  ylab("Proportion of tumours with co-deletion")
  
##Save plot
ggsave("BRCA_co-deletion_distance_BW.tiff")

## Plot all data in coloured by chromosome:
ggplot(co_deletions_removed_zeros_plot_table2, aes(x = inter_gene_distance.proportion, 
                                                   y = co_deletions,
                                                   colour = inter_gene_distance.chromosome)) +
  geom_point(size = 1, shape = 1) +
  geom_smooth() +
  xlab("Inter gene distance") +
  ylab("Proportion of tumours with co-deletion") +
  scale_color_discrete()+
  labs(colour ="Chromosome")

##Save plot
ggsave("BRCA_co-deletion_distance_colour.tiff")


#axis! legend?
# geom_point(colour = "red", size = 3)
# geom_smooth()...data processing intensive

##Comment: do not plot any with a score of 0 for x or y?
#colour and shape by chromosome or cancer, ?

##Facet wrap plot:


## Plot all data in coloured by chromosome:

head(co_deletions_removed_zeros_plot_table2)

ggplot(co_deletions_removed_zeros_plot_table2, aes(x = inter_gene_distance.proportion, 
                                                   y = co_deletions,
                                                   colour = inter_gene_distance.chromosome)) +
  geom_point(size = 1, shape = 1) +
  geom_smooth() +
  xlab("Inter gene distance") +
  ylab("Proportion of tumours with co-deletion") +
  scale_color_discrete( guide = FALSE)+
  labs(colour ="Chromosome") +
  facet_wrap(~inter_gene_distance.chromosome) +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))

##Save plot
ggsave("BRCA_co-deletion_distance_colour_wrap.tiff")

###########
## Find outlying genes proportion > 0.015 and distance > 5.0e+07

outlying_genes<-co_deletions_removed_zeros_plot_table2 %>%
  dplyr::filter(inter_gene_distance.proportion > 5.0e+07) %>%
  dplyr::filter(co_deletions > 0.015)
dim(outlying_genes)
head(outlying_genes)
tail(outlying_genes)

outlying_genes %>% 
  dplyr::filter(inter_gene_distance.chromosome == 8) %>%
  dim()

##########################
### Repeat scatter plots but using distance from gene of interest 
###and investigate co-deletions surrounding gene of interest only
##########################


##########
##Find average cytoband sizes:

# cnv.table<- chromosomal_location(short.cnv.list[[1]])
# cnv.table[1:2, 1:12]
# 
# cnv.table %>%
#   dplyr::group_by(Cytoband) %>%
#   dplyr::summarise(max())

##Create input file of cytoband.cordinates
cytoband.cordinates<- read.delim("../../Input data/Annotations/Cytoband_start_end_hgTables", header = TRUE, stringsAsFactors = FALSE )

cytoband.cordinates.X.chrom<-cytoband.cordinates$X.chrom
cytoband.cordinates.X.chrom<- gsub(".*chr","",cytoband.cordinates.X.chrom)
cytoband.cordinates.X.chrom<- paste(cytoband.cordinates.X.chrom, cytoband.cordinates$name, sep = "")
cytoband.cordinates<-cbind(cytoband_name = cytoband.cordinates.X.chrom, cytoband.cordinates)
head(cytoband.cordinates)

cytoband.cordinates<- cytoband.cordinates %>% 
  dplyr::mutate(cytoband_size = chromEnd-chromStart)

head(cytoband.cordinates)

##Average cytoband distances:
mean(cytoband.cordinates$cytoband_size)
sd(cytoband.cordinates$cytoband_size)
range(cytoband.cordinates$cytoband_size)

##Comment: mean cytoband size = 2482046, sd = 2430150, range = 970 to 30627415
##Comment:What is the average size of focal deletions?


############
### Create a long datafame of co-deletions 2.5MB upstream and downstream of gene of interest.


##Use co.deletion_co.amplification_matrix() function to calculate co-deletion of genes 
#2MB upstream and downstream of gene of interest 

# matrix<- co.deletion_co.amplification_matrix(cnv.table, column_start =  11, threshold = threshold, selection_criteria = cytoband, Cytoband = TRUE, deletion = TRUE, normalisation = "tumours.with.event")
# test5<- co.deletion_co.amplification_matrix(acc.cnv.chr.location, column_start = 11, threshold = -1, start = TRUE, selection_criteria = c(100000, 250000), deletion = TRUE)


##Find start and end of gene: BRCA
# cnv.table[1:2, 1:12]
# 
# cnv.table %>% 
#   dplyr::filter(Gene.Symbol == "CDKN2A") %>%
#   dplyr::select(Gene.Symbol, CHR, Cytoband, start, end)

# CHR Cytoband    start      end
# 1   9   9p21.3 21967751 21974827


##take all whole genes 2.5MB down stream of start of gene:
# 
# cnv.table %>% 
#   dplyr::filter(CHR == 9) %>%
#   dplyr::filter(start <= 21967751) %>%
#   dplyr::filter(start >= 21967751 - 2.5e+6) %>%
#   dplyr::select(Gene.Symbol, CHR, Cytoband, start, end)

##take all whole genes 2.5MB up stream of start of gene:

# cnv.table %>% 
#   dplyr::filter(CHR == 9) %>%
#   dplyr::filter(start >= 21974827) %>%
#   dplyr::filter(start <= 21974827 + 2.5e+6) %>%
#   dplyr::select(Gene.Symbol, CHR, Cytoband, start, end)

## Calculate proportion of co-deletions

target_genes<- c("MET", "CDKN2A", "RB1", "WWOX", 
      "LRP1B", "PDE4D", "CCNE1", "TP53",
      "FGFR1", "MYC", "EGFR","WHSC1L1",
      "ERBB2", "MCL1", "MDM2", "CCND1", "ATM",
      "NOTCH1", "PPP2R2A", "BRD4", "ARID1A",
      "STK11", "PARK2")

# x = list(list("CDKN2A",9, 21967751, 21974827))
# x
# 
# y<- c(x, list(list("CDKN2A",9, 21967751, 21974827)))
# y
# class(x[[1]][[2]])
# distance<- 2.5e+06
# 
# co.deletion.per.target.gene<- lapply(x, function(x) co.deletion_co.amplification_matrix(cnv.table, column_start = 11, threshold = -1, start = TRUE, Chromosome = x[[2]], selection_criteria = c(x[[3]] - distance, x[[4]] + distance), deletion = TRUE, normalisation = "total.tumour.number"))


## times two nested apply functions. Inside = target genes co-deletions, outside = different cancers.
#co.deletion.per.chromosome<- lapply(x, function(x) co.deletion_co.amplification_matrix(cancer.table, column_start = 11, threshold = -1, start = TRUE, selection_criteria = c(downstream_limit, upstream_limit), deletion = TRUE, normalisation = "total.tumour.number"))
##Comment: try different normalisation method e.g. number of individuals with deletion or number of 
#individual with deletion in gene x.


#########
###Create list of target genes and their chromosome and start and end

## Create an empty list to store gene information
gene_information_list<- vector("list", length(target_genes)) 
gene_information_list

##loop to create list of genes and their start and stop locations
for (i in 1: length(target_genes)){
  
  gene<- target_genes[i]
  
  gene_information<- cnv.table %>% 
    dplyr::filter(Gene.Symbol == gene) %>%
    dplyr::select(Gene.Symbol, CHR, Cytoband, start, end)
  
  gene_information<- as.list(gene_information)
  gene_information_list[[i]]<- gene_information
  
}

gene_information_list
gene_information_list[[4]][[4]]

##save gene_information.list
saveRDS(gene_information_list, "/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/target_gene_information_list.rds")


##Create co-deletion matricies for each target gene
co.deletion.per.target.gene<- lapply(gene_information_list, function(x) co.deletion_co.amplification_matrix(cnv.table, column_start = 11, threshold = -1, start = TRUE, Chromosome = x[[2]], selection_criteria = c(x[[4]] - distance, x[[5]] + distance), deletion = TRUE, normalisation = "total.tumour.number"))
co.deletion.per.target.gene[[5]]
dim(co.deletion.per.target.gene[[5]])

##Add gene name to each column to be used with gather function later
co.deletion.per.target.gene<- lapply(co.deletion.per.target.gene, function(x) as.data.frame(cbind(Gene.Symbol.row = rownames(x), x)))
co.deletion.per.target.gene[[2]]
dim(co.deletion.per.target.gene[[2]])


##Create a long 3 column wide table with pair-wise proportion of pair wise deletions
gathered<- lapply(co.deletion.per.target.gene, function(x) tidyr::gather(x, Gene.Symbol.col,proportion, 2:ncol(x)))
glimpse(gathered)

##Keep rows relating to MET v's all genes only and not all genes v's all genes
gathered[[1]]

gathered_target_genes<- vector("list", length(target_genes))

for(i in 1: length(target_genes)){
  
  gene<- target_genes[[i]][1]
  
  gathered_target_genes[[i]]<- gathered[[i]] %>% 
    dplyr::filter(Gene.Symbol.col == gene)
}

gathered_target_genes

##Bind all dataframes in list together
gathered.co.deletion.per.target.gene<- do.call(rbind, gathered_target_genes)
dim(gathered.co.deletion.per.target.gene)

##Check the correct number of gene:gene pairwise co-deletions
# number.of.genes<- sapply(co.deletion.per.target.gene, function(x) ncol(x)-1)
# number.of.genes
# sum(number.of.genes^2)

sum(sapply(co.deletion.per.target.gene, function(x) ncol(x)-1))

## comment: gathered.co.deletion.per.target.gene has correct number of dimensions, LRP1B has only 3 genes
#maybe increase the distance from the target gene?


############
### Create a dataframe of gene distances from gene of interest.



target_genes<- c("MET", "CDKN2A", "RB1", "WWOX", 
                 "LRP1B", "PDE4D", "CCNE1", "TP53",
                 "FGFR1", "MYC", "EGFR","WHSC1L1",
                 "ERBB2", "MCL1", "MDM2", "CCND1", "ATM",
                 "NOTCH1", "PPP2R2A", "BRD4", "ARID1A",
                 "STK11", "PARK2")

gene_information_list
gene_information_list[[4]][[4]]
x<- gene_information_list[[5]]
x

distance_from_target_gene_function<- function(cnv.table, x, distance){

##Get start of target gene:
# cnv.table[1:2, 1:12]
# 
# x<- gene_information_list[[1]]
# 
# x[[4]]



##Get end site of genes 2.5MB away of 5' start site of target gene

end_sites_5prime_genes<- cnv.table %>%
  dplyr::filter(CHR == x[[2]]) %>%
  dplyr::filter(start <= x[[4]], start >=  x[[4]] - distance) %>%
  dplyr::select(end)

# end_sites_5prime_genes

##Calculate the distances between genes
distance_5prime_genes<- x[[4]] - end_sites_5prime_genes
# distance_5prime_genes
##Any value <0 = 0 i.e. the gene of interest and any overlapping genes
distance_5prime_genes[distance_5prime_genes < 0]<- 0

##Get the end of the gene

# x[[5]]

##Get start site of genes 2.5MB away of 3' end of end of target gene
start_sites_3prime_genes<-cnv.table %>%
  dplyr::filter(CHR == x[[2]]) %>%
  dplyr::filter(start > x[[4]], end <= x[[5]] + distance ) %>%
  dplyr::select(start)

# start_sites_3prime_genes

##calculate the distance to the end of the genes 5' of the start of the target gene.
distance_3prime_genes<- start_sites_3prime_genes - x[[5]]
# distance_3prime_genes

##Any value <0 = 0 i.e. the gene of interest and any overlapping genes
distance_3prime_genes[distance_3prime_genes < 0]<- 0

###combine start and end distance lists and add 0 in place of MET
## Add zero:
# distance_5prime_genes<-rbind(distance_5prime_genes, 0)
# distance_5prime_genes

##Make sure both tables have same colnames
colnames(distance_5prime_genes)<- "start"

## Join data together

distance_from_target_gene<- rbind(distance_5prime_genes, distance_3prime_genes)
# distance_from_target_gene

return(distance_from_target_gene)

}

##Use function:

distance_from_target_gene_table<- lapply(gene_information_list, function(x) distance_from_target_gene_function(cnv.table = cnv.table, x = x, distance = 2.5e+06))
distance_from_target_gene_table[[2]]

identical(distance_from_target_gene, distance_from_target_gene_table[[2]])


##check distance data is correct

# cnv.table %>%
#   dplyr::filter(CHR == x[[2]]) %>%
#   dplyr::filter(end <= x[[4]]) %>%
#   dplyr::filter(end >= x[[4]] - 2.5e+6) %>%
#   dplyr::select(Gene.Symbol, start, end)
# x[[4]]
# distance_from_target_gene
# 
# cnv.table %>%
#   dplyr::filter(CHR == x[[2]]) %>%
#   dplyr::filter(start >= x[[5]]) %>%
#   dplyr::filter(start <= x[[5]] + 2.5e+6) %>%
#    dplyr::select(Gene.Symbol, start, end)
# x[[5]]
# distance_from_target_gene

########
##Bind all dataframes in list together
dim(distance_from_target_gene_table)
distance_from_target_gene_table<- do.call(rbind, distance_from_target_gene_table)
dim(distance_from_target_gene_table)

dim(gathered.co.deletion.per.target.gene)


sapply(co.deletion.per.target.gene, function(x) nrow(x))
# sapply(distance_from_target_gene_table, function(x) nrow(x))


#########
### join pair-wise distance table to pair-wise proportion of co-deletions table

co_deletions_distance_from_target_gene_plot_table<- cbind(proportion_of_co_deletion = gathered.co.deletion.per.target.gene, distance_from_target_gene = distance_from_target_gene_table)
head(co_deletions_distance_from_target_gene_plot_table)
tail(co_deletions_distance_from_target_gene_plot_table)
dim(co_deletions_distance_from_target_gene_plot_table)

colnames(co_deletions_distance_from_target_gene_plot_table)<- c("Comparison_gene", "Target_gene", "proportion_co_del_amp", "distance_from_target_genes")

##Save data

saveRDS(co_deletions_distance_from_target_gene_plot_table, "/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/BRCA_co_deletion_distance_from_target_gene_plot_table.rds")



##############
###Plot data

ggplot(co_deletions_distance_from_target_gene_plot_table, aes(x = distance_from_target_genes, 
                                                   y = as.numeric(proportion_co_del_amp))) +
  geom_point(size = 1, shape = 1) +
  scale_x_continuous(breaks = c(0, 0.5e+6, 1.0e+6, 1.5e+6, 2.0e+6, 2.5e+6), 
                     labels = c(0, 0.5, 1.0, 1.5, 2.0, 2.5)) +
  scale_y_continuous(breaks = c(0, 0.005, 0.01, 0.015, 0.02,0.025, 0.03)) +
  xlab("Distance from target gene (MB)") +
  ylab("Proportion of tumours with co-deletion") +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))

##Save plot
ggsave("BRCA_co-deletion_2.5MB_distance_from_target_gene_BW.tiff")

## Plot all data coloured by target gene:
ggplot(co_deletions_distance_from_target_gene_plot_table, aes(x = distance_from_target_genes, 
                                                   y = as.numeric(proportion_co_del_amp),
                                                   colour = Target_gene)) +
  geom_point(size = 1, shape = 1) +
  scale_x_continuous(breaks = c(0, 0.5e+6, 1.0e+6, 1.5e+6, 2.0e+6, 2.5e+6), 
                     labels = c(0, 0.5, 1.0, 1.5, 2.0, 2.5)) +
  scale_y_continuous(breaks = c(0, 0.005, 0.01, 0.015, 0.02,0.025, 0.03)) +
  xlab("Distance from target gene (MB)") +
  ylab("Proportion of tumours with co-deletion") +
  scale_color_discrete()+
  labs(colour ="Gene") +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))

##Save plot
ggsave("BRCA_co-deletion_2.5MB_distance_from_target_gene_colour.tiff")

##Facet wrap plot coloured by gene:


head(co_deletions_distance_from_target_gene_plot_table)

ggplot(co_deletions_distance_from_target_gene_plot_table, aes(x = distance_from_target_genes, 
                                                   y = as.numeric(proportion_co_del_amp),
                                                   colour = Target_gene)) +
  geom_point(size = 1, shape = 1) +
  scale_x_continuous(breaks = c(0, 0.5e+6, 1.0e+6, 1.5e+6, 2.0e+6, 2.5e+6), 
                     labels = c(0, 0.5, 1.0, 1.5, 2.0, 2.5)) +
  scale_y_continuous(breaks = c(0, 0.005, 0.01, 0.015, 0.02,0.025, 0.03)) +
  xlab("Distance from target gene (MB)") +
  ylab("Proportion of tumours with co-deletion") +
  scale_color_discrete(guide = FALSE)+
  labs(colour ="Gene") +
  facet_wrap(~Target_gene) +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))

##Save plot
ggsave("BRCA_co-deletion_2.5MB_distance_from_target_gene_colour_wrap.tiff")

##Comment: could add column saying if gene is 5' of 3' of gene? Maybe give it a different shape?




#####################
### Repeat function and plot but with greater distance from gene of interest (5MB up and downstream)



######
### Repeat using lapply for each cancer type. Maybe make one huge table with cancer type as a new column.







