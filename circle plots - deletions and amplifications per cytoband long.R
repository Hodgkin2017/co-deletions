################
### Circle plots of deletions and amplifications per cytoband
###############

library(ggplot2)
library(scales)
library(mgcv)

### Trial circle plots:

####################
###Dummy data

cytoband.name<- c("1q31.2", "2q31.2", "3q31.2", "4q31.2", "5q31.2", "6q31.2", "7q31.2", "8q31.2", "9q31.2", "10q31.2")
cytoband<- rep(seq(1:10), 4)
tumour.type<- rep(seq(1:4), each = 10)
proportion.of.deletions.per.cytoband<- runif(40)
df<- cbind(x = cytoband, y = tumour.type, proportion_deletions = proportion.of.deletions.per.cytoband)
df<- as.data.frame(df)
df

df3<-cbind(x = cytoband.name, y = letters[1:10], proportion_deletions = proportion.of.deletions.per.cytoband)
df3<- as.data.frame(df3)


################
###Trial plotting dummy data

ggplot(df, aes(x,y)) +
  geom_point(aes(size = proportion_deletions)) +
  scale_x_continuous(name="Cytoband", breaks=pretty_breaks(n=10), labels = c("crap", cytoband.name[1:10], "crap")) +
  scale_y_continuous(name="Cancer", breaks=pretty_breaks(n=4), labels = cytoband.name[1:6]) +
  scale_size_area("proportion of\n tumours with\n deletions", breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))


ggplot(df3, aes(x,y)) +
  geom_point(aes(size = proportion_deletions)) +
  scale_x_discrete(name="Cytoband") +
  scale_y_continuous(name="Cancer")

#############
##Another test

set.seed (1234)
rectheat = sample(c(rnorm (10, 5,1), NA, NA), 7*14, replace = T)

dataf <- data.frame (rowv = rep (1:7, 14), columnv = rep(1:14, each = 7),
                     rectheat, circlesize = rectheat*1.5,
                     circlefill =  rectheat*10 )
dataf

##plot another test

ggplot(dataf, aes(y = factor(rowv),
                  x = factor(columnv))) +        ## global aes
  geom_tile(aes(fill = rectheat)) +         ## to get the rect filled
  geom_point(aes(colour = circlefill, 
                 size =circlesize))  +    ## geom_point for circle illusion
  scale_color_gradient(low = "yellow",  
                       high = "red")+       ## color of the corresponding aes
  scale_size(range = c(1, 20))+             ## to tune the size of circles
  theme_bw()

ggplot(dataf, aes(y = factor(rowv),
                  x = factor(columnv))) +
  geom_point(aes(size = circlesize))

dataf <- data.frame (rowv = rep (LETTERS[1:7], 14), columnv = rep(letters[1:14], each = 7),
                     rectheat, circlesize = rectheat*1.5,
                     circlefill =  rectheat*10 )
dataf

############
### Get deletions per cytoband data using events.per.cytoband function

##Target genes:
target.genes<- c("MET", "CDKN2A", "RB1", "WWOX", 
      "LRP1B", "PDE4D", "CCNE1", "TP53",
      "FGFR1", "MYC", "EGFR","WHSC1L1",
      "ERBB2", "MCL1", "MDM2", "CCND1", "ATM",
      "NOTCH1", "PPP2R2A", "BRD4", "ARID1A",
      "STK11", "PARK2")

##New cnv list of cancer types we are most interested in:
short.cnv.list<- cnv.list[c(3, 7, 9, 12, 20, 21, 23, 24, 26, 29, 30)]
length(short.cnv.list)
names(short.cnv.list)

cytoband.table<- unique(short.cnv.list[[1]]$Cytoband)
length(cytoband.table)

length(unique(short.cnv.list$BRCA$Cytoband))

deletions.per.cytoband.circle.plots.table<- data.frame(cytoband = NA, cancer = NA, proportion_deletions = NA)


for (i in 1:length(short.cnv.list)){
  
  cancer.type<- names(short.cnv.list)[i]
  cancer.table<- chromosomal_location(short.cnv.list[[i]])
  cancer.list<- events.per.cytoband(cancer.table,
                              threshold = -1,
                              cytoband_column = 10,
                              column_data_start = 11,
                              chromosome_interval = 0,
                              deletion = TRUE)
  cancer.deletions.per.cytoband<- cbind(cytoband = cytoband.table, cancer = rep(cancer.type, 806), proportion_deletions = cancer.list[[2]]$proportion.of.deletions)
  
  deletions.per.cytoband.circle.plots.table<-rbind(deletions.per.cytoband.circle.plots.table, cancer.deletions.per.cytoband) 
  
  
}

## delete first row of deletions.per.cytoband.circle.plots.table

deletions.per.cytoband.circle.plots.table<-deletions.per.cytoband.circle.plots.table[-1,]
head(deletions.per.cytoband.circle.plots.table)
tail(deletions.per.cytoband.circle.plots.table)
nrow(deletions.per.cytoband.circle.plots.table)

## Plot circle plot

df2<- deletions.per.cytoband.circle.plots.table[c(1:10, 807:817),]

ggplot(deletions.per.cytoband.circle.plots.table, aes(y = factor(cancer),
                  x = factor(cytoband))) +
  xlab("Cytoband") +
  ylab("Cancer") +
  geom_point(aes(size = as.numeric(proportion_deletions))) + 
  scale_size_area("Proportion of\ntumours with\ndeletions")

##Save plot
ggsave("deletions_per_cytoband_circle_plots.tiff")

#############
### Repeat for certain cytobands only:


# short.cnv.list[[1]] %>% 
#   dplyr::filter(Gene.Symbol == "CDKN2A") %>%
#   dplyr::select(Cytoband)

##Get cytobands of target genes
target.genes.cytoband <- short.cnv.list[[1]] %>% 
  dplyr::filter(Gene.Symbol %in% target.genes) %>%
  dplyr::select(Cytoband) %>%
  t() %>% #required to convert output from data.frame to vector
  as.character()

##create table of proportion of tumours with deletions per cytoband for target genes
target.genes.deletions.per.cytoband.circle.plots.table<- deletions.per.cytoband.circle.plots.table %>% 
  dplyr::filter(cytoband %in% target.genes.cytoband)

##Plot:
##Comment: in order to get the cytobands to be in numerical order need to use levels =
ggplot(target.genes.deletions.per.cytoband.circle.plots.table, aes(y = factor(cancer),
                                                      x = factor(cytoband, levels = unique(target.genes.deletions.per.cytoband.circle.plots.table$cytoband)))) +
  xlab("Cytoband") +
  ylab("Cancer") +
  geom_point(aes(size = as.numeric(proportion_deletions))) + 
  scale_size_area("Proportion of\nevents per\ncytoband")+
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
        ) + #theme_bw()
##Save plot
ggsave("target_genes_deletions_per_cytoband_circle_plots.tiff")



##############
### Get deletions per genomic region data using events.per.cytoband function?



################
### Circle plot: Proportion of co-deletions per cytoband
###############

##Think about how to do it: 
##Use my co-deletion function..remind myself how it works

cytoband.table
co.deletions.per.cytoband.circle.plots.table<- data.frame(cytoband = NA, cancer = NA, proportion_deletions = NA)
co.deletions.per.cytoband.circle.plots.table

for (i in 1:length(short.cnv.list)){

cancer.type<- names(short.cnv.list)[i]
cancer.table<- chromosomal_location(short.cnv.list[[i]])
cancer.list<- sapply(cytoband.table, function(x) co.deletion_co.amplification_matrix(cancer.table, column_start = 11, threshold = -1, Cytoband = TRUE, selection_criteria = x, deletion = TRUE, normalisation = "frequency.for.whole.sample", remove_NA = FALSE))
#cancer.list2<-unlist(cancer.list)
cancer.deletions.per.cytoband<- cbind(cytoband = cytoband.table, cancer = rep(cancer.type, 806), proportion_deletions = cancer.list)
co.deletions.per.cytoband.circle.plots.table<-rbind(co.deletions.per.cytoband.circle.plots.table, cancer.deletions.per.cytoband) 
print(cancer.type)
}

##Problem solving empty cytoband entries:
# which(!(names(cancer.list) %in% names(cancer.list2)), arr.ind = TRUE)
# cancer.list[82:85]
# test<- sapply("2p11.1", function(x) co.deletion_co.amplification_matrix(cancer.table, column_start = 11, threshold = -1, Cytoband = TRUE, selection_criteria = x, deletion = TRUE, normalisation = "frequency.for.whole.sample"))
# test
# co.deletion_co.amplification_matrix(cancer.table, column_start = 11, threshold = -1, Cytoband = TRUE, selection_criteria = "2p11.1", deletion = TRUE, normalisation = "none")
# names(cancer.table)
# cancer.table %>% dplyr::filter(Cytoband == "2p11.1") %>% dplyr::select(Gene.Symbol, CHRLOC, start)
# co.deletion_co.amplification_matrix(cancer.table, column_start = 11, threshold = -1, Cytoband = TRUE, selection_criteria = "2p11.1", deletion = TRUE, normalisation = "frequency.for.whole.sample", remove_NA = FALSE)


## delete first row of deletions.per.cytoband.circle.plots.table
dim(co.deletions.per.cytoband.circle.plots.table)
co.deletions.per.cytoband.circle.plots.table<-co.deletions.per.cytoband.circle.plots.table[-1,]
head(co.deletions.per.cytoband.circle.plots.table)
tail(co.deletions.per.cytoband.circle.plots.table)
nrow(co.deletions.per.cytoband.circle.plots.table)

## Plot circle plot

df2<- co.deletions.per.cytoband.circle.plots.table[c(1:10, 807:817),]

ggplot(co.deletions.per.cytoband.circle.plots.table, aes(y = factor(cancer),
                                                      x = factor(cytoband))) +
  xlab("Cytoband") +
  ylab("Cancer") +
  geom_point(aes(size = as.numeric(proportion_deletions))) + 
  scale_size_area("Proportion of\nco-deletion\nevents") +
  
##Save plot
ggsave("co_deletions_per_cytoband_circle_plots.tiff")

#############
### Repeat co-deltion plot for certain cytobands only:

##Get cytobands of target genes
target.genes.cytoband <- short.cnv.list[[1]] %>% 
  dplyr::filter(Gene.Symbol %in% target.genes) %>%
  dplyr::select(Cytoband) %>%
  t() %>% #required to convert output from data.frame to vector
  as.character()

##create table of proportion of tumours with deletions per cytoband for target genes
target.genes.co.deletions.per.cytoband.circle.plots.table<- co.deletions.per.cytoband.circle.plots.table %>% 
  dplyr::filter(cytoband %in% target.genes.cytoband)

##Plot:
##Comment: in order to get the cytobands to be in numerical order need to use levels =
ggplot(target.genes.co.deletions.per.cytoband.circle.plots.table, aes(y = factor(cancer),
                                                                   x = factor(cytoband, levels = unique(target.genes.deletions.per.cytoband.circle.plots.table$cytoband)))) +
  xlab("Cytoband") +
  ylab("Cancer") +
  geom_point(aes(size = as.numeric(proportion_deletions))) + 
  scale_size_area("Proportion of\nco-deletion\nevents")+
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
  ) + #theme_bw()

##Save plot
ggsave("target_genes_co_deletions_per_cytoband_circle_plots.tiff")

##############
###save objects:

##Keeps object name:
# save(mod, file = "mymodel.rda")
# load(file = "mymodel.rda")

##Allows you to re-assign objects name:
saveRDS(co.deletions.per.cytoband.circle.plots.table, "/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/total.co-deletion.events.per.cytoband.rds")
saveRDS(deletions.per.cytoband.circle.plots.table, "/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/total.deletion.events.per.cytoband.rds")
saveRDS(short.cnv.list, "/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/target.cancer.list.rds")


target.genes.co.deletions.per.cytoband.circle.plots.table2 <- readRDS("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/total.co-deletion.events.per.cytoband.rds")
identical(target.genes.co.deletions.per.cytoband.circle.plots.table, target.genes.co.deletions.per.cytoband.circle.plots.table2, ignore.environment = TRUE)
identical(target.genes.co.deletions.per.cytoband.circle.plots.table, target.genes.co.deletions.per.cytoband.circle.plots.table2)
