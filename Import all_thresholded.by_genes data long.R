#####################
### Importing all_thresholded.by_genes.txt and creating new cnv.list and short.cnv.list
#####################

##Use import.files.from.directories function to import data:

path<-"/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/unzipped original broad TCGA CNV data"
file.to.import<-"all_thresholded.by_genes.txt"

threshold_cnv_list<- import.files.from.directories(path, file.to.import)
length(threshold_cnv_list)
dim(threshold_cnv_list[[1]])
acc.cnv.list.loc<- chromosomal_location(threshold_cnv_list[[1]])

##Save object
saveRDS(threshold_cnv_list, "/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/threshold_cnv_list.rds")



############
### Create a table containing all tumours 
threshold_CNV_all_table<-join.cnv.datasets(threshold_cnv_list, column = 4)
dim(threshold_CNV_all_table)
which(is.na(threshold_CNV_all_table), arr.ind = T)
threshold_CNV_all_table[1:2, 1:4]

##Save object
saveRDS(threshold_CNV_all_table, "/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/threshold_CNV_all_table.rds")


#############
### Append chromosomal location to threshold_CNV_all_table
threshold_CNV_all_table_loc<- chromosomal_location(threshold_CNV_all_table)

##Save object
saveRDS(threshold_CNV_all_table_loc, "/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/threshold_CNV_all_table_loc.rds")

#############
### Append chromosomal location to each CNV table in list: threshold_cnv_list

threshold_cnv_list_loc<- lapply(threshold_cnv_list, function(x) chromosomal_location(x))
identical(acc.cnv.list.loc, threshold_cnv_list_loc[[1]])
threshold_cnv_list_loc[[1]][1:2, 1:12]

##Save object
saveRDS(threshold_cnv_list_loc, "/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/threshold_cnv_list_loc.rds")

############
### Create shorter CNV.list with cancers that are particularly interesting to AZ:

##Target genes:
x<- c("MET", "CDKN2A", "RB1", "WWOX", 
      "LRP1B", "PDE4D", "CCNE1", "TP53",
      "FGFR1", "MYC", "EGFR","WHSC1L1",
      "ERBB2", "MCL1", "MDM2", "CCND1", "ATM",
      "NOTCH1", "PPP2R2A", "BRD4", "ARID1A",
      "STK11", "PARK2")

##New cnv list of cancer types we are most interested in:
threshold_short_cnv_list<- threshold_cnv_list[c(3, 7, 9, 12, 20, 21, 23, 24, 26, 29, 30)]
length(threshold_short_cnv_list)
names(threshold_short_cnv_list)
##Save object
saveRDS(threshold_short_cnv_list, "/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/threshold_short_cnv_list.rds")


threshold_short_cnv_list_loc<- threshold_cnv_list_loc[c(3, 7, 9, 12, 20, 21, 23, 24, 26, 29, 30)]
length(threshold_short_cnv_list_loc)
names(threshold_short_cnv_list_loc)
##Save object
saveRDS(threshold_short_cnv_list_loc, "/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/threshold_short_cnv_list_loc.rds")


# -        Breast invasive carcinoma (BRCA)
# -        Esophageal cancer (ESCA)
# -        Head and neck squamous cell carcinoma (HNSC)
# -        Lung adenocarcinoma (LUAD)
# -        Lung squamous cell carcinoma (LUSC)
# -        Ovarian serous cystadenocarcinoma (OV)
# -        Pancreatic ductal adenocarcinoma (PAAD)
# -        Stomach adenocarcinoma (STAD)
# -        Skin cutaneous melanoma (SKCM)
# -        Prostate adenocarcinoma (PRAD)
# -        Colorectal adenocarcinoma (COADREAD)

