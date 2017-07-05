################
### Import and test clinical data
#################

#Summary: All CNV data sets (threshold_short_cnv_list) have clinical data for all tumours excepts for BRCA which is missing sample
# 518 and OV which is missing 6 samples: 183 184 304 382 393 468.



#####################
##Import BRCA clinical data
BRCA_clinical<- read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/clinical/BRCA_clinical.tsv", header = TRUE, stringsAsFactors = FALSE)
BRCA_clinical[1:5, 1:10]
colnames(threshold_short_cnv_list[[1]])
BRCA_cnv_names<- threshold_short_cnv_list$BRCA %>% 
  dplyr::select(-Gene.Symbol, -Locus.ID, -Cytoband) %>% 
  colnames()
length(BRCA_cnv_names)

length(BRCA_clinical$tcga_participant_barcode)
class(BRCA_cnv_names)
class(BRCA_clinical$tcga_participant_barcode)
patients<- BRCA_clinical$tcga_participant_barcode
BRCA_cnv_names[1:5]
BRCA_cnv_names<- substr(BRCA_cnv_names, 0, 12)
BRCA_cnv_names<- gsub("\\.", "-", BRCA_cnv_names)

# a<- c(1,2,3,4)
# b<- c(2,3,5,6)
# which(a %in% b)

##Test BRCA clinical data

#which(BRCA_cnv_names %in% BRCA_cnv_names)
which(BRCA_cnv_names %in% patients)
length(which(BRCA_cnv_names %in% BRCA_clinical$tcga_participant_barcode))
which(!(BRCA_cnv_names %in% BRCA_clinical$tcga_participant_barcode))
length(which(!(BRCA_cnv_names %in% BRCA_clinical$tcga_participant_barcode)))
BRCA_cnv_names[518]
##Comment: only BRCA CNV sample 518 does not have clinical data


######################
##Import COADREAD clinical data
COADREAD_clinical<- read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/clinical/COADREAD_clinical.tsv", header = TRUE, stringsAsFactors = FALSE)
# COADREAD_clinical[1:5, 1:10]
# colnames(threshold_short_cnv_list$COADREAD)

clinical_tumour_ID<-COADREAD_clinical$tcga_participant_barcode 

## Extract patient tumour ID from CNV data and convert to same format as clinical data
cnv_tumour_ID<- threshold_short_cnv_list$COADREAD %>% 
  dplyr::select(-Gene.Symbol, -Locus.ID, -Cytoband) %>% 
  colnames() %>% 
  substr(0, 12) %>%
  gsub("\\.", "-", .)

cnv_tumour_ID[1:5]
clinical_tumour_ID[1:5]

length(cnv_tumour_ID)
length(clinical_tumour_ID)

##Test number of CNV tumour IDs have clinical data:
which(cnv_tumour_ID %in% clinical_tumour_ID)
length(which(cnv_tumour_ID %in% clinical_tumour_ID))
##Test number of CNV tumour IDs have clinical data:
which(!(cnv_tumour_ID %in% clinical_tumour_ID))
length(which(!(cnv_tumour_ID %in% clinical_tumour_ID)))
#cnv_tumour_ID[518]

##Comment: ALl COADREAD CNV samples have clinical data

######################
##Import ESCA clinical data
ESCA_clinical<- read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/clinical/ESCA_clinical.tsv", header = TRUE, stringsAsFactors = FALSE)

clinical_tumour_ID<-ESCA_clinical$tcga_participant_barcode 

## Extract patient tumour ID from CNV data and convert to same format as clinical data
cnv_tumour_ID<- threshold_short_cnv_list$ESCA %>% 
  dplyr::select(-Gene.Symbol, -Locus.ID, -Cytoband) %>% 
  colnames() %>% 
  substr(0, 12) %>%
  gsub("\\.", "-", .)

cnv_tumour_ID[1:5]
clinical_tumour_ID[1:5]

length(cnv_tumour_ID)
length(clinical_tumour_ID)

##Test number of CNV tumour IDs have clinical data:
which(cnv_tumour_ID %in% clinical_tumour_ID)
length(which(cnv_tumour_ID %in% clinical_tumour_ID))
##Test number of CNV tumour IDs have clinical data:
which(!(cnv_tumour_ID %in% clinical_tumour_ID))
length(which(!(cnv_tumour_ID %in% clinical_tumour_ID)))
#cnv_tumour_ID[518]

##Comment: ALl ESCA CNV samples have clinical data

######################
##Import HNSC clinical data
HNSC_clinical<- read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/clinical/HNSC_clinical.tsv", header = TRUE, stringsAsFactors = FALSE)

clinical_tumour_ID<-HNSC_clinical$tcga_participant_barcode 

## Extract patient tumour ID from CNV data and convert to same format as clinical data
cnv_tumour_ID<- threshold_short_cnv_list$HNSC %>% 
  dplyr::select(-Gene.Symbol, -Locus.ID, -Cytoband) %>% 
  colnames() %>% 
  substr(0, 12) %>%
  gsub("\\.", "-", .)

cnv_tumour_ID[1:5]
clinical_tumour_ID[1:5]

length(cnv_tumour_ID)
length(clinical_tumour_ID)

##Test number of CNV tumour IDs have clinical data:
which(cnv_tumour_ID %in% clinical_tumour_ID)
length(which(cnv_tumour_ID %in% clinical_tumour_ID))
##Test number of CNV tumour IDs have clinical data:
which(!(cnv_tumour_ID %in% clinical_tumour_ID))
length(which(!(cnv_tumour_ID %in% clinical_tumour_ID)))
#cnv_tumour_ID[518]

##Comment: ALl HNSC CNV samples have clinical data


######################
##Import LUAD clinical data
LUAD_clinical<- read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/clinical/LUAD_clinical.tsv", header = TRUE, stringsAsFactors = FALSE)

clinical_tumour_ID<-LUAD_clinical$tcga_participant_barcode 

## Extract patient tumour ID from CNV data and convert to same format as clinical data
cnv_tumour_ID<- threshold_short_cnv_list$LUAD %>% 
  dplyr::select(-Gene.Symbol, -Locus.ID, -Cytoband) %>% 
  colnames() %>% 
  substr(0, 12) %>%
  gsub("\\.", "-", .)

cnv_tumour_ID[1:5]
clinical_tumour_ID[1:5]

length(cnv_tumour_ID)
length(clinical_tumour_ID)

##Test number of CNV tumour IDs have clinical data:
which(cnv_tumour_ID %in% clinical_tumour_ID)
length(which(cnv_tumour_ID %in% clinical_tumour_ID))
##Test number of CNV tumour IDs have clinical data:
which(!(cnv_tumour_ID %in% clinical_tumour_ID))
length(which(!(cnv_tumour_ID %in% clinical_tumour_ID)))

##Comment: ALl LUAD CNV samples have clinical data


######################
##Import LUSC clinical data
LUSC_clinical<- read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/clinical/LUSC_clinical.tsv", header = TRUE, stringsAsFactors = FALSE)

clinical_tumour_ID<-LUSC_clinical$tcga_participant_barcode 

## Extract patient tumour ID from CNV data and convert to same format as clinical data
cnv_tumour_ID<- threshold_short_cnv_list$LUSC %>% 
  dplyr::select(-Gene.Symbol, -Locus.ID, -Cytoband) %>% 
  colnames() %>% 
  substr(0, 12) %>%
  gsub("\\.", "-", .)

cnv_tumour_ID[1:5]
clinical_tumour_ID[1:5]

length(cnv_tumour_ID)
length(clinical_tumour_ID)

##Test number of CNV tumour IDs have clinical data:
which(cnv_tumour_ID %in% clinical_tumour_ID)
length(which(cnv_tumour_ID %in% clinical_tumour_ID))
##Test number of CNV tumour IDs have clinical data:
which(!(cnv_tumour_ID %in% clinical_tumour_ID))
length(which(!(cnv_tumour_ID %in% clinical_tumour_ID)))

##Comment: ALl LUSC CNV samples have clinical data


######################
##Import OV clinical data
OV_clinical<- read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/clinical/OV_clinical.tsv", header = TRUE, stringsAsFactors = FALSE)

clinical_tumour_ID<-OV_clinical$tcga_participant_barcode 

## Extract patient tumour ID from CNV data and convert to same format as clinical data
cnv_tumour_ID<- threshold_short_cnv_list$OV %>% 
  dplyr::select(-Gene.Symbol, -Locus.ID, -Cytoband) %>% 
  colnames() %>% 
  substr(0, 12) %>%
  gsub("\\.", "-", .)

cnv_tumour_ID[1:5]
clinical_tumour_ID[1:5]

length(cnv_tumour_ID)
length(clinical_tumour_ID)

##Test number of CNV tumour IDs have clinical data:
which(cnv_tumour_ID %in% clinical_tumour_ID)
length(which(cnv_tumour_ID %in% clinical_tumour_ID))
##Test number of CNV tumour IDs have clinical data:
which(!(cnv_tumour_ID %in% clinical_tumour_ID))
length(which(!(cnv_tumour_ID %in% clinical_tumour_ID)))

##Comment: 6 OV CNV samples do not have clinical data (samples: 183 184 304 382 393 468)


######################
##Import PAAD clinical data
PAAD_clinical<- read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/clinical/PAAD_clinical.tsv", header = TRUE, stringsAsFactors = FALSE)

clinical_tumour_ID<-PAAD_clinical$tcga_participant_barcode 

## Extract patient tumour ID from CNV data and convert to same format as clinical data
cnv_tumour_ID<- threshold_short_cnv_list$PAAD %>% 
  dplyr::select(-Gene.Symbol, -Locus.ID, -Cytoband) %>% 
  colnames() %>% 
  substr(0, 12) %>%
  gsub("\\.", "-", .)

cnv_tumour_ID[1:5]
clinical_tumour_ID[1:5]

length(cnv_tumour_ID)
length(clinical_tumour_ID)

##Test number of CNV tumour IDs have clinical data:
which(cnv_tumour_ID %in% clinical_tumour_ID)
length(which(cnv_tumour_ID %in% clinical_tumour_ID))
##Test number of CNV tumour IDs have clinical data:
which(!(cnv_tumour_ID %in% clinical_tumour_ID))
length(which(!(cnv_tumour_ID %in% clinical_tumour_ID)))

##Comment: ALl PAAD CNV samples have clinical data


######################
##Import PRAD clinical data
PRAD_clinical<- read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/clinical/PRAD_clinical.tsv", header = TRUE, stringsAsFactors = FALSE)

clinical_tumour_ID<-PRAD_clinical$tcga_participant_barcode 

## Extract patient tumour ID from CNV data and convert to same format as clinical data
cnv_tumour_ID<- threshold_short_cnv_list$PRAD %>% 
  dplyr::select(-Gene.Symbol, -Locus.ID, -Cytoband) %>% 
  colnames() %>% 
  substr(0, 12) %>%
  gsub("\\.", "-", .)

cnv_tumour_ID[1:5]
clinical_tumour_ID[1:5]

length(cnv_tumour_ID)
length(clinical_tumour_ID)

##Test number of CNV tumour IDs have clinical data:
which(cnv_tumour_ID %in% clinical_tumour_ID)
length(which(cnv_tumour_ID %in% clinical_tumour_ID))
##Test number of CNV tumour IDs have clinical data:
which(!(cnv_tumour_ID %in% clinical_tumour_ID))
length(which(!(cnv_tumour_ID %in% clinical_tumour_ID)))

##Comment: ALl PRAD CNV samples have clinical data



######################
##Import SKCM clinical data
SKCM_clinical<- read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/clinical/SKCM_clinical.tsv", header = TRUE, stringsAsFactors = FALSE)

clinical_tumour_ID<-SKCM_clinical$tcga_participant_barcode 

## Extract patient tumour ID from CNV data and convert to same format as clinical data
cnv_tumour_ID<- threshold_short_cnv_list$SKCM %>% 
  dplyr::select(-Gene.Symbol, -Locus.ID, -Cytoband) %>% 
  colnames() %>% 
  substr(0, 12) %>%
  gsub("\\.", "-", .)

cnv_tumour_ID[1:5]
clinical_tumour_ID[1:5]

length(cnv_tumour_ID)
length(clinical_tumour_ID)

##Test number of CNV tumour IDs have clinical data:
which(cnv_tumour_ID %in% clinical_tumour_ID)
length(which(cnv_tumour_ID %in% clinical_tumour_ID))
##Test number of CNV tumour IDs have clinical data:
which(!(cnv_tumour_ID %in% clinical_tumour_ID))
length(which(!(cnv_tumour_ID %in% clinical_tumour_ID)))

##Comment: ALl SKCM CNV samples have clinical data


######################
##Import STAD clinical data
STAD_clinical<- read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/clinical/STAD_clinical.tsv", header = TRUE, stringsAsFactors = FALSE)

clinical_tumour_ID<-STAD_clinical$tcga_participant_barcode 

## Extract patient tumour ID from CNV data and convert to same format as clinical data
cnv_tumour_ID<- threshold_short_cnv_list$STAD %>% 
  dplyr::select(-Gene.Symbol, -Locus.ID, -Cytoband) %>% 
  colnames() %>% 
  substr(0, 12) %>%
  gsub("\\.", "-", .)

cnv_tumour_ID[1:5]
clinical_tumour_ID[1:5]

length(cnv_tumour_ID)
length(clinical_tumour_ID)

##Test number of CNV tumour IDs have clinical data:
which(cnv_tumour_ID %in% clinical_tumour_ID)
length(which(cnv_tumour_ID %in% clinical_tumour_ID))
##Test number of CNV tumour IDs have clinical data:
which(!(cnv_tumour_ID %in% clinical_tumour_ID))
length(which(!(cnv_tumour_ID %in% clinical_tumour_ID)))

##Comment: ALl STAD CNV samples have clinical data


#############################
### Tidy data to keep only interesting columns and then combine datasets together into 
#one long table with a new column for cancer type
#############################



  

###############
###Join tables together to make one long table?

dim(STAD_clinical)
dim(SKCM_clinical)
##Not all tables have the same number of columns...can dplyr help?

ALL_clinical_table<- rbind()


##Save data


#######################
### Create a list of Clinical data using fbget
######################

setwd("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/clinical/fbget/")

file_names<-dir()
file_names
clinical_fbget_list <- vector("list", length(file_names))
for (i in 1: length(file_names)){
  clinical_fbget_list[[i]]<-read.delim(file_names[i], header = TRUE, stringsAsFactors = FALSE)
  print(file_names[i])
}
names(clinical_fbget_list)<- file_names
clinical_fbget_list[[1]][1:2, 1:4]

##Save data
saveRDS(clinical_fbget_list, file = "/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/clinical_fbget_list.rds")

#######################
### Create a list of Clinical data obtained from CBioPortal
######################

setwd("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/clinical/cbioportal/")

##Create list:
file_names<-dir()
file_names
clinical_cbioportal_list <- vector("list", length(file_names))
for (i in 1: length(file_names)){
  clinical_cbioportal_list[[i]]<-read.delim(file_names[i], header = TRUE, stringsAsFactors = FALSE)
  print(file_names[i])
}
names(clinical_cbioportal_list)<- file_names

clinical_cbioportal_list[[1]][1:2, 1:5]

lapply(clinical_cbioportal_list, function(x) dim(x))

##Loop to Remove duplicate tumour entries

for (i in 1:length(clinical_cbioportal_list)) {
##Get sample IDs
patient_IDs<- clinical_cbioportal_list[[i]]$Patient.ID

##Find duplicate patients IDS
indicies<- which(duplicated(patient_IDs))
indicies

#clinical_cbioportal_list[[i]][indicies, 4:6]

if (any(indicies)) {

##Find original and repeated patient IDS
repeated_samples<- patient_IDs[indicies]
indicies<- which(patient_IDs %in% repeated_samples)
#indicies

###Keep column with less NA's
sample_pairs<- c()
##find sample pairs:
for (j in 1: length(indicies)){
  sampleID<- patient_IDs[indicies[j]]
  sampleID<- grep(sampleID,patient_IDs) 
  sample_pairs<- cbind(sample_pairs, sampleID)
}
#sample_pairs


##keep unique sample_pairs columns only
sample_pairs<- sample_pairs[,!duplicated(sample_pairs, MARGIN = 2)]
sample_pairs<- as.data.frame(sample_pairs)
#sample_pairs


## Identify which column has the lowest number of NAs
n<- nrow(sample_pairs)
df<- matrix(NA, nrow = 1, ncol = n)
#df
which_row_to_remove<- c()

for (a in 1:ncol(sample_pairs)){
  
  for (b in 1:n){
    df[1,b]<- sum(is.na(clinical_cbioportal_list[[i]][sample_pairs[b,a],]))
  }
  which_row<- which(max(df[1,]) == df[1,])
  which_row_to_remove<-c(which_row_to_remove, sample_pairs[which_row, a])
}
# sample_pairs
# which_row_to_remove
which_row_to_remove<- as.vector(which_row_to_remove)

clinical_cbioportal_list[[i]]<- clinical_cbioportal_list[[i]][-which_row_to_remove,]
} else {}

print(names(clinical_cbioportal_list[i]))
}

sapply(clinical_cbioportal_list, function(x) dim(x))
sapply(clinical_fbget_list, function(x) dim(x))
sapply(threshold_CNV_all_table_loc ,function(x) dim(x)-11)

###############
###Test I have removed duplicate samples from clinical_cbioportal_list:

##Get sample IDs
for (i in 1:length(clinical_cbioportal_list)){
patient_IDs<- clinical_cbioportal_list[[i]]$Patient.ID

indicies<- which(duplicated(patient_IDs))
print(names(clinical_cbioportal_list[i]))
print(indicies)
}

#############
###Test I have removed duplicate samples from clinical_fbget_list:

##Get sample IDs
for (i in 1:length(clinical_fbget_list)){
  patient_IDs<- clinical_fbget_list[[i]]$Patient.ID
  
  indicies<- which(duplicated(patient_IDs))
  print(names(clinical_fbget_list[i]))
  print(indicies)
}

#################
###Compare number of 
sapply(clinical_cbioportal_list, function(x) nrow(x))
sapply(clinical_fbget_list, function(x) nrow(x))
sapply(threshold_short_cnv_list_loc, function(x) ncol(x)-11)
