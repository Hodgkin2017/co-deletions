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

setwd("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions")
saveRDS(clinical_cbioportal_list, file = "./R workspaces/clinical_cbioportal_list")


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


################
##Create a list for each cancer type containing all survival data: Overall Survival data from fbget 
#and disease free survival from cbioportal 

survival_list

clinical_cbioportal_list[[1]]$Patient.ID

x<- clinical_cbioportal_list[[1]]
x$Patient.ID

  
###################
### Create table using BRCA cbioportal dataset to plot survival data:

##Get patient IDs:
clinical_IDS<- x$Patient.ID

# get the columns that contain data we can use: days to death, new tumor event, last day contact to....
new_tum_collapsed<- as.numeric(x$Disease.Free..Months.)

#Convert months to days:
new_tum_collapsed<- new_tum_collapsed*(365.25/12)
new_tum_collapsed

# do the same for death
death_collapsed<- as.numeric(x$Death.from.Initial.Pathologic.Diagnosis.Date)
death_collapsed

# and days last follow up here we take the most recent which is the max number
followup_collapsed<- as.numeric(x$Days.to.Last.Followup)
followup_collapsed

# and put everything together
all_clin <- data.frame(new_tum_collapsed, death_collapsed,followup_collapsed)
colnames(all_clin) <- c('new_tumor_days', 'death_days', 'followUp_days')
head(all_clin, 20)
tail(all_clin, 20)

# create vector with time to new tumor containing data to censor for new_tumor
all_clin$new_time <- c()
for (i in 1:nrow(all_clin)){
  ## Combine new_tumor days and followUp_days
  all_clin$new_time[i] <- ifelse (is.na(as.numeric(as.character(all_clin$new_tumor_days)))[i],
                                  as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$new_tumor_days))[i])
  ## If no followUp_days or new_tumor_days then take death_days
  all_clin$new_time[i]<- ifelse (is.na(as.numeric(as.character(all_clin$new_time)))[i],
                                 as.numeric(as.character(all_clin$death_days))[i], as.numeric(as.character(all_clin$new_time))[i])
  
}
all_clin$new_time

# create vector time to death containing values to censor for death
all_clin$new_death <- c()
for (i in 1:nrow(all_clin)){
  all_clin$new_death[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                    as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
}

head(all_clin, 20)
tail(all_clin, 20)

# create vector for death censoring
table(x$Patient.s.Vital.Status)
# alive dead
# 944   152

all_clin$death_event <- ifelse(x$Patient.s.Vital.Status == 'Alive', 0,1)

#finally add row.names to clinical
length(clinical_IDS)
nrow(all_clin)
rownames(all_clin) <- clinical_IDS

head(all_clin, 40)


all_clin_cbio<- all_clin

################
###fbget clinical data
x<- clinical_fbget_list[[1]]

clinical_IDS<- x$tcga_participant_barcode
clinical_IDS

### get the columns that contain data we can use: days to death, new tumor event, last day contact to....

## do the same for death
death_collapsed<- as.numeric(x$days_to_death)
death_collapsed

## and days last follow up here we take the most recent which is the max number
followup_collapsed<- as.numeric(x$days_to_last_followup)
followup_collapsed

## and put everything together
all_clin <- data.frame(death_collapsed,followup_collapsed)
colnames(all_clin) <- c('death_days', 'followUp_days')
head(all_clin, 20)
tail(all_clin, 20)

# create vector time to death containing values to censor for death
all_clin$new_death <- c()
for (i in 1:nrow(all_clin)){
  all_clin$new_death[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                    as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
}

head(all_clin, 20)
tail(all_clin, 20)

# create vector for death censoring
table(x$vital_status)
# alive dead
# 944   152

all_clin$death_event <- ifelse(x$vital_status == 'alive', 0,1)

#finally add row.names to clinical
length(clinical_IDS)
nrow(all_clin)
rownames(all_clin) <- clinical_IDS

head(all_clin, 40)

############
###Compare datasets:

all_clin<-cbind(patient_ID = rownames(all_clin), all_clin)
head(all_clin, 40)
dim(all_clin)
all_clin_cbio<-cbind(patient_ID = rownames(all_clin_cbio), all_clin_cbio)
head(all_clin_cbio, 40)
dim(all_clin_cbio)
all_clin_joined<- dplyr::full_join(all_clin, all_clin_cbio, by = "patient_ID")
dim(all_clin_joined)
head(all_clin_joined, 40)

all_clin_joined2 <- all_clin_joined %>% 
  dplyr::filter(!is.na(death_event.x)) %>%
  dplyr::filter(!is.na(death_event.y))

dim(all_clin_joined)
dim(all_clin_joined2)
identical(all_clin_joined2$death_event.x, all_clin_joined2$death_event.y)
all.equal(all_clin_joined2$death_event.x, all_clin_joined2$death_event.y)
all.equal(all_clin_joined2$death_days.x, all_clin_joined2$death_days.y)
all.equal(all_clin_joined2$new_death.x, all_clin_joined2$new_death.y)
which(is.na(all_clin_joined2$new_death.y))
identical(all_clin_joined2$followUp_days.x, all_clin_joined2$followUp_days.y)
which(is.na(all_clin_joined2$followUp_days.x))
which(is.na(all_clin_joined2$followUp_days.y))


#################
###Create list of all survival data:

##empty list

clinical_survival<- vector("list", length(clinical_fbget_list))
clinical_survival

##Loop:
for ( j in 1: length(clinical_fbget_list)){
  
  ##Get fbget clinical data
  x<- clinical_fbget_list[[j]]
  
  clinical_IDS<- x$tcga_participant_barcode
  
  ## Get days to death data:
  death_collapsed<- as.numeric(x$days_to_death)
  
  ## Get days to follow up data
  followup_collapsed<- as.numeric(x$days_to_last_followup)
  
  ## combine data:
  all_clin <- data.frame(death_collapsed,followup_collapsed)
  colnames(all_clin) <- c('death_days', 'followUp_days')
  
  ## create new column with last follow up or days to death data:
  all_clin$new_death <- c()
  for (i in 1:nrow(all_clin)){
    all_clin$new_death[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                      as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
  }
  
  ## Create new column with 1= patient died and 0 - patient survived
  all_clin$death_event <- ifelse(x$vital_status == 'alive', 0,1)
  
  ##Add row names to all_clin and a new column with clinical IDs
  rownames(all_clin) <- clinical_IDS
  all_clin<- cbind(patient_IDs = rownames(all_clin), all_clin)
  all_clin$patient_IDs<- as.character(all_clin$patient_IDs)
  
  ###############
  ##Get disease free survival days from cbioportal clinical data:
  
  x<- clinical_cbioportal_list[[j]]
  
  ##Get patient IDs:
  clinical_IDS<- x$Patient.ID
  
  # get the number of month to new tumour
  new_tum_collapsed<- as.numeric(x$Disease.Free..Months.)
  
  #Convert number of month to new tumour to days:
  new_tum_collapsed<- new_tum_collapsed*(365.25/12)
  
  new_tum_collapsed<- cbind(patient_IDs = clinical_IDS, new_tumor_days = new_tum_collapsed)
  new_tum_collapsed<- as.data.frame(new_tum_collapsed)
  new_tum_collapsed$patient_IDs<- as.character(new_tum_collapsed$patient_IDs)
  
  ##Join days to new tumour and overall survival together:
  all_clin_joined<- dplyr::full_join(new_tum_collapsed, all_clin, by = "patient_IDs")
  
  ##Calculate disease free survival time:
  all_clin_joined$disease_free_survival<- c()
  for (i in 1:nrow(all_clin_joined)){
      ## Combine new_tumor days and followUp_days
    all_clin_joined$disease_free_survival[i] <- ifelse (is.na(as.numeric(as.character(all_clin_joined$new_tumor_days)))[i],
                                      as.numeric(as.character(all_clin_joined$followUp_days))[i],as.numeric(as.character(all_clin_joined$new_tumor_days))[i])
      ## If no followUp_days or new_tumor_days then take death_days
    all_clin_joined$disease_free_survival[i]<- ifelse (is.na(as.numeric(as.character(all_clin_joined$disease_free_survival)))[i],
                                     as.numeric(as.character(all_clin_joined$death_days))[i], as.numeric(as.character(all_clin_joined$disease_free_survival))[i])

    }
  
  clinical_survival[[j]]<- all_clin_joined
  
}

names(clinical_survival)<- names(clinical_fbget_list)
head(clinical_survival[[1]], 40)
head(clinical_survival[[2]], 40)
head(clinical_survival[[11]], 40)

##save data:
clinical_survival_list<- clinical_survival
#saveRDS(clinical_survival_list, file = "./R workspaces/clinical_survival_list")

###########################################################

#################
### Create new longer lists of survival data:
################

#######################
### Create a list of Clinical data using fbget
######################

setwd("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/clinical/fbget/")

file_names<-dir()
file_names
clinical_fbget_long_list <- vector("list", length(file_names))
for (i in 1: length(file_names)){
  clinical_fbget_long_list[[i]]<-read.delim(file_names[i], header = TRUE, stringsAsFactors = FALSE)
  print(file_names[i])
}
names(clinical_fbget_long_list)<- file_names
clinical_fbget_long_list[[1]][1:2, 1:4]
length(clinical_fbget_long_list)

##Save data
saveRDS(clinical_fbget_long_list, file = "/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions/R workspaces/clinical_fbget_long_list.rds")

#######################
### Create a list of Clinical data obtained from CBioPortal
######################

setwd("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/clinical/cbioportal/")

##Create list:
file_names<-dir()
file_names
clinical_cbioportal_long_list <- vector("list", length(file_names))
for (i in 1: length(file_names)){
  clinical_cbioportal_long_list[[i]]<-read.delim(file_names[i], header = TRUE, stringsAsFactors = FALSE)
  print(file_names[i])
}
names(clinical_cbioportal_long_list)<- file_names

clinical_cbioportal_long_list[[1]][1:2, 1:5]
length(clinical_cbioportal_long_list)

lapply(clinical_cbioportal_long_list, function(x) dim(x))

##Loop to Remove duplicate tumour entries

for (i in 1:length(clinical_cbioportal_long_list)) {
  ##Get sample IDs
  patient_IDs<- clinical_cbioportal_long_list[[i]]$Patient.ID
  
  ##Find duplicate patients IDS
  indicies<- which(duplicated(patient_IDs))
  indicies
  
  #clinical_cbioportal_long_list[[i]][indicies, 4:6]
  
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
        df[1,b]<- sum(is.na(clinical_cbioportal_long_list[[i]][sample_pairs[b,a],]))
      }
      which_row<- which(max(df[1,]) == df[1,])
      which_row_to_remove<-c(which_row_to_remove, sample_pairs[which_row, a])
    }
    # sample_pairs
    # which_row_to_remove
    which_row_to_remove<- as.vector(which_row_to_remove)
    
    clinical_cbioportal_long_list[[i]]<- clinical_cbioportal_long_list[[i]][-which_row_to_remove,]
  } else {}
  
  print(names(clinical_cbioportal_long_list[i]))
}

setwd("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Code/co-deletions")
saveRDS(clinical_cbioportal_long_list, file = "./R workspaces/clinical_cbioportal_long_list")


sapply(clinical_cbioportal_long_list, function(x) dim(x))
sapply(clinical_fbget_long_list, function(x) dim(x))
sapply(threshold_CNV_all_table_loc ,function(x) ncol(x)-11)

###############
###Test I have removed duplicate samples from clinical_cbioportal_long_list:

##Get sample IDs
for (i in 1:length(clinical_cbioportal_long_list)){
  patient_IDs<- clinical_cbioportal_long_list[[i]]$Patient.ID
  
  indicies<- which(duplicated(patient_IDs))
  print(names(clinical_cbioportal_long_list[i]))
  print(indicies)
}

#############
###Test I have removed duplicate samples from clinical_fbget_long_list:

##Get sample IDs
for (i in 1:length(clinical_fbget_long_list)){
  patient_IDs<- clinical_fbget_long_list[[i]]$Patient.ID
  
  indicies<- which(duplicated(patient_IDs))
  print(names(clinical_fbget_long_list[i]))
  print(indicies)
}

#################
###Compare number of 
sapply(clinical_cbioportal_long_list, function(x) nrow(x))
sapply(clinical_fbget_long_list, function(x) nrow(x))
sapply(threshold_short_cnv_list_loc, function(x) ncol(x)-11)


##################
###Create clinical_cbioportal_long_list and clinical_fbget_long_list containing identical cancer types
##################

length(clinical_cbioportal_long_list)
length(clinical_fbget_long_list)
names(clinical_cbioportal_long_list)
names(clinical_fbget_long_list)

clinical_fbget_selected_long_list<-clinical_fbget_long_list[c(1,2,3,4,5,7,8,9,10,12,13,15,16,17,18,19,20,
                                                              21,22,23,24,25,26,28,29,30,32,33,34,35,36,37)]
names(clinical_cbioportal_long_list)
names(clinical_fbget_selected_long_list)
length(clinical_cbioportal_long_list)
length(clinical_fbget_selected_long_list)

#################
###Create list of all survival data:
##################

##empty list

clinical_long_survival<- vector("list", length(clinical_fbget_selected_long_list))
clinical_long_survival

##Loop:
for ( j in 1: length(clinical_fbget_selected_long_list)){
  
  ##Get fbget clinical data
  x<- clinical_fbget_selected_long_list[[j]]
  
  clinical_IDS<- x$tcga_participant_barcode
  
  ## Get days to death data:
  death_collapsed<- as.numeric(x$days_to_death)
  
  ## Get days to follow up data
  followup_collapsed<- as.numeric(x$days_to_last_followup)
  
  ## combine data:
  all_clin <- data.frame(death_collapsed,followup_collapsed)
  colnames(all_clin) <- c('death_days', 'followUp_days')
  
  ## create new column with last follow up or days to death data:
  all_clin$new_death <- c()
  for (i in 1:nrow(all_clin)){
    all_clin$new_death[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                      as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
  }
  
  ## Create new column with 1= patient died and 0 - patient survived
  all_clin$death_event <- ifelse(x$vital_status == 'alive', 0,1)
  
  ##Add row names to all_clin and a new column with clinical IDs
  rownames(all_clin) <- clinical_IDS
  all_clin<- cbind(patient_IDs = rownames(all_clin), all_clin)
  all_clin$patient_IDs<- as.character(all_clin$patient_IDs)
  
  ###############
  ##Get disease free survival days from cbioportal clinical data:
  
  x<- clinical_cbioportal_long_list[[j]]
  
  ##Get patient IDs:
  clinical_IDS<- x$Patient.ID
  
  # get the number of month to new tumour
  new_tum_collapsed<- as.numeric(x$Disease.Free..Months.)
  
  #Convert number of month to new tumour to days:
  new_tum_collapsed<- new_tum_collapsed*(365.25/12)
  
  new_tum_collapsed<- cbind(patient_IDs = clinical_IDS, new_tumor_days = new_tum_collapsed)
  new_tum_collapsed<- as.data.frame(new_tum_collapsed)
  new_tum_collapsed$patient_IDs<- as.character(new_tum_collapsed$patient_IDs)
  
  ##Join days to new tumour and overall survival together:
  all_clin_joined<- dplyr::full_join(new_tum_collapsed, all_clin, by = "patient_IDs")
  
  ##Calculate disease free survival time:
  all_clin_joined$disease_free_survival<- c()
  for (i in 1:nrow(all_clin_joined)){
    ## Combine new_tumor days and followUp_days
    all_clin_joined$disease_free_survival[i] <- ifelse (is.na(as.numeric(as.character(all_clin_joined$new_tumor_days)))[i],
                                                        as.numeric(as.character(all_clin_joined$followUp_days))[i],as.numeric(as.character(all_clin_joined$new_tumor_days))[i])
    ## If no followUp_days or new_tumor_days then take death_days
    all_clin_joined$disease_free_survival[i]<- ifelse (is.na(as.numeric(as.character(all_clin_joined$disease_free_survival)))[i],
                                                       as.numeric(as.character(all_clin_joined$death_days))[i], as.numeric(as.character(all_clin_joined$disease_free_survival))[i])
    
  }
  
  clinical_long_survival[[j]]<- all_clin_joined
  
}

names(clinical_long_survival)<- names(clinical_fbget_selected_long_list)
names(clinical_long_survival)
length(clinical_long_survival)
head(clinical_long_survival[[1]], 40)
head(clinical_long_survival[[2]], 40)
head(clinical_long_survival[[32]], 40)

##save data:
clinical_survival_long_list<- clinical_long_survival
saveRDS(clinical_survival_long_list, file = "./R workspaces/clinical_survival_long_list")




