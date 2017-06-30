####################
### Survival analaysis practice
####################

library(survival)
library(ISwR)

###################
### Chapter 14: Introductory Statistics with R (Peter Dalgaard)
####################

attach(melanom)
names(melanom)
melanom$status
melanom$days

##Create Surv object
Surv(days, status == 1)

##Perform Kaplan-Meier estimates requires  ~1 to work 

survfit(Surv(days, status == 1) ~1)

surv.all<- survfit(Surv(days, status == 1)~1)
surv.all
summary(surv.all)


##Plot survival curve
plot(surv.all)
plot(surv.all, conf.int = F)

##Plot surival curve by sex
surv_by_sex<- survfit(Surv(days, status==1)~sex)
surv_by_sex
plot(surv_by_sex)
plot(surv_by_sex, conf.int=T, col=c("black", "grey"))

surv_by_sex<- survfit(Surv(days, status==1)~sex, conf.int=0.99)
surv_by_sex
plot(surv_by_sex, conf.int=T, col=c("black", "grey"))

##log-ranks test
survdiff(Surv(days, status==1)~sex)

##Stratify by ulcer
survdiff(Surv(days, status==1)~sex + strata(ulc))

##Cox proportional hazards model
summary(coxph(Surv(days, status==1)~sex))

##Cox proportional hazards model stratified by ulcer
summary(coxph(Surv(days, status==1)~sex+log(thick)+strata(ulc)))

##Plot stratified cox model
#Taking sex and thickness of tumour into account the two lines represent 
#ulcerated and non-ulcerated tumours
plot(survfit(coxph(Surv(days, status ==1)~log(thick)+sex+strata(ulc))))
plot(survfit(coxph(Surv(days, status ==1)~log(thick)+sex)))
plot(survfit(coxph(Surv(days, status ==1)~sex)))
plot(survfit(coxph(Surv(days, status ==1)~sex+strata(ulc))))



# KM Survival Analysis cannot use multiple predictors, whereas Cox Regression can.
# KM Survival Analysis can run only on a single binary predictor, whereas Cox Regression can use both continuous and binary predictors.
# KM is a non-parametric procedure, whereas Cox Regression is a semi-parametric procedure.
# Cox model gives you HR showing % increse or decrease of hazard
# Cox proportional hazards makes fewer assumptions than other forms of survival analysis; that's one reason it is dominant in the literature

###############
### Load BRCA data and perform survival analysis
###############

#############
###From BRCA_clinical.tsv:

BRCA_clinical<- read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/clinical/BRCA_clinical.tsv", header = TRUE, stringsAsFactors = FALSE)
BRCA_clinical$tcga_participant_barcode
colnames(BRCA_clinical)
BRCA_clinical[1,1]

#############
###From : ....merged_only_clinical_clin_format.txt

# clinical <- t(read.table("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/Survival analysis and gene expression tutorial/BRCA/gdac.broadinstitute.org_BRCA.Merge_Clinical.Level_1.2016012800.0.0/BRCA.merged_only_clinical_clin_format.txt",header=T, row.names=1, sep='\t'))
# clinical
clinical <- t(read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/Survival analysis and gene expression tutorial/BRCA/gdac.broadinstitute.org_BRCA.Merge_Clinical.Level_1.2016012800.0.0/BRCA.merged_only_clinical_clin_format.txt",header=T, row.names=1))
# clinical <- read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/Survival analysis and gene expression tutorial/BRCA/gdac.broadinstitute.org_BRCA.Merge_Clinical.Level_1.2016012800.0.0/BRCA.merged_only_clinical_clin_format.txt",header=T, row.names=1)
# t_clinical<- t(clinical)
# clinical <- t(read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/Survival analysis and gene expression tutorial/BRCA/gdac.broadinstitute.org_BRCA.Merge_Clinical.Level_1.2016012800.0.0/BRCA.merged_only_clinical_clin_format.txt",header=T))
# clinical_col_names<-clinical[1,] 
# clinical_col_names
# colnames(clinical)<- clinical_col_names

dim(clinical)
clinical[1:2, 1:10]
#clinical_names<- colnames(clinical)
colnames(clinical)
clinical[1,1]

##Get patient IDs:
ind_keep <- grep("patient.bcr_patient_barcode",colnames(clinical))
clinical_IDS<- clinical[,ind_keep]
clinical_IDS




##########
### Create data for survival analysis using: merged_only_clinical_clin_format file
##Follow instrictions from: https://www.biostars.org/p/153013/

# get the columns that contain data we can use: days to death, new tumor event, last day contact to....
ind_keep <- grep('days_to_new_tumor_event_after_initial_treatment',colnames(clinical))
ind_keep

# this is a bit tedious, since there are numerous follow ups, let's collapse them together and keep the first value (the higher one) if more than one is available
new_tum <- as.matrix(clinical[,ind_keep])
dim(new_tum)
new_tum_collapsed <- c()
for (i in 1:dim(new_tum)[1]){
  if ( sum ( is.na(new_tum[i,])) < dim(new_tum)[2]){
    m <- min(new_tum[i,],na.rm=T)
    new_tum_collapsed <- c(new_tum_collapsed,m)
  } else {
    new_tum_collapsed <- c(new_tum_collapsed,'NA')
  }
}

new_tum_collapsed

## Get tumour neoplasm status
ind_keep <- grep('person_neoplasm_cancer_status',colnames(clinical))
ind_keep
cancer_status <- as.matrix(clinical[,ind_keep])
cancer_status
cancer_status_collapsed <- c()
for (i in 1:dim(cancer_status)[1]){
  if ( sum ( is.na(cancer_status[i,])) < dim(cancer_status)[2]){
    m <- max(cancer_status[i,],na.rm=T)
    cancer_status_collapsed <- c(cancer_status_collapsed,m)
  } else {
    cancer_status_collapsed <- c(cancer_status_collapsed,'NA')
  }
}

cancer_status_collapsed

## comment: Why more than one value for death? Use min or max values?

# do the same to death
ind_keep <- grep('days_to_death',colnames(clinical))
death <- as.matrix(clinical[,ind_keep])
death_collapsed <- c()
for (i in 1:dim(death)[1]){
  if ( sum ( is.na(death[i,])) < dim(death)[2]){
    m <- max(death[i,],na.rm=T)
    death_collapsed <- c(death_collapsed,m)
  } else {
    death_collapsed <- c(death_collapsed,'NA')
  }
}

death_collapsed


# and days last follow up here we take the most recent which is the max number
ind_keep <- grep('days_to_last_followup',colnames(clinical))
fl <- as.matrix(clinical[,ind_keep])
fl_collapsed <- c()
for (i in 1:dim(fl)[1]){
  if ( sum (is.na(fl[i,])) < dim(fl)[2]){
    m <- max(fl[i,],na.rm=T)
    fl_collapsed <- c(fl_collapsed,m)
  } else {
    fl_collapsed <- c(fl_collapsed,'NA')
  }
}

fl_collapsed

# and put everything together
all_clin <- data.frame(new_tum_collapsed,cancer_status_collapsed, death_collapsed,fl_collapsed)
colnames(all_clin) <- c('new_tumor_days', "cancer_status", 'death_days', 'followUp_days')
head(all_clin, 20)
tail(all_clin, 20)

# create vector with time to new tumor containing data to censor for new_tumor
all_clin$new_time <- c()
for (i in 1:length(as.numeric(as.character(all_clin$new_tumor_days)))){
  ## Combine new_tumor days and followUp_days
  all_clin$new_time[i] <- ifelse (is.na(as.numeric(as.character(all_clin$new_tumor_days)))[i],
                                   as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$new_tumor_days))[i])
  ## If no followUp_days or new_tumor_days then take death_days
  all_clin$new_time[i]<- ifelse (is.na(as.numeric(as.character(all_clin$new_time)))[i],
                                 as.numeric(as.character(all_clin$death_days))[i], as.numeric(as.character(all_clin$new_time))[i])

  }


# create vector time to death containing values to censor for death
all_clin$new_death <- c()
for (i in 1:length(as.numeric(as.character(all_clin$death_days)))){
  all_clin$new_death[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                    as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
}

head(all_clin, 20)

# create vector for death censoring
ind_keep <- grep('patient.vital_status',colnames(clinical))
clinical[,ind_keep]
table(clinical[,ind_keep])
# alive dead
# 993   104

all_clin$death_event <- ifelse(clinical[,ind_keep] == 'alive', 0,1)

#finally add row.names to clinical
length(clinical_IDS)
nrow(all_clin)
rownames(all_clin) <- clinical_IDS

head(all_clin, 40)

# run survival analysis
# s <- survfit(Surv(as.numeric(as.character(all_clin$new_death)),all_clin$death_event)~1)
# s1 <- tryCatch(survdiff(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum]), error = function(e) return(NA))

##Overall survival:
s<- survfit(Surv(new_death, event = death_event == 1)~1, data = all_clin)
s
summary(s)
plot(s)

##Disease free:?
##Create a new column for 1 = had new tumour event?
s1<- survfit(Surv(new_time, event = death_event == 1)~1, data = all_clin)
s1
summary(s1)
plot(s1)
plot(s1, mark.time = T)



################
### Create data for survival analysis using :BRCA_clinical file

##As far as I am aware there is no tumour free time in the data

## Get tumour neoplasm status
ind_keep <- grep('person_neoplasm_cancer_status',colnames(BRCA_clinical))
ind_keep
cancer_status2 <- as.matrix(BRCA_clinical[,ind_keep])
cancer_status2
cancer_status_collapsed2 <- c()
for (i in 1:dim(cancer_status2)[1]){
  if ( sum ( is.na(cancer_status2[i,])) < dim(cancer_status2)[2]){
    m <- max(cancer_status2[i,],na.rm=T)
    cancer_status_collapsed2 <- c(cancer_status_collapsed2,m)
  } else {
    cancer_status_collapsed2 <- c(cancer_status_collapsed2,'NA')
  }
}

cancer_status_collapsed2


# do the same to death
ind_keep <- grep('days_to_death',colnames(BRCA_clinical))
ind_keep
death2 <- as.matrix(BRCA_clinical[,ind_keep])
death2
death_collapsed2 <- c()
for (i in 1:dim(death2)[1]){
  if ( sum ( is.na(death2[i,])) < dim(death2)[2]){
    m <- max(death2[i,],na.rm=T)
    death_collapsed2 <- c(death_collapsed2,m)
  } else {
    death_collapsed2 <- c(death_collapsed2,'NA')
  }
}

death_collapsed2


# and days last follow up here we take the most recent which is the max number
ind_keep <- grep('days_to_last_followup',colnames(BRCA_clinical))
fl2 <- as.matrix(BRCA_clinical[,ind_keep])
fl_collapsed2 <- c()
for (i in 1:dim(fl2)[1]){
  if ( sum (is.na(fl2[i,])) < dim(fl2)[2]){
    m <- max(fl2[i,],na.rm=T)
    fl_collapsed2 <- c(fl_collapsed2,m)
  } else {
    fl_collapsed2 <- c(fl_collapsed2,'NA')
  }
}

fl_collapsed2

# and put everything together
all_clin2 <- data.frame(death_collapsed2,fl_collapsed2, cancer_status_collapsed2)
colnames(all_clin2) <- c('death_days', 'followUp_days', "cancer_status")
head(all_clin2, 20)
tail(all_clin2, 20)


# create vector time to death containing values to censor for death
all_clin2$new_death <- c()
for (i in 1:length(as.numeric(as.character(all_clin2$death_days)))){
  all_clin2$new_death[i] <- ifelse ( is.na(as.numeric(as.character(all_clin2$death_days))[i]),
                                    as.numeric(as.character(all_clin2$followUp_days))[i],as.numeric(as.character(all_clin2$death_days))[i])
}

head(all_clin2, 20)

# create vector for death censoring
ind_keep <- grep('vital_status',colnames(BRCA_clinical))
BRCA_clinical[,ind_keep]
table(BRCA_clinical[,ind_keep])
# alive dead
# 945   152 (this analysis)
# 993   104 (previous analysis)

all_clin2$death_event <- ifelse(BRCA_clinical[,ind_keep] == 'alive', 0,1)

##finally add row.names to clinical
clinical_IDS2<- BRCA_clinical$tcga_participant_barcode
clinical_IDS2
length(clinical_IDS2)
nrow(all_clin2)
rownames(all_clin2) <- clinical_IDS2

head(all_clin2, 20)

# run survival analysis
# s <- survfit(Surv(as.numeric(as.character(all_clin$new_death)),all_clin$death_event)~1)
# s1 <- tryCatch(survdiff(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum]), error = function(e) return(NA))

##Overall survival:
s<- survfit(Surv(new_death, death_event == 1)~1, data = all_clin2)
s
summary(s)
plot(s)

##Disease free:?
# s1<- survfit(Surv(new_time, death_event == 1)~1, data = all_clin2)
# s1
# summary(s1)
# plot(s1)

##Comment: On first glance the survival analysis using the two datasets is quite different. Need to
#Find out why and how different they are.


##############
### Compare patient IDs between two datasets:

clinical_ID_dataset_fget<- sort(clinical_IDS2, decreasing = FALSE)
clinical_ID_dataset1_firehose<- as.vector(toupper(clinical_IDS))
clinical_ID_dataset1_firehose<- sort(clinical_ID_dataset1_firehose, decreasing = FALSE)

clinical_ID_dataset_fget[1:10]
clinical_ID_dataset1_firehose[1:10]
identical(clinical_ID_dataset_fget, clinical_ID_dataset1_firehose)

##Comment: Samples are the same between datasets


############
### Compare vital status between datasets

all_clin$ID<- as.vector(toupper(clinical_IDS))
head(all_clin, 20)

all_clin2$ID<- clinical_IDS2
head(all_clin2, 20)

##Order dataframes based on patient barcode:
all_clin<- all_clin[order(all_clin$ID),]
all_clin2<- all_clin2[order(all_clin2$ID),]

identical(all_clin$death_days, all_clin2$death_days)
identical(all_clin$death_event, all_clin2$death_event)

indicies<- which(!(all_clin$death_event == all_clin2$death_event))
indicies
# a<- c(0, 0, 1, 0, 1, 1)
# b<- c(0, 0, 1, 1, 0, 0)
# which(a %in% b)
# match(a,b)
# which(!(a==b))

all_clin_combined<- cbind(all_clin, all_clin2)
head(all_clin_combined)

all_clin_combined[indicies,]

##Comment: Some sample say the patient is alive when they have a days to death number and BRCA_clinical and
#'picked data' files say they should be dead. This suggests the death event data is incorrect for the firehose dataset.
#However, there is no info for BRCA of what they died from while for some cancers there is information. Also, for 
# example TCGA-A2-A0CO was tumour free when she died but is classed as a death event for the fget data but not for the 
#firebrowse data so perhaps this individual died from other causes?

##Comment, repeat analysis with PRAD/STAD data that has patient_death_reason column and see what happens.
#Looking at CBio survival curve (http://www.cbioportal.org/study?id=brca_tcga#summary) my survival curve using the (fget) BRCA_clinical dataset is most similar.







############
##Use disease free survival with (new_time) with death_event from fbget data table (BRCA_clinical)

all_clin2_with_DFS<- cbind(all_clin2, ID2 = all_clin$ID, new_time = all_clin$new_time)
head(all_clin2_with_DFS, 20)

identical(all_clin2_with_DFS$ID, as.character(all_clin2_with_DFS$ID2))

s<- survfit(Surv(new_death, death_event == 1)~1, data = all_clin2_with_DFS)
s
summary(s)
plot(s)
plot(s, mark.time = T)

##This individual does not have death time:
all_clin2_with_DFS[932,]

##Disease free survival
s1<- survfit(Surv(new_time, death_event == 1)~1, data = all_clin2_with_DFS)
s1
summary(s1)
plot(s1)
plot(s1, mark.time = T)

##################
### Repeat Survival analysis with STAD data to see if patient_reason_death explains difference in vitality status
#################


#############
###From STAD_clinical.tsv:

STAD_clinical<- read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/clinical/STAD_clinical.tsv", header = TRUE, stringsAsFactors = FALSE)
STAD_clinical$tcga_participant_barcode
colnames(STAD_clinical)
STAD_clinical[1,1]

#############
###From : ....merged_only_clinical_clin_format.txt


STAD_clinical_firebrowse <- t(read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/Survival analysis and gene expression tutorial/STAD/gdac.broadinstitute.org_STAD.Merge_Clinical.Level_1.2016012800.0.0/STAD.merged_only_clinical_clin_format.txt",header=T, row.names=1))

dim(STAD_clinical_firebrowse)
STAD_clinical_firebrowse[1:2, 1:10]
colnames(STAD_clinical_firebrowse)
STAD_clinical_firebrowse[1,1]

##Get patient IDs:
ind_keep <- grep("patient.bcr_patient_barcode",colnames(STAD_clinical_firebrowse))
clinical_IDS<- STAD_clinical_firebrowse[,ind_keep]
clinical_IDS




##########
### Create data for survival analysis using: merged_only_clinical_clin_format file
##Follow instrictions from: https://www.biostars.org/p/153013/

# get the columns that contain data we can use: days to death, new tumor event, last day contact to....
ind_keep <- grep('days_to_new_tumor_event_after_initial_treatment',colnames(STAD_clinical_firebrowse))
ind_keep

# this is a bit tedious, since there are numerous follow ups, let's collapse them together and keep the first value (the higher one) if more than one is available
new_tum <- as.matrix(STAD_clinical_firebrowse[,ind_keep])
dim(new_tum)
new_tum_collapsed <- c()
for (i in 1:dim(new_tum)[1]){
  if ( sum ( is.na(new_tum[i,])) < dim(new_tum)[2]){
    m <- min(new_tum[i,],na.rm=T)
    new_tum_collapsed <- c(new_tum_collapsed,m)
  } else {
    new_tum_collapsed <- c(new_tum_collapsed,'NA')
  }
}

new_tum_collapsed

## comment: Why more than one value for death? Use min or max values?

# do the same to death
ind_keep <- grep('days_to_death',colnames(STAD_clinical_firebrowse))
death <- as.matrix(STAD_clinical_firebrowse[,ind_keep])
death_collapsed <- c()
for (i in 1:dim(death)[1]){
  if ( sum ( is.na(death[i,])) < dim(death)[2]){
    m <- max(death[i,],na.rm=T)
    death_collapsed <- c(death_collapsed,m)
  } else {
    death_collapsed <- c(death_collapsed,'NA')
  }
}

death_collapsed


# and days last follow up here we take the most recent which is the max number
ind_keep <- grep('days_to_last_followup',colnames(STAD_clinical_firebrowse))
fl <- as.matrix(STAD_clinical_firebrowse[,ind_keep])
fl_collapsed <- c()
for (i in 1:dim(fl)[1]){
  if ( sum (is.na(fl[i,])) < dim(fl)[2]){
    m <- max(fl[i,],na.rm=T)
    fl_collapsed <- c(fl_collapsed,m)
  } else {
    fl_collapsed <- c(fl_collapsed,'NA')
  }
}

fl_collapsed

# and put everything together
STAD_clinical_firebrowse_all_clin <- data.frame(new_tum_collapsed,death_collapsed,fl_collapsed)
colnames(STAD_clinical_firebrowse_all_clin) <- c('new_tumor_days', 'death_days', 'followUp_days')
head(STAD_clinical_firebrowse_all_clin, 20)
tail(STAD_clinical_firebrowse_all_clin, 20)

# create vector with time to new tumor containing data to censor for new_tumor
STAD_clinical_firebrowse_all_clin$new_time <- c()
for (i in 1:length(as.numeric(as.character(STAD_clinical_firebrowse_all_clin$new_tumor_days)))){
  STAD_clinical_firebrowse_all_clin$new_time[i] <- ifelse ( is.na(as.numeric(as.character(STAD_clinical_firebrowse_all_clin$new_tumor_days))[i]),
                                   as.numeric(as.character(STAD_clinical_firebrowse_all_clin$followUp_days))[i],as.numeric(as.character(STAD_clinical_firebrowse_all_clin$new_tumor_days))[i])
}


# create vector time to death containing values to censor for death
STAD_clinical_firebrowse_all_clin$new_death <- c()
for (i in 1:length(as.numeric(as.character(STAD_clinical_firebrowse_all_clin$death_days)))){
  STAD_clinical_firebrowse_all_clin$new_death[i] <- ifelse ( is.na(as.numeric(as.character(STAD_clinical_firebrowse_all_clin$death_days))[i]),
                                    as.numeric(as.character(STAD_clinical_firebrowse_all_clin$followUp_days))[i],as.numeric(as.character(STAD_clinical_firebrowse_all_clin$death_days))[i])
}

head(STAD_clinical_firebrowse_all_clin, 20)

# create vector for death censoring
ind_keep <- grep('patient.vital_status',colnames(STAD_clinical_firebrowse))
ind_keep
STAD_clinical_firebrowse[,ind_keep]
table(STAD_clinical_firebrowse[,ind_keep])
# alive dead
# 356    87 

STAD_clinical_firebrowse_all_clin$death_event <- ifelse(STAD_clinical_firebrowse[,ind_keep] == 'alive', 0,1)

#finally add row.names to clinical
length(clinical_IDS)
nrow(STAD_clinical_firebrowse_all_clin)
rownames(STAD_clinical_firebrowse_all_clin) <- clinical_IDS

head(STAD_clinical_firebrowse_all_clin, 20)

# run survival analysis
# s <- survfit(Surv(as.numeric(as.character(all_clin$new_death)),all_clin$death_event)~1)
# s1 <- tryCatch(survdiff(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum]), error = function(e) return(NA))

##Overall survival:
s<- survfit(Surv(new_death, death_event == 1)~1, data = STAD_clinical_firebrowse_all_clin)
s
summary(s)
plot(s)

##Disease free:?
##Create a new column for 1 = had new tumour event?
s1<- survfit(Surv(new_time, death_event == 1)~1, data = STAD_clinical_firebrowse_all_clin)
s1
summary(s1)
plot(s1)





################
### Create data for survival analysis using :STAD_clinical file

##As far as I am awars there is no tumour free time in the data

# do the same to death
ind_keep <- grep('days_to_death',colnames(STAD_clinical))
ind_keep
death2 <- as.matrix(STAD_clinical[,ind_keep])
death2
death_collapsed2 <- c()
for (i in 1:dim(death2)[1]){
  if ( sum ( is.na(death2[i,])) < dim(death2)[2]){
    m <- max(death2[i,],na.rm=T)
    death_collapsed2 <- c(death_collapsed2,m)
  } else {
    death_collapsed2 <- c(death_collapsed2,'NA')
  }
}

death_collapsed2


# and days last follow up here we take the most recent which is the max number
ind_keep <- grep('days_to_last_followup',colnames(STAD_clinical))
fl2 <- as.matrix(STAD_clinical[,ind_keep])
fl_collapsed2 <- c()
for (i in 1:dim(fl2)[1]){
  if ( sum (is.na(fl2[i,])) < dim(fl2)[2]){
    m <- max(fl2[i,],na.rm=T)
    fl_collapsed2 <- c(fl_collapsed2,m)
  } else {
    fl_collapsed2 <- c(fl_collapsed2,'NA')
  }
}

fl_collapsed2

# and put everything together
STAD_clinical_all_clin2 <- data.frame(death_collapsed2,fl_collapsed2)
colnames(STAD_clinical_all_clin2) <- c('death_days', 'followUp_days')
head(STAD_clinical_all_clin2, 20)
tail(STAD_clinical_all_clin2, 20)


# create vector time to death containing values to censor for death
STAD_clinical_all_clin2$new_death <- c()
for (i in 1:length(as.numeric(as.character(STAD_clinical_all_clin2$death_days)))){
  STAD_clinical_all_clin2$new_death[i] <- ifelse ( is.na(as.numeric(as.character(STAD_clinical_all_clin2$death_days))[i]),
                                     as.numeric(as.character(STAD_clinical_all_clin2$followUp_days))[i],as.numeric(as.character(STAD_clinical_all_clin2$death_days))[i])
}

head(STAD_clinical_all_clin2, 20)

# create vector for death censoring
ind_keep <- grep('vital_status',colnames(STAD_clinical))
STAD_clinical[,ind_keep]
table(STAD_clinical[,ind_keep])
# alive dead
# 268   175 (this analysis)
# 356    87 (previous analysis)

STAD_clinical_all_clin2$death_event <- ifelse(STAD_clinical[,ind_keep] == 'alive', 0,1)

##finally add row.names to clinical
clinical_IDS2<- STAD_clinical$tcga_participant_barcode
clinical_IDS2
length(clinical_IDS2)
nrow(STAD_clinical_all_clin2)
rownames(STAD_clinical_all_clin2) <- clinical_IDS2

head(STAD_clinical_all_clin2, 20)

# run survival analysis
# s <- survfit(Surv(as.numeric(as.character(all_clin$new_death)),all_clin$death_event)~1)
# s1 <- tryCatch(survdiff(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum]), error = function(e) return(NA))

##Overall survival:
s<- survfit(Surv(new_death, death_event == 1)~1, data = STAD_clinical_all_clin2)
s
summary(s)
plot(s)

##Disease free:?
# s1<- survfit(Surv(new_time, death_event == 1)~1, data = all_clin2)
# s1
# summary(s1)
# plot(s1)

##Comment: On first glance the survival analysis using the two datasets is quite different. Need to
#Find out why and how different they are.


##############
### Compare patient IDs between two datasets:

STAD_clinical_ID_dataset_fget<- sort(clinical_IDS2, decreasing = FALSE)
STAD_clinical_ID_dataset1_firehose<- as.vector(toupper(clinical_IDS))
STAD_clinical_ID_dataset1_firehose<- sort(STAD_clinical_ID_dataset1_firehose, decreasing = FALSE)

STAD_clinical_ID_dataset_fget[1:10]
STAD_clinical_ID_dataset1_firehose[1:10]
identical(STAD_clinical_ID_dataset_fget, STAD_clinical_ID_dataset1_firehose)

##Comment: Samples are the same between datasets


############
### Compare vital status between datasets

STAD_clinical_firebrowse_all_clin$ID<- as.vector(toupper(clinical_IDS))
head(STAD_clinical_firebrowse_all_clin, 20)

STAD_clinical_all_clin2$ID<- clinical_IDS2
head(STAD_clinical_all_clin2, 20)

##Order dataframes based on patient barcode:
STAD_clinical_firebrowse_all_clin<- STAD_clinical_firebrowse_all_clin[order(STAD_clinical_firebrowse_all_clin$ID),]
STAD_clinical_all_clin2<- STAD_clinical_all_clin2[order(STAD_clinical_all_clin2$ID),]

identical(STAD_clinical_firebrowse_all_clin$death_days, STAD_clinical_all_clin2$death_days)
identical(STAD_clinical_firebrowse_all_clin$death_event, STAD_clinical_all_clin2$death_event)

indicies<- which(!(STAD_clinical_firebrowse_all_clin$death_event == STAD_clinical_all_clin2$death_event))
indicies
# a<- c(0, 0, 1, 0, 1, 1)
# b<- c(0, 0, 1, 1, 0, 0)
# which(a %in% b)
# match(a,b)
# which(!(a==b))

all_clin_combined<- cbind(STAD_clinical_firebrowse_all_clin, STAD_clinical_all_clin2)
head(all_clin_combined)

all_clin_combined[indicies,]

##Comment: Differences in datasets not related to cause of death. Therefore, maybe firebrowse dataset is older?
##From now on use survival data from fget downloaded data but get days_to_tumor_event_after_initial_treatment
#to determine disease free status?









################
### Use BRCA data from cbioportal to plot survival plots
###############

BRCA_cbioportal_clinical <- t(read.delim("/Users/Matt/Documents/Masters_Bioinformatics/Internships/Input data/Survival analysis and gene expression tutorial/BRCA/BRCA_cbioportal_clinical_data.txt",header=T, row.names=5))
dim(BRCA_cbioportal_clinical)
BRCA_cbioportal_clinical[1:2, 1:10]
colnames(BRCA_cbioportal_clinical)
rownames(BRCA_cbioportal_clinical)
BRCA_cbioportal_clinical[1,1]

###Remove duplicate samples?
##Get sample IDs
patient_IDs<- colnames(BRCA_cbioportal_clinical)

#Keep patient ID only
patient_IDs<- substr(patient_IDs, 0, 12)
ncol(BRCA_cbioportal_clinical)
length(unique(patient_IDs))

indicies<- which(duplicated(patient_IDs))
indicies

BRCA_cbioportal_clinical[1:2, indicies]

repeated_samples<- patient_IDs[indicies]
indicies<- which(patient_IDs %in% repeated_samples)
indicies
patient_IDs[indicies]
BRCA_cbioportal_clinical[1:2,indicies]

###Keep column with less NA's
sample_pairs<- c()
##find sample pairs:
for (i in 1: length(indicies)){
  sampleID<- patient_IDs[indicies[i]]
  sampleID<- grep(sampleID,patient_IDs) 
  sample_pairs<- cbind(sample_pairs, sampleID)
}
sample_pairs


##keep unique sample_pairs columns only
sample_pairs<- sample_pairs[,!duplicated(sample_pairs, MARGIN = 2)]
sample_pairs

## Identify which column has the lowest number of NAs
n<- nrow(sample_pairs)
df<- matrix(NA, nrow = 1, ncol = n)
df
which_column_to_remove<- c()

for (i in 1:ncol(sample_pairs)){

  for (j in 1:n){
  df[1,j]<- sum(is.na(BRCA_cbioportal_clinical[,sample_pairs[j,i]]))
  }
  which_row<- which(max(df[1,]) == df[1,])
  which_column_to_remove<-c(which_column_to_remove, sample_pairs[which_row, i])
}
sample_pairs
which_column_to_remove
which_column_to_remove<- as.vector(which_column_to_remove)

BRCA_cbioportal_clinical_unique<- BRCA_cbioportal_clinical[,-which_column_to_remove]

###Test I have removed duplicate samples:

##Get sample IDs
patient_IDs<- colnames(BRCA_cbioportal_clinical_unique)

##Keep patient ID only
patient_IDs<- substr(patient_IDs, 0, 12)
ncol(BRCA_cbioportal_clinical_unique)
length(unique(patient_IDs))

indicies<- which(duplicated(patient_IDs))
indicies

BRCA_cbioportal_clinical_unique[1:2, indicies]
##Comment: Worked! Only unique columns now! Turn above into function?






###################
### Create table using BRCA cbioportal dataset to plot survival data:

##Get patient IDs:
ind_keep <- grep("Patient ID",colnames(BRCA_cbioportal_clinical_unique))
ind_keep
clinical_IDS<- colnames(BRCA_cbioportal_clinical_unique)
clinical_IDS<- substr(clinical_IDS, 0, 12)
clinical_IDS

# get the columns that contain data we can use: days to death, new tumor event, last day contact to....
ind_keep <- grep("Disease.Free..Months.",rownames(BRCA_cbioportal_clinical_unique))
ind_keep
new_tum_collapsed<- as.numeric(BRCA_cbioportal_clinical_unique[ind_keep,])
new_tum_collapsed

# do the same for death
ind_keep <- grep("Death.from.Initial.Pathologic.Diagnosis.Date",rownames(BRCA_cbioportal_clinical_unique))
ind_keep
death_collapsed<- as.numeric(BRCA_cbioportal_clinical_unique[ind_keep,])
death_collapsed


# and days last follow up here we take the most recent which is the max number
ind_keep <- grep('days_to_last_followup',colnames(clinical))
fl <- as.matrix(clinical[,ind_keep])
fl_collapsed<- as.numeric(BRCA_cbioportal_clinical_unique[ind_keep,])
fl_collapsed

# and put everything together
all_clin <- data.frame(new_tum_collapsed,cancer_status_collapsed, death_collapsed,fl_collapsed)
colnames(all_clin) <- c('new_tumor_days', "cancer_status", 'death_days', 'followUp_days')
head(all_clin, 20)
tail(all_clin, 20)

# create vector with time to new tumor containing data to censor for new_tumor
all_clin$new_time <- c()
for (i in 1:length(as.numeric(as.character(all_clin$new_tumor_days)))){
  ## Combine new_tumor days and followUp_days
  all_clin$new_time[i] <- ifelse (is.na(as.numeric(as.character(all_clin$new_tumor_days)))[i],
                                  as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$new_tumor_days))[i])
  ## If no followUp_days or new_tumor_days then take death_days
  all_clin$new_time[i]<- ifelse (is.na(as.numeric(as.character(all_clin$new_time)))[i],
                                 as.numeric(as.character(all_clin$death_days))[i], as.numeric(as.character(all_clin$new_time))[i])
  
}


# create vector time to death containing values to censor for death
all_clin$new_death <- c()
for (i in 1:length(as.numeric(as.character(all_clin$death_days)))){
  all_clin$new_death[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                    as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
}

head(all_clin, 20)

# create vector for death censoring
ind_keep <- grep('patient.vital_status',colnames(clinical))
clinical[,ind_keep]
table(clinical[,ind_keep])
# alive dead
# 993   104

all_clin$death_event <- ifelse(clinical[,ind_keep] == 'alive', 0,1)

#finally add row.names to clinical
length(clinical_IDS)
nrow(all_clin)
rownames(all_clin) <- clinical_IDS

head(all_clin, 40)

# run survival analysis
# s <- survfit(Surv(as.numeric(as.character(all_clin$new_death)),all_clin$death_event)~1)
# s1 <- tryCatch(survdiff(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum]), error = function(e) return(NA))

##Overall survival:
s<- survfit(Surv(new_death, event = death_event == 1)~1, data = all_clin)
s
summary(s)
plot(s)

##Disease free:?
##Create a new column for 1 = had new tumour event?
s1<- survfit(Surv(new_time, event = death_event == 1)~1, data = all_clin)
s1
summary(s1)
plot(s1)
plot(s1, mark.time = T)




