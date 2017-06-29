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
###From 



#############
###From : ....merged_only_clinical_clin_format.txt


