##########
### Sample script for performing survival analysis

library(survival)
library(RColorBrewer)
library(reshape)

### Function that performs survival analysis:
performSurvivalAnalysis <- function(surv,dfCov,plotTitle="",ylabel="Overall survival") {
  
  if (is.null(ncol(dfCov))) {
    fittedSurv <- survfit(surv~dfCov, na.action = na.exclude)
    df.categ <- cbind(sapply(names(fittedSurv$strata), function(x) strsplit(x,"=")[[1]][2]),
                      fittedSurv$n)
    df.categ <- data.frame(df.categ)
  } else {
    fittedSurv <- survfit(surv~dfCov[,1]+dfCov[,2], na.action = na.exclude)
    df.categ <- melt(fittedSurv$strata)
    df.categ$name <- rownames(df.categ)
  }
  categNames <- apply(df.categ, 1, function(x) paste0(x[1]," (",x[2],")"))
  if (!is.null(ncol(dfCov))) {
    categNames <- sapply(categNames, function(x) gsub("dfCov\\[, 1\\]",colnames(dfCov)[1],x))
    categNames <- sapply(categNames, function(x) gsub("dfCov\\[, 2\\]",colnames(dfCov)[2],x))
  }
  
  plot(fittedSurv, main=plotTitle,
       xlab="Time (days)", ylab=ylabel, 
       col=brewer.pal(9,"Set1"), mark.time=T)
  legend("topright", legend=categNames, 
         col=brewer.pal(9,"Set1"), 
         lwd=2, cex=0.9)
  
  print("Chi-sq test:")
  if (is.null(ncol(dfCov))) {
    print(survdiff(surv~dfCov,rho = 0))
  } else {
    print(survdiff(surv~dfCov[,1]+dfCov[,2],rho = 0))
  }
  
  print("Cox PH test:")
  if (is.null(ncol(dfCov))) {
    print(summary(coxph(surv~dfCov)))
    coxfit <- coxph(surv~dfCov)
  } else {
    print(summary(coxph(surv~dfCov[,1]+dfCov[,2])))
    coxfit <- coxph(surv~dfCov[,1]+dfCov[,2])
  }
  
  text(1000,0,labels=paste0("HR=",round(exp(summary(coxfit)$coefficients[1]),2),"; p=",
                            round(summary(coxfit)$logtest[3],3)))
}


#### Load sample data:
load("../../Input data/From Maria/survival example/dat.example.RData")

# remove NA values:
dat.example <- dat.example[which(!is.na(dat.example$OS.days)),]

# construct survival object:
dat.example$OSsurvObj <- with(dat.example, Surv(OS.days, AliveDead==1))

# survival by sex:
performSurvivalAnalysis(dat.example$OSsurvObj,dat.example$Sex, 
                        plotTitle = "Survival by sex")
# survival by sex and T stage:
performSurvivalAnalysis(dat.example$OSsurvObj,dat.example[,c("Sex","Tstage.EarlyLate")], 
                        plotTitle = "Survival by sex and T stage")
