################
### Circle plots of deletions and amplifications per cytoband
###############

library(ggplot2)


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
################
###Trial plotting dummy data

ggplot(df, aes(x,y)) +
  geom_point(aes(size = proportion_deletions)) +
  xlab("Cytoband") +
  ylab("Cancer") +
  scale_size_area("proportion of\n tumours with\n deletions", breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))


############













