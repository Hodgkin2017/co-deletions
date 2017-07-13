################
### Fishers exact test practice
################

##Create a data table

df<- data.frame(col1=c(10, 10), col2=c(0,1000))
df

##Perform fishers exact test
fisher.test(df)
ft<- fisher.test(df)
ft

## Extract important information:
str(ft)
ft$p.value
ft$estimate
