md5.cat.files<- read.table(file = "/Users/Matt/Documents/co-deletions/md5_cat_Level_4.txt")
md5.cat.files

md5.checksum<- read.table(file = "/Users/Matt/Documents/co-deletions/md5_checksum.txt")
md5.checksum

dim(md5.cat.files)
dim(md5.checksum)

class(md5.cat.files[,2])
class(md5.checksum[,2])

md5.cat.files[,3]<-  as.character(md5.cat.files[,2])
md5.checksum[,3]<-  as.character(md5.checksum[,2])
dim(md5.cat.files)
dim(md5.checksum)


ordered.md5.cat.files<- md5.cat.files[order(md5.cat.files[,3]),]
ordered.md5.checksum<- md5.checksum[order(md5.checksum[,3]),]


ordered.md5.cat.files[1,3] == ordered.md5.checksum[1,3]



identical.checksum<- rep(FALSE, nrow(md5.cat.files))
identical.checksum

for ( i in 1: nrow(md5.cat.files)) {
  
  identical.checksum[i]<- ordered.md5.cat.files[i,1] == ordered.md5.checksum[i,1]
  
}


identical.checksum

nrow(md5.checksum)


new.table<- cbind( ordered.md5.cat.files[,3], ordered.md5.checksum[,3]) 
new.table

write.csv(new.table, file = "table.csv")



library(dplyr)
library(tidyr)

md5.cat.files<- read.table(file = "/Users/Matt/Documents/co-deletions/md5_cat_Level_4.txt", stringsAsFactors = FALSE, header = FALSE)
md5.cat.files

md5.checksum<- read.table(file = "/Users/Matt/Documents/co-deletions/md5_checksum.txt", stringsAsFactors = FALSE, header = FALSE)
md5.checksum

colnames(md5.cat.files)<- c("checksum", "file_name")
colnames(md5.checksum)<- c("checksum_download", "file_name")

head()

dplyr.checksum<- dplyr::full_join(md5.cat.files, md5.checksum, by = "file_name")

head(dplyr.checksum)

help(mutate)

dplyr.checksum2<- dplyr::full_join(md5.cat.files, md5.checksum, by = "file_name") %>% 
  dplyr::mutate(identical_checksum = checksum == checksum_download)

sum(dplyr.checksum2$identical_checksum)
sum(!dplyr.checksum2$identical_checksum)

dplyr.checksum2
dim(dplyr.checksum2)

##




