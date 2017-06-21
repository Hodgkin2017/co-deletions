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

which(BRCA_cnv_names %in% BRCA_cnv_names)
which(BRCA_cnv_names %in% patients)

a<- c(1,2,3,4)
b<- c(2,3,5,6)
which(a %in% b)
