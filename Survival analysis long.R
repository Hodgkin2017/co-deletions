#################
### Create a list of tables for each cancer containing overall survival and disease free survival
#################








###############
### Characterise Tunours based on their type of deletion and perform survival analysis:
##############

##Parameters:
target_gene<- target_genes[2]
gene_information_list

##########
###For loop to calculate p-vale and
for (i in 1:1){}
#############
## Get genes surrounding target gene

if (remove_NA == TRUE){
  
  cnv.table<- cnv.table %>%
    dplyr::filter(!is.na(start))
}

############
## Select first gene and catagorise tumours 1,2,3 or 4 depending on what type of deletion each tumour has
# Make a vector for each tumour in BRCA if it contains a deletion in the target gene only (1), 
#target gene and proximal gene (2), proximal gene only (3) or deletion in neither gene (4)

###########
## Perform survival analysis and record p-value and mean/median survival
#print(survfit(my.surv ~ 1), print.rmean=TRUE)




##########







#################
###
################





















