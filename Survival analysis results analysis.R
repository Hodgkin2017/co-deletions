####################
###
####################

###################
### Maria's and Christophe's gene list
#################

#########################
### BRCA survival results analysis

##Import gene survival analysis results:

BRCA_co_deletion_overall_survival_stats_long<- read.csv("../../Output/Survival analysis/tables/Maria's-Christophe gene list/BRCA_co-deletion_overall_survival_stats_long.csv", header = TRUE, stringsAsFactors = FALSE)
dim(BRCA_co_deletion_overall_survival_stats_long)

BRCA_co_deletion_overall_survival_stats_long_significant<- BRCA_co_deletion_overall_survival_stats_long %>% dplyr::filter(p.value_logrank_test <= 0.05)
dim(BRCA_co_deletion_overall_survival_stats_long_significant)


BRCA_co_deletion_overall_survival_stats_long_significant_n10<- BRCA_co_deletion_overall_survival_stats_long_significant %>% dplyr::filter(number_of_samples_cat1 >= 10) %>%
  dplyr::filter(number_of_samples_cat2 >= 10)

dim(BRCA_co_deletion_overall_survival_stats_long_significant_n10)

write.csv(BRCA_co_deletion_overall_survival_stats_long_significant_n10, file = "BRCA_co_deletion_overall_survival_stats_long_significant_n10.csv")




#######################



