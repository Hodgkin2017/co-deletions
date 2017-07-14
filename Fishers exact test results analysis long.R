##################
###
###################

####################
###
#####################

fishers_co_deletion_test_results_list<- readRDS(file = "./R workspaces/fishers_co_deletion_test_results_list")
fishers_co_deletion_test_results_list

length(fishers_co_deletion_test_results_list)
names(fishers_co_deletion_test_results_list)<- names(threshold_selected_cnv_list_plus_all_loc)
names(fishers_co_deletion_test_results_list)

fishers_co_deletion_test_results_list[[1]]$`co-deletion_p-value`<= 0.05

significant_hits<- lapply(fishers_co_deletion_test_results_list, function(x) x$`co-deletion_p-value`<= 0.05)

lapply(fishers_co_deletion_test_results_list, function(x) x$`co-deletion_p-value`<= 0.05)

###Loop to create new list with co-deletions with p-value <= 0.5
fishers_co_deletion_test_results_list_significant<- vector("list", length(fishers_co_deletion_test_results_list))

names(fishers_co_deletion_test_results_list)<- names(threshold_selected_cnv_list_plus_all_loc)

for (i in 1: length(fishers_co_deletion_test_results_list)) {
  
  df<- fishers_co_deletion_test_results_list[[i]]
  
  fishers_co_deletion_test_results_list_significant[[i]]<- df %>% dplyr::filter(`co-deletion_p-value`<= 0.05)
  
}

names(fishers_co_deletion_test_results_list_significant)<- names(threshold_selected_cnv_list_plus_all_loc)

sapply(fishers_co_deletion_test_results_list, function(x) nrow(x))
sapply(fishers_co_deletion_test_results_list_significant, function(x) nrow(x))

BRCA_fishers_exact_test<- fishers_co_deletion_test_results_list_significant[[3]]
BRCA_fishers_exact_test<- BRCA_fishers_exact_test[order(BRCA_fishers_exact_test$`co-deletion_p-value`),]
head(BRCA_fishers_exact_test, 40)
BRCA_fishers_exact_test[41:80,]

