# CS 598 JP 
# ER project multiclass 

# RESULT OF EVALUATION OF THE LINEAR MODELS: Sim_Ann_weighted_148_restart_all_models_eval
aa_all_loss <- rowSums(do.call(rbind, lapply(Sim_Ann_weighted_148_restart_all_models_eval, "[[", 3)))
aa_GR_loss_sort_ind <- sort(aa_all_loss, decreasing = F, index.return=T)$ix
aa_GR_enh <- do.call(rbind, lapply(Sim_Ann_weighted_148_restart_all_models_eval, "[[", 2))
aaa <- unlist(lapply(ER_52_opt_simAnn_2_input$Affinity_scores, nrow))
par(mfrow=c(1,1))
aa <- Enhancer_vote_barplot(aa_GR_enh[aa_GR_loss_sort_ind[1:10], ], aaa/15, .ylim = c(-10, 10))
sum(aa[1,] >= 3)

# index of the genes to continue with
CS598_gene_index <- which(aa[1,] >= 3)
CS598_gene_names <- names(ER_52_opt_simAnn_2_input$Affinity_scores)[CS598_gene_index]

#index of enhacners per gene
CS598_enhancer_index <- integer(length(CS598_gene_index))
names(CS598_enhancer_index) <- CS598_gene_names
for(i in 1:length(CS598_gene_index)){
  aa <- table(aa_GR_enh[aa_GR_loss_sort_ind[1:10], CS598_gene_index[i]])
  CS598_enhancer_index[i] <- as.integer(names(aa)[which.max(aa)])
}
table(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52[union(CS598_gene_index, 
                                                              my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52_gene_gte3common_index), ])
##################################################################################################################################
##################################################################################################################################
# affinity of each TF for this enhancer
CS_598_enhancer_affinity <- matrix(nrow = length(CS598_enhancer_index), 
                                   ncol = ncol(ER_52_opt_simAnn_2_input$Affinity_scores[[1]]))
rownames(CS_598_enhancer_affinity) <- character(nrow(CS_598_enhancer_affinity))
colnames(CS_598_enhancer_affinity) <- colnames(ER_52_opt_simAnn_2_input$Affinity_scores[[1]])
for(i in 1:nrow(CS_598_enhancer_affinity)){
  rownames(CS_598_enhancer_affinity)[i] <- paste(CS598_gene_names[i], CS598_enhancer_index[i], sep = "_")
  CS_598_enhancer_affinity[i, ] <- ER_52_opt_simAnn_2_input$Affinity_scores[[CS598_gene_index[i]]][CS598_enhancer_index[i], ]
}

aa_enhancer_score_all <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = lapply(Sim_Ann_weighted_148_restart_all_models_eval, "[[", 4),
                                          enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                                          model_index = aa_GR_loss_sort_ind[1:10],
                                          filename = "whatever.png",
                                          loss = T,
                                          real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                                          draw_plot = T)
############################################################################################################################################################
############################################################################################################################################################
# sequence of the enhancer
CS_598_enhancer_sequence <- character(length(CS598_enhancer_index))
names(CS_598_enhancer_sequence) <- rownames(CS_598_enhancer_affinity)
for(i in 1:length(CS_598_enhancer_sequence)){
  CS_598_enhancer_sequence[i] <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Seq[[CS598_gene_index[i]]][[CS598_enhancer_index[i]]]
}

nchar(CS_598_enhancer_sequence)

############################################################################################################################################################
############################################################################################################################################################
# Kmer analysis
# list all possible gapped 6 mers.
# count the number of each kmer in the sequence
# feed the counts as features
library(tcR)
CS_598_all_5mers <- generate.kmers(.k = 5)
CS_598_enhancer_5mers <- list()
for(i in 1:length(CS_598_enhancer_sequence)){
  CS_598_enhancer_5mers[[i]] <- get.kmers(CS_598_enhancer_sequence[i], .k = 5)
}
names(CS_598_enhancer_5mers) <- names(CS_598_enhancer_sequence)

CS_598_enhancer_5mers_full <- list()
for(i in 1:length(CS_598_enhancer_sequence)){
  CS_598_enhancer_5mers_full[[i]] <- cbind(CS_598_all_5mers, integer(length(CS_598_all_5mers)))
  CS_598_enhancer_5mers_full[[i]][, 2] <- as.integer(CS_598_enhancer_5mers[[i]]$Count[match(CS_598_enhancer_5mers_full[[i]][, 1], 
                                                                                 CS_598_enhancer_5mers[[i]]$Kmers)])
  CS_598_enhancer_5mers_full[[i]][is.na(CS_598_enhancer_5mers_full[[i]])] <- 0
  
}
names(CS_598_enhancer_5mers_full) <- names(CS_598_enhancer_sequence)
# create a matrix where each row corresponds to one enhancer, each column corresponds to one kmer. 
# each entry indicates the count of the kmer in that sequence
CS_598_enhancer_5mers_full_mat <- t(do.call(cbind, CS_598_enhancer_5mers_full))
CS_598_enhancer_5mers_full_mat <- CS_598_enhancer_5mers_full_mat[seq(2, nrow(CS_598_enhancer_5mers_full_mat), 2), ]
CS_598_enhancer_5mers_full_mat <- apply(CS_598_enhancer_5mers_full_mat, 2, as.integer)
rownames(CS_598_enhancer_5mers_full_mat) <- names(CS_598_enhancer_sequence)
colnames(CS_598_enhancer_5mers_full_mat) <- CS_598_enhancer_5mers_full[[1]][, 1]

hist(colSums(CS_598_enhancer_5mers_full_mat), breaks = 100)

sum(colSums(CS_598_enhancer_5mers_full_mat) == 0)
############################################################################################################################################################
############################################################################################################################################################
# set up svm model
# use this:
CS_598_dataset <- createOptimDataset(GeneExpMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52[CS598_gene_index, ],
                                     TFExpMat=TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01,
                                     ChoppedEnhancerScores = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score[CS598_gene_index])

aa_gene_names_dataset <- unlist(lapply(strsplit(rownames(CS_598_dataset$gene_conditon_mat),
                                                split = "\\."), "[[", 1))
aa_gene_names_enhancers <- unlist(lapply(strsplit(rownames(CS_598_enhancer_5mers_full_mat),
                                                  split = "_"), "[[", 1))
aa_match <- match(aa_gene_names_dataset, aa_gene_names_enhancers)

CS_598_enh_feature <- cbind(CS_598_enhancer_affinity[aa_match, ],
                            CS_598_enhancer_5mers_full_mat[aa_match, ])
colnames(CS_598_enh_feature) <- c(paste(colnames(CS_598_enhancer_affinity), "Affinity", sep = "_"),
                                  colnames(CS_598_enhancer_5mers_full_mat))
CS_598_all_feature <- cbind(CS_598_dataset$gene_conditon_mat[, 1:(ncol(CS_598_dataset$gene_conditon_mat)-1)],
                            CS_598_enh_feature)

CS_598_all_label <- CS_598_dataset$gene_conditon_mat[, ncol(CS_598_dataset$gene_conditon_mat)]

install.packages("e1071")
library(e1071)
aa_test_ind <- sample(x = c(1:nrow(CS_598_all_feature)),
                      size = floor(nrow(CS_598_all_feature)/ 4),
                      replace = F) 
aa_test_set <- CS_598_all_feature

CS_598_all_feature_2 <- CS_598_all_feature[, -which(colnames(CS_598_all_feature) %in% c('ESR1_C', 'ESR1_T'))]
aa_svm_model2 <- svm(x = CS_598_all_feature_2, y = factor(as.character(CS_598_all_label)), cross = 5, type = "C-classification")
aa_svm_model_polynomialkernel <- svm(x = CS_598_all_feature_2, y = factor(as.character(CS_598_all_label)), cross = 5, type = "C-classification", kernel = "polynomial")
aa_svm_model_radial <- svm(x = CS_598_all_feature_2, y = factor(as.character(CS_598_all_label)), cross = 5, type = "C-classification", kernel = "radial")
aa_svm_model_sigmoid<- svm(x = CS_598_all_feature_2, y = factor(as.character(CS_598_all_label)), cross = 5, type = "C-classification", kernel = "sigmoid")
aa_svm_model_linear <- svm(x = CS_598_all_feature_2, y = factor(as.character(CS_598_all_label)), cross = 5, type = "C-classification", kernel = "linear")
aa_svm_model_linear2 <- svm(x = CS_598_all_feature_2, y = factor(as.character(CS_598_all_label)), cross = 5, type = "C-classification", kernel = "linear", cost = 1000)
aa_svm_model_linear2 <- svm(x = CS_598_all_feature_2, y = factor(as.character(CS_598_all_label)), cross = 5, type = "C-classification", kernel = "linear",  cost = 100)
aa_cls_wght <- c(length(CS_598_all_label)/ sum(CS_598_all_label == 1),
                 length(CS_598_all_label)/ sum(CS_598_all_label == 0), 
                 length(CS_598_all_label)/ sum(CS_598_all_label == -1))
names(aa_cls_wght) <- c("1", "0", "-1")
aa_svm_model_linear_wt <- svm(x = CS_598_all_feature_2,
                              y = factor(as.character(CS_598_all_label)),
                              cross = 5, type = "C-classification",
                              kernel = "linear", class.weights = aa_cls_wght)

summary(aa_svm_model_linear_wt)
aa_svm_model2$accuracies
aa_svm_model2$index
nrow(aa_svm_model2$coefs)
aa_tuned_svm <- tune.svm(x = CS_598_all_feature_2, y = factor(as.character(CS_598_all_label)))
aa_tuned_svm$best.parameters
aa_svm_model_2
aa_svm_model_polynomialkernel$accuracies
aa_svm_model_radial$accuracies
aa_svm_model_sigmoid$accuracies

aa_svm_tune <- tune(svm, train.x = CS_598_all_feature_2, train.y = factor(as.character(CS_598_all_label)), kernel = "linear",
                 ranges = list(epsilon = seq(0,1,0.1), cost = 2^(0:5)))

summary(aa_svm_model_linear)
############################################################################################################################################################
############################################################################################################################################################
# random forest model
library(randomForest)
aa_test_ind <- sample(x = c(1:nrow(CS_598_all_feature_2)),
                      size = floor(nrow(CS_598_all_feature_2)/ 4),
                      replace = F) 
aa_test_set <- CS_598_all_feature_2[aa_test_ind,]
aa_test_lab <-CS_598_all_label[aa_test_ind]

aa_train_set <- CS_598_all_feature_2[-aa_test_ind,]
aa_train_lab <-CS_598_all_label[-aa_test_ind]


aa_randomForest_wt2 <- randomForest(x = aa_train_set,
                                y = factor(as.character(aa_train_lab)),
                                xtest = aa_test_set, ytest =factor(as.character(aa_test_lab)) ,
                                importance = T, do.trace = T, 
                                classwt = c("1"= length(aa_train_lab/sum(aa_train_lab == 1)),
                                            "0"=length(aa_train_lab)/sum(aa_train_lab == 0), 
                                            "-1"=length(aa_train_lab)/sum(aa_train_lab == -1)) )

rownames(aa_randomForest_wt$importance)[sort(aa_randomForest_wt$importance[,4], decreasing = T, index.return=T)$ix[1:20]]
aa_randomForest$test$confusion
summary(aa_randomForest)
aa_improtance <- importance(aa_randomForest)

aa_randomForest_wt$confusion
aa_randomForest_wt$test$confusion
aa_randomForest$confusion
aa_randomForest$test$confusion

aa_randomForest_ss <- randomForest(x = aa_train_set,
                                   y = factor(as.character(aa_train_lab)),
                                   xtest = aa_test_set, ytest =factor(as.character(aa_test_lab)) ,
                                   importance = T, do.trace = T, 
                                   sampsize = c(174, 174, 174) )
aa_randomForest_ss$confusion
aa_randomForest_wt$confusion
aa_randomForest_ss$test$confusion
aa_randomForest_wt$test$confusion
aa_randomForest_ss$test$err.rate
# 6 Fold cross validation
aa_shuffle <- sample(c(1:nrow(CS_598_all_feature_2)), size = nrow(CS_598_all_feature_2), replace = F)
aa_nfolds <- 6
aa_z <- floor(length(aa_shuffle)/aa_nfolds)
aa_rf_Res <- list()
for(i in 1:6){
  aa_test_ind  <- aa_shuffle[((i-1)*aa_z + 1):(i*aa_z)]
  aa_train_ind <- setdiff(c(1:length(aa_shuffle)), aa_test_ind)
  aa_train_added <- add_train_example( CS_598_all_feature_2[aa_train_ind,], 
                                       CS_598_all_label[aa_train_ind])
  aa_rf_Res[[i]]  <- randomForest(x = aa_train_added$new_train_set,
                                  y = factor(as.character(aa_train_added$new_train_lab)),
                                  xtest = CS_598_all_feature_2[aa_test_ind,], 
                                  ytest =factor(as.character(CS_598_all_label[aa_test_ind])) ,
                                  importance = T, do.trace = T, ntree = 100)
  
}
plot(aa_rf_Res[[1]]$test$err.rate[,1], type= "l", col = col_vector[1], ylim = c(0.2, 0.5))
for(i in 2:6){
  lines(aa_rf_Res[[i]]$test$err.rate[,1], type= "l", col = col_vector[i+10])
}
aa_rf_Res[[6]]

plot(aa_rf_Res[[1]], main = "OOB error", lwd = 3)
legend("topright", legend = c("total", "-1", "1", "0"), col = c(1, 4,2, 3), fill=c(1, 4,2, 3))

aa_rf_Res[[1]]$test$confusion
plot(aa_rf_Res[[2]])
varImpPlot(aa_rf_Res[[1]], main = "Feature Importance")
varImpPlot(aa_rf_Res[[2]], main = "Feature Importance")


aa_mean <- integer(length(aa_rf_Res))
for(i in 1:length(aa_rf_Res)){
  aa_mean[i] <- min(aa_rf_Res[[i]]$test$err.rate[, 1])
}
mean(aa_mean)
sd(aa_mean)

plot(aa_rf_Res[[1]]$test$err.rate[, 1])

aa_imp <- integer(0)
for(i in 1:length(aa_rf_Res)){
  aa_srt <- sort(aa_rf_Res[[1]]$importance[, 5], 
                 decreasing = T, index.return = T)
  aa_imp <- c(aa_imp, aa_srt$ix[aa_srt$x > 2])
}
aa_imp <- unique(aa_imp)
colnames(CS_598_all_feature_2)[aa_imp]

CS_598_all_feature_repeat <- CS_598_all_feature_2
CS_598_all_label_repeat <- CS_598_all_label


add_train_example <- function(train_set, train_lab){
  my_table <- table(train_lab)
  lab_max<- names(my_table)[which.max(my_table)]
  other_lab = setdiff(names(my_table), lab_max)
  my_added_index <- list()
  for(i in 1:length(other_lab)){
    my_added_index[[i]] <- sample(x = which(train_lab==other_lab[i]), 
                             size = (my_table[which.max(my_table)]-my_table[other_lab[i]]), 
                             replace = T)
    
  }
  new_train_ind <- c(c(1:length(train_lab)), unlist(my_added_index))
  return(list(new_train_set = train_set[new_train_ind,],
              new_train_lab = train_lab[new_train_ind]))
}

aa_Added <- add_train_example(CS_598_all_feature_2, CS_598_all_label)
aa_randomForest_added <- randomForest(x = aa_Added$new_train_set,
                                   y = factor(as.character(aa_Added$new_train_lab)),
                                 #  xtest = aa_test_set, ytest =factor(as.character(aa_test_lab)) ,
                                   importance = T, do.trace = T )

## cross validation
# keeping enhancers out
aa_ind_pergene <- list()
aa_name_split <- unlist(lapply(strsplit(names(CS_598_all_label), split = "\\."), "[[", 1))
for(i in 1:length(CS598_gene_names)){
  aa_ind_pergene[[i]] <- which(aa_name_split %in% CS598_gene_names[i])
}
names(aa_ind_pergene) <- CS598_gene_names
aa_nfolds <- 6
aa_z <- floor(length(aa_ind_pergene)/aa_nfolds)
aa_rf_Res <- list()
aa_shuffle <- sample(c(1:length(aa_ind_pergene)), size = length(aa_ind_pergene), replace = F)

for(i in 1:6){
  aa_test_ind  <- unlist(aa_ind_pergene[aa_shuffle[((i-1)*aa_z + 1):(i*aa_z)]])
  aa_train_ind <- setdiff(c(1:length(CS_598_all_label)), aa_test_ind)
  aa_train_added <- add_train_example( CS_598_all_feature_2[aa_train_ind,], 
                                       CS_598_all_label[aa_train_ind])
  aa_rf_Res[[i]]  <- randomForest(x = aa_train_added$new_train_set,
                                  y = factor(as.character(aa_train_added$new_train_lab)),
                                  xtest = CS_598_all_feature_2[aa_test_ind,], 
                                  ytest =factor(as.character(CS_598_all_label[aa_test_ind])) ,
                                  importance = T, do.trace = T, ntree = 500)
  
}

aa_mean_old <- integer(length(aa_rf_Res))
for(i in 1:length(aa_rf_Res)){
  aa_mean_old[i] <- min(aa_rf_Res[[i]]$test$err.rate[, 1])
}
mean(aa_mean_old)
sd(aa_mean_old)

# building new feature set
CS_598_all_feature_2
aa_tf_dif <- CS_598_all_feature[, 1:19] - CS_598_all_feature[, 20:38]
aa_tf_dif_aff <- aa_tf_dif * CS_598_all_feature[, 39:57]

aa_rf_Res_nf <- list()
for(i in 1:6){
  aa_test_ind  <- unlist(aa_ind_pergene[aa_shuffle[((i-1)*aa_z + 1):(i*aa_z)]])
  aa_train_ind <- setdiff(c(1:length(CS_598_all_label)), aa_test_ind)
  aa_train_added <- add_train_example( aa_tf_dif_aff[aa_train_ind,], 
                                       CS_598_all_label[aa_train_ind])
  aa_rf_Res_nf[[i]]  <- randomForest(x = aa_train_added$new_train_set,
                                  y = factor(as.character(aa_train_added$new_train_lab)),
                                  xtest = aa_tf_dif_aff[aa_test_ind,], 
                                  ytest =factor(as.character(CS_598_all_label[aa_test_ind])) ,
                                  importance = T, do.trace = T, ntree = 500)
  
}

aa_mean <- integer(length(aa_rf_Res_nf))
for(i in 1:length(aa_rf_Res_nf)){
  aa_mean[i] <- min(aa_rf_Res_nf[[i]]$test$err.rate[, 1])
}
mean(aa_mean)
sd(aa_mean)
varImpPlot(aa_rf_Res_nf[[1]])


aa_tf_dif_aff_4kmer <- cbind(aa_tf_dif_aff, CS_598_all_feature[,which(colnames(CS_598_all_feature) %in% c("TAACA", "CAAGA", "TCCCG", "TTGGG"))])

aa_rf_Res_nf_4k <- list()
for(i in 1:6){
  aa_test_ind  <- unlist(aa_ind_pergene[aa_shuffle[((i-1)*aa_z + 1):(i*aa_z)]])
  aa_train_ind <- setdiff(c(1:length(CS_598_all_label)), aa_test_ind)
  aa_train_added <- add_train_example( aa_tf_dif_aff_4kmer[aa_train_ind,], 
                                       CS_598_all_label[aa_train_ind])
  aa_rf_Res_nf_4k[[i]]  <- randomForest(x = aa_train_added$new_train_set,
                                     y = factor(as.character(aa_train_added$new_train_lab)),
                                     xtest = aa_tf_dif_aff_4kmer[aa_test_ind,], 
                                     ytest =factor(as.character(CS_598_all_label[aa_test_ind])) ,
                                     importance = T, do.trace = T, ntree = 500)
  
}


aa_mean_2 <- integer(length(aa_rf_Res_nf))
for(i in 1:length(aa_rf_Res_nf_4k)){
  aa_mean_2[i] <- min(aa_rf_Res_nf_4k[[i]]$test$err.rate[, 1])
}
mean(aa_mean_2)
sd(aa_mean_2)
varImpPlot(aa_rf_Res_nf_4k[[1]])
aa_rf_Res_nf_4k[[6]]


varImpPlot(aa_rf_Res_nf_4k[[1]])



