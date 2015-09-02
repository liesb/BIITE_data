library(BIITE)
library(ggplot2)
library(reshape2)


beta_cutoff <- get_KL_vs_unif(rbeta(100000,1,2))

# 1. Get the results from the models from Pseudomonas
pseu_res <- read.table("../../Pseudomonas/unif_prior/results_pseudomonas_unif.txt",
                       header=T, sep="\t")
# 2. Get the predictions from NetMHCPan for pseu
predictions <- read.table("predictions_pseudomonas.txt", 
                          header=T, sep="\t")

pseu_pred <- melt(predictions)
colnames(pseu_pred) <- c("Pep", "Molec", "Value")
pseu_pred$tg <- NA
pseu_pred$ba <- NA
pseu_pred$DKL <- NA
for ( i in 1:dim(pseu_pred)[1] ){
  pseu_pred$tg[i] <- pseu_res[pseu_res$Pep==pseu_pred$Pep[i] & 
                                pseu_res$Molec==pseu_pred$Molec[i] &
                                pseu_res$What=="Mode",]$tg
  
  pseu_pred$ba[i] <- pseu_res[pseu_res$Pep==pseu_pred$Pep[i] & 
                                pseu_res$Molec==pseu_pred$Molec[i] &
                                pseu_res$What=="Mode",]$ba
  pseu_pred$DKL[i] <- pseu_res[pseu_res$Pep==pseu_pred$Pep[i] & 
                                 pseu_res$Molec==pseu_pred$Molec[i] &
                                 pseu_res$What=="Mode",]$DKL
}


auc_pred_cutoff.df <- data.frame(matrix(nrow=1000, ncol=4))
colnames(auc_pred_cutoff.df) <- c("What", "Against", "AUC", "Comment")

auc_pred_cutoff_line <- 1

auc_pred_cutoff_list <- list()

res_short <- pseu_pred[pseu_pred$DKL>beta_cutoff,]
df <- res_short[order(res_short$Value, decreasing=T),]
## Set decreasing to T because binding affinities are 'inverse measured' or so
for ( against in c("tg") ){
  tpr_stuff <- get_tpr_fpr(df, against, dont_order=T)
  
  auc_pred_cutoff.df$AUC[auc_pred_cutoff_line] <- get_auc(tpr_stuff)
  auc_pred_cutoff.df$What[auc_pred_cutoff_line] <- "pred"
  auc_pred_cutoff.df$Against[auc_pred_cutoff_line] <- against
  auc_pred_cutoff.df$Comment[auc_pred_cutoff_line] <- "cutoff"
  
  tpr_stuff$Against <- against
  tpr_stuff$What <- "pred"
  
  auc_pred_cutoff_list[[auc_pred_cutoff_line]] <- tpr_stuff
  
  auc_pred_cutoff_line <- auc_pred_cutoff_line + 1  
}

auc_pred_cutoff.df <- auc_pred_cutoff.df[!is.na(auc_pred_cutoff.df$What),]
# > auc_pred_cutoff.df
# What Against       AUC Comment
# 1 pred      tg 0.8134021  cutoff
