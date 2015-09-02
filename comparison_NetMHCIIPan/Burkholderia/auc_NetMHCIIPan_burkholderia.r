###############
# 21-Nov-2014 #
###############

library(ggplot2)
library(BIITE)
library(reshape2)

beta_cutoff <- get_KL_vs_unif(rbeta(100000,1,2))

# 1. Get the results from the models from Burkholderia
burk_res <- read.table("../../Burkholderia/unif_prior/results_Burkholderia_unif.txt",
                       header=T, sep="\t")
# 2. Get the predictions from NetMHCPan for Burk
predictions <- read.table("predictions_burkholderia.csv",
                        header=T, sep="\t")
burk_pred <- melt(predictions)
colnames(burk_pred) <- c("Pep", "Molec", "Value")
burk_pred$tg <- NA
burk_pred$ba <- NA
burk_pred$DKL <- NA
for ( i in 1:dim(burk_pred)[1] ){
  burk_pred$tg[i] <- burk_res[burk_res$Pep==burk_pred$Pep[i] & 
                                burk_res$Molec==burk_pred$Molec[i] &
                                burk_res$What=="Mode",]$tg
  
  burk_pred$ba[i] <- burk_res[burk_res$Pep==burk_pred$Pep[i] & 
                                burk_res$Molec==burk_pred$Molec[i] &
                                burk_res$What=="Mode",]$ba
  burk_pred$DKL[i] <- burk_res[burk_res$Pep==burk_pred$Pep[i] & 
                                burk_res$Molec==burk_pred$Molec[i] &
                                burk_res$What=="Mode",]$DKL
}

##########################
# burk preds with cutoff #
##########################
auc_pred_cutoff.df <- data.frame(matrix(nrow=1000, ncol=4))
colnames(auc_pred_cutoff.df) <- c("What", "Against", "AUC", "Comment")

auc_pred_cutoff_line <- 1

auc_pred_cutoff_list <- list()

res_short <- burk_pred[burk_pred$DKL>beta_cutoff,]
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
# What Against       AUC Comment
# 1 pred      tg 0.6765957  cutoff

