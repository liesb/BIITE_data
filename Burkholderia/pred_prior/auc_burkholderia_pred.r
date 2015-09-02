library(BIITE)
library(entropy)
library(ggplot2)
library(reshape2)


res <- read.table("results_Burkholderia_pred.txt", header=T, sep="\t")


auc3.df <- data.frame(matrix(nrow=1000, ncol=4))
colnames(auc3.df) <- c("What", "Against", "AUC", "Comment")

auc_line <- 1

what_tpr_list3 <- list()

beta_cutoff <- get_KL_vs_unif(rbeta(100000,1,2))

for ( what in c("Mode") ){
  res_short <- res[res$What==what & res$DKL>beta_cutoff,]
  df <- res_short[order(res_short$Value, decreasing=F),]
  for ( against in c("tg") ){
    tpr_stuff <- get_tpr_fpr(df, against)
    
    auc3.df$AUC[auc_line] <- get_auc(tpr_stuff)
    auc3.df$What[auc_line] <- what
    auc3.df$Against[auc_line] <- against
    auc3.df$Comment[auc_line] <- "Cutoff Beta"
    
    tpr_stuff$Against <- against
    tpr_stuff$What <- what
    
    what_tpr_list3[[auc_line]] <- tpr_stuff
    
    auc_line <- auc_line + 1  
    
  }
  rm(df)
}

auc3.df <- auc3.df[!is.na(auc3.df$What),]

# > auc3.df
# What Against       AUC     Comment
# 1 Mode      tg 0.6539773 Cutoff Beta



g <- ggplot(do.call(rbind,what_tpr_list3))
g <- g + geom_line(aes(x=FPR, y=TPR, colour=What, linetype=Against))
g <- g + geom_abline(intercept=0, slope=1)
g <- g + scale_x_continuous(limits=c(0,1))
g <- g + scale_y_continuous(limits=c(0,1))
g

