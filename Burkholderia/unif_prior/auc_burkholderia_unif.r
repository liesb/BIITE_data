
library(BIITE)
library(ggplot2)
library(reshape2)

res <- read.table("results_Burkholderia_unif.txt", header=T, sep="\t")


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
#   What Against      AUC     Comment
# 1 Mode      tg 0.712766 Cutoff Beta

###################
# Only mode vs tg # # with beta cutoff
###################
res_short <- res[res$What=="Mode" & res$DKL>beta_cutoff,]
df <- res_short[order(res_short$Value, decreasing=F),]
tpr_stuff <- get_tpr_fpr(df, against="tg")

get_auc(tpr_stuff)
tpr_stuff <- rbind(tpr_stuff, c(1,1))
  
g <- ggplot(tpr_stuff)
g <- g + geom_line(aes(x=FPR, y=TPR), colour="darkorange", size=1.25)
g <- g + geom_abline(intercept=0, slope=1)
g <- g + scale_x_continuous("False Positive Rate",limits=c(0,1))
g <- g + scale_y_continuous("True Positive Rate", limits=c(0,1))
g <- g + theme(axis.text=element_text(size=25), axis.title=element_text(size=25))
g  


  



