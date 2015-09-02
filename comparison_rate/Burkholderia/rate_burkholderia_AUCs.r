
library(BIITE)
library(ggplot2)

out <- read.table("rate_burkholderia.csv")



beta_cutoff <- get_KL_vs_unif(rbeta(100000,1,2))


# Our AUC

res_short2 <- out[out$DKL>beta_cutoff,]
df <- res_short2[order(res_short2$Value, decreasing=F),]
against <- "tg"
tpr_stuff <- get_tpr_fpr(df, against)

get_auc(tpr_stuff) 


# Their AUC

df <- res_short2[order(res_short2$Relative_freq, decreasing=F),]
against <- "tg"
tpr_stuff <- get_tpr_fpr(df, against, dont_order = T)

get_auc(tpr_stuff) 
