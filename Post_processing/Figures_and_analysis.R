#################################################
##
##
## We have so far worked out prevalence
## we compare like for like with prevalence

## Efficacy: here we use the control arm, and compare
##           performance over time with this to treated arms

## Now we need to compare effiacy to the baseline of the
## respective arm instead as sometimes the baseline is mismatched

##########################
##
## Figures
##
##



setwd("C:/Users/esherrar/Documents/Rprojects/ibm_rct_prediction")


effi = read.csv("Post_processing/summary_data/DATA_RESOURCE_EFFICACY_global_for_fig1_drop_targetedFINAL.csv",header=TRUE)

effi$col_efficacy = ifelse(effi$Obs_mean_eff < 0,"grey",
                           ifelse(effi$Obs_mean_eff > 0 & effi$Obs_mean_eff < 0.1, "darkred",
                                  ifelse(effi$Obs_mean_eff > 0.1 & effi$Obs_mean_eff < 0.2, "red",
                                         ifelse(effi$Obs_mean_eff > 0.2 & effi$Obs_mean_eff < 0.3, "orange",
                                                ifelse(effi$Obs_mean_eff > 0.3 & effi$Obs_mean_eff < 0.4, "yellow",
                                                       ifelse(effi$Obs_mean_eff > 0.4 & effi$Obs_mean_eff < 0.5, "green","darkgreen"))))))

effi$col_efficacy_MOD = ifelse(effi$Pred_1_eff < 0,"grey",
                               ifelse(effi$Pred_1_eff > 0 & effi$Pred_1_eff < 0.1, "darkred",
                                      ifelse(effi$Pred_1_eff > 0.1 & effi$Pred_1_eff < 0.2, "red",
                                             ifelse(effi$Pred_1_eff > 0.2 & effi$Pred_1_eff < 0.3, "orange",
                                                    ifelse(effi$Pred_1_eff > 0.3 & effi$Pred_1_eff < 0.4, "yellow",
                                                           ifelse(effi$Pred_1_eff > 0.4 & effi$Pred_1_eff < 0.5, "green","darkgreen"))))))



##############################################################
##
## Figure 1
## to show the efficacy over time across trials
## Months tested in the trials

## TSC feedback - wants to drop targeted ITNs
##              - and add in reference numbers for studies
##              - country as column


timeline = 0:max(effi$Month,na.rm=TRUE)

## Set up the image as wanted
plot.new()
par(new = "TRUE",
    plt = c(0.03,0.15,0.1,0.9),
    las = 1,
    cex.axis = 1)

# lis = 56:1
# test_col1 = c(rep("black",31),rep("blue",6),rep("darkred",19))


effi2 = effi[!duplicated(effi[,c('reps')]),]
# effi2 = effi2[c(3,5:21),]
lis = 20:1
test_col2 = c(#"gold4","gold4",
  rep("black",6),
  rep("grey25",5),
  rep("grey5",2),
  rep("darkred",7),"darkgreen","blue")

plot(1:length(unique(effi2$reps)),
     main = "Control arm",
     xaxt="n",yaxt="n",bty="n",xlab="",ylab="",pch="")
for(i in 1:nrow(effi2)){
  text(12,i,print(effi2$control_arm_intervention[lis[i]]),
       col=rev(test_col2)[i])
}

par(new = "TRUE",
    plt = c(0.16,0.30,0.1,0.9),
    las = 1,
    cex.axis = 1)

test_col3 = as.character(effi2$intn_col)
plot(1:length(effi2$reps),
     main = "Intervention",
     xaxt="n",yaxt="n",bty="n",xlab="",ylab="",pch="")
for(i in 1:nrow(effi2)){
  text(12,i,print(effi2$Intervention[lis[i]]),
       col=rev(test_col3)[i])
}


par(new = "TRUE",
    plt = c(0.32,0.48,0.1,0.9),
    las = 1,
    cex.axis = 1)
plot(1:length(effi2$reps),
     main = "Location      ",
     xaxt="n",yaxt="n",bty="n",xlab="",ylab="",pch="")
for(i in 1:nrow(effi2)){
  text(8,i,print(effi2$country[lis[i]]))
}

par(new = "TRUE",
    plt = c(0.42,0.88,0.1,0.9),
    las = 1,
    cex.axis = 1)

effi2$intervention_colour = ifelse(effi2$Intervention == "ITN","darkred",
                                   ifelse(effi2$Intervention == "ITN_IRS","aquamarine3",
                                          ifelse(effi2$Intervention == "targetedITN","orange",
                                                 ifelse(effi2$Intervention == "targetedITN_IRS","black",
                                                        ifelse(effi2$Intervention == "PBO_ITN","blue",
                                                               ifelse(effi2$Intervention == "PBO_ITN_IRS","purple","grey"))))))

plot(lis~effi2$Month[1:20],#effi2$Month[1:21],
     main = "Observed efficacy against prevalence (%)",
     xlim=c(-8,50),xlab="Months after deployment",yaxt="n",
     ylab="",bty="n",pch="",cex=1.6,col=effi$col_efficacy,xaxt="n")
axis(1,seq(0,50,10))
lis1 = rev(1:20)#+0.2
effi$reps
for(i in 1:20){#1:21){
  points(c(rep(lis1[i],length(effi$Month[effi$reps == i])))~c(effi$Month[effi$reps == i]),
         pch=15,cex=2,col=effi$col_efficacy[effi$reps == i])
}


#circle24, point-up triangle12, cross25, x26, diamond27, point-down triangle28, square with x29, asterisk2, diamond with cross30, circle with cross31, star23, square with cross13, circle with x32. 
ref_paper = c(9,#curtis
              65,#henry
              67,#marbiah
              68,#phil-how
              10,#nevill
              9,
              rep(64,5),#D'Ales
              66,#kafy
              70,#west
              2,#proto
              63,#chaccour
              2,
              69,#staedke
              2,
              62,#bradley
              2)
points(lis~c(rep(-1.5,20)),pch=effi2$study_symb,cex=1.6)
for(i in 1:20){
  text(-4,lis[i],paste0("(",ref_paper[i],")"),cex=0.8)
  
}
# 
# legend("bottomright",legend=c("Observed efficacy","Model predicted efficacy"),
#        col=c("grey","black"),pch=c(15,0),cex=1)

par(new = "TRUE",
    plt = c(0.9,0.95,0.3,0.7),
    las = 1,
    cex.axis = 1)

unique(effi$col_efficacy)
col_test2 = c("grey","darkred","red","orange","yellow","green","darkgreen")
levels = seq(-10,60,10)
plot.window(xlim = c(0.2, 0.8), ylim = range(levels))
rect(0, levels[-length(levels)], 1, levels[-1L], col = col_test2)
axis(2,las=2,at=c(0,10,20,30,40,50,60),cex.axis=1)
rect(0, levels[-length(levels)], 1, levels[-1L], col = col_test2)


##############################
##
##
## Figure 2 we can add in this one


########################################
##
##
## Figure 2
par(mfrow=c(1,3))
par(mar=c(5,5,2,2))

##histogram of obersved efficacy against prevalence
##Best prevalence estimates (log-logistic, all data)
##Best efficacy estiamtes

effi = read.csv("Post_processing/summary_data/DATA_RESOURCE_EFFICACY_global_for_fig1_drop_targetedFINAL.csv",header=TRUE)

prev = read.csv("Post_processing/summary_data/DATA_RESOURCE_PREVALENCE_global_update_droptargetedFINAL.csv",header=TRUE)

# effi = read.csv("Post_processing/data/DATA_RESOURCE_EFFICACY_global.csv",header=TRUE)
# effi = read.csv("Post_processing/data/DATA_RESOURCE_EFFICACY_global_for_fig1.csv",header=TRUE)

## a) hist of observed efficacy by intervention type
hist(effi$Obs_mean_eff,breaks=20,xlim=c(-0.6,1),
     ylab="Number of trial data observations",
     xlab = "Efficacy against prevalence (%)",
     cex.axis=1.6,cex.lab=1.6,
     yaxt="n",
     xaxt="n",bty="n",main="",col = "darkred")
axis(1,at=seq(-0.6,1,0.2),labels=seq(-60,100,20))
axis(2,las=2,at=seq(0,8,2),cex.axis=1.6,cex.lab=1.6)
hist(effi$Obs_mean_eff[effi$Intervention == "ITN" |
                         effi$Intervention == "pyr-IRS" |
                         effi$Intervention == "ITN + pyr/ben-IRS" |
                         effi$Intervention == "ITN + op-IRS" |
                         effi$Intervention == "PBO-ITN" |
                         effi$Intervention == "PBO-ITN + op-IRS"],add=TRUE,breaks=10,xlim=c(-0.6,1),
     col="red",border=NA)

hist(effi$Obs_mean_eff[effi$Intervention == "pyr-IRS" |
                         effi$Intervention == "ITN + pyr/ben-IRS" |
                         effi$Intervention == "ITN + op-IRS" |
                         effi$Intervention == "PBO-ITN" |
                         effi$Intervention == "PBO-ITN + op-IRS"],add=TRUE,breaks=10,xlim=c(-0.6,1),
     col="purple",border=NA)

hist(effi$Obs_mean_eff[effi$Intervention == "pyr-IRS" |
                         effi$Intervention == "ITN + pyr/ben-IRS" |
                         effi$Intervention == "ITN + op-IRS" |
                         effi$Intervention == "PBO-ITN"],add=TRUE,breaks=10,xlim=c(-0.6,1),
     col="blue",border=NA)

hist(effi$Obs_mean_eff[effi$Intervention == "pyr-IRS" |
                         effi$Intervention == "ITN + pyr/ben-IRS" |
                         effi$Intervention == "ITN + op-IRS"],add=TRUE,breaks=10,xlim=c(-0.6,1),
     col="aquamarine3",border=NA)
hist(effi$Obs_mean_eff[effi$Intervention == "pyr-IRS"] ,add=TRUE,breaks=3,xlim=c(-0.6,1),
     col="orange",border=NA)

legend("topleft",legend=c("Pyr-ITN",
                          "Pyr-ITN + IRS",
                          "Pyr-PBO ITN",
                          "Pyr-PBO ITN + IRS",
                          "IRS"),
       col=c("darkred","aquamarine3",
             "blue","purple","orange"),
       pch=15,cex=1.2,title="Intervention",
       bty="n")

## b) prevalence from log-logistic
## Investigate effi1

plot(prev$prev_obs~prev$prev_pred_4,ylim=c(0,1),xlim=c(0,1),
     main = "",
     xlab = "Transmission model prevalence estimate (%)", ##model etsimts
     ylab = "RCT observed prevalence (%)",##predict the data
     bty="n",
     xaxt = "n", yaxt = "n",cex.lab=1.6,cex.axis=1.6,pch="")
axis(1,at=seq(0,1,0.2),labels = seq(0,100,20),cex.lab=1.6,cex.axis=1.4)
axis(2,las=2,at=seq(0,1,0.2),labels = seq(0,100,20),cex.lab=1.4,cex.axis=1.4)
Y = seq(-10,10,length=20)
Yu = seq(-10,10,length=20)+0.1
Yl = seq(-10,10,length=20)-0.1
X = seq(-10,10,length=20)
lines(Y~X,lty=2)

points(prev$prev_obs~prev$prev_pred_4,
       col=as.character(prev$intvn_colour),pch=prev$study_symbol,cex=1.8)  



segments(x0=prev$prev_pred_4,x1=prev$prev_pred_4,
         y0=prev$prev_obs_max,y1=prev$prev_obs_min,
         col=as.character(prev$intvn_colour),lwd=1)
segments(x0=prev$prev_pred_4_max,x1=prev$prev_pred_4_min,
         y0=prev$prev_obs,y1=prev$prev_obs,
         col=as.character(prev$intvn_colour),lwd=1)

tapply(prev$study_symbol,prev$PI,mean)
legend("topleft",
       legend = c("Marbiah (Jun 1992)",
                  "D'Alessandro (Jul 1992)",
                  "Nevill (Jul 1993)",
                  "Curtis (Dec 1995)",
                  "Phillips-Howard (Jan 1997)",
                  "Henry (Jun 1999)",#,
                  "Corbel (Jun 2008)"),
       title = "RCT (start month, year)",
       col="black",
       pch=c(6,3,12,2,
             7,4,1),
       cex=1.2,bty="n",pt.cex = 1.4)

legend("bottomright",
       legend = c("Kafy (Apr 2011)",
                 # "Abuaku (May 2011)",
                  "West (Dec 2011)",
                  "Protopopoff (Mar 2014)",
                  "Bradley (Mar-Apr 2014)",
                  # "Loha (Sep 2014)",
                  "Chaccour (Oct 2016)",
                  "Staedke (Jul 2017)"),
       col="black",
       pch=c(5,#14,
             9,8,13,10,11),
       cex=1.2,bty="n",pt.cex = 1.4)





## c) Efficacy logistic
# setwd("C:/Users/esherrar/Documents/Rprojects/Model_validation_ITN_IRS")
# effi = read.csv("Post_processing/data/DATA_RESOURCE_EFFICACY_global.csv",header=TRUE)
# effi = read.csv("Post_processing/data/DATA_RESOURCE_EFFICACY_global_for_fig1_drop_targeted.csv",header=TRUE)
effi$col_time = ifelse(effi$Month < 0,"black",
                       ifelse(effi$Month > 0 & effi$Month < 6, "grey",
                              ifelse(effi$Month > 6 & effi$Month < 12, "aquamarine3",
                                     ifelse(effi$Month > 12 & effi$Month < 18, "blue",
                                            ifelse(effi$Month > 18 & effi$Month < 24, "darkblue","purple")))))






plot(effi$Obs_mean_eff ~ effi$Pred_4_eff,ylim=c(-0.6,1),xlim=c(-0.6,1),
     xlab = "Model predicted efficacy against prevalence (%)",
     ylab = "RCT measured efficacy against prevalence (%)",
     bty="n",
     xaxt = "n", yaxt = "n",cex.lab=1.6,cex.axis=1.6,pch="")
axis(1,at=seq(-0.6,1,0.2),labels = seq(-60,100,20),cex.lab=1.6,cex.axis=1.6)
axis(2,las=2,at=seq(-0.4,1,0.2),labels = seq(-40,100,20),cex.lab=1.6,cex.axis=1.6)
Y = seq(-10,10,length=20)
Yu = seq(-10,10,length=20)+0.1
Yl = seq(-10,10,length=20)-0.1
X = seq(-10,10,length=20)
lines(Y~X,lty=2)

points(effi$Obs_mean_eff ~ effi$Pred_4_eff,
       col=as.character(effi$col_time),pch=effi$study_symb,cex=1.8)

segments(x0=effi$Pred_4_eff,x1=effi$Pred_4_eff,
         y0=effi$eff_obs_min,y1=effi$eff_obs_max,
         col=as.character(effi$col_time),lwd=1)
segments(x0=effi$Pred_4_eff_min,x1=effi$Pred_4_eff_max,
         y0=effi$Obs_mean_eff,y1=effi$Obs_mean_eff,
         col=as.character(effi$col_time),lwd=1)


legend("topleft",legend=c("within 6-months",
                          "6 - 12-months",
                          "12 - 18-months",
                          "18 - 24-months",
                          "later than 24-months"),
       col=c("grey","lightgreen","darkgreen","blue","purple"),
       pch=15,cex=1.2,title="Observation time",
       bty="n")



par(xpd=NA,cex = 1.1)
text(x = -5, y = 1.1,"(A)")
text(x = -2.9, y = 1.1,"(B)")
text(x = -0.8, y = 1.1,"(C)")


#################################
##
##
## TABLE 1 RESULTS
##
##

## Table 1
length(prev$prev_obs) ## n = 73
summary.lm(lm(prev$prev_obs~prev$prev_pred_1+0)) ## 93.22 ## 0.90 
summary.lm(lm(prev$prev_obs~prev$prev_pred_2+0)) ## 92.89 ## 0.89
summary.lm(lm(prev$prev_obs~prev$prev_pred_3+0)) ## 93.61 ## 0.91
summary.lm(lm(prev$prev_obs~prev$prev_pred_4+0)) ## 95.31 ## 0.96
summary.lm(lm(prev$prev_obs~prev$prev_pred_5+0)) ## 94.41 ## 0.93
summary.lm(lm(prev$prev_obs~prev$prev_pred_6+0)) ## 95.35 ## 0.96

og1 = lm(prev$prev_pred_4~prev$prev_obs+0)




###################################

## Explore differences in efficacy
## for the different interventions

## ITNs or IRS
unique(prev$Intervention)
unique(effi$Intervention)

## Any net without IRS
effi$INTN2 = ifelse(effi$Intervention == "pyr-IRS","IRS",
                    ifelse(effi$Intervention == "ITN + pyr/ben-IRS","IRS",
                           ifelse(effi$Intervention == "ITN + op-IRS","IRS",
                                  ifelse(effi$Intervention == "PBO-ITN + op-IRS","IRS","ITN"))))
prev$INTN2 = ifelse(prev$Intervention == "None","None",
                    ifelse(prev$Intervention == "ITN_IRS","IRS",
                           ifelse(prev$Intervention == "IRS","IRS",
                           ifelse(prev$Intervention == "PBO_ITN_IRS","IRS","ITN"))))
net_effi= subset(effi,effi$INTN2 == "ITN")
net_prev= subset(prev,prev$INTN2 == "ITN")

itn_prev = subset(prev,prev$Intervention == "ITN")


## Any IRS              
effi$INTN3 = ifelse(effi$Intervention == "CTN","NO",
                    ifelse(effi$Intervention == "PBO-ITN","No","YES-IRS"))
prev$INTN3 = ifelse(prev$Intervention == "None","No",
                    ifelse(prev$Intervention == "ITN","No",
                           ifelse(prev$Intervention == "CTN","No",
                                  ifelse(prev$Intervention == "PBO_ITN","No","Yes"))))

SPRAY_effi= subset(effi,effi$INTN3 == "YES-IRS")
SPRAY_prev= subset(prev,prev$INTN3 == "Yes")


## PBO_nets
pbo_effi = subset(effi,effi$Intervention == "PBO-ITN")
pbo_prev = subset(prev,prev$Intervention == "PBO_ITN")

## PBO_nets
pbo_irs_effi = subset(effi,effi$Intervention == "PBO-ITN + op-IRS")
pbo_irs_prev = subset(prev,prev$Intervention == "PBO-ITN + op-IRS")


## Studies using any net without IRS
##Table 1 
##Table prevalence ## Any net (no IRS)
summary.lm(lm(net_prev$prev_obs~net_prev$prev_pred_1+0))
summary.lm(lm(net_prev$prev_obs~net_prev$prev_pred_2+0)) ## BEST 
summary.lm(lm(net_prev$prev_obs~net_prev$prev_pred_3+0)) 
summary.lm(lm(net_prev$prev_obs~net_prev$prev_pred_4+0))
summary.lm(lm(net_prev$prev_obs~net_prev$prev_pred_5+0))
summary.lm(lm(net_prev$prev_obs~net_prev$prev_pred_6+0))

## insecticide nets - pyr-itn
summary.lm(lm(itn_prev$prev_obs~itn_prev$prev_pred_1+0))
summary.lm(lm(itn_prev$prev_obs~itn_prev$prev_pred_2+0)) ## BEST 
summary.lm(lm(itn_prev$prev_obs~itn_prev$prev_pred_3+0)) 
summary.lm(lm(itn_prev$prev_obs~itn_prev$prev_pred_4+0))
summary.lm(lm(itn_prev$prev_obs~itn_prev$prev_pred_5+0))
summary.lm(lm(itn_prev$prev_obs~itn_prev$prev_pred_6+0))

## pyr-PBO
summary.lm(lm(pbo_prev$prev_obs~pbo_prev$prev_pred_1+0))
summary.lm(lm(pbo_prev$prev_obs~pbo_prev$prev_pred_2+0)) ## BEST 
summary.lm(lm(pbo_prev$prev_obs~pbo_prev$prev_pred_3+0)) 
summary.lm(lm(pbo_prev$prev_obs~pbo_prev$prev_pred_4+0))
summary.lm(lm(pbo_prev$prev_obs~pbo_prev$prev_pred_5+0))
summary.lm(lm(pbo_prev$prev_obs~pbo_prev$prev_pred_6+0))

## spray
summary.lm(lm(SPRAY_prev$prev_obs~SPRAY_prev$prev_pred_1+0))
summary.lm(lm(SPRAY_prev$prev_obs~SPRAY_prev$prev_pred_2+0)) ## BEST 
summary.lm(lm(SPRAY_prev$prev_obs~SPRAY_prev$prev_pred_3+0)) 
summary.lm(lm(SPRAY_prev$prev_obs~SPRAY_prev$prev_pred_4+0))
summary.lm(lm(SPRAY_prev$prev_obs~SPRAY_prev$prev_pred_5+0))
summary.lm(lm(SPRAY_prev$prev_obs~SPRAY_prev$prev_pred_6+0))


## Table 1
length(effi$Obs_mean_eff) ## N = 46
summary.lm(lm(effi$Obs_mean_eff~effi$Pred_1_eff+0))
summary.lm(lm(effi$Obs_mean_eff~effi$Pred_2_eff+0)) ## BEST 
summary.lm(lm(effi$Obs_mean_eff~effi$Pred_3_eff+0)) 
summary.lm(lm(effi$Obs_mean_eff~effi$Pred_4_eff+0))
summary.lm(lm(effi$Obs_mean_eff~effi$Pred_5_eff+0))
summary.lm(lm(effi$Obs_mean_eff~effi$Pred_6_eff+0))

## Table Efficacy Any net (no IRS)
summary.lm(lm(net_effi$Obs_mean_eff~net_effi$Pred_1_eff+0))
summary.lm(lm(net_effi$Obs_mean_eff~net_effi$Pred_2_eff+0)) 
summary.lm(lm(net_effi$Obs_mean_eff~net_effi$Pred_3_eff+0)) 
summary.lm(lm(net_effi$Obs_mean_eff~net_effi$Pred_4_eff+0))
summary.lm(lm(net_effi$Obs_mean_eff~net_effi$Pred_5_eff+0))
summary.lm(lm(net_effi$Obs_mean_eff~net_effi$Pred_6_eff+0))

## pbo 
summary.lm(lm(pbo_effi$Obs_mean_eff~pbo_effi$Pred_1_eff+0))
summary.lm(lm(pbo_effi$Obs_mean_eff~pbo_effi$Pred_2_eff+0)) ## BEST 
summary.lm(lm(pbo_effi$Obs_mean_eff~pbo_effi$Pred_3_eff+0)) 
summary.lm(lm(pbo_effi$Obs_mean_eff~pbo_effi$Pred_4_eff+0))
summary.lm(lm(pbo_effi$Obs_mean_eff~pbo_effi$Pred_5_eff+0))
summary.lm(lm(pbo_effi$Obs_mean_eff~pbo_effi$Pred_6_eff+0))

## any spray
summary.lm(lm(SPRAY_effi$Obs_mean_eff~SPRAY_effi$Pred_1_eff+0))
summary.lm(lm(SPRAY_effi$Obs_mean_eff~SPRAY_effi$Pred_2_eff+0)) ## BEST 
summary.lm(lm(SPRAY_effi$Obs_mean_eff~SPRAY_effi$Pred_3_eff+0)) 
summary.lm(lm(SPRAY_effi$Obs_mean_eff~SPRAY_effi$Pred_4_eff+0))
summary.lm(lm(SPRAY_effi$Obs_mean_eff~SPRAY_effi$Pred_5_eff+0))
summary.lm(lm(SPRAY_effi$Obs_mean_eff~SPRAY_effi$Pred_6_eff+0))




####################
##
## SUPPLEMENTS

## Table S2 Params

setwd("C:/Users/esherrar/Documents/Rprojects/Model_validation_ITN_IRS")

# Global set
ll_1 = readRDS("general_functions/Entomological model fits/rds files/Combined/Bioassay/fit_ew_comb_log_logistic.rds")
l_1 = readRDS("general_functions/Entomological model fits/rds files/Combined/Bioassay/fit_ew_comb_logistic.rds")

LL_fit <- extract(ll_1, permuted = TRUE)
L_fit <- extract(l_1, permuted = TRUE)

quantile(LL_fit$a,c(0.05,0.5,0.950))
quantile(LL_fit$c,c(0.05,0.5,0.950))

quantile(L_fit$a,c(0.05,0.5,0.950))
quantile(L_fit$c,c(0.05,0.5,0.950))

#West
ll_1w = readRDS("general_functions/Entomological model fits/rds files/Huts separate/Bioassay/fit_west_log_logistic.rds")
l_1w = readRDS("general_functions/Entomological model fits/rds files/Huts separate/Bioassay/fit_west_logistic.rds")

LL_fitw <- extract(ll_1w, permuted = TRUE)
L_fitw <- extract(l_1w, permuted = TRUE)

quantile(LL_fitw$a,c(0.05,0.5,0.950))
quantile(LL_fitw$c,c(0.05,0.5,0.950))

quantile(L_fitw$a,c(0.05,0.5,0.950))
quantile(L_fitw$c,c(0.05,0.5,0.950))


#East
ll_1e = readRDS("general_functions/Entomological model fits/rds files/Huts separate/Bioassay/fit_east_log_logistic.rds")
l_1e = readRDS("general_functions/Entomological model fits/rds files/Huts separate/Bioassay/fit_east_logistic.rds")

LL_fite <- extract(ll_1e, permuted = TRUE)
L_fite <- extract(l_1e, permuted = TRUE)

quantile(LL_fite$a,c(0.05,0.5,0.950))
quantile(LL_fite$c,c(0.05,0.5,0.950))

quantile(L_fite$a,c(0.05,0.5,0.950))
quantile(L_fite$c,c(0.05,0.5,0.950))



## Benefit of PBO
benefitall <- readRDS("general_functions/Entomological model fits/rds files/ento_pbo_benefit.rds")
pbo_bene <- extract(benefitall, permuted = TRUE)


quantile(pbo_bene$alpha1,c(0.05,0.5,0.95))
quantile(pbo_bene$alpha2,c(0.05,0.5,0.95))


## Deterrence and mortality
fit1_a <- readRDS("general_functions/Entomological model fits/rds files/Combined/Feeding attempt/deterrence_ew_fit.rds")
fit1_a_w <- readRDS("general_functions/Entomological model fits/rds files/Huts separate/Feeding attempt/deterrence_west_fit_corrected.rds")
fit1_a_e <- readRDS("general_functions/Entomological model fits/rds files/Huts separate/Feeding attempt/deterrence_east_fit_corrected.rds")

fit1_a_fit <- extract(fit1_a, permuted = TRUE)
fit1_a_fitW <- extract(fit1_a_w, permuted = TRUE)
fit1_a_fitE <- extract(fit1_a_e, permuted = TRUE)

quantile(fit1_a_fit$c,c(0.05,0.5,0.95))
quantile(fit1_a_fitW$c,c(0.05,0.5,0.95))
quantile(fit1_a_fitE$c,c(0.05,0.5,0.95))

quantile(fit1_a_fit$d,c(0.05,0.5,0.95))
quantile(fit1_a_fitW$d,c(0.05,0.5,0.95))
quantile(fit1_a_fitE$d,c(0.05,0.5,0.95))

quantile(fit1_a_fit$e,c(0.05,0.5,0.95))
quantile(fit1_a_fitW$e,c(0.05,0.5,0.95))
quantile(fit1_a_fitE$e,c(0.05,0.5,0.95))


## Success and mortality
fit3_a <- readRDS("general_functions/Entomological model fits/rds files/Combined/Feeding attempt/Succfed_ew_fit.rds")
fit3_aW <- readRDS("general_functions/Entomological model fits/rds files/Huts separate/Feeding attempt/Succfed_west_fit_corrected.rds")
fit3_aE <- readRDS("general_functions/Entomological model fits/rds files/Huts separate/Feeding attempt/Succfed_east_fit_corrected.rds")

fit3_a_fit <- extract(fit3_a, permuted = TRUE)
fit3_a_fitW <- extract(fit3_aW, permuted = TRUE)
fit3_a_fitE <- extract(fit3_aE, permuted = TRUE)

quantile(fit3_a_fit$a,c(0.05,0.5,0.95))
quantile(fit3_a_fitW$a,c(0.05,0.5,0.95))
quantile(fit3_a_fitE$a,c(0.05,0.5,0.95))

quantile(fit3_a_fit$b,c(0.05,0.5,0.95))
quantile(fit3_a_fitW$b,c(0.05,0.5,0.95))
quantile(fit3_a_fitE$b,c(0.05,0.5,0.95))

####################
##
## SUPPLEMENTS
## Not including 
## Fig S1 - S16
##
##

setwd("C:/Users/esherrar/Documents/Rprojects/Model_validation_ITN_IRS")

effi = read.csv("Post_processing/data/efficacy_relative_baselineFINAL.csv",header=TRUE)

effi$col_efficacy = ifelse(effi$Obs_mean_eff < 0,"grey",
                           ifelse(effi$Obs_mean_eff > 0 & effi$Obs_mean_eff < 0.1, "darkred",
                                  ifelse(effi$Obs_mean_eff > 0.1 & effi$Obs_mean_eff < 0.2, "red",
                                         ifelse(effi$Obs_mean_eff > 0.2 & effi$Obs_mean_eff < 0.3, "orange",
                                                ifelse(effi$Obs_mean_eff > 0.3 & effi$Obs_mean_eff < 0.4, "yellow",
                                                       ifelse(effi$Obs_mean_eff > 0.4 & effi$Obs_mean_eff < 0.5, "green","darkgreen"))))))

effi$col_time = ifelse(effi$time_months < 0,"black",
                               ifelse(effi$time_months > 0 & effi$time_months < 6, "grey",
                                      ifelse(effi$time_months > 6 & effi$time_months < 12, "aquamarine3",
                                             ifelse(effi$time_months > 12 & effi$time_months < 18, "blue",
                                                    ifelse(effi$time_months > 18 & effi$time_months < 24, "darkblue","purple")))))




plot(effi$Obs_mean_eff ~ effi$Pred_2_eff,ylim=c(-0.6,1),xlim=c(-0.6,1),
     xlab = "Model predicted efficacy against prevalence (%)",
     ylab = "RCT measured efficacy against prevalence (%)",
     bty="n",
     xaxt = "n", yaxt = "n",cex.lab=1.6,cex.axis=1.6,pch="")
axis(1,at=seq(-0.6,1,0.2),labels = seq(-60,100,20),cex.lab=1.6,cex.axis=1.6)
axis(2,las=2,at=seq(-0.4,1,0.2),labels = seq(-40,100,20),cex.lab=1.6,cex.axis=1.6)
Y = seq(-10,10,length=20)
Yu = seq(-10,10,length=20)+0.1
Yl = seq(-10,10,length=20)-0.1
X = seq(-10,10,length=20)
lines(Y~X,lty=2)

points(effi$Obs_mean_eff ~ effi$Pred_2_eff,
       col=as.character(effi$col_time),
       pch=effi$study_symb,cex=1.8)

segments(x0=effi$Pred_2_eff,x1=effi$Pred_2_eff,
         y0=effi$eff_obs_min,y1=effi$eff_obs_max,
         col=as.character(effi$col_time),
         lwd=1)
segments(x0=effi$Pred_2_eff_min,x1=effi$Pred_2_eff_max,
         y0=effi$Obs_mean_eff,y1=effi$Obs_mean_eff,
         col=as.character(effi$col_time),
         lwd=1)


legend("topleft",legend=c("within 6-months",
                          "6 - 12-months",
                          "12 - 18-months",
                          "18 - 24-months",
                          "later than 24-months"),
       col=c("grey","lightgreen","darkgreen","blue","purple"),
       pch=15,cex=1.2,title="Observation time",
       bty="n")



par(xpd=NA,cex = 1.1)
text(x = -0.80, y = 1.3,"(D)")
text(x = -2.9, y = 1.3,"(C)")
text(x = -0.8, y = 3.4,"(B)")
text(x = -2.9, y = 3.4,"(A)")


########################################
##
## Figure S18

## create a plotting function for the repeated image

plot_prev_f = function(data_variable,
                       data_varU,data_varL,
                       main_header){
  
  plot(prev$prev_obs~prev[,data_variable],ylim=c(0,1),xlim=c(0,1),
       main = paste(main_header),
       xlab = "Transmission model prevalence estimate (%)", ##model etsimts
       ylab = "RCT observed prevalence (%)",##predict the data
       bty="n",
       xaxt = "n", yaxt = "n",cex.lab=1.6,cex.axis=1.6,pch="")
  axis(1,at=seq(0,1,0.2),labels = seq(0,100,20),cex.lab=1.6,cex.axis=1.4)
  axis(2,las=2,at=seq(0,1,0.2),labels = seq(0,100,20),cex.lab=1.4,cex.axis=1.4)
  Y = seq(-10,10,length=20)
  Yu = seq(-10,10,length=20)+0.1
  Yl = seq(-10,10,length=20)-0.1
  X = seq(-10,10,length=20)
  lines(Y~X,lty=2)
  
  points(prev$prev_obs~prev[,data_variable],
         col=as.character(prev$intvn_colour),pch=prev$study_symbol,cex=1.8)  
  
  
  
  segments(x0=prev[,data_variable],x1=prev[,data_variable],
           y0=prev$prev_obs_max,y1=prev$prev_obs_min,
           col=as.character(prev$intvn_colour),lwd=1)
  segments(x0=prev[,data_varL],x1=prev[,data_varU],
           y0=prev$prev_obs,y1=prev$prev_obs,
           col=as.character(prev$intvn_colour),lwd=1)
  
  
}

names(prev)
par(mfrow=c(2,3))
plot_prev_f(data_variable = 7,
                       data_varU = 8,data_varL = 9,
                       main_header = "Logistic function and all EHT data")
plot_prev_f(data_variable = 10,
            data_varU = 11,data_varL = 12,
            main_header = "Logistic function and West African EHT data")
plot_prev_f(data_variable = 13,
            data_varU = 14,data_varL = 15,
            main_header = "Logistic function and East African EHT data")

tapply(prev$study_symbol,prev$PI,mean)
legend("topleft",
       legend = c("Marbiah (Jun 1992)",
                  "D'Alessandro (Jul 1992)",
                  "Nevill (Jul 1993)"),
       title = "RCT (start month, year)",
       col="black",
       pch=c(6,3,12),
       cex=1.2,bty="n",pt.cex = 1.4)

legend("bottomright",
       legend = c("Kafy (Apr 2011)",
                  # "Abuaku (May 2011)",
                  "West (Dec 2011)",
                  "Protopopoff (Mar 2014)"),
       col="black",
       pch=c(5,#14,
             9,8),
       cex=1.2,bty="n",pt.cex = 1.4)

plot_prev_f(data_variable = 16,
            data_varU = 17,data_varL = 18,
            main_header = "Log-logistic function and all EHT data")
plot_prev_f(data_variable = 19,
            data_varU = 20,data_varL = 21,
            main_header = "Log-logistic function and West African EHT data")
plot_prev_f(data_variable = 22,
            data_varU = 23,data_varL = 24,
            main_header = "Log-logistic function and East African EHT data")

legend("topleft",
       legend = c("Curtis (Dec 1995)",
                  "Phillips-Howard (Jan 1997)",
                  "Henry (Jun 1999)",#,
                  "Corbel (Jun 2008)"),
       title = "",
       col="black",
       pch=c(2,
             7,4,1),
       cex=1.2,bty="n",pt.cex = 1.4)

legend("bottomright",
       legend = c("Bradley (Mar-Apr 2014)",
                  # "Loha (Sep 2014)",
                  "Chaccour (Oct 2016)",
                  "Staedke (Jul 2017)"),
       col="black",
       pch=c(13,10,11),
       cex=1.2,bty="n",pt.cex = 1.4)


##################################
##
## Absolute difference
setwd("C:/Users/esherrar/Documents/Rprojects/ibm_rct_prediction")
abs_dat = read.csv("Post_processing/summary_data/absolute_difference.csv",header=TRUE)
par(mfrow=c(1,2))
par(mar=c(4,4,2,2))
plot(abs_dat$absolute_diff_OX_OBS~abs_dat$absolute_diff_OX_MOD,ylim=c(-0.2,0.6),xlim=c(-0.2,0.6),
     main = "",
     xlab = "Difference to baseline prevalence, model estimate", ##model etsimts
     ylab = "Difference to baseline prevalence, RCT measure",##predict the data
     bty="n",
     xaxt = "n", yaxt = "n",cex.lab=1,cex.axis=1,pch="")
axis(1,at=seq(-0.2,0.6,0.2),labels = seq(-20,60,20),cex.lab=1,cex.axis=1)
axis(2,las=2,at=seq(-0.2,0.6,0.2),labels = seq(-20,60,20),cex.lab=1,cex.axis=1)
Y = seq(-10,10,length=20)
Yu = seq(-10,10,length=20)+0.1
Yl = seq(-10,10,length=20)-0.1
X = seq(-10,10,length=20)
lines(Y~X,lty=2)

points(abs_dat$absolute_diff_OX_OBS~abs_dat$absolute_diff_OX_MOD,
       col=as.character(abs_dat$intvn_colour),pch=abs_dat$study_symbol,cex=1.8)  

tapply(abs_dat$study_symbol,abs_dat$PI,mean)
legend("topleft",
       legend = c("Marbiah (Jun 1992)",
                  "D'Alessandro (Jul 1992)",
                  "Nevill (Jul 1993)",
                  "Curtis (Dec 1995)",
                  "Phillips-Howard (Jan 1997)",
                  "Henry (Jun 1999)",#,
                  "Corbel (Jun 2008)"),
       title = "RCT (start month, year)",
       col="black",
       pch=c(6,3,12,2,
             7,4,1),
       cex=1,bty="n",pt.cex = 1.2)

legend("bottomright",
       legend = c("Kafy (Apr 2011)",
                  # "Abuaku (May 2011)",
                  "West (Dec 2011)",
                  "Protopopoff (Mar 2014)",
                  "Bradley (Mar-Apr 2014)",
                  # "Loha (Sep 2014)",
                  "Chaccour (Oct 2016)",
                  "Staedke (Jul 2017)"),
       col="black",
       pch=c(5,#14,
             9,8,13,10,11),
       cex=1,bty="n",pt.cex = 1.2)


abs_dat$Calc_diff = abs_dat$absolute_diff_OX_OBS-abs_dat$absolute_diff_OX_MOD
par(mar=c(8,7,6,2))
boxplot(abs_dat$Calc_diff ~ abs_dat$intvn_order,
        ylab = "Observed - modelled difference to baseline prevalence",
        xlab = "",ylim=c(-0.3,0.3),
        xaxt="n",col=c("grey","darkred","red","blue","aquamarine3","purple","orange"))
##1 None, 2 CTN, 3 ITN, 4 PBO, 5 ITN + IRS, 6 PBO + IRS, 7 IRS only
axis(1, las=2,at=1:7, labels=c("None (16)","CTN (16)","LLIN (14)","PBO-ITN (7)","LLIN + IRS (12)","PBO-ITN + IRS (4)","IRS (4)"))
abline(h=0,lty=2)
abline(h=-0.1,lty=4,col="grey30")
abline(h=0.1,lty=4,col="grey30")
boxplot(abs_dat$Calc_diff ~ abs_dat$intvn_order,
        ylab = "Observed - modelled difference to baseline prevalence",
        xlab = "",ylim=c(-0.3,0.3),
        xaxt="n",col=c("grey","darkred","red","blue","aquamarine3","purple","orange"),add=TRUE)


par(xpd=NA,cex = 1.1)
text(x = -0.8, y = 0.4,"(B)")
text(x = -12, y = 0.45,"(A)")

tapply(abs_dat$Intervention,
       abs_dat$intvn_order,length)

summary.lm(aov(abs_dat$Calc_diff ~ as.factor(abs_dat$intvn_order)))


################################################
##
## Table S4 East African EHT to predict studies in the East, and West equivalent

effi = read.csv("Post_processing/summary_data/DATA_RESOURCE_EFFICACY_global_for_fig1_drop_targetedFINAL.csv",header=TRUE)

prev = read.csv("Post_processing/summary_data/DATA_RESOURCE_PREVALENCE_global_update_droptargetedFINAL.csv",header=TRUE)

## create west and east data
west_prev = subset(prev,prev$location == "West")
East_prev = subset(prev,prev$location == "East")

## Table 4
length(west_prev$prev_obs) ## n = 17
length(East_prev$prev_obs) ## n = 52
## just log-logistic
summary.lm(lm(west_prev$prev_obs~west_prev$prev_pred_4+0)) ## 95.31 ## 0.96
summary.lm(lm(west_prev$prev_obs~west_prev$prev_pred_5+0)) ## 94.41 ## 0.93
summary.lm(lm(west_prev$prev_obs~west_prev$prev_pred_6+0)) ## 95.35 ## 0.96

summary.lm(lm(East_prev$prev_obs~East_prev$prev_pred_4+0)) ## 95.31 ## 0.96
summary.lm(lm(East_prev$prev_obs~East_prev$prev_pred_5+0)) ## 94.41 ## 0.93
summary.lm(lm(East_prev$prev_obs~East_prev$prev_pred_6+0)) ## 95.35 ## 0.96

## Table 4
west_eff = subset(effi,effi$location == "West")
East_eff = subset(effi,effi$location == "East")

length(west_eff$Obs_mean_eff) ## N = 8
length(East_eff$Obs_mean_eff) ## N = 36

summary.lm(lm(west_eff$Obs_mean_eff~west_eff$Pred_4_eff+0))
summary.lm(lm(west_eff$Obs_mean_eff~west_eff$Pred_5_eff+0))
summary.lm(lm(west_eff$Obs_mean_eff~west_eff$Pred_6_eff+0))

summary.lm(lm(East_eff$Obs_mean_eff~East_eff$Pred_4_eff+0))
summary.lm(lm(East_eff$Obs_mean_eff~East_eff$Pred_5_eff+0))
summary.lm(lm(East_eff$Obs_mean_eff~East_eff$Pred_6_eff+0))
