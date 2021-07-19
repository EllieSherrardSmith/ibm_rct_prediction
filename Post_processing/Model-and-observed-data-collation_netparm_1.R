##########################################################
##
## Post-processing
##
## Extracting each study data for prevalence estimates
## Matching the time for observed prevalence data to the model predictions

## Drawing Figures S3 - S17 in Supplementary Material

## Model outputs not provided - files too large - can be made available on request

## Define the 'bootstrapped' runs for uncertainty
run = 1:1000
 

##########################
##
## ID 1 Bradley
run=1:1000
Tup_arm1a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/bradley/arm1/VALIDATION_trial_arm_1_net_1_0_1.txt",header=TRUE)  ## filepath
Tup_arm2a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/bradley/arm2/VALIDATION_trial_arm_2_net_1_0_1.txt",header=TRUE)  ## filepath

Tup_arm1 = array(dim=c(nrow(Tup_arm1a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:1000){
  Tup_arm1[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/bradley/arm1/VALIDATION_trial_arm_1_net_1_0_",run[i],".txt"),header=TRUE)$prev_2_14
}
Tup_arm2 = array(dim=c(nrow(Tup_arm2a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:1000){
  Tup_arm2[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/bradley/arm2/VALIDATION_trial_arm_2_net_1_0_",run[i],".txt"),header=TRUE)$prev_2_14
}

plot(Tup_arm1[,1] ~ Tup_arm1a[,1],ylim = c(0,0.4),xlim=c(0.5,2),ylab="Prevalence in children 2 to 14 years (%)",
     xlab = "Time in years",pch="",yaxt="n",main="Bradley",xaxt="n")
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
axis(1,at=c(1:2),labels=c("Jan 2014","Jan 2015"))

abline(v=1,lty=2)
arrows(x0 = 1.25, y0 = 0.05, x1 = 1.25, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)

legend("topright",legend = c("Pyrethroid ITN + Pyrethroid IRS",
                             "Pyrethroid ITN + Carbamate IRS"),
       col = c("blue","darkgreen"), pch = 15,bty="n")

upp_Tu = low_Tu = med_Tu = array(dim=c(1301,2))
for(t in 1:1301){
  upp_Tu[t,1] = quantile(Tup_arm1[t,],0.975,na.rm=TRUE)
  low_Tu[t,1] = quantile(Tup_arm1[t,],0.025,na.rm=TRUE)
  med_Tu[t,1] = quantile(Tup_arm1[t,],0.5,na.rm=TRUE)
  
  upp_Tu[t,2] = quantile(Tup_arm2[t,],0.975,na.rm=TRUE)
  low_Tu[t,2] = quantile(Tup_arm2[t,],0.025,na.rm=TRUE)
  med_Tu[t,2] = quantile(Tup_arm2[t,],0.5,na.rm=TRUE)
  
}
cols_Tu = adegenet::transp(c("darkgreen","blue"),0.4)
cols_Tu2 = c("darkgreen","blue")
lwd_Tu = c(2,2)
lty_Tu=c(1,4)
for(i in 1:2){
  polygon(c(Tup_arm1a[,1],rev(Tup_arm1a[,1])),c(upp_Tu[,i],rev(low_Tu[,i])),col = cols_Tu[i],border=NA)
  lines(med_Tu[,i] ~ Tup_arm1a[,1],lwd=lwd_Tu[i],lty=lty_Tu[i],col=cols_Tu2[i])
}
Tup_prev_recorded_t1 = c(0.175,0.168)
Tup_prev_recorded_t2 = c(0.184, 0.232)

Tup_time_match = c((1/12) * c(15,18))

effect_TU = (Tup_prev_recorded_t2 - Tup_prev_recorded_t1)/Tup_prev_recorded_t2



points(Tup_prev_recorded_t1 ~ Tup_time_match, col="darkgreen",pch=19,cex=1.5)
points(Tup_prev_recorded_t2 ~ c(Tup_time_match+0.04), col="blue",pch=19,cex=1.5)

## Add in any uncertainty bars from the trial (need to acess these if there are any)
Tup_prev_recorded_t1 = c(0.175,0.168)
Tup_prev_recorded_t2 = c(0.184, 0.232)
Tup_prev_record_t1_min <- c(0.123,0.111)
Tup_prev_record_t1_max <- c(0.243,0.247)
Tup_prev_record_t2_min <- c(0.127,0.16)
Tup_prev_record_t2_max <- c(0.259,0.323)

segments(x0 = Tup_time_match, x1 = Tup_time_match, y0 = Tup_prev_record_t1_min, y1 = Tup_prev_record_t1_max, col="darkgreen")
segments(x0 = c(Tup_time_match+0.04), x1 = c(Tup_time_match+0.04), y0 = Tup_prev_record_t2_min, y1 = Tup_prev_record_t2_max, col="blue")

VALES = 261+52*Tup_time_match
Tup_arm1a$year[261+52*Tup_time_match]

model_estimates_TU1 = c(med_Tu[326,1],med_Tu[339,1])
model_estimates_TU2 = c(med_Tu[326,2],med_Tu[339,2])

model_Effect_TU =  (model_estimates_TU2 - model_estimates_TU1)/model_estimates_TU2

model_estimates_TU1_min = c(low_Tu[326,1],low_Tu[339,1])
model_estimates_TU1_max = c(upp_Tu[326,1],upp_Tu[339,1])
model_estimates_TU2_min = c(low_Tu[326,2],low_Tu[339,2])
model_estimates_TU2_max = c(upp_Tu[326,2],upp_Tu[339,2])

model_Effect_TU_min =  (model_estimates_TU2 - model_estimates_TU1_min)/model_estimates_TU2
model_Effect_TU_max =  (model_estimates_TU2 - model_estimates_TU1_max)/model_estimates_TU2

effect_TU_min = (Tup_prev_recorded_t2 - Tup_prev_record_t1_min)/Tup_prev_recorded_t2
effect_TU_max = (Tup_prev_recorded_t2 - Tup_prev_record_t1_max)/Tup_prev_recorded_t2




######################################
##
## ID 2 Chaccour et al 2020

Tpp_arm1a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/chaccour/arm1/VALIDATION_trial_arm_1_net_1_0_1.txt",header=TRUE)
Tpp_arm1 = array(dim=c(nrow(Tpp_arm1a),1000)) ## this will be 1000 once all the runs are done
for(i in c(1:667,669:1000)){
  Tpp_arm1[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/chaccour/arm1/VALIDATION_trial_arm_1_net_1_0_",run[i],".txt"),header=TRUE)$prev_all
}

Tpp_arm2a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/chaccour/arm2/VALIDATION_trial_arm_2_net_1_0_1.txt",header=TRUE)
Tpp_arm2 = array(dim=c(nrow(Tpp_arm2a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Tpp_arm2)){
  Tpp_arm2[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/chaccour/arm2/VALIDATION_trial_arm_2_net_1_0_",run[i],".txt"),header=TRUE)$prev_all
}

plot(Tpp_arm1[,1] ~ Tpp_arm1a[,1],ylim = c(0,1),xlim=c(-0.5,2.5),ylab="All-age prevalence (%)",
     xlab = "Time in years",pch="",yaxt="n",xaxt="n",main = "Chaccour")
axis(1,at = 0:3,labels=c("Jun 2016","Jun 2017","Jun 2018","Jun 2019"))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))


abline(v=0.6666667,lty=2)## November 2016 distribution of IRS
arrows(x0 = 1, y0 = 0.1, x1 = 1, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 1)##nets
arrows(x0 = 0.4166667, y0 = 0.1, x1 = 0.4166667, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)
arrows(x0 = 1.4166667, y0 = 0.1, x1 = 1.4166667, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)

legend("topright",legend = c("Pyrethroid ITN",
                             "Pyrethroid ITN + Organophosphate IRS"),
       col = c("darkred","gold2"), pch = 15,bty="n")


upp_Tp = low_Tp = med_Tp = array(dim=c(1301,2))
for(t in 1:1301){
  upp_Tp[t,1] = quantile(Tpp_arm1[t,],0.975,na.rm=TRUE)
  low_Tp[t,1] = quantile(Tpp_arm1[t,],0.025,na.rm=TRUE)
  med_Tp[t,1] = quantile(Tpp_arm1[t,],0.5,na.rm=TRUE)
  
  upp_Tp[t,2] = quantile(Tpp_arm2[t,],0.975,na.rm=TRUE)
  low_Tp[t,2] = quantile(Tpp_arm2[t,],0.025,na.rm=TRUE)
  med_Tp[t,2] = quantile(Tpp_arm2[t,],0.5,na.rm=TRUE)
  
}
cols_Tp = adegenet::transp(c("darkred","gold2"),0.4)
cols_Tp2 = c("darkred","gold2")
lwd_Tp = c(2,2)
lty_Tp = c(1,4)

for(i in 1:2){
  polygon(c(Tpp_arm1a[,1],rev(Tpp_arm1a[,1])),c(upp_Tp[,i],rev(low_Tp[,i])),col = cols_Tp[i],border=NA)
  lines(med_Tp[,i] ~ Tpp_arm1a[,1],col=cols_Tp[i],lwd=lwd_Tp[i],lty=lty_Tp[i])
}

# All age
# model_Effect_TP = (Tpp_arm1[325,]-Tpp_arm2[325,])/Tpp_arm1[325,]


Tpp_prev_recorded_t1 = c(0.44,0.43)
Tpp_prev_recorded_t2 = c(0.43,0.34)

Tpp_time_match = c(0.6666667,1.6666667)

effect_TP = (Tpp_prev_recorded_t1[2] - Tpp_prev_recorded_t2[2])/Tpp_prev_recorded_t1[2]

points(Tpp_prev_recorded_t1 ~ Tpp_time_match, col="red",pch=18,cex=1.5)
points(Tpp_prev_recorded_t2 ~ c(Tpp_time_match+0.04), col="gold2",pch=18,cex=1.5)

model_estimates_TP1 = c(med_Tp[348,1])
model_estimates_TP2 = c(med_Tp[348,2])

model_Effect_TP2 =  (model_estimates_TP1 - model_estimates_TP2)/model_estimates_TP1

model_estimates_TP1_min = c(low_Tp[348,1])
model_estimates_TP1_max = c(upp_Tp[348,1])
model_estimates_TP2_min = c(low_Tp[348,2])
model_estimates_TP2_max = c(upp_Tp[348,2])

model_Effect_TP_min =  (model_estimates_TP1 - model_estimates_TP2_min)/model_estimates_TP1
model_Effect_TP_max =  (model_estimates_TP1 - model_estimates_TP2_max)/model_estimates_TP1

# effect_TP_min = (Tpp_prev_record_t1_min - Tpp_prev_record_t2_min)/model_estimates_TP1
# effect_TP_max = (Tpp_prev_record_t1_max - Tpp_prev_record_t2_max)/model_estimates_TP1





################################################################
## A 
## ID 3 Corbel et al. 2012
## Data TRIAL
# Tap_arm1a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/corbel/arm1/VALIDATION_trial_arm_1_net_1_0_1.txt",header=TRUE)
# Tap_arm1 = array(dim=c(nrow(Tap_arm1a),1000)) ## this will be 1000 once all the runs are done
# for(i in 1:ncol(Tap_arm1)){
#   Tap_arm1[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/corbel/arm1/VALIDATION_trial_arm_1_net_1_0_",run[i],".txt"),header=TRUE)$prev_0_6
# }
Tap_arm2a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/corbel/arm2/VALIDATION_trial_arm_2_net_1_0_1.txt",header=TRUE)
Tap_arm2 = array(dim=c(nrow(Tap_arm2a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Tap_arm2)){
  Tap_arm2[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/corbel/arm2/VALIDATION_trial_arm_2_net_1_0_",run[i],".txt"),header=TRUE)$prev_0_6
}
# Tap_arm3a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/corbel/arm3/VALIDATION_trial_arm_3_net_1_0_1.txt",header=TRUE)
# Tap_arm3 = array(dim=c(nrow(Tap_arm3a),1000)) ## this will be 1000 once all the runs are done
# for(i in 1:ncol(Tap_arm3)){
#   Tap_arm3[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/corbel/arm3/VALIDATION_trial_arm_3_net_1_0_",run[i],".txt"),header=TRUE)$prev_0_6
# }


plot(Tap_arm2[,1] ~ Tap_arm2a[,1],ylim = c(0,0.4),xlim=c(-1,2),ylab="Prevalence in children Under 6 years (%)",
     xlab = "Time in years",pch="",yaxt="n",main="Corbel",xaxt="n")
axis(2,las=2,at=seq(0,0.4,0.2),labels=seq(0,40,20))
axis(1,at=seq(-1,2,1),labels=c("Jun 2007","Jun 2008","Jun 2009","Jun 2010"))

abline(v=0.0630137,lty=2) ## June 23rd
arrows(x0 = 0.0630137, y0 = 0.05, x1 = 0.0630137, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)
arrows(x0 = 0.7296804, y0 = 0.05, x1 = 0.7296804, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)
arrows(x0 = 1.396347, y0 = 0.05, x1 = 1.396347, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)



legend("topright",legend = c("Pyrethroid ITN (targeted)",
                             "Pyrethroid ITN (universal)",
                             "Pyrethroid ITN (universal) + Carbamate IRS"),
       col = c("blue","darkred","darkgreen"), pch = 15,bty="n")



upp_Ta = low_Ta = med_Ta = array(dim=c(1301,4))
for(t in 1:1301){
  upp_Ta[t,1] = quantile(Tap_arm2[t,],0.975,na.rm=TRUE)
  low_Ta[t,1] = quantile(Tap_arm2[t,],0.025,na.rm=TRUE)
  med_Ta[t,1] = quantile(Tap_arm2[t,],0.5,na.rm=TRUE)
  
}
## Interventions
## 1 is [X1] LLIN targeted to pregnant women and children <6 years (TLLIN)
## 2 is [X2] LLIN universal coverage (ULLIN)
## 3 is [X3] TLLIN + full coverage of carbamate IRS (TLLIN + IRS)

cols_Ta = adegenet::transp(c("blue","darkred","darkgreen","darkgreen"),0.4)
cols_Ta2 = c("blue","darkred","darkgreen","darkgreen")
lty_Ta = c(1,1,1,6)
lwd_Ta = c(2,2,2,2)

polygon(c(Tap_arm2a[,1],rev(Tap_arm2a[,1])),c(upp_Ta[,1],rev(low_Ta[,1])),col = cols_Ta[2],border=NA)
lines(med_Ta[,1] ~ Tap_arm2a[,1],lty=lty_Ta[2],lwd=lwd_Ta[2],col=cols_Ta2[2])

Tao_arm2 = c(0.301,0.200)
Tao_arm2_max = c(0.353,0.24)
Tao_arm2_min = c(0.248,0.16)

time_matchA = (1/12) * c(0,17)##because the surveys were done 2 months into the year

points(Tao_arm2 ~ c(time_matchA), col="red",pch=15,cex=1.5)

segments(x0=c(time_matchA),x1=c(time_matchA),y0=Tao_arm2_min,y1=Tao_arm2_max,lty=1,lwd=1.5,col="red")

# model_estimates_TA1 = mean(Tap_arm1[335,])
model_estimates_TA2 = mean(Tap_arm2[335,])
# model_estimates_TA3 = mean(Tap_arm3[335,])
# # model_estimates_TA4 = mean(Tap_arm4[343,])
# model_Effect_TA = c((med_Ta[335,1] - med_Ta[335,2])/med_Ta[335,1],
#                     (med_Ta[335,1] - med_Ta[335,3])/med_Ta[335,1])
# 
# model_Effect_TA_min = c((med_Ta[335,1] - low_Ta[335,2])/med_Ta[335,1],
#                         (med_Ta[335,1] - low_Ta[335,3])/med_Ta[335,1])
# 
# model_Effect_TA_max = c((med_Ta[335,1] - upp_Ta[335,2])/med_Ta[335,1],
#                         (med_Ta[335,1] - upp_Ta[335,3])/med_Ta[335,1])
# 


################################################################
## B 
## ID 4 Curtis et al. 1998
## Data TRIAL

## Identify the variable of interest
## i.e. prevalence in under 6 year olds
Tbp_arm1a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/curtis/arm1/VALIDATION_trial_arm_1_net_1_0_1.txt",header=TRUE)
Tbp_arm1 = Bclin_inc_U5_arm1 = array(dim=c(nrow(Tbp_arm1a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Tbp_arm1)){
  Tbp_arm1[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/curtis/arm1/VALIDATION_trial_arm_1_net_1_0_",run[i],".txt"),header=TRUE)$prev_0_6
  # Bclin_inc_U5_arm1[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL/curtis/arm1/VALIDATION_trial_arm_1_net_1_0_",run[i],".txt"),header=TRUE)$clin_inc_0_5
}
Tbp_arm2a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/curtis/arm2/VALIDATION_trial_arm_2_net_1_0_1.txt",header=TRUE)
Tbp_arm2 = Bclin_inc_U5_arm2 = array(dim=c(nrow(Tbp_arm2a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Tbp_arm2)){
  Tbp_arm2[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/curtis/arm2/VALIDATION_trial_arm_2_net_1_0_",run[i],".txt"),header=TRUE)$prev_0_6
  Bclin_inc_U5_arm2[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/curtis/arm2/VALIDATION_trial_arm_2_net_1_0_",run[i],".txt"),header=TRUE)$clin_inc_0_5
}
Tbp_arm3a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/curtis/arm3/VALIDATION_trial_arm_3_net_1_0_1.txt",header=TRUE)
Tbp_arm3 = Bclin_inc_U5_arm3 = array(dim=c(nrow(Tbp_arm3a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Tbp_arm3)){
  Tbp_arm3[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/curtis/arm3/VALIDATION_trial_arm_3_net_1_0_",run[i],".txt"),header=TRUE)$prev_0_6
  Bclin_inc_U5_arm3[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/curtis/arm3/VALIDATION_trial_arm_3_net_1_0_",run[i],".txt"),header=TRUE)$clin_inc_0_5
}

plot(Tbp_arm1[,1] ~ Tbp_arm1a[,1],ylim = c(0,1),xlim=c(-1,2),ylab="Prevalence in children Under 6 years (%)",
     xlab = "Time in years",pch="",yaxt="n",main="Curtis",xaxt="n")
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
axis(1,at=seq(-1,2,1),labels=c("Jan 1995","Jan 1996","Jan 1997","Jan 1998"))

abline(v=-0.0833,lty=2)
arrows(x0 = 0, y0 = 0.05, x1 = 0, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)
arrows(x0 = 0.6666667, y0 = 0.05, x1 = 0.6666667, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)

legend("topright",legend = c("No intervention",
                             "Pyrethroid ITN",
                             "Pyrethroid IRS"),
       col = c("grey","darkred","purple"), pch = 15,bty="n")
upp_Tb = low_Tb = med_Tb = array(dim=c(1301,3))
for(t in 1:1301){
  upp_Tb[t,1] = quantile(Tbp_arm1[t,],0.975,na.rm=TRUE)
  low_Tb[t,1] = quantile(Tbp_arm1[t,],0.025,na.rm=TRUE)
  med_Tb[t,1] = quantile(Tbp_arm1[t,],0.5,na.rm=TRUE)
  
  upp_Tb[t,2] = quantile(Tbp_arm2[t,],0.975,na.rm=TRUE)
  low_Tb[t,2] = quantile(Tbp_arm2[t,],0.025,na.rm=TRUE)
  med_Tb[t,2] = quantile(Tbp_arm2[t,],0.5,na.rm=TRUE)
  
  upp_Tb[t,3] = quantile(Tbp_arm3[t,],0.975,na.rm=TRUE)
  low_Tb[t,3] = quantile(Tbp_arm3[t,],0.025,na.rm=TRUE)
  med_Tb[t,3] = quantile(Tbp_arm3[t,],0.5,na.rm=TRUE)
  
}
## Interventions
## 1 is [X1] LLIN targeted to pregnant women and children <6 years (TLLIN)
## 2 is [X2] LLIN universal coverage (ULLIN)
## 3 is [X3] TLLIN + full coverage of carbamate IRS (TLLIN + IRS)

cols_Tb = adegenet::transp(c("grey","darkred","purple"),0.4)
lty_Tb = c(1,1,4)
lwd_Tb = c(1,2,2)

for(i in 1:3){
  polygon(c(Tbp_arm1a[,1],rev(Tbp_arm1a[,1])),c(upp_Tb[,i],rev(low_Tb[,i])),col = cols_Tb[i],border=NA)
  lines(med_Tb[,i] ~ Tbp_arm1a[,1],lwd=lwd_Tb[i],lty=lty_Tb[i],col=cols_Tb[i])
}

Tbo_arm1 = c(0.947,0.834,0.853,0.817,0.675)
Tbo_arm2 = c(0.957,0.841,0.792,0.734,0.603)
Tbo_arm3 = c(0.884,0.762,0.801,0.730,0.637)

(1/12) * c(-3,2,5,8,11) ## months of testing prevalence Table 5 Curtis 1998
time_matchB = (1/12) * c(-3,2,5,8,11)#

points(Tbo_arm1 ~ time_matchB, col="grey",pch=4,cex=1.5)
points(Tbo_arm2 ~ c(time_matchB+0.01), col="red",pch=1,cex=1.5)
points(Tbo_arm3 ~ c(time_matchB+0.02), col="purple",pch=1,cex=1.5)

effect_TB = c((Tbo_arm1[2:5] - Tbo_arm2[2:5])/Tbo_arm1[2:5],
              (Tbo_arm1[2:5] - Tbo_arm3[2:5])/Tbo_arm1[2:5])

Tbp_arm1a$year[c(270,283,296,309)]

model_estimates_TB1 = c(med_Tb[270,1],med_Tb[283,1],med_Tb[296,1],med_Tb[309,1])
model_estimates_TB2 = c(med_Tb[270,2],med_Tb[283,2],med_Tb[296,2],med_Tb[309,2])
model_estimates_TB3 = c(med_Tb[270,3],med_Tb[283,3],med_Tb[296,3],med_Tb[309,3])

model_estimates_TB1u = c(upp_Tb[270,1],upp_Tb[283,1],upp_Tb[296,1],upp_Tb[309,1])
model_estimates_TB2u = c(upp_Tb[270,2],upp_Tb[283,2],upp_Tb[296,2],upp_Tb[309,2])
model_estimates_TB3u = c(upp_Tb[270,3],upp_Tb[283,3],upp_Tb[296,3],upp_Tb[309,3])

model_estimates_TB1l = c(low_Tb[270,1],low_Tb[283,1],low_Tb[296,1],low_Tb[309,1])
model_estimates_TB2l = c(low_Tb[270,2],low_Tb[283,2],low_Tb[296,2],low_Tb[309,2])
model_estimates_TB3l = c(low_Tb[270,3],low_Tb[283,3],low_Tb[296,3],low_Tb[309,3])



model_Effect_TB = c((model_estimates_TB1 - model_estimates_TB2)/model_estimates_TB1,
                    (model_estimates_TB1 - model_estimates_TB3)/model_estimates_TB1)

model_Effect_TB_min = c((model_estimates_TB1 - model_estimates_TB2u)/model_estimates_TB1,
                        (model_estimates_TB1 - model_estimates_TB3u)/model_estimates_TB1)

model_Effect_TB_max = c((model_estimates_TB1 - model_estimates_TB2l)/model_estimates_TB1,
                        (model_estimates_TB1 - model_estimates_TB3l)/model_estimates_TB1)


################################################################
## C 
## ID 5 D'Alessandro et al. 1992
## Data TRIAL



run=1:1000
site_num = rep(1:5,each=1000)
Tcp_arm1a = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/dalessandro/arm1/VALIDATION_trial_arm_1_1_net_1_0_1.txt"),header=TRUE)
Tcp_arm2a = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/dalessandro/arm2/VALIDATION_trial_arm_2_1_net_1_0_1.txt"),header=TRUE)
Tcp_arm1_1 = Tcp_arm1_2 = Tcp_arm1_3 = Tcp_arm1_4 = Tcp_arm1_5 = array(dim=c(nrow(Tcp_arm1a),1000)) ## this will be 1000 for each 5 sites once all the runs are done
for(i in 1:ncol(Tcp_arm1_1)){
  
  Tcp_arm1_1[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/dalessandro/arm1/VALIDATION_trial_arm_1_1_net_1_0_",run[i],".txt"),header=TRUE)$prev_1_4
  Tcp_arm1_2[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/dalessandro/arm1/VALIDATION_trial_arm_1_2_net_1_0_",run[i],".txt"),header=TRUE)$prev_1_4
  Tcp_arm1_3[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/dalessandro/arm1/VALIDATION_trial_arm_1_3_net_1_0_",run[i],".txt"),header=TRUE)$prev_1_4
  Tcp_arm1_4[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/dalessandro/arm1/VALIDATION_trial_arm_1_4_net_1_0_",run[i],".txt"),header=TRUE)$prev_1_4
  Tcp_arm1_5[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/dalessandro/arm1/VALIDATION_trial_arm_1_5_net_1_0_",run[i],".txt"),header=TRUE)$prev_1_4
}
Tcp_arm2_1 = Tcp_arm2_2 = Tcp_arm2_3 = Tcp_arm2_4 = Tcp_arm2_5 = array(dim=c(nrow(Tcp_arm1a),1000)) ## this will be 1000 for each 5 sites once all the runs are done
for(i in 1:ncol(Tcp_arm2_1)){
  Tcp_arm2_1[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/dalessandro/arm2/VALIDATION_trial_arm_2_1_net_1_0_",run[i],".txt"),header=TRUE)$prev_1_4
  Tcp_arm2_2[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/dalessandro/arm2/VALIDATION_trial_arm_2_2_net_1_0_",run[i],".txt"),header=TRUE)$prev_1_4
  Tcp_arm2_3[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/dalessandro/arm2/VALIDATION_trial_arm_2_3_net_1_0_",run[i],".txt"),header=TRUE)$prev_1_4
  Tcp_arm2_4[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/dalessandro/arm2/VALIDATION_trial_arm_2_4_net_1_0_",run[i],".txt"),header=TRUE)$prev_1_4
  Tcp_arm2_5[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/dalessandro/arm2/VALIDATION_trial_arm_2_5_net_1_0_",run[i],".txt"),header=TRUE)$prev_1_4
}

plot(Tcp_arm1_1[,1] ~ Tcp_arm1a[,1],ylim = c(0,1),xlim=c(-1,2),ylab="Prevalence in children 1 to 4 years (%)",
     xlab = "Time in years",pch="",yaxt="n",main="D'Alessandro",xaxt="n")
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
axis(1,at=seq(-1,2,1),labels=c("Jul 1991","Jul 1992","Jul 1993","Jul 1994"))

abline(v=0.083331,lty=2)
arrows(x0 = 0, y0 = 0.05, x1 = 0, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)



upp_Tc = low_Tc = med_Tc = array(dim=c(1301,2,5))
for(t in 1:1301){
  upp_Tc[t,1,1] = quantile(Tcp_arm1_1[t,],0.975,na.rm=TRUE)
  low_Tc[t,1,1] = quantile(Tcp_arm1_1[t,],0.025,na.rm=TRUE)
  med_Tc[t,1,1] = quantile(Tcp_arm1_1[t,],0.5,na.rm=TRUE)
  
  upp_Tc[t,1,2] = quantile(Tcp_arm1_2[t,],0.975,na.rm=TRUE)
  low_Tc[t,1,2] = quantile(Tcp_arm1_2[t,],0.025,na.rm=TRUE)
  med_Tc[t,1,2] = quantile(Tcp_arm1_2[t,],0.5,na.rm=TRUE)
  
  upp_Tc[t,1,3] = quantile(Tcp_arm1_3[t,],0.975,na.rm=TRUE)
  low_Tc[t,1,3] = quantile(Tcp_arm1_3[t,],0.025,na.rm=TRUE)
  med_Tc[t,1,3] = quantile(Tcp_arm1_3[t,],0.5,na.rm=TRUE)
  
  upp_Tc[t,1,4] = quantile(Tcp_arm1_4[t,],0.975,na.rm=TRUE)
  low_Tc[t,1,4] = quantile(Tcp_arm1_4[t,],0.025,na.rm=TRUE)
  med_Tc[t,1,4] = quantile(Tcp_arm1_4[t,],0.5,na.rm=TRUE)
  
  upp_Tc[t,1,5] = quantile(Tcp_arm1_5[t,],0.975,na.rm=TRUE)
  low_Tc[t,1,5] = quantile(Tcp_arm1_5[t,],0.025,na.rm=TRUE)
  med_Tc[t,1,5] = quantile(Tcp_arm1_5[t,],0.5,na.rm=TRUE)
  
  upp_Tc[t,2,1] = quantile(Tcp_arm2_1[t,],0.975,na.rm=TRUE)
  low_Tc[t,2,1] = quantile(Tcp_arm2_1[t,],0.025,na.rm=TRUE)
  med_Tc[t,2,1] = quantile(Tcp_arm2_1[t,],0.5,na.rm=TRUE)
  
  upp_Tc[t,2,2] = quantile(Tcp_arm2_2[t,],0.975,na.rm=TRUE)
  low_Tc[t,2,2] = quantile(Tcp_arm2_2[t,],0.025,na.rm=TRUE)
  med_Tc[t,2,2] = quantile(Tcp_arm2_2[t,],0.5,na.rm=TRUE)
  
  upp_Tc[t,2,3] = quantile(Tcp_arm2_3[t,],0.975,na.rm=TRUE)
  low_Tc[t,2,3] = quantile(Tcp_arm2_3[t,],0.025,na.rm=TRUE)
  med_Tc[t,2,3] = quantile(Tcp_arm2_3[t,],0.5,na.rm=TRUE)
  
  upp_Tc[t,2,4] = quantile(Tcp_arm2_4[t,],0.975,na.rm=TRUE)
  low_Tc[t,2,4] = quantile(Tcp_arm2_4[t,],0.025,na.rm=TRUE)
  med_Tc[t,2,4] = quantile(Tcp_arm2_4[t,],0.5,na.rm=TRUE)
  
  upp_Tc[t,2,5] = quantile(Tcp_arm2_5[t,],0.975,na.rm=TRUE)
  low_Tc[t,2,5] = quantile(Tcp_arm2_5[t,],0.025,na.rm=TRUE)
  med_Tc[t,2,5] = quantile(Tcp_arm2_5[t,],0.5,na.rm=TRUE)
  
}
cols_Tc = adegenet::transp(c("grey","darkred"),0.4)
lty_Tc = c(1,1,1,1,1)
lwd_Tc = c(1,1,1,1)

cols2_Tc = c("grey","darkred")


for(j in 1:5){
  for(i in 1:2){
    polygon(c(Tcp_arm1a[,1],rev(Tcp_arm1a[,1])),c(upp_Tc[,i,j],rev(low_Tc[,i,j])),col = cols_Tc[i],border=NA)
    lines(med_Tc[,i,j] ~ Tcp_arm1a[,1,1],col=cols2_Tc[i],lty=1,lwd=i)
  }  
}

Tco_arm1 = c(0.395,	0.287,	0.342,	0.573,	0.712)
Tco_arm2 = c(0.395,	0.287,	0.342,	0.573,	0.712)

Tco_arm1b = c(0.367,	0.330,	0.258,	0.533, 0.447)
Tco_arm2b = c(0.282,	0.228,	0.159,	0.431,	0.710)

(1/12) * c(-0.416667,12) ## months of testing prevalence in the trial (from Mark Rowland presentation, email 19/12/2018)
time_matchC = rep(c(-0.416667,1.583333),each=5)##because the surveys were done 2 months into the year

points(Tco_arm1 ~ time_matchC[1:5], col="grey20",pch=4,cex=1.5)
points(Tco_arm1b ~ time_matchC[6:10], col="grey20",pch=4,cex=1.5)
points(Tco_arm2b ~ c(time_matchC[6:10]+0.02), col="red",pch=1,cex=1.5)

effect_TC = c((Tco_arm2 - Tco_arm2b)/Tco_arm2)
effect_TC = c((Tco_arm1b - Tco_arm2b)/Tco_arm1b)

Tcp_arm1a$year[343]

model_estimates_TC1 = c(med_Tc[343,1,1],med_Tc[343,1,2],med_Tc[343,1,3],med_Tc[343,1,4],med_Tc[343,1,5])
model_estimates_TC1 = c(med_Tc[343,2,1],med_Tc[343,2,2],med_Tc[343,2,3],med_Tc[343,2,4],med_Tc[343,2,5])

model_Effect2_TC = c((med_Tc[343,1,1] - med_Tc[343,2,1])/med_Tc[343,1,1],
                     (med_Tc[343,1,2] - med_Tc[343,2,2])/med_Tc[343,1,2],
                     (med_Tc[343,1,3] - med_Tc[343,2,3])/med_Tc[343,1,3],
                     (med_Tc[343,1,4] - med_Tc[343,2,4])/med_Tc[343,1,4],
                     (med_Tc[343,1,5] - med_Tc[343,2,5])/med_Tc[343,1,5])

model_Effect2_TC_max = c((med_Tc[343,1,1] - upp_Tc[343,2,1])/med_Tc[343,1,1],
                         (med_Tc[343,1,2] - upp_Tc[343,2,2])/med_Tc[343,1,2],
                         (med_Tc[343,1,3] - upp_Tc[343,2,3])/med_Tc[343,1,3],
                         (med_Tc[343,1,4] - upp_Tc[343,2,4])/med_Tc[343,1,4],
                         (med_Tc[343,1,5] - upp_Tc[343,2,5])/med_Tc[343,1,5])

model_Effect2_TC_min = c((med_Tc[343,1,1] - low_Tc[343,2,1])/med_Tc[343,1,1],
                         (med_Tc[343,1,2] - low_Tc[343,2,2])/med_Tc[343,1,2],
                         (med_Tc[343,1,3] - low_Tc[343,2,3])/med_Tc[343,1,3],
                         (med_Tc[343,1,4] - low_Tc[343,2,4])/med_Tc[343,1,4],
                         (med_Tc[343,1,5] - low_Tc[343,2,5])/med_Tc[343,1,5])

legend("topright",legend=c("No intervention","Pyrethroid ITN"),
       pch=15,lty=1,col=c("grey","darkred"),bty="n")


## D NOT APPROPRIATE STUDY

#########################################
##
## E
## ID 6 Henry et al 2005
##
output1 = read.table("P:/validating_modelling/output_files/FINAL_netparm1/henry/arm1/VALIDATION_trial_arm_1_net_1_0_1.txt",header=TRUE)
Tep_arm1 = clin_inc_U5_arm1 = array(dim=c(nrow(output1),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Tep_arm1)){
  Tep_arm1[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/henry/arm1/VALIDATION_trial_arm_1_net_1_0_",run[i],".txt"),header=TRUE)$prev_0_5
}
Tep_arm2 = array(dim=c(nrow(output1),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Tep_arm2)){
  Tep_arm2[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/henry/arm2/VALIDATION_trial_arm_2_net_1_0_",run[i],".txt"),header=TRUE)$prev_0_5
}

plot(Tep_arm1[,1] ~ output1[,1],ylim = c(0,1),xlim=c(-1,2),ylab="Prevalence in children Under 5 years (%)",
     xlab = "Time in years",pch="",yaxt="n",main="Henry",xaxt="n")
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
axis(1,at=c(-1,0,1,2),labels=c("Jun 1998","Jun 1999","Jun 2000","Jun 2001"))

abline(v=0,lty=2)
arrows(x0 = 0.083333, y0 = 0.05, x1 = 0.083333, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)


upp_Te = low_Te = med_Te = array(dim=c(1301,2))
for(t in 1:1301){
  upp_Te[t,1] = quantile(Tep_arm1[t,],0.975,na.rm=TRUE)
  low_Te[t,1] = quantile(Tep_arm1[t,],0.025,na.rm=TRUE)
  med_Te[t,1] = quantile(Tep_arm1[t,],0.5,na.rm=TRUE)
  
  upp_Te[t,2] = quantile(Tep_arm2[t,],0.975,na.rm=TRUE)
  low_Te[t,2] = quantile(Tep_arm2[t,],0.025,na.rm=TRUE)
  med_Te[t,2] = quantile(Tep_arm2[t,],0.5,na.rm=TRUE)
  
}
cols_Te = adegenet::transp(c("grey","darkred"),0.4)
lwd_Te = c(1,2)
lty_Te = c(1,1)
for(i in 1:2){
  polygon(c(output1[,1],rev(output1[,1])),c(upp_Te[,i],rev(low_Te[,i])),col = cols_Te[i],border=NA)
  lines(med_Te[,i] ~ output1[,1],col=cols_Te[i],lwd = lwd_Te[i],lty=lty_Te[i])
}


Teo_arm1 = c(0.84,0.685)
Teo_arm2 = c(0.80,0.566)

Teo_arm1U = c(0.903,0.721)
Teo_arm2U = c(0.872,0.602)

Teo_arm1L = c(0.77,0.649)
Teo_arm2L = c(0.728,0.53)

effect_TE  = (Teo_arm1 - Teo_arm2)/Teo_arm1

effect_TE_min = (Teo_arm1L - Teo_arm2)/Teo_arm1
effect_TE_max = (Teo_arm1U - Teo_arm2)/Teo_arm1


time_matchE = c(-0.083333,0.91667)#
cols_Te2 = c("grey20","red")
points(Teo_arm1 ~ time_matchE, col="grey",pch=4,cex=1.5)
points(Teo_arm2 ~ c(time_matchE+0.01), col="red",pch=1,cex=1.5)
segments(x0=time_matchE[1],x1=time_matchE[1],y0=Teo_arm1U[1],Teo_arm1L[1],col=cols_Te2[1]) 
segments(x0=time_matchE[2],x1=time_matchE[2],y0=Teo_arm1U[2],Teo_arm1L[2],col=cols_Te2[1]) 

segments(x0=c(time_matchE[1]+0.01),x1=c(time_matchE[1]+0.01),y0=Teo_arm2U[1],Teo_arm2L[1],col=cols_Te2[2]) 
segments(x0=c(time_matchE[2]+0.01),x1=c(time_matchE[2]+0.01),y0=Teo_arm2U[2],Teo_arm2L[2],col=cols_Te2[2]) 

output1$year[c(257,309)]

model_estimates_TE1 = c(med_Te[257,1],med_Te[309,1])
model_estimates_TE2 = c(med_Te[257,2],med_Te[309,2])
model_Effect_TE = (model_estimates_TE1 - model_estimates_TE2)/model_estimates_TE1

model_estimates_TE1_min = c(low_Te[257,1],low_Te[309,1])
model_estimates_TE2_min = c(low_Te[257,2],low_Te[309,2])

model_estimates_TE1_max = c(upp_Te[257,1],upp_Te[309,1])
model_estimates_TE2_max = c(upp_Te[257,2],upp_Te[309,2])

model_Effect_TE_min = (model_estimates_TE1 - model_estimates_TE2_min)/model_estimates_TE1
model_Effect_TE_max = (model_estimates_TE1 - model_estimates_TE2_max)/model_estimates_TE1



################################################################
## F 
## ID 7 Kafy et al. 2017
## Data TRIAL
Tfo_arm1 = c(0.07,0.05,0.05)
Tfo_arm1_min = c(0.14,0.10,0.09)
Tfo_arm1_max = c(0.03,0.02,0.03)

Tfo_arm2 = c(0.1,0.04,0.03)
Tfo_arm2_min = c(0.16,0.07,0.05)
Tfo_arm2_max = c(0.06,0.02,0.02)

time_match_F1 = (1/12) * c(-1,12,24) ##march 2011,apr 2012, 2013
time_match_F2 = (1/12) * c(-1,12,24) + 1/24

effect_TF = (Tfo_arm1 - Tfo_arm2)/Tfo_arm1

effect_TF_min = (Tfo_arm1_min - Tfo_arm2_min)/Tfo_arm1
effect_TF_max = (Tfo_arm1_max - Tfo_arm2_max)/Tfo_arm1

## Data MODEL
Tfp_arm1a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/kafy/arm1/VALIDATION_trial_arm_1_net_1_0_1.txt",header=TRUE)
Tfp_arm1 = array(dim=c(nrow(Tfp_arm1a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:1000){
  Tfp_arm1[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/kafy/arm1/VALIDATION_trial_arm_1_net_1_0_",run[i],".txt"),header=TRUE)$prev_4_6
}

Tfp_arm2a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/kafy/arm2/VALIDATION_trial_arm_2_net_1_0_1.txt",header=TRUE)
Tfp_arm2 = array(dim=c(nrow(Tfp_arm2a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:1000){
  Tfp_arm2[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/kafy/arm2/VALIDATION_trial_arm_2_net_1_0_",run[i],".txt"),header=TRUE)$prev_4_6
}

plot(Tfp_arm1[,1] ~ Tfp_arm1a[,1],ylim = c(0,0.4),xlim=c(-1,3),ylab="Prevalence in children 4 - 6 years (%)",
     xlab = "Time in years",pch="",yaxt="n",xaxt = "n",main="Kafy")
axis(1,at=c(-1,0,1,2,3),labels = c("Apr 2010","Apr 2011","Apr 2012","Apr 2013",""))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))

abline(v=0,lty=2)
arrows(x0 = 0, y0 = 0.05, x1 = 0, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)
arrows(x0 = 0.3333, y0 = 0.05, x1 = 0.333, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)

legend("topright",legend = c("Pyrethroid ITN",
                             "Pyrethroid ITN + Pyrethroid then Carbamate IRS"),
       col = c("darkred","darkgreen"), pch = 15,bty="n")

upp_Tf = low_Tf = med_Tf = array(dim=c(1301,2))
for(t in 1:1301){
  upp_Tf[t,1] = quantile(Tfp_arm1[t,],0.975,na.rm=TRUE)
  low_Tf[t,1] = quantile(Tfp_arm1[t,],0.025,na.rm=TRUE)
  med_Tf[t,1] = quantile(Tfp_arm1[t,],0.5,na.rm=TRUE)
  
  upp_Tf[t,2] = quantile(Tfp_arm2[t,],0.975,na.rm=TRUE)
  low_Tf[t,2] = quantile(Tfp_arm2[t,],0.025,na.rm=TRUE)
  med_Tf[t,2] = quantile(Tfp_arm2[t,],0.5,na.rm=TRUE)
  
}
cols_Tf = adegenet::transp(c("darkred","darkgreen"),0.4)
cols_Tf2 = c("darkred","darkgreen")
lwd_Tf = c(2,2)
lty_Tf = c(1,4)
for(i in 1:2){
  polygon(c(Tfp_arm1a[,1],rev(Tfp_arm1a[,1])),c(upp_Tf[,i],rev(low_Tf[,i])),col = cols_Tf[i],border=NA)
  lines(med_Tf[,i] ~ Tfp_arm1a[,1], col=cols_Tf2[i],lty=lty_Tf[i],lwd=lwd_Tf[i])
}

time_match_F1
Tfp_arm1a$year[c(257,313,365)]

model_estimates_TF1 = c(med_Tf[257,1],med_Tf[313,1],med_Tf[365,1])
model_estimates_TF2 = c(med_Tf[257,2],med_Tf[313,2],med_Tf[365,2])
model_Effect_TF = (model_estimates_TF1 - model_estimates_TF2)/model_estimates_TF1

model_estimates_TF1_min = c(low_Tf[257,1],low_Tf[313,1],low_Tf[365,1])
model_estimates_TF2_min = c(low_Tf[257,2],low_Tf[313,2],low_Tf[365,2])

model_estimates_TF1_max = c(upp_Tf[257,1],upp_Tf[313,1],upp_Tf[365,1])
model_estimates_TF2_max = c(upp_Tf[257,2],upp_Tf[313,2],upp_Tf[365,2])

model_Effect_TF_min = (model_estimates_TF1 - model_estimates_TF2_min)/model_estimates_TF1
model_Effect_TF_max = (model_estimates_TF1 - model_estimates_TF2_max)/model_estimates_TF1

points(Tfo_arm1 ~ time_match_F1,col="red",pch=19,cex=1.5)
points(Tfo_arm2 ~ time_match_F2,col="darkgreen",pch=19,cex=1.5)
segments(x0=time_match_F1,x1=time_match_F1,y0=Tfo_arm1_min,y1=Tfo_arm1_max,lty=1,lwd=1.5,col="red")
segments(x0=time_match_F2,x1=time_match_F2,y0=Tfo_arm2_min,y1=Tfo_arm2_max,lty=1,lwd=1.5,col="darkgreen")


##mass distribution in aug 2012 == 261+35weeks
##mass distribution in dec 2012 == 261+50weeks
##mass distribution in aug 2013 == 261+35+52weeks
##mass distribution in dec 2013 == 261+50+52weeks
##mass distribution in aug 2013 == 261+35+104weeks
##mass distribution in dec 2013 == 261+50+104weeks



################################################################
## H 
## ID 9 MARBIAH 1992 - 1994

#  ## Bo - Southern Province, Sierre Leone

Thp_arm1a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/marbiah/arm1/VALIDATION_trial_arm_1_net_1_0_1.txt",header=TRUE)
Thp_arm1 = array(dim=c(nrow(Thp_arm1a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Thp_arm1)){
  Thp_arm1[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/marbiah/arm1/VALIDATION_trial_arm_1_net_1_0_",run[i],".txt"),header=TRUE)$prev_0_7
} ## Non-interventions
Thp_arm2a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/marbiah/arm2/VALIDATION_trial_arm_2_net_1_0_1.txt",header=TRUE)
Thp_arm2 = array(dim=c(nrow(Thp_arm2a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Thp_arm2)){
  Thp_arm2[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/marbiah/arm2/VALIDATION_trial_arm_2_net_1_0_",run[i],".txt"),header=TRUE)$prev_0_7
} ## Bo - Southern Province, Sierre Leone

plot(Thp_arm1[,1] ~ Thp_arm1a[,1],ylim = c(0,1),xlim=c(-1,2),ylab="Prevalence in children Under 7 years (%)",
     xlab = "Time in years",pch="",yaxt="n",main="Marbiah",xaxt="n")
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
axis(1,at=seq(-1,2,1),labels=c("Jun 1991","Jun 1992","Jun 1993","Jun 1994"))


abline(v=0,lty=2)
arrows(x0 = 0, y0 = 0.05, x1 = 0, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)

legend("topright",legend = c("No intervention",
                             "Pyrethroid ITN"),
       col = c("grey","darkred"), pch = 15,bty="n")
upp_Th = low_Th = med_Th = array(dim=c(1301,2))
for(t in 1:1301){
  upp_Th[t,1] = quantile(Thp_arm1[t,],0.975,na.rm=TRUE)
  low_Th[t,1] = quantile(Thp_arm1[t,],0.025,na.rm=TRUE)
  med_Th[t,1] = quantile(Thp_arm1[t,],0.5,na.rm=TRUE)
  
  upp_Th[t,2] = quantile(Thp_arm2[t,],0.975,na.rm=TRUE)
  low_Th[t,2] = quantile(Thp_arm2[t,],0.025,na.rm=TRUE)
  med_Th[t,2] = quantile(Thp_arm2[t,],0.5,na.rm=TRUE)
}
cols_Th = adegenet::transp(c("grey","darkred"),0.4)
col_Th2 = c("grey","darkred")
lwd_Th=c(1,2)
lty_Th=c(1,1)
for(i in 1:2){
  polygon(c(Thp_arm1a[,1],rev(Thp_arm1a[,1])),c(upp_Th[,i],rev(low_Th[,i])),col = cols_Th[i],border=NA)
  lines(med_Th[,i] ~ Thp_arm1a[,1],lwd=lwd_Th[i],lty=lty_Th[i],col=col_Th2[i])
}

Tho_arm1 = c(0.61,0.476)
Tho_arm2 = c(0.61,0.345)

effect_TH = c((Tho_arm1[2] - Tho_arm2[2])/Tho_arm1[2])

time_matchH = (1/12) * c(-4.5,9)

points(Tho_arm1 ~ time_matchH, col="grey",pch=4,cex=1.5)
points(Tho_arm2 ~ c(time_matchH+0.01), col="red",pch=17,cex=1.5)

Thp_arm1a$year[300]
model_estimates_TH1 = mean(Thp_arm1[300,])
model_estimates_TH2 = mean(Thp_arm2[300,])

model_Effect_TH = c((model_estimates_TH1 - model_estimates_TH2)/model_estimates_TH1)

model_estimates_TH1_min = low_Th[300,1]
model_estimates_TH1_max = upp_Th[300,1]
model_estimates_TH2_min = low_Th[300,2]
model_estimates_TH2_max = upp_Th[300,2]

model_Effect_TH_min =  (model_estimates_TH1 - model_estimates_TH2_min)/model_estimates_TH1
model_Effect_TH_max =  (model_estimates_TH1 - model_estimates_TH2_max)/model_estimates_TH1


#############################
##
## S 
## ID 10 Nevill et al 1996


Tsp_arm1a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/nevill/arm1/VALIDATION_trial_arm_1_net_1_0_1.txt",header=TRUE)
Tsp_arm1 = array(dim=c(nrow(Tsp_arm1a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Tsp_arm1)){
  Tsp_arm1[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/nevill/arm1/VALIDATION_trial_arm_1_net_1_0_",run[i],".txt"),header=TRUE)$prev_0_1
}

Tsp_arm2a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/nevill/arm2/VALIDATION_trial_arm_2_net_1_0_1.txt",header=TRUE)
Tsp_arm2 = array(dim=c(nrow(Tsp_arm2a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Tsp_arm2)){
  Tsp_arm2[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/nevill/arm2/VALIDATION_trial_arm_2_net_1_0_",run[i],".txt"),header=TRUE)$prev_0_1
}


plot(Tsp_arm1[,1] ~ Tsp_arm1a[,1],ylim = c(0,0.4),xlim=c(-1,3),ylab="Point slide-prevalence under 12-months (%)",
     xlab = "Time in years",pch="",yaxt="n",xaxt="n",main="Nevill")
axis(1,at = -1:2,labels=c("Jul 1992","Jul 1993","Jul 1994","Jul 1995"))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))


abline(v=0,lty=2)
arrows(x0 = 0, y0 = 0.05, x1 = 0, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)
arrows(x0 = 0.5, y0 = 0.05, x1 = 0.5, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)
arrows(x0 = 1, y0 = 0.05, x1 = 1, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)
arrows(x0 = 1.5, y0 = 0.05, x1 = 1.5, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)

legend("topright",legend = c("No intervention",
                             "Pyrethroid ITN"),
       col = c("grey","darkred"), pch = 15,bty="n")

upp_Ts = low_Ts = med_Ts = array(dim=c(1301,2))
for(t in 1:1301){
  upp_Ts[t,1] = quantile(Tsp_arm1[t,],0.975,na.rm=TRUE)
  low_Ts[t,1] = quantile(Tsp_arm1[t,],0.025,na.rm=TRUE)
  med_Ts[t,1] = quantile(Tsp_arm1[t,],0.5,na.rm=TRUE)
  
  upp_Ts[t,2] = quantile(Tsp_arm2[t,],0.975,na.rm=TRUE)
  low_Ts[t,2] = quantile(Tsp_arm2[t,],0.025,na.rm=TRUE)
  med_Ts[t,2] = quantile(Tsp_arm2[t,],0.5,na.rm=TRUE)
  
}
cols_Ts = adegenet::transp(c("grey","darkred"),0.4)
cols_Ts2 = c("grey","darkred")
lwd_Ts = c(1,2)
lty_Ts = c(1,1)

for(i in 1:2){
  polygon(c(Tsp_arm1a[,1],rev(Tsp_arm1a[,1])),c(upp_Ts[,i],rev(low_Ts[,i])),col = cols_Ts[i],border=NA)
  lines(med_Ts[,i] ~ Tsp_arm1a[,1],col=cols_Ts2[i],lwd=lwd_Ts[i],lty=lty_Ts[i])
}


# UNDER 1s
##from Snow et al. 1996 Table 1
Tsp_prev_recorded_t1 = c(0.392,0.254,0.115)
Tsp_prev_recorded_t2 = c(0.168,0.12,0.078)

Tsp_time_match = c(20/12,24/12,27/12)

effect_TS = c((Tsp_prev_recorded_t1[1] - Tsp_prev_recorded_t2[1])/Tsp_prev_recorded_t1[1],
              (Tsp_prev_recorded_t1[2] - Tsp_prev_recorded_t2[2])/Tsp_prev_recorded_t1[2],
              (Tsp_prev_recorded_t1[3] - Tsp_prev_recorded_t2[3])/Tsp_prev_recorded_t1[3])

points(Tsp_prev_recorded_t1 ~ Tsp_time_match, col="grey",pch=1,cex=1.5)
points(Tsp_prev_recorded_t2 ~ c(Tsp_time_match+0.04), col="red",pch=1,cex=1.5)

nevil_base = c(0.251,0.251)
nevil_time_base = c(-3/12,-3.05/12)
points(nevil_base ~ nevil_time_base,col=c("grey","red"),pch=19)

#1.666667 2.000000 2.250000
261+Tsp_time_match*52
## 347.6667 365.0000 378.0000
model_estimates_TS1 = c(med_Tp[348,1],med_Tp[365,1],med_Tp[378,1])
model_estimates_TS2 = c(med_Tp[348,2],med_Tp[365,2],med_Tp[378,2])

model_Effect_TS2 =  (model_estimates_TS1 - model_estimates_TS2)/model_estimates_TS1

model_estimates_TS1_min = c(low_Tp[348,1],low_Tp[365,1],low_Tp[378,1])
model_estimates_TS2_min = c(low_Tp[348,2],low_Tp[365,2],low_Tp[378,2])
model_estimates_TS1_max = c(upp_Tp[348,1],upp_Tp[365,1],upp_Tp[378,1])
model_estimates_TS2_max = c(upp_Tp[348,2],upp_Tp[365,2],upp_Tp[378,2])

model_Effect_TS_min =  (model_estimates_TS1 - model_estimates_TS2_min)/model_estimates_TS1
model_Effect_TS_max =  (model_estimates_TS1 - model_estimates_TS2_max)/model_estimates_TS1



################################################################
## J 
## ID 11 Philips Howard et al 2003

## Data TRIAL
Tjp_arm1a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/philipshoward/arm1/VALIDATION_trial_arm_1_1_net_1_0_1.txt",header=TRUE)
Tjp_arm1 = array(dim=c(nrow(Tjp_arm1a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Tjp_arm1)){
  Tjp_arm1[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/philipshoward/arm1/VALIDATION_trial_arm_1_1_net_1_0_",run[i],".txt"),header=TRUE)$prev_0_5
} ## Non-interventions
Tjp_arm2a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/philipshoward/arm2/VALIDATION_trial_arm_2_1_net_1_0_1.txt",header=TRUE)
Tjp_arm2 = array(dim=c(nrow(Tjp_arm2a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Tjp_arm2)){
  Tjp_arm2[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/philipshoward/arm2/VALIDATION_trial_arm_2_1_net_1_0_",run[i],".txt"),header=TRUE)$prev_0_5
} ## Asembo
Tjp_arm3a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/philipshoward/arm3/VALIDATION_trial_arm_3_2_net_1_0_1.txt",header=TRUE)
Tjp_arm3 = array(dim=c(nrow(Tjp_arm3a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Tjp_arm3)){
  Tjp_arm3[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/philipshoward/arm3/VALIDATION_trial_arm_3_2_net_1_0_",run[i],".txt"),header=TRUE)$prev_0_5
} ## Gem

plot(Tjp_arm1[,1] ~ Tjp_arm1a[,1],ylim = c(0,1),xlim=c(-1,4),ylab="Prevalence in children Under 5 years (%)",
     xlab = "Time in years",pch="",yaxt="n",main="Philips-Howard",xaxt="n")
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
axis(1,at=seq(-1,4,1),labels=c("Jan 1996","Jan 1997","Jan 1998","Jan 1999","Jan 2000","Jan 2001"))

abline(v=-0.2,lty=2)
arrows(x0 = 0, y0 = 0.05, x1 = 0, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)

legend("topright",legend = c("Do nothing then Pyrethroid ITN",
                             "Pyrethroid ITN (Asembo)",
                             "Pyrethroid ITN (Gem)"),
       col = c("grey","darkred","darkred"), pch = 15,bty="n")

upp_Tj = low_Tj = med_Tj = array(dim=c(1301,3))
for(t in 1:1301){
  upp_Tj[t,1] = quantile(Tjp_arm1[t,],0.975,na.rm=TRUE)
  low_Tj[t,1] = quantile(Tjp_arm1[t,],0.025,na.rm=TRUE)
  med_Tj[t,1] = quantile(Tjp_arm1[t,],0.5,na.rm=TRUE)
  
  upp_Tj[t,2] = quantile(Tjp_arm2[t,],0.975,na.rm=TRUE)
  low_Tj[t,2] = quantile(Tjp_arm2[t,],0.025,na.rm=TRUE)
  med_Tj[t,2] = quantile(Tjp_arm2[t,],0.5,na.rm=TRUE)
  
  upp_Tj[t,3] = quantile(Tjp_arm3[t,],0.975,na.rm=TRUE)
  low_Tj[t,3] = quantile(Tjp_arm3[t,],0.025,na.rm=TRUE)
  med_Tj[t,3] = quantile(Tjp_arm3[t,],0.5,na.rm=TRUE)
}
cols_Tj = adegenet::transp(c("grey","darkred","darkred"),0.4)
cols_Tj2 = c("grey","darkred","darkred")
lwd_Tj =c(1,2,2)
lty_Tj =c(1,1,1)
for(i in 1:3){
  polygon(c(Tjp_arm1a[,1],rev(Tjp_arm1a[,1])),c(upp_Tj[,i],rev(low_Tj[,i])),col = cols_Tj[i],border=NA)
  lines(med_Tj[,i] ~ Tjp_arm1a[,1],col=cols_Tj2[i],lwd=lwd_Tj[i],lty=lty_Tj[i])
}

Tjo_arm1 = c(0.75,0.514,0.391)
Tjo_arm2 = c(0.75,0.325,0.338)

effect_TJ = c((Tjo_arm1[2] - Tjo_arm2[2])/Tjo_arm1[2],
              (Tjo_arm1[3] - Tjo_arm2[3])/Tjo_arm1[3])

(1/12) * c(0,18) ## months of testing prevalence in the trial (from Mark Rowland presentation, email 19/12/2018)
time_matchJ = (1/12) * c(0,22,40)##because the surveys were done 2 months into the year

points(Tjo_arm1 ~ time_matchJ, col="grey",pch=19,cex=1.5)
points(Tjo_arm2 ~ c(time_matchJ+0.01), col="red",pch=19,cex=1.5)


model_estimates_TJ1 = c(mean(Tjp_arm1[356,]),mean(Tjp_arm1[434,]))
model_estimates_TJ2 = c(mean(Tjp_arm2[356,]),mean(Tjp_arm2[434,]))
model_estimates_TJ3 = c(mean(Tjp_arm3[356,]),mean(Tjp_arm3[434,]))

model_Effect_TJ = c((model_estimates_TJ1 - model_estimates_TJ3)/model_estimates_TJ1)

model_estimates_TJ1min = c(low_Tj[356,1],low_Tj[434,1])
model_estimates_TJ2min =c(low_Tj[356,2],low_Tj[434,2])
model_estimates_TJ3min = c(low_Tj[356,3],low_Tj[434,3])
model_estimates_TJ1max = c(upp_Tj[356,1],upp_Tj[434,1])
model_estimates_TJ2max =c(upp_Tj[356,2],upp_Tj[434,2])
model_estimates_TJ3max = c(upp_Tj[356,3],upp_Tj[434,3])

model_Effect_TJ_max = c((model_estimates_TJ1 - model_estimates_TJ3min)/model_estimates_TJ1)
model_Effect_TJ_min = c((model_estimates_TJ1 - model_estimates_TJ3max)/model_estimates_TJ1)



################################################################
## L
## ID 12 Protopopoff et al. 2018
## Data TRIAL
Tlo_arm1 = c(0.555, 0.553, 0.530, 0.68,0.827,0.648)
Tlo_arm2 = c(0.458, 0.311, 0.348, 0.459,0.6835,0.490)
Tlo_arm3 = c(0.386, 0.264, 0.399, 0.551,NA,NA)
Tlo_arm4= c(0.375, 0.286, 0.291, 0.437,NA,NA)
time_match_L1 = c(0.33,0.75,1.33,1.75,2.33,2.75) 

effect_TL2 = (Tlo_arm1 - Tlo_arm2)/Tlo_arm1
effect_TL3 = (Tlo_arm1 - Tlo_arm3)/Tlo_arm1
effect_TL4 = (Tlo_arm1 - Tlo_arm4)/Tlo_arm1

EFF_Pbo_irs = (Tlo_arm2 - Tlo_arm4)/Tlo_arm2

## Data MODEL
Tlp_arm1a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/protopopoff/arm1/VALIDATION_trial_arm_1_net_1_0_1.txt",header=TRUE)
Tlp_arm1 = array(dim=c(nrow(Tlp_arm1a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Tlp_arm1)){
  Tlp_arm1[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/protopopoff/arm1/VALIDATION_trial_arm_1_net_1_0_",run[i],".txt"),header=TRUE)$prev_0_14
}
Tlp_arm2a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/protopopoff/arm2/VALIDATION_trial_arm_2_net_2_0_1.txt",header=TRUE)
Tlp_arm2 = array(dim=c(nrow(Tlp_arm2a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Tlp_arm2)){
  Tlp_arm2[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/protopopoff/arm2/VALIDATION_trial_arm_2_net_2_0_",run[i],".txt"),header=TRUE)$prev_0_14
}
Tlp_arm3a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/protopopoff/arm3/VALIDATION_trial_arm_3_net_1_0_1.txt",header=TRUE)
Tlp_arm3 = array(dim=c(nrow(Tlp_arm3a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Tlp_arm3)){
  Tlp_arm3[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/protopopoff/arm3/VALIDATION_trial_arm_3_net_1_0_",run[i],".txt"),header=TRUE)$prev_0_14
}
Tlp_arm4a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/protopopoff/arm4/VALIDATION_trial_arm_4_net_2_0_1.txt",header=TRUE)
Tlp_arm4 = array(dim=c(nrow(Tlp_arm4a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Tlp_arm4)){
  Tlp_arm4[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/protopopoff/arm4/VALIDATION_trial_arm_4_net_2_0_",run[i],".txt"),header=TRUE)$prev_0_14
}

VALES = round(261+52*time_match_L1,0)
Tlp_arm1a$year[261+52*time_match_L1]

model_estimates_TL1 = c(mean(Tlp_arm1[VALES[1],],na.rm=TRUE),mean(Tlp_arm1[VALES[2],],na.rm=TRUE),mean(Tlp_arm1[VALES[3],],na.rm=TRUE),mean(Tlp_arm1[VALES[4],],na.rm=TRUE),mean(Tlp_arm1[VALES[5],],na.rm=TRUE),mean(Tlp_arm1[VALES[6],],na.rm=TRUE))
model_estimates_TL2 = c(mean(Tlp_arm2[VALES[1],]),mean(Tlp_arm2[VALES[2],]),mean(Tlp_arm2[330,]),mean(Tlp_arm2[VALES[4],]),mean(Tlp_arm2[VALES[5],]),mean(Tlp_arm2[VALES[6],]))
model_estimates_TL3 = c(mean(Tlp_arm3[VALES[1],]),mean(Tlp_arm3[VALES[2],]),mean(Tlp_arm3[330,]),mean(Tlp_arm3[VALES[4],]),mean(Tlp_arm3[VALES[5],]),mean(Tlp_arm3[VALES[6],]))
model_estimates_TL4 = c(mean(Tlp_arm4[VALES[1],]),mean(Tlp_arm4[VALES[2],]),mean(Tlp_arm4[330,]),mean(Tlp_arm4[VALES[4],]),mean(Tlp_arm4[VALES[5],]),mean(Tlp_arm4[VALES[6],]))

model_Effect_TL2 = (model_estimates_TL1 - model_estimates_TL2)/model_estimates_TL1
model_Effect_TL3 = (model_estimates_TL1 - model_estimates_TL3)/model_estimates_TL1
model_Effect_TL4 = (model_estimates_TL1 - model_estimates_TL4)/model_estimates_TL1


## Add in any uncertainty bars from the trial (need to acess these if there are any)
plot(Tlp_arm1[,1] ~ Tlp_arm1a[,1],ylim = c(0,1),xlim=c(-0.5,3),ylab="Prevalence in children 6 months to 14 years (%)",
     xlab = "Time in years",pch="",xaxt="n",yaxt="n",main="Protopopoff")
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
axis(1, at=c(0:3),labels=c("Jan 2015","Jan 2016","Jan 2017","Jan 2018"))

abline(v=-0.1666667,lty=2)
arrows(x0 = 0, y0 = 0.05, x1 = 0, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)

legend("bottomright",legend = c("Pyrethroid ITN",
                                "Pyrethroid + PBO ITN",
                                "Pyrethroid ITN + IRS",
                                "Pyrethroid + PBO ITN + IRS"),
       col = c("darkred","aquamarine3","gold2","blue"), pch = 15,bty="n")
upp_Tl = low_Tl = med_Tl = array(dim=c(1301,4))
for(t in 1:1301){
  upp_Tl[t,1] = quantile(Tlp_arm1[t,],0.975,na.rm=TRUE)
  low_Tl[t,1] = quantile(Tlp_arm1[t,],0.025,na.rm=TRUE)
  med_Tl[t,1] = quantile(Tlp_arm1[t,],0.5,na.rm=TRUE)
  
  upp_Tl[t,2] = quantile(Tlp_arm2[t,],0.975,na.rm=TRUE)
  low_Tl[t,2] = quantile(Tlp_arm2[t,],0.025,na.rm=TRUE)
  med_Tl[t,2] = quantile(Tlp_arm2[t,],0.5,na.rm=TRUE)
  
  upp_Tl[t,3] = quantile(Tlp_arm3[t,],0.975,na.rm=TRUE)
  low_Tl[t,3] = quantile(Tlp_arm3[t,],0.025,na.rm=TRUE)
  med_Tl[t,3] = quantile(Tlp_arm3[t,],0.5,na.rm=TRUE)
  
  upp_Tl[t,4] = quantile(Tlp_arm4[t,],0.975,na.rm=TRUE)
  low_Tl[t,4] = quantile(Tlp_arm4[t,],0.025,na.rm=TRUE)
  med_Tl[t,4] = quantile(Tlp_arm4[t,],0.5,na.rm=TRUE)
  
  
}

model_estimates_TL1_min = c(low_Tl[VALES[1],1],low_Tl[VALES[2],1],low_Tl[VALES[3],1],low_Tl[VALES[4],1],low_Tl[VALES[5],1],low_Tl[VALES[6],1])
model_estimates_TL2_min = c(low_Tl[VALES[1],2],low_Tl[VALES[2],2],low_Tl[VALES[3],2],low_Tl[VALES[4],2],low_Tl[VALES[5],2],low_Tl[VALES[6],2])
model_estimates_TL3_min = c(low_Tl[VALES[1],3],low_Tl[VALES[2],3],low_Tl[VALES[3],3],low_Tl[VALES[4],3],low_Tl[VALES[5],3],low_Tl[VALES[6],3])
model_estimates_TL4_min = c(low_Tl[VALES[1],4],low_Tl[VALES[2],4],low_Tl[VALES[3],4],low_Tl[VALES[4],4],low_Tl[VALES[5],4],low_Tl[VALES[6],4])

model_estimates_TL1_max = c(upp_Tl[VALES[1],1],upp_Tl[VALES[2],1],upp_Tl[VALES[3],1],upp_Tl[VALES[4],1],upp_Tl[VALES[5],1],upp_Tl[VALES[6],1])
model_estimates_TL2_max = c(upp_Tl[VALES[1],2],upp_Tl[VALES[2],2],upp_Tl[VALES[3],2],upp_Tl[VALES[4],2],upp_Tl[VALES[5],2],upp_Tl[VALES[6],2])
model_estimates_TL3_max = c(upp_Tl[VALES[1],3],upp_Tl[VALES[2],3],upp_Tl[VALES[3],3],upp_Tl[VALES[4],3],upp_Tl[VALES[5],3],upp_Tl[VALES[6],3])
model_estimates_TL4_max = c(upp_Tl[VALES[1],4],upp_Tl[VALES[2],4],upp_Tl[VALES[3],4],upp_Tl[VALES[4],4],upp_Tl[VALES[5],4],upp_Tl[VALES[6],4])

model_Effect_TL2_min = (model_estimates_TL1 - model_estimates_TL2_min)/model_estimates_TL1
model_Effect_TL3_min = (model_estimates_TL1 - model_estimates_TL3_min)/model_estimates_TL1
model_Effect_TL4_min = (model_estimates_TL1 - model_estimates_TL4_min)/model_estimates_TL1

model_Effect_TL2_max = (model_estimates_TL1 - model_estimates_TL2_max)/model_estimates_TL1
model_Effect_TL3_max = (model_estimates_TL1 - model_estimates_TL3_max)/model_estimates_TL1
model_Effect_TL4_max = (model_estimates_TL1 - model_estimates_TL4_max)/model_estimates_TL1


mod_pbo_vs_pbo_irs = (model_estimates_TL2 - model_estimates_TL4)/model_estimates_TL2
mod_pbo_vs_pbo_irsU = (model_estimates_TL2 - model_estimates_TL4_max)/model_estimates_TL2
mod_pbo_vs_pbo_irsL = (model_estimates_TL2 - model_estimates_TL4_min)/model_estimates_TL2


cols_Tl = adegenet::transp(c("darkred","aquamarine2","gold2","darkblue"),0.4)
cols_Tl2 = c("darkred","aquamarine2","gold2","darkblue")
lty_Tl = c(1,1,4,4)
lwd_Tl=c(2,2,2,2)
for(i in 1:4){
  polygon(c(Tlp_arm1a[,1],rev(Tlp_arm1a[,1])),c(upp_Tl[,i],rev(low_Tl[,i])),col = cols_Tl[i],border=NA)
  lines(med_Tl[,i] ~ Tlp_arm1a[,1],lty=lty_Tl[i],lwd=lwd_Tl[i],col = cols_Tl2[i])
}


points(Tlo_arm1 ~ time_match_L1,col="red",pch=19,cex=1.5)
points(Tlo_arm2 ~ c(time_match_L1+1/48),col="aquamarine3",pch=8,cex=1.5)
points(Tlo_arm3 ~ c(time_match_L1+2/48),col="gold2",pch=19,cex=1.5)
points(Tlo_arm4 ~ c(time_match_L1+3/48),col="darkblue",pch=8,cex=1.5)


## add baselines
points(c(0.68) ~ c(-0.1666667), col="red",pch=19,cex=1.5)
points(c(0.61) ~ c(-0.1666667), col="aquamarine3",pch=8,cex=1.5)
points(c(0.67) ~ c(-0.1666667), col="gold2",pch=19,cex=1.5)
points(c(0.64) ~ c(-0.1666667), col="darkblue",pch=8,cex=1.5)




#############################
##
## Q 
## ID 13 Staedke et al 2020
run = 1:1000

Tqo_arm1 = c(0.193,0.145,0.13,0.14)
Tqo_arm2 = c(0.192,0.107,0.106,0.118)
time_match_Q1 = c(4/12,12/12,18/12,24/12)
time_match_Q2 = time_match_Q1 +0.02

effect_TQ2 = (Tqo_arm1 - Tqo_arm2)/Tqo_arm1

## Data MODEL
Tqp_arm1a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/staedke/arm1/VALIDATION_trial_arm_1_net_1_0_1.txt",header=TRUE)
Tqp_arm1 = array(dim=c(nrow(Tqp_arm1a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Tqp_arm1)){
  Tqp_arm1[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/staedke/arm1/VALIDATION_trial_arm_1_net_1_0_",run[i],".txt"),header=TRUE)$prev_2_10
}

Tqp_arm2a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/staedke/arm2/VALIDATION_trial_arm_2_net_2_0_1.txt",header=TRUE)
Tqp_arm2 = array(dim=c(nrow(Tqp_arm2a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Tqp_arm2)){
  Tqp_arm2[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/staedke/arm2/VALIDATION_trial_arm_2_net_2_0_",run[i],".txt"),header=TRUE)$prev_2_10
}

upp_Tq = low_Tq = med_Tq = array(dim=c(1301,2))
for(t in 1:1301){
  upp_Tq[t,1] = quantile(Tqp_arm1[t,],0.975,na.rm=TRUE)
  low_Tq[t,1] = quantile(Tqp_arm1[t,],0.025,na.rm=TRUE)
  med_Tq[t,1] = quantile(Tqp_arm1[t,],0.5,na.rm=TRUE)
  
  upp_Tq[t,2] = quantile(Tqp_arm2[t,],0.975,na.rm=TRUE)
  low_Tq[t,2] = quantile(Tqp_arm2[t,],0.025,na.rm=TRUE)
  med_Tq[t,2] = quantile(Tqp_arm2[t,],0.5,na.rm=TRUE)
  
}
VALES2 = round(261+52*time_match_Q1,0)
Tqp_arm1a$year[261+52*time_match_Q1]

model_estimates_TQ1 = c(mean(Tqp_arm1[VALES2[1],]),
                        mean(Tqp_arm1[VALES2[2],]),
                        mean(Tqp_arm1[VALES2[3],]),
                        mean(Tqp_arm1[VALES2[4],]))
model_estimates_TQ2 = c(mean(Tqp_arm2[VALES2[1],]),
                        mean(Tqp_arm2[VALES2[2],]),
                        mean(Tqp_arm2[VALES2[3],]),
                        mean(Tqp_arm2[VALES2[4],]))

model_Effect_TQ2 = (model_estimates_TQ1 - model_estimates_TQ2)/model_estimates_TQ1

model_estimates_TQ1_min = c(low_Tq[VALES2[1],1],low_Tq[VALES2[2],1],low_Tq[VALES2[3],1],low_Tq[VALES2[4],1])
model_estimates_TQ2_min = c(low_Tq[VALES2[1],2],low_Tq[VALES2[2],2],low_Tq[VALES2[3],2],low_Tq[VALES2[4],2])

model_estimates_TQ1_max = c(upp_Tq[VALES2[1],1],upp_Tq[VALES2[2],1],upp_Tq[VALES2[3],1],upp_Tq[VALES2[4],1])
model_estimates_TQ2_max = c(upp_Tq[VALES2[1],2],upp_Tq[VALES2[2],2],upp_Tq[VALES2[3],2],upp_Tq[VALES2[4],2])

model_Effect_TQ2_min = (model_estimates_TQ1 - model_estimates_TQ2_min)/model_estimates_TQ1
model_Effect_TQ2_max = (model_estimates_TQ1 - model_estimates_TQ2_max)/model_estimates_TQ1


plot(Tqp_arm1[,1] ~ Tqp_arm1a[,1],ylim = c(0,0.4),xlim=c(0,3),ylab="Point slide-prevalence 2-10 years (%)",
     xlab = "",pch="",yaxt="n",xaxt="n",main = "Staedke")
axis(1,at = 0:2,labels=c("Jan 2017","Jan 2018","Jan 2019"))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))

abline(v=0,lty=2)
arrows(x0 = 0.677, y0 = 0.05, x1 = 0.677, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)

cols_Tq = adegenet::transp(c("darkred","aquamarine3"),0.6)
for(i in 1:2){
  polygon(c(Tqp_arm1a[,1],rev(Tqp_arm1a[,1])),c(upp_Tq[,i],rev(low_Tq[,i])),col = cols_Tq[i],border=NA)
  lines(med_Tq[,i] ~ Tqp_arm1a[,1],col = cols_Tq[i],lwd=2)
}
Tqo_arm_min = c(0.099,0.035)
Tqo_arm_max = c(0.416,0.401)
segments(x0=c(time_match_Q1[1],time_match_Q2[1]),
         x1=c(time_match_Q1[1],time_match_Q2[1]),
         y0 = Tqo_arm_min,
         y1 = Tqo_arm_max,col=c("darkred","aquamarine3"),lwd=2)

## Baseline all 4 arms
## 23.7% (12.8 ? 42)	25.1% (7.2 ? 41.7)	15.9% (1.5 ? 37.5)	8.2% (2.8 ? 34.1)

# abline(v=0.677,lty=2)

points(Tqo_arm1 ~ time_match_Q1,pch=15,cex=2,col="darkred")
points(Tqo_arm2 ~ time_match_Q2,pch=17,cex=2,col="aquamarine3")

legend("topright",legend = c("Pyrethroid ITN","Pyrethroid + PBO ITN"),
       col = c("darkred","aquamarine3"), pch = c(15,17),bty="n")


####################################################################
## O
## ID 14  West 2014
Top_arm1a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/west/arm1/VALIDATION_trial_arm_1_net_0_0_1.txt",header=TRUE)  ## filepath
Top_arm2a = read.table("P:/validating_modelling/output_files/FINAL_netparm1/west/arm2/VALIDATION_trial_arm_2_net_0_0_1.txt",header=TRUE)  ## filepath

Top_arm1 = array(dim=c(nrow(Top_arm1a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Top_arm1)){
  Top_arm1[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/west/arm1/VALIDATION_trial_arm_1_net_0_0_",run[i],".txt"),header=TRUE)$prev_0_14
}
Top_arm2 = array(dim=c(nrow(Top_arm2a),1000)) ## this will be 1000 once all the runs are done
for(i in 1:ncol(Top_arm2)){
  Top_arm2[,i] = read.table(paste0("P:/validating_modelling/output_files/FINAL_netparm1/west/arm2/VALIDATION_trial_arm_2_net_0_0_",run[i],".txt"),header=TRUE)$prev_0_14
}

plot(Top_arm1[,1] ~ Top_arm1a[,1],ylim = c(0,0.6),xlim=c(-0.5,2.5),ylab="Prevalence in children 6 months to 14 years (%)",
     xlab = "Time in years",pch="",yaxt="n",main="West",xaxt="n")
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
axis(1,at=c(0:2),labels=c("Jan 2011","Jan 2012","Jan 2013"))

abline(v=0.25,lty=2)
abline(v=0.5833,lty=2)
arrows(x0 = 0.9583333, y0 = 0.05, x1 = 0.9583333, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)
arrows(x0 = 1.2916667, y0 = 0.05, x1 = 1.2916667, y1 = 0, length = 0.1, angle = 30,
       code = 2, col = "black", lty = 1,
       lwd = 2)

legend("topright",legend = c("Pyrethroid ITN",
                             "Pyrethroid ITN + Pyrethroid then Carbamate IRS"),
       col = c("darkred","darkgreen"), pch = 15,bty="n")



upp_To = low_To = med_To = array(dim=c(1301,2))
for(t in 1:1301){
  upp_To[t,1] = quantile(Top_arm1[t,],0.975,na.rm=TRUE)
  low_To[t,1] = quantile(Top_arm1[t,],0.025,na.rm=TRUE)
  med_To[t,1] = quantile(Top_arm1[t,],0.5,na.rm=TRUE)
  
  upp_To[t,2] = quantile(Top_arm2[t,],0.975,na.rm=TRUE)
  low_To[t,2] = quantile(Top_arm2[t,],0.025,na.rm=TRUE)
  med_To[t,2] = quantile(Top_arm2[t,],0.5,na.rm=TRUE)
  
}
cols_To = adegenet::transp(c("darkred","darkgreen"),0.4)
cols_To2 = c("darkred","darkgreen")
lwd_To = c(2,2)
lty_To=c(1,4)
for(i in 1:2){
  polygon(c(Top_arm1a[,1],rev(Top_arm1a[,1])),c(upp_To[,i],rev(low_To[,i])),col = cols_To[i],border=NA)
  lines(med_To[,i] ~ Top_arm1a[,1],lwd=lwd_To[i],lty=lty_To[i],col=cols_To2[i])
}
Top_prev_recorded_t1 = c(0.246,0.23,0.3,0.24)
Top_prev_recorded_t2 = c(0.21, 0.136, 0.127, 0.134)

Top_time_match = c((1/12) * c(6,15,18,23))

effect_T0 = (Top_prev_recorded_t1 - Top_prev_recorded_t2)/Top_prev_recorded_t1



points(Top_prev_recorded_t1 ~ Top_time_match, col="red",pch=19,cex=1.5)
points(Top_prev_recorded_t2 ~ c(Top_time_match+0.04), col="darkgreen",pch=19,cex=1.5)

## Add in any uncertainty bars from the trial (need to acess these if there are any)
Top_prev_record_t1_min <- c(0.17,0.154, 0.202, 0.142)
Top_prev_record_t1_max <- c(0.343,0.342, 0.434, 0.389)
Top_prev_record_t2_min <- c(0.138,0.083, 0.074, 0.073)
Top_prev_record_t2_max <- c(0.305,0.214, 0.21, 0.234)

segments(x0 = Top_time_match, x1 = Top_time_match, y0 = Top_prev_record_t1_min, y1 = Top_prev_record_t1_max, col="red")
segments(x0 = c(Top_time_match+0.04), x1 = c(Top_time_match+0.04), y0 = Top_prev_record_t2_min, y1 = Top_prev_record_t2_max, col="darkgreen")

VALES = 261+52*Top_time_match
Top_arm1a$year[261+52*Top_time_match]

model_estimates_TO1 = c(med_To[287,1],med_To[326,1],med_To[339,1],med_To[361,1])
model_estimates_TO2 = c(med_To[287,2],med_To[326,2],med_To[339,2],med_To[361,2])

model_Effect_TO2 =  (model_estimates_TO1 - model_estimates_TO2)/model_estimates_TO1

model_estimates_TO1_min = c(low_To[287,1],low_To[326,1],low_To[339,1],low_To[361,1])
model_estimates_TO1_max = c(upp_To[287,1],upp_To[326,1],upp_To[339,1],upp_To[361,1])
model_estimates_TO2_min = c(low_To[287,2],low_To[326,2],low_To[339,2],low_To[361,2])
model_estimates_TO2_max = c(upp_To[287,2],upp_To[326,2],upp_To[339,2],upp_To[361,2])

model_Effect_TO_min =  (model_estimates_TO1 - model_estimates_TO2_min)/model_estimates_TO1
model_Effect_TO_max =  (model_estimates_TO1 - model_estimates_TO2_max)/model_estimates_TO1

effect_TO_min = (Top_prev_recorded_t1 - Top_prev_record_t2_min)/model_estimates_TO1
effect_TO_max = (Top_prev_recorded_t1 - Top_prev_record_t2_max)/model_estimates_TO1






#########################################################
##
## Efficacy figure
DATA_RESOURCE_Ef = expand.grid(study = rep(c("Curtis",#B
                                             "D'Alessandro",#C
                                             "Henry",#E
                                             "Kafy",#F
                                             "Marbiah",#H
                                             "Philips-Howard",#J
                                             "Protopopoff",#L
                                             "West",#O
                                             "Chaccour",#P
                                             "Staedke", #Q
                                             "Nevill",#S
                                             "Bradley"),#x
                                           c(8,5,2,2,#4,
                                             1,2,18,3,1,3,3,2)) )
#type of interventions
# "Corbel",#A         ITN,      
# "Curtis",#B         None, ITN, IRS
# "D'Alessandro",#C   None, None, None, None, None, ITN, ITN, ITN, ITN, ITN 
# "Henry",#E          None, ITN
# "Kafy",#F           ITN, ITN+IRS
# "Marbiah",#H        None, ITN
# "Philips-Howard",#J None, ITN, ITN
# "Protopopoff",#L    ITN, PBO_ITN, ITN + IRS, PBO_ITN + IRS
# "West",#O           ITN, ITN+IRS
# "Chaccour",#P       ITN, ITN + IRS
# "Staedke", #Q       ITN, PBO_ITN
# "Nevill",#S         None, ITN
# "Bradley"           IRS, IRS
#

DATA_RESOURCE_Ef[,2] = c(rep(c("ITN","IRS"),each=4),
                         rep("ITN",5),
                         c("ITN","ITN"),
                         c("ITN_IRS","ITN_IRS"),
                         #"ITN","ITN","ITN","ITN",
                         "ITN",
                         "ITN","ITN",
                         rep("PBO_ITN",6),rep(c("ITN_IRS","PBO_ITN_IRS"),each=6),
                         c("ITN_IRS","ITN_IRS","ITN_IRS"),
                         "ITN_IRS",
                         rep("PBO_ITN",3),
                         "ITN","ITN","ITN",
                         "IRS","IRS")


DATA_RESOURCE_Ef[,3] = c(2,5,8,11,2,5,8,11,
                         1,1,1,1,1,
                         48,48,
                         21,33,
                         #9,15,21,24,
                         15,
                         22,40,
                         4,9,16,21,28,33,4,9,16,21,NA,NA,4,9,16,21,NA,NA,
                         11,15,19,#2, 6, 10
                         17,
                         12,18,24,
                         20,24,27,
                         3,3)

OBSERVED_Effect_size = c(effect_TB,
                         effect_TC,
                         effect_TE,
                         effect_TF[2:3],
                         #effect_TG,# effect_TG2[2:4],effect_TG3[2:4],effect_TG4[2:4],
                         effect_TH,
                         effect_TJ,
                         effect_TL2,effect_TL3,effect_TL4,
                         # effect_TK[2:3],
                         effect_T0[2:4],
                         effect_TP,
                         effect_TQ2[2:4],
                         effect_TS,
                         effect_TU)
DATA_RESOURCE_Ef[,4] = OBSERVED_Effect_size
PREDICTED_Effect_size= c(model_Effect_TB,
                         model_Effect2_TC,
                         model_Effect_TE,
                         model_Effect_TF[2:3],
                         #model_Effect_TG,# model_Effect_TG2[2:4],model_Effect_TG3[2:4],model_Effect_TG4[2:4],
                         model_Effect_TH, 
                         model_Effect_TJ,
                         model_Effect_TL2,model_Effect_TL3,model_Effect_TL4,
                         # model_Effect_TK,
                         model_Effect_TO2[2:4],
                         model_Effect_TP2,
                         model_Effect_TQ2[2:4],
                         model_Effect_TS2,
                         model_Effect_TU)
DATA_RESOURCE_Ef[,5] = PREDICTED_Effect_size
summary.lm(glm(PREDICTED_Effect_size~OBSERVED_Effect_size+0))
og2 = glm(PREDICTED_Effect_size~OBSERVED_Effect_size+0)
quantile(resid(og2),c(0.2,0.8))

OBSERVED_Effect_size_upper = c(effect_TB,
                               effect_TC,
                               effect_TE_max,
                               effect_TF_max[2:3],
                               #effect_TG_max,# effect_TG2_max[2:4],effect_TG3_max[2:4],effect_TG4_max[2:4],
                               effect_TH,
                               effect_TJ,# effect_TK,
                               effect_TL2,effect_TL3,effect_TL4,
                               effect_TO_max[2:4],
                               effect_TP,
                               effect_TQ2[2:4],
                               effect_TS,
                               effect_TU_max)
DATA_RESOURCE_Ef[,6] = OBSERVED_Effect_size_upper

OBSERVED_Effect_size_lower = c(effect_TB,
                               effect_TC,
                               effect_TE_min,
                               effect_TF_min[2:3],
                               #effect_TG_min,# effect_TG2_min[2:4],effect_TG3_min[2:4],effect_TG4_min[2:4],
                               effect_TH,
                               effect_TJ,# effect_TK,
                               effect_TL2,effect_TL3,effect_TL4,
                               effect_TO_min[2:4],
                               effect_TP,
                               effect_TQ2[2:4],
                               effect_TS,
                               effect_TU_min)
DATA_RESOURCE_Ef[,7] = OBSERVED_Effect_size_lower


PREDICTED_Effect_size_upper = c(model_Effect_TB_max,
                                model_Effect2_TC_max,
                                model_Effect_TE_max,
                                model_Effect_TF_max[2:3],
                                #model_Effect_TG_max,# model_Effect_TG2_max[2:4],model_Effect_TG3_max[2:4],model_Effect_TG4_max[2:4],
                                model_Effect_TH_max,
                                model_Effect_TJ_max,# model_Effect_TK,
                                model_Effect_TL2_max,model_Effect_TL3_max,model_Effect_TL4_max,
                                model_Effect_TO_max[2:4],
                                model_Effect_TP_max,
                                model_Effect_TQ2_max[2:4],
                                model_Effect_TS_max,
                                model_Effect_TU_max)

PREDICTED_Effect_size_lower = c(model_Effect_TB_min,
                                model_Effect2_TC_min,
                                model_Effect_TE_min,
                                model_Effect_TF_min[2:3],
                                #model_Effect_TG_min,# model_Effect_TG2_min[2:4],model_Effect_TG3_min[2:4],model_Effect_TG4_min[2:4],
                                model_Effect_TH_min,
                                mean(c(model_Effect_TJ[1],model_Effect_TJ[3])),mean(c(model_Effect_TJ[2],model_Effect_TJ[4])),# model_Effect_TK,
                                model_Effect_TL2_min,model_Effect_TL3_min,model_Effect_TL4_min,
                                model_Effect_TO_min[2:4],
                                model_Effect_TP_min,
                                model_Effect_TQ2_min[2:4],
                                model_Effect_TS_min,
                                model_Effect_TU_min)
DATA_RESOURCE_Ef[,8] = PREDICTED_Effect_size_upper

DATA_RESOURCE_Ef[,9] = PREDICTED_Effect_size_lower

# colnames(DATA_RESOURCE_Ef) = c("Study","Intervention","Month","Obs_mean_eff","Pred_mean_eff","Obs_upp_eff","Obs_low_eff","Pred_upp_eff","Pred_low_eff")

# write.csv(DATA_RESOURCE_Ef,"Post_processing/DATA_RESOURCE_EFFICACY_1.csv")
###################################
###
### Actual estimates

# How well can we match the observed prevalence?

##Corbel
Tao_prev_measured = c(0.200)
Tao_prev_measured_max = c(0.24)
Tao_prev_measured_min = c(0.16)
model_TA_prev_measured = c(med_Ta[335,1])
model_TA_prev_measured_min = c(low_Ta[335,1])
model_TA_prev_measured_max = c(upp_Ta[335,1])

##Curtis
Tbo_prev_measured = c(0.834,0.853,0.817,0.675,
                      0.841,0.792,0.734,0.603,
                      0.762,0.801,0.730,0.637)
model_TB_prev_measured = c(med_Tb[270,1],med_Tb[283,1],med_Tb[296,1],med_Tb[309,1],
                           med_Tb[270,2],med_Tb[283,2],med_Tb[296,2],med_Tb[309,2],
                           med_Tb[270,3],med_Tb[283,3],med_Tb[296,3],med_Tb[309,3])
model_TB_prev_measured_min = c(upp_Tb[270,1],upp_Tb[283,1],upp_Tb[296,1],upp_Tb[309,1],
                               upp_Tb[270,2],upp_Tb[283,2],upp_Tb[296,2],upp_Tb[309,2],
                               upp_Tb[270,3],upp_Tb[283,3],upp_Tb[296,3],upp_Tb[309,3])
model_TB_prev_measured_max = c(low_Tb[270,1],low_Tb[283,1],low_Tb[296,1],low_Tb[309,1],
                               low_Tb[270,2],low_Tb[283,2],low_Tb[296,2],low_Tb[309,2],
                               low_Tb[270,3],low_Tb[283,3],low_Tb[296,3],low_Tb[309,3])

##D'Alessandro
Tco_prev_measured = c(0.367,	0.330,	0.258,	0.533, 0.447,0.282,	0.228,	0.159,	0.431,	0.710)
model_TC_prev_measured = c(med_Tc[343,1,1],med_Tc[343,1,2],med_Tc[343,1,3],med_Tc[343,1,4],med_Tc[343,1,5],
                           med_Tc[343,2,1],med_Tc[343,2,2],med_Tc[343,2,3],med_Tc[343,2,4],med_Tc[343,2,5])
model_TC_prev_measured_min = c(low_Tc[343,1,1],low_Tc[343,1,2],low_Tc[343,1,3],low_Tc[343,1,4],low_Tc[343,1,5],
                               low_Tc[343,2,1],low_Tc[343,2,2],low_Tc[343,2,3],low_Tc[343,2,4],low_Tc[343,2,5])
model_TC_prev_measured_max = c(upp_Tc[343,1,1],upp_Tc[343,1,2],upp_Tc[343,1,3],upp_Tc[343,1,4],upp_Tc[343,1,5],
                               upp_Tc[343,2,1],upp_Tc[343,2,2],upp_Tc[343,2,3],upp_Tc[343,2,4],upp_Tc[343,2,5])

##Henry
Teo_prev_measured = c(0.685,0.566)
Teo_prev_measured_max = c(0.721,0.602)
Teo_prev_measured_min = c(0.649,0.53)
model_TE_prev_measured = c(med_Te[309,1],med_Te[309,2])
model_TE_prev_measured_min = c(low_Te[309,1],low_Te[309,2])
model_TE_prev_measured_max = c(upp_Te[309,1],upp_Te[309,2])


##Kafy
Tfo_prev_measured = c(0.05,0.05,0.04,0.03)
Tfo_prev_measured_min = c(0.10,0.09,0.07,0.05)
Tfo_prev_measured_max = c(0.02,0.03,0.02,0.02)
model_TF_prev_measured = c(med_Tf[313,1],med_Tf[365,1],med_Tf[313,2],med_Tf[365,2])
model_TF_prev_measured_min = c(low_Tf[313,1],low_Tf[365,1],low_Tf[313,2],low_Tf[365,2])
model_TF_prev_measured_max = c(upp_Tf[313,1],upp_Tf[365,1],upp_Tf[313,2],upp_Tf[365,2])



##Marbiah
Tho_prev_measured = c(0.476,0.345)
model_TH_prev_measured = c(med_Th[300,1],med_Th[300,2])
model_TH_prev_measured_min = c(low_Th[300,1],low_Th[300,2])
model_TH_prev_measured_max = c(upp_Th[300,1],upp_Th[300,2])


##Philips-Howard
Tjo_prev_measured = c(0.514,0.391,0.325,0.338,0.325,0.338)
model_TJ_prev_measured = c(med_Tj[356,1],med_Tj[434,1],med_Tj[356,2],med_Tj[434,2],med_Tj[356,3],med_Tj[434,3])
model_TJ_prev_measured_min = c(low_Tj[356,1],low_Tj[434,1],low_Tj[356,2],low_Tj[434,2],low_Tj[356,3],low_Tj[434,3])
model_TJ_prev_measured_max = c(upp_Tj[356,1],upp_Tj[434,1],upp_Tj[356,2],upp_Tj[434,2],upp_Tj[356,3],upp_Tj[434,3])


# ##Pinder
# Tko_prev_measured = c(0.15,0.17)
# model_TK_prev_measured = c(med_Tk[356,1],med_Tk[408,1])
# model_TK_prev_measured_min = c(low_Tk[356,1],low_Tk[408,1])
# model_TK_prev_measured_max = c(upp_Tk[317,1],upp_Tk[356,1])


##Protopopoff
# VALES = c(278, 300, 330, 352, 382, 404)
VALES = round(261+52*time_match_L1,0)
Tlp_arm1a$year[261+52*time_match_L1]
Tlo_prev_measured = c(0.555, 0.553, 0.530, 0.68,0.827,0.648,
                      0.458, 0.311, 0.348, 0.459,0.6835,0.490,
                      0.386, 0.264, 0.399, 0.551,NA,NA,
                      0.375, 0.286, 0.291, 0.437,NA,NA)
model_TL_prev_measured = c(med_Tl[VALES[1],1],med_Tl[VALES[2],1],med_Tl[VALES[3],1],med_Tl[VALES[4],1],med_Tl[VALES[5],1],med_Tl[VALES[6],1],
                           med_Tl[VALES[1],2],med_Tl[VALES[2],2],med_Tl[VALES[3],2],med_Tl[VALES[4],2],med_Tl[VALES[5],2],med_Tl[VALES[6],2],
                           med_Tl[VALES[1],3],med_Tl[VALES[2],3],med_Tl[VALES[3],3],med_Tl[VALES[4],3],med_Tl[VALES[5],3],med_Tl[VALES[6],3],
                           med_Tl[VALES[1],4],med_Tl[VALES[2],4],med_Tl[VALES[3],4],med_Tl[VALES[4],4],med_Tl[VALES[5],4],med_Tl[VALES[6],4])
model_TL_prev_measured_min = c(low_Tl[VALES[1],1],low_Tl[VALES[2],1],low_Tl[VALES[3],1],low_Tl[VALES[4],1],low_Tl[VALES[5],1],low_Tl[VALES[6],1],
                               low_Tl[VALES[1],2],low_Tl[VALES[2],2],low_Tl[VALES[3],2],low_Tl[VALES[4],2],low_Tl[VALES[5],2],low_Tl[VALES[6],2],
                               low_Tl[VALES[1],3],low_Tl[VALES[2],3],low_Tl[VALES[3],3],low_Tl[VALES[4],3],low_Tl[VALES[5],3],low_Tl[VALES[6],3],
                               low_Tl[VALES[1],4],low_Tl[VALES[2],4],low_Tl[VALES[3],4],low_Tl[VALES[4],4],low_Tl[VALES[5],4],low_Tl[VALES[6],4])
model_TL_prev_measured_max = c(upp_Tl[VALES[1],1],upp_Tl[VALES[2],1],upp_Tl[VALES[3],1],upp_Tl[VALES[4],1],upp_Tl[VALES[5],1],upp_Tl[VALES[6],1],
                               upp_Tl[VALES[1],2],upp_Tl[VALES[2],2],upp_Tl[VALES[3],2],upp_Tl[VALES[4],2],upp_Tl[VALES[5],2],upp_Tl[VALES[6],2],
                               upp_Tl[VALES[1],3],upp_Tl[VALES[2],3],upp_Tl[VALES[3],3],upp_Tl[VALES[4],3],upp_Tl[VALES[5],3],upp_Tl[VALES[6],3],
                               upp_Tl[VALES[1],4],upp_Tl[VALES[2],4],upp_Tl[VALES[3],4],upp_Tl[VALES[4],4],upp_Tl[VALES[5],4],upp_Tl[VALES[6],4])

##West
Too_prev_measured = c(0.23,0.3,0.24,0.136, 0.127, 0.134)
Too_prev_measured_min <- c(0.154, 0.202, 0.142,0.083, 0.074, 0.073)
Too_prev_measured_max <- c(0.342, 0.434, 0.389,0.214, 0.21, 0.234)
model_TO_prev_measured = c(med_To[326,1],med_To[339,1],med_To[361,1],med_To[326,2],med_To[339,2],med_To[361,2])
model_TO_prev_measured_min = c(low_To[326,1],low_To[339,1],low_To[361,1],low_To[326,2],low_To[339,2],low_To[361,2])
model_TO_prev_measured_max = c(upp_To[326,1],upp_To[339,1],upp_To[361,1],upp_To[326,2],upp_To[339,2],upp_To[361,2])

##Chaccour
Tpo_prev_measured = c(0.43,0.34)
model_TP_prev_measured = c(med_Tp[348,1],med_Tp[348,2])
model_TP_prev_measured_min = c(low_Tp[348,1],low_Tp[348,2])
model_TP_prev_measured_max = c(upp_Tp[348,1],upp_Tp[348,2])

##Staedke
Tqo_prev_measured = c(0.145,0.13,0.14,0.107,0.106,0.118)
Tqo_arm1 = c(0.145,0.13,0.14)
Tqo_arm2 = c(0.107,0.106,0.118)
time_Q = 261+52*c(4/12,12/12,18/12,24/12)
time_match_Q2 = time_match_Q1 +0.02


model_TQ_prev_measured = c(med_Tq[time_Q[2],1],med_Tq[time_Q[3],1],med_Tq[time_Q[4],1],
                           med_Tq[time_Q[2],2],med_Tq[time_Q[3],2],med_Tq[time_Q[4],2])
model_TQ_prev_measured_min = c(low_Tq[time_Q[2],1],low_Tq[time_Q[3],1],low_Tq[time_Q[4],1],
                               low_Tq[time_Q[2],2],low_Tq[time_Q[3],2],low_Tq[time_Q[4],2])
model_TQ_prev_measured_max = c(upp_Tq[time_Q[2],1],upp_Tq[time_Q[3],1],upp_Tq[time_Q[4],1],
                               upp_Tq[time_Q[2],2],upp_Tq[time_Q[3],2],upp_Tq[time_Q[4],2])



##Nevill
Tso_prev_measured = c(0.392,0.254,0.115,0.168,0.12,0.078)
model_TS_prev_measured = c(med_Tp[348,1],med_Tp[365,1],med_Tp[378,1],med_Tp[348,2],med_Tp[365,2],med_Tp[378,2])
model_TS_prev_measured_min = c(low_Tp[348,1],low_Tp[365,1],low_Tp[378,1],low_Tp[348,2],low_Tp[365,2],low_Tp[378,2])
model_TS_prev_measured_max = c(upp_Tp[348,1],upp_Tp[365,1],upp_Tp[378,1],upp_Tp[348,2],upp_Tp[365,2],upp_Tp[378,2])

## Bradley
Tuo_prev_measured = c(0.168, 0.232)
Tuo_prev_measured_min <- c(0.111,0.16)
Tuo_prev_measured_max <- c(0.247,0.323)

model_TU_prev_measured = c(med_Tu[339,1],med_Tu[339,2])
model_TU_prev_measured_min = c(low_Tu[339,1],low_Tu[339,2])
model_TU_prev_measured_max = c(upp_Tu[339,1],upp_Tu[339,2])


## prevalence
Prevalence_observed_RCT = c(Tao_prev_measured,Tbo_prev_measured,Tco_prev_measured,
                            Teo_prev_measured,Tfo_prev_measured,#Tgo_prev_measured,
                            Tho_prev_measured,Tjo_prev_measured,
                            Tlo_prev_measured,Too_prev_measured,Tpo_prev_measured,
                            Tqo_prev_measured,Tso_prev_measured,Tuo_prev_measured)#,Txp_prev_measured)
Prevalence_observed_RCT_min = c(Tao_prev_measured_min,Tbo_prev_measured,Tco_prev_measured,
                                Teo_prev_measured_min,Tfo_prev_measured_min,#Tgo_prev_measured_min,
                                Tho_prev_measured,Tjo_prev_measured,
                                Tlo_prev_measured,Too_prev_measured_min,Tpo_prev_measured,
                                Tqo_prev_measured,Tso_prev_measured,Tuo_prev_measured_min)#,Txp_prev_measured)
Prevalence_observed_RCT_max = c(Tao_prev_measured_max,Tbo_prev_measured,Tco_prev_measured,
                                Teo_prev_measured_max,Tfo_prev_measured_max,#Tgo_prev_measured_max,
                                Tho_prev_measured,Tjo_prev_measured,#Tko_prev_measured,
                                Tlo_prev_measured,Too_prev_measured_max,Tpo_prev_measured,
                                Tqo_prev_measured,Tso_prev_measured,Tuo_prev_measured_max)#,Txp_prev_measured)

Prevalence_modelled_RCT =c(model_TA_prev_measured,model_TB_prev_measured,model_TC_prev_measured,
                           model_TE_prev_measured,model_TF_prev_measured,#model_TG_prev_measured,
                           model_TH_prev_measured,model_TJ_prev_measured,#model_TK_prev_measured,
                           model_TL_prev_measured,model_TO_prev_measured,model_TP_prev_measured,
                           model_TQ_prev_measured,model_TS_prev_measured,model_TU_prev_measured)#,model_TX_prev_measured)
Prevalence_modelled_RCT_min =c(model_TA_prev_measured_min,model_TB_prev_measured_min,model_TC_prev_measured_min,
                               model_TE_prev_measured_min,model_TF_prev_measured_min,#model_TG_prev_measured_min,
                               model_TH_prev_measured_min,model_TJ_prev_measured_min,#model_TK_prev_measured_min,
                               model_TL_prev_measured_min,model_TO_prev_measured_min,model_TP_prev_measured_min,
                               model_TQ_prev_measured_min,model_TS_prev_measured_min,model_TU_prev_measured_min)#,model_TX_prev_measured_min)
Prevalence_modelled_RCT_max =c(model_TA_prev_measured_max,model_TB_prev_measured_max,model_TC_prev_measured_max,
                               model_TE_prev_measured_max,model_TF_prev_measured_max,#model_TG_prev_measured_max,
                               model_TH_prev_measured_max,model_TJ_prev_measured_max,#model_TK_prev_measured_max,
                               model_TL_prev_measured_max,model_TO_prev_measured_max,model_TP_prev_measured_max,
                               model_TQ_prev_measured_max,model_TS_prev_measured_max,model_TU_prev_measured_max)#,model_TX_prev_measured_max)


#type of interventions
# "Corbel",#A         ITN, ITN, ITN+IRS     
# "Curtis",#B         None, ITN, IRS
# "D'Alessandro",#C   None, None, None, None, None, ITN, ITN, ITN, ITN, ITN 
# "Henry",#E          None, ITN
# "Kafy",#F           ITN, ITN+IRS
## "Loha",#G           None, ITN
# "Marbiah",#H        None, ITN
# "Philips-Howard",#J None, ITN, ITN
# "Protopopoff",#L    ITN, PBO_ITN, ITN + IRS, PBO_ITN + IRS
# "West",#O           ITN, ITN+IRS
# "Chaccour",#P       ITN, ITN + IRS
# "Staedke", #Q       ITN, PBO_ITN
# "Nevill",#S         None, ITN
# "Bradley"           IRS, IRS
DATA_RESOURCE_Pr = expand.grid(studdy_col2 = rep(c("Corbel",#A
                                                   "Curtis",#B
                                                   "D'Alessandro",#C
                                                   "Henry",#E
                                                   "Kafy",#F
                                                   #"Loha",#G
                                                   "Marbiah",#H
                                                   "Philips-Howard",#J
                                                   ##    "pink",#K eXTRA 
                                                   "Protopopoff",#L
                                                   "West",#O
                                                   "Chaccour",#P
                                                   "Staedke", #Q
                                                   "Nevill",#S
                                                   "Bradley"),#u
                                                 ## "Hamainza"),#X
                                                 c(1,12,10,2,4,#8,
                                                   2,6,#2,
                                                   24,6,2,6,6,2)))##,9
DATA_RESOURCE_Pr[,2] = c("ITN",#A Corbel; Delt itn targ vs uni, vs targ+IRS, vs uni+IRS
                         rep(c("ITN_IRS","IRS"),each=6),#B Curtis: dip nets vs dip+IRS, vs IRS 
                         rep(c("None","ITN"),each=5),#C D'Alessandro undipped vs dipped permethrin
                         "None","ITN",#E Henry undipped vs dipped lambda
                         "ITN","ITN","ITN_IRS","ITN_IRS",#F Kafy
                         #"None","None","None","None","ITN","ITN","ITN","ITN",#G Loha do nothing then ITN delta
                         "None","ITN",#H Marbiah
                         "None","None","ITN","ITN","ITN","ITN",#J Phillips-Howard
                         #  "ITN","ITN_IRS",#K Pinder
                         rep(c("ITN","PBO_ITN","ITN_IRS","PBO_ITN_IRS"),each=6),#L Proto
                         "ITN","ITN","ITN","ITN_IRS","ITN_IRS","ITN_IRS",#O West
                         "ITN","ITN_IRS",#P chaccour
                         "ITN","ITN","ITN","PBO_ITN","PBO_ITN","PBO_ITN",#Q Staedke
                         "None","None","None","ITN","ITN","ITN",#S Nevill
                         "ITN","ITN")#U Bradley
                       #  "ITN_IRS","ITN_IRS","ITN_IRS","ITN_IRS","ITN_IRS")#,#W Abuaku ##first 2 points alpha-IRS, last 3 Actellic
# "None","None","None","ITN_IRS","ITN_IRS","ITN_IRS","ITN_IRS","ITN_IRS","ITN_IRS")

DATA_RESOURCE_Pr[,3] = Prevalence_observed_RCT
DATA_RESOURCE_Pr[,4] = Prevalence_observed_RCT_max
DATA_RESOURCE_Pr[,5] = Prevalence_observed_RCT_min
DATA_RESOURCE_Pr[,6] = Prevalence_modelled_RCT
DATA_RESOURCE_Pr[,7] = Prevalence_modelled_RCT_max
DATA_RESOURCE_Pr[,8] = Prevalence_modelled_RCT_min

DATA_RESOURCE_Pr[,9] = ifelse(DATA_RESOURCE_Pr[,2] == "None","grey",
                              ifelse(DATA_RESOURCE_Pr[,2] == "ITN", "darkred",
                                     ifelse(DATA_RESOURCE_Pr[,2] == "IRS","purple",
                                            ifelse(DATA_RESOURCE_Pr[,2] == "ITN_IRS", "aquamarine3",
                                                   ifelse(DATA_RESOURCE_Pr[,2] == "PBO_ITN", "blue","orange")))))
DATA_RESOURCE_Pr[,10] = rep(1:13,c(1,12,10,2,4,#8,
                                   2,6,#2,
                                   24,6,2,6,6,2))#,9

# studdy2 = rep(c("Corbel",#A
#                "Curtis",#B
#                "D'Alessandro",#C
#                "Henry",#E
#                "Kafy",#F
#                "Loha",#G
#                "Marbiah",#H
#                "Philips-Howard",#J
#     ##out           "Pinder", ##K EXTRA! 
#                "Protopopoff",#L
#                "West",#O
#                "Chaccour",#p
#                "Staedke", #q
#                "Nevill"),
#              c(3,12,10,2,4,8,2,6,2,24,6,2,))
studdy_col2 = rep(c("purple",#A
                    "darkred",#B
                    "red",#C
                    "orange",#E
                    "blue",#F
                    # "lightblue",#G
                    "darkgreen",#H
                    "darkblue",#J
                    ##    "pink",#K eXTRA 
                    "green",#L
                    "gold2",#O
                    "grey",#P
                    "aquamarine1", #Q
                    "black",#S
                    "pink"#u
                   ),#W
                  # "aquamarine4"),#X
                  c(1,12,10,2,4,#8,
                    2,6,#2,
                    24,6,2,6,6,2))#,9
DATA_RESOURCE_Pr[,11] = studdy_col2
intvn_symb2 = rep(c(15,#A
                    1,#B
                    1,#C
                    1,#E
                    19,#F
                    # 19,#G
                    17,#H
                    19,#J
                    # 19,#K
                    8,#L
                    19,#O
                    18,#P
                    8,#Q
                    1,#S
                    19),#u
                  # 19),#X
                  c(1,12,10,2,4,#8,
                    2,6,#2,
                    24,6,2,6,6,2))#,9
DATA_RESOURCE_Pr[,12] = intvn_symb2
## Type of intervention
intvn_symb3 = c(0,#A Corbel; Delt itn targ vs uni, vs targ+IRS, vs uni+IRS
                rep(c(0,5),each=6),#B Curtis: dip nets vs dip+IRS, vs IRS 
                rep(c(1,17),each=5),#C D'Alessandro undipped vs dipped permethrin
                1,19,#E Henry undipped vs dipped lambda
                15,15,0,0,#F Kafy
                # 1,1,1,1,15,15,15,15,#G Loha do nothing then ITN delta
                1,17,#H Marbiah
                1,1,19,19,19,19,#J Phillips-Howard
                19,0,#K Pinder
                rep(c(19,8,0,11),each=6),#L Proto
                19,19,19,0,0,0,#O West
                18,0,#P chaccour
                15,15,15,8,8,8,#Q Staedke
                1,1,1,19,19,19,#S Nevill
                19,19)#,#U Bradley
                # 0,0,0,0,0)#,#W Abuaku ##first 2 points alpha-IRS, last 3 Actellic
# 1,1,1,0,0,0,0,0,0)
head(DATA_RESOURCE_Pr)
colnames(DATA_RESOURCE_Pr) = c("PI","Intervention",
                               "Prev_obs","prev_obs_max","prev_obs_min",
                               "Prev_model","prev_model_max","prev_model_min",
                               "intvn_colour","study_order",
                               "study_colour","study_symbol")

write.csv(DATA_RESOURCE_Pr,"Post_processing/DATA_RESOURCE_PREVALENCE_1.csv")




par(mar=c(5,5,2,2))
plot(Prevalence_modelled_RCT~Prevalence_observed_RCT,ylim=c(0,1),xlim=c(0,1),
     ylab = "Transmission model prevalence estimate (%)",
     xlab = "Randomised control trial oberved prevalence (%)",
     xaxt = "n", yaxt = "n",cex.lab=1.4,cex.axis=1.4,pch="")
axis(1,at=seq(0,1,0.2),labels = seq(0,100,20),cex.lab=1.4,cex.axis=1.4)
axis(2,las=2,at=seq(0,1,0.2),labels = seq(0,100,20),cex.lab=1.4,cex.axis=1.4)
Y = seq(-10,10,length=20)
Yu = seq(-10,10,length=20)+0.1
Yl = seq(-10,10,length=20)-0.1
X = seq(-10,10,length=20)
lines(Y~X,lty=2)
# abline(lm(PREDICTED_Effect_size ~ OBSERVED_Effect_size+0),lty=2,lwd=2)
# polygon(c(X,rev(X)),c(Yu,rev(Yl)),border=NA,col=adegenet::transp("grey",0.4))
for(i in 1:length(unique(DATA_RESOURCE_Pr[,1]))){
  points(Prevalence_modelled_RCT~Prevalence_observed_RCT,
         col=DATA_RESOURCE_Pr[,9],pch=DATA_RESOURCE_Pr[,10],cex=1.8)  
}


segments(x0=Prevalence_observed_RCT,x1=Prevalence_observed_RCT,y0=Prevalence_modelled_RCT_min,y1=Prevalence_modelled_RCT_max,
         col=DATA_RESOURCE_Pr[,9],lwd=1)
segments(x0=Prevalence_observed_RCT_min,x1=Prevalence_observed_RCT_max,y0=Prevalence_modelled_RCT,y1=Prevalence_modelled_RCT,
         col=DATA_RESOURCE_Pr[,9],lwd=1)

legend("topleft",
       legend = c("Marbiah (Jun 1992)",
                  "D'Alessandro (Jul 1992)",
                  "Nevill (Jul 1993)",
                  "Curtis (Dec 1995)",
                  "Phillips-Howard (Jan 1997)",
                  "Henry (Jun 1999)",
                  "Corbel (Jun 2008)"),
       title = "RCT Principle Investigator (start month, year)",
       
       col="black",ncol=2,
       pch=c(7,3,13,2,8,4,1),
       cex=0.9,bty="n")

legend("bottomright",
       legend = c("Kafy (Apr 2011)",
                  "Abuaku (May 2011)",
                  "West (Dec 2011)",
                  "Protopopoff (Mar 2014)",
                  "Bradley (Mar-Apr 2014)",
                  "Loha (Sep 2014)",
                  "Chaccour (Oct 2016)",
                  "Staedke (Jul 2017)"),
       # col = c("darkgreen",
       #         "red",
       #         "black",
       #         "darkred",
       #         "darkblue",
       #         "orange",
       #         "purple",
       #         "aquamarine4",
       #         "blue",
       #         "yellow",
       #         "gold2",
       #         "green",
       #         "pink",
       #         "lightblue",
       #         "grey",
       #         "aquamarine1"),
       col="black",
       pch=c(5,15,10,9,14,6,11,12),
       cex=0.9,bty="n")

study_list = rep(c("corbel",#A
                   "curtis",#B
                   "dalessandro",#C
                   "henry",#E
                   "kafy",#F
                   #"loha",#G
                   "marbiah",#H
                   "philsHoward",#J
                   ##    "pink",#K eXTRA 
                   "protopopoff",#L
                   "west",#O
                   "chaccour",#P
                   "staedke", #Q
                   "nevil",#S
                   "bradley",#u
                   "abuaku"),#,#W
                 # "hamainza"),#X
                 c(3,12,10,2,4,#8,
                   2,6,#2,
                   24,6,2,6,6,2,5))#,9

summary.lm(lm(Prevalence_modelled_RCT~Prevalence_observed_RCT+0))
og1 = lm(Prevalence_modelled_RCT~Prevalence_observed_RCT+0)
resid(og1)
quantile(resid(og1),c(0.1,0.9))
dat_resid_prev_T = data.frame(Prevalence_modelled_RCT,
                              Prevalence_observed_RCT,study_list)
dat_resid_prev = dat_resid_prev_T[complete.cases(dat_resid_prev_T), ]
dat_resid_prev$RESIDUALS = resid(og1)
# write.csv(dat_resid_prev,"H:\\Ellie\\Vector interventions Model validation\\Submission v1\\analysis_data_prevalence.csv")
# abline(lm(Prevalence_modelled_RCT~Prevalence_observed_RCT+0),lwd=2,lty=4)


studdy = rep(c("Corbel",#A
               "Curtis",#B
               "D'Alessandro",#C
               "Henry",#E
               "Kafy",#F
               #"Loha",#G
               "Marbiah",#H
               "Philips-Howard",#J
               "Protopopoff",#L
               "West",#O
               "Chaccour",#P
               "Staedke", #Q
               "Nevill",#S
               "Bradley"),#x
             c(2,8,5,2,2,#4,
               1,2,18,3,1,3,3,2))
studdy_col = rep(c("purple",#A
                   "darkred",#B
                   "red",#C
                   "orange",#E
                   "blue",#F
                   #"lightblue",#G
                   "darkgreen",#H
                   "darkblue",#J
                   "green",#L
                   "gold2",#O
                   "grey",#P
                   "aquamarine1",#Q 
                   "black",#S
                   "pink"),#W
                 c(2,8,5,2,2,#4,
                   1,2,18,3,1,3,3,2))
intvn_symb = rep(c(15,#A
                   1,#B
                   1,#C
                   1,#E
                   19,#F
                   #19,#G
                   17,#H
                   19,#J
                   8,#L
                   19,#O
                   18,#P
                   8,#Q
                   1,#s
                   19),#u
                 c(2,8,5,2,2,#4,
                   1,2,18,3,1,3,3,2))
study_symb = c(1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,4,4,5,5,#6,6,6,6,
               7,8,8,
               9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,10,10,10,11,12,12,12,
               13,13,13,14,14)
month_prev_COL = c("darkgreen","darkgreen",#A
                   "grey","grey","lightgreen","lightgreen","grey","grey","lightgreen","lightgreen",#B
                   "grey","grey","grey","grey","grey",#C
                   "black","black",#E
                   "darkblue","black",#F
                   # "lightgreen","darkgreen","darkblue","darkblue",#G
                   "darkblue",#H
                   "darkblue","black",#J
                   "grey","lightgreen","darkgreen","darkblue","black","black","grey","lightgreen","darkgreen","darkblue","black","black","grey","lightgreen","darkgreen","darkblue","black","black",#L
                   "lightgreen","darkgreen","darkblue",#O
                   "darkgreen",#P
                   "lightgreen","darkgreen","darkblue",#Q
                   "darkblue","darkblue","black",#s
                   "grey","grey")#u

DATA_RESOURCE_Ef[,10] = studdy
DATA_RESOURCE_Ef[,11] = studdy_col
DATA_RESOURCE_Ef[,12] = intvn_symb
DATA_RESOURCE_Ef[,13] = month_prev_COL
DATA_RESOURCE_Ef[,14] = study_symb
colnames(DATA_RESOURCE_Ef) = c("PI","Intervention","MONTH_prev_observed",
                               "Obs_mean_eff","Pred_eff","Obs_eff_max","Obs_eff_min",
                               "Eff_model_min","Eff_model_max",
                               "study","study_col","intervention_symbol",
                               "col_rel_eff_month","study_symb")
# write.csv(DATA_RESOURCE_Ef,"Q:/RProjects/Model_validation_ITN_IRS/Post_processing/data/DATA_RESOURCE_EFFICACY_net_parms1.csv")
