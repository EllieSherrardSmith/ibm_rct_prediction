################################
##
## Analysis

## 1 We focus on predicting prevalence given that there are far more data for 
##   this epidemiological outcome (Data S1.2). The model is calibrated individually 
##   to the baseline prevalence in a defined age cohort for a total of 33 trial 
##   arms from 15 RCTs (Table S1). The absolute mean estimates for prevalence at 
##   multiple timepoints through each trial follow up period (where available) are 
##   linearly associated with the mean prediction from the model outputs matching 
##   the age-cohort and time point of each trial arm. There are xx cross-sectional
##   surveys in total.

##   We calculated the absolute difference between the prevalence observed and 
##   predicted from the 1:1 line (values of 0 indicate the model is perfectly 
##   predicting the trial result). The accuracy of the model at predicting ITNs and 
##   IRS results were compared (Fig. S19) using analysis of variance. 

##   The analyses were repeated for the mean efficacy of each trial, calculated as 
##   the average efficacy given all time points reported in the follow up process. 
##   These ranged from 1 to 6 cross sectional prevalence surveys across the trials. 


# prev = read.csv("Post_processing/data/DATA_RESOURCE_PREVALENCE_global.csv",header=TRUE)
# effi = read.csv("Post_processing/data/DATA_RESOURCE_EFFICACY_global.csv",header=TRUE)

effi = read.csv("Post_processing/data/DATA_RESOURCE_EFFICACY_global_for_fig1_drop_targetedFINAL.csv",header=TRUE)
prev = read.csv("Post_processing/data/DATA_RESOURCE_PREVALENCE_global_update_droptargetedFINAL.csv",header=TRUE)

##############################################################
##
## Analysing the residuals

## can the model predictions predict the observed data?

## 1 is logistic All data
## 2 is log-logistic All data
## 3 is log-logistic West
## 4 is log-logistic East
## 5 is logistic West
## 6 in logistic East

## Can the model estimate predict the trial obseration?
## Table 1
length(prev$prev_obs)
summary.lm(lm(prev$prev_obs~prev$prev_pred_1+0)) ## 93.68
summary.lm(lm(prev$prev_obs~prev$prev_pred_2+0)) ##BEST 95.51, F=1832
summary.lm(lm(prev$prev_obs~prev$prev_pred_3+0)) ## 94.95, F=1619
summary.lm(lm(prev$prev_obs~prev$prev_pred_4+0)) ## 95.49, F=1822
summary.lm(lm(prev$prev_obs~prev$prev_pred_5+0)) ## 93.36, F=1210, 
summary.lm(lm(prev$prev_obs~prev$prev_pred_6+0)) ## 94.24, F=1408, df = 85, p < 0.0001

og1 = lm(prev$prev_pred_2~prev$prev_obs+0)

## Table 1
summary.lm(lm(effi$Obs_mean_eff~effi$Pred_1_eff+0))
summary.lm(lm(effi$Obs_mean_eff~effi$Pred_2_eff+0)) ## BEST 
summary.lm(lm(effi$Obs_mean_eff~effi$Pred_3_eff+0)) 
summary.lm(lm(effi$Obs_mean_eff~effi$Pred_4_eff+0))
summary.lm(lm(effi$Obs_mean_eff~effi$Pred_5_eff+0))
summary.lm(lm(effi$Obs_mean_eff~effi$Pred_6_eff+0))


## Supplementary plot for prevalence
pred_plot_f = function(data_source,
                       data_source_min,
                       data_source_max,
                       main_title){
  
  plot(prev$prev_obs~data_source,
       ylab="Trial observed prevalence (%)",
       xlab="Model predicted prevalence (%)",
       main=main_title,
       cex.lab = 1.6,cex.axis=1.6,cex=1.5,
       ylim=c(0,1),yaxt="n",xlim=c(0,1),xaxt="n",pch=prev$study_symbol,
       col=as.character(prev$intvn_colour))
  axis(1,at=seq(0,1,0.2),labels=seq(0,100,20),cex.axis=1.6)
  axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.axis=1.6)
  for(i in 1:length(prev$prev_obs)){
    segments(x0=data_source_max[i],x1=data_source_min[i],
             y0=prev$prev_obs[i],y1=prev$prev_obs[i],
             col=as.character(prev$intvn_colour)[i])
    segments(x0=data_source[i],x1=data_source[i],
             y0=prev$prev_obs_max[i],y1=prev$prev_obs_min[i],
             col=as.character(prev$intvn_colour)[i])
    
  }
  abline(a = 0,b = 1,lty=2)
  
}

## Figure S4
##width 1300, height 700
par(mfrow=c(2,3))

pred_plot_f(prev$prev_pred_1,
            prev$prev_pred_1_min,
            prev$prev_pred_1_max,
            "Logistic association, all data")
text(0.3,1,"Adj R-squared: 0.9387",cex = 1.6)

pred_plot_f(prev$prev_pred_5,
            prev$prev_pred_5_min,
            prev$prev_pred_5_max,
            "Logistic association, West data")
text(0.3,1,"Adj R-squared: 0.9359",cex = 1.6)

pred_plot_f(prev$prev_pred_6,
            prev$prev_pred_6_min,
            prev$prev_pred_6_max,
            "Logistic association, East data")
text(0.3,1,"Adj R-squared: 0.9424",cex = 1.6)

pred_plot_f(prev$prev_pred_2,
            prev$prev_pred_2_min,
            prev$prev_pred_2_max,
            "Log-logistic association, all data")
text(0.3,1,"Adj R-squared: 0.9566",cex = 1.6)

pred_plot_f(prev$prev_pred_3,
            prev$prev_pred_3_min,
            prev$prev_pred_3_max,
            "Log-logistic association, West data")
text(0.3,1,"Adj R-squared: 0.9493",cex = 1.6)

pred_plot_f(prev$prev_pred_4,
            prev$prev_pred_4_min,
            prev$prev_pred_4_max,
            "Log-logistic association, East data")
text(0.3,1,"Adj R-squared: 0.9564",cex = 1.6)

par(xpd=NA,cex = 1.11)
# 
# text(x = -585, y = 2.75,"(A)")
# text(x = -315, y = 2.75,"(B)")
# text(x = -50, y = 2.75,"(C)")
# 
# text(x = -585, y = 1.1,"(D)")
# text(x = -315, y = 1.1,"(E)")
# text(x = -50, y = 1.1,"(F)")
# 

text(x = -2.77, y = 2.8,"(A)")
text(x = -1.51, y = 2.8,"(B)")
text(x = -0.15, y = 2.8,"(C)")

text(x = -2.77, y = 1.27,"(D)")
text(x = -1.51, y = 1.27,"(E)")
text(x = -0.15, y = 1.27,"(F)")


# resid(og1)
quantile(resid(og1),c(0.1,0.9))

##best performing
prev$Prev_model = prev$prev_pred_2


summary.lm(lm(prev$Prev_model~prev$prev_obs))
og1a = lm(prev$Prev_model~prev$prev_obs)
# resid(og1)
quantile(resid(og1a),c(0.1,0.9))

dat_resid_prev_T = data.frame(Prevalence_modelled_RCT = prev$Prev_model,
                              Prevalence_observed_RCT = prev$prev_obs,
                              study_list = prev$study_colour,
                              intervention = prev$Intervention,
                              col_intn = prev$intvn_colour,
                              symbol_study = prev$study_symbol,
                              time_months = prev$time_months)
dat_resid_prev = dat_resid_prev_T[complete.cases(dat_resid_prev_T), ]
dat_resid_prev$RESIDUALS = resid(og1)

effi2$means_obse = as.numeric(tapply(effi$Obs_mean_eff,effi$reps,mean))
effi2$means_pred1 = as.numeric(tapply(effi$Pred_1_eff,effi$reps,mean))
effi2$means_pred2 = as.numeric(tapply(effi$Pred_2_eff,effi$reps,mean))
effi2$means_pred3 = as.numeric(tapply(effi$Pred_3_eff,effi$reps,mean))
effi2$means_pred4 = as.numeric(tapply(effi$Pred_4_eff,effi$reps,mean))
effi2$means_pred5 = as.numeric(tapply(effi$Pred_5_eff,effi$reps,mean))
effi2$means_pred6 = as.numeric(tapply(effi$Pred_6_eff,effi$reps,mean))


summary.lm(lm(effi2$means_pred1~effi2$means_obse+0))
summary.lm(lm(effi2$means_pred2~effi2$means_obse+0))
summary.lm(lm(effi2$means_pred3~effi2$means_obse+0))
summary.lm(lm(effi2$means_pred4~effi2$means_obse+0))
summary.lm(lm(effi2$means_pred5~effi2$means_obse+0))
summary.lm(lm(effi2$means_pred6~effi2$means_obse+0))


og1 = lm(effi2$means_pred2~effi2$means_obse+0)
# resid(og1)
quantile(resid(og1),c(0.1,0.9))

dat_resid_eff_T = data.frame(Mean_efficacy_pred = effi2$means_pred2,
                             Mean_efficacy_obs = effi2$means_obse,
                             study_list = effi2$PI,
                             intervention = effi2$Intervention,
                             col_intn = effi2$intervention_colour,
                             symbol_study = effi2$study_symb,
                             time_months = effi2$Month)
dat_resid_eff = dat_resid_eff_T[complete.cases(dat_resid_eff_T), ]
dat_resid_eff$RESIDUALS = resid(og1)
# write.csv(dat_resid_prev,"H:\\Ellie\\Vector interventions Model validation\\Submission v1\\analysis_data_mean_efficacy.csv")

##############################################
##
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
                    ifelse(prev$Intervention == "PBO_ITN_IRS","IRS","ITN")))

## Any net with or without IRS              
effi$INTN3 = ifelse(effi$Intervention == "targetITN_IRS","IRS","ITN")
prev$INTN3 = ifelse(prev$Intervention == "None","None","ITN")


## Studies using any net without IRS
##Table 1
length(effi$Pred_2_eff[effi$INTN2 == "ITN"])
summary.lm(lm(effi$Obs_mean_eff[effi$INTN2 == "ITN"]~
                effi$Pred_1_eff[effi$INTN2 == "ITN"]+0))
summary.lm(lm(effi$Obs_mean_eff[effi$INTN2 == "ITN"]~
                effi$Pred_2_eff[effi$INTN2 == "ITN"]+0))
summary.lm(lm(effi$Obs_mean_eff[effi$INTN2 == "ITN"]~
                effi$Pred_3_eff[effi$INTN2 == "ITN"]+0))
summary.lm(lm(effi$Obs_mean_eff[effi$INTN2 == "ITN"]~
                effi$Pred_4_eff[effi$INTN2 == "ITN"]+0))
summary.lm(lm(effi$Obs_mean_eff[effi$INTN2 == "ITN"]~
                effi$Pred_5_eff[effi$INTN2 == "ITN"]+0))
summary.lm(lm(effi$Obs_mean_eff[effi$INTN2 == "ITN"]~
                effi$Pred_6_eff[effi$INTN2 == "ITN"]+0))

##Table 1
length(prev$prev_obs[prev$INTN2 == "ITN"])
summary.lm(lm(prev$prev_obs[prev$INTN2 == "ITN"]~
                prev$prev_pred_1[prev$INTN2 == "ITN"]+0))
summary.lm(lm(prev$prev_obs[prev$INTN2 == "ITN"]~
                prev$prev_pred_2[prev$INTN2 == "ITN"]+0))
summary.lm(lm(prev$prev_obs[prev$INTN2 == "ITN"]~
                prev$prev_pred_3[prev$INTN2 == "ITN"]+0))
summary.lm(lm(prev$prev_obs[prev$INTN2 == "ITN"]~
                prev$prev_pred_4[prev$INTN2 == "ITN"]+0))
summary.lm(lm(prev$prev_obs[prev$INTN2 == "ITN"]~
                prev$prev_pred_5[prev$INTN2 == "ITN"]+0))
summary.lm(lm(prev$prev_obs[prev$INTN2 == "ITN"]~
                prev$prev_pred_6[prev$INTN2 == "ITN"]+0))


## Studies using any net with or without IRS
length(effi$Pred_mean_eff[effi$INTN3 == "ITN"])
summary.lm(lm(effi$Pred_2_eff[effi$INTN3 == "ITN"]~
                effi$Obs_mean_eff[effi$INTN3 == "ITN"]+0))

length(prev$prev_obs[prev$INTN3 == "ITN"])
summary.lm(lm(prev$prev_obs[prev$INTN3 == "ITN"]~
                prev$prev_pred_2[prev$INTN3 == "ITN"]+0))


## Studies any IRS
length(effi$Obs_mean_eff[effi$INTN2 == "IRS"])
summary.lm(lm(effi$Obs_mean_eff[effi$INTN2 == "IRS"]~effi$Pred_1_eff[effi$INTN2 == "IRS"]+0))
summary.lm(lm(effi$Obs_mean_eff[effi$INTN2 == "IRS"]~effi$Pred_2_eff[effi$INTN2 == "IRS"]+0))
summary.lm(lm(effi$Obs_mean_eff[effi$INTN2 == "IRS"]~effi$Pred_3_eff[effi$INTN2 == "IRS"]+0))
summary.lm(lm(effi$Obs_mean_eff[effi$INTN2 == "IRS"]~effi$Pred_4_eff[effi$INTN2 == "IRS"]+0))
summary.lm(lm(effi$Obs_mean_eff[effi$INTN2 == "IRS"]~effi$Pred_5_eff[effi$INTN2 == "IRS"]+0))
summary.lm(lm(effi$Obs_mean_eff[effi$INTN2 == "IRS"]~effi$Pred_6_eff[effi$INTN2 == "IRS"]+0))

## Table 1
length(prev$prev_obs[prev$INTN2 == "IRS"])
summary.lm(lm(prev$prev_obs[prev$INTN2 == "IRS"]~prev$prev_pred_1[prev$INTN2 == "IRS"]+0))
summary.lm(lm(prev$prev_obs[prev$INTN2 == "IRS"]~prev$prev_pred_2[prev$INTN2 == "IRS"]+0))
summary.lm(lm(prev$prev_obs[prev$INTN2 == "IRS"]~prev$prev_pred_3[prev$INTN2 == "IRS"]+0))
summary.lm(lm(prev$prev_obs[prev$INTN2 == "IRS"]~prev$prev_pred_4[prev$INTN2 == "IRS"]+0))
summary.lm(lm(prev$prev_obs[prev$INTN2 == "IRS"]~prev$prev_pred_5[prev$INTN2 == "IRS"]+0))
summary.lm(lm(prev$prev_obs[prev$INTN2 == "IRS"]~prev$prev_pred_6[prev$INTN2 == "IRS"]+0))


## Studies using Pyrethroid ITNs
length(effi$Pred_mean_eff[effi$Intervention == "ITN"])
summary.lm(lm(effi$Pred_mean_eff[effi$Intervention == "ITN"]~effi$Obs_mean_eff[effi$Intervention == "ITN"]+0))

length(prev$Prev_model[prev$Intervention == "ITN"])
summary.lm(lm(prev$Prev_model[prev$Intervention == "ITN"]~prev$prev_obs[prev$Intervention == "ITN"]+0))

## Studies using Pyrethroid PBO ITNs
length(effi$Pred_mean_eff[effi$Intervention == "PBO_ITN"])
summary.lm(lm(effi$Pred_mean_eff[effi$Intervention == "PBO_ITN"]~effi$Obs_mean_eff[effi$Intervention == "PBO_ITN"]+0))

## tABLE 1
length(prev$prev_obs[prev$Intervention == "PBO_ITN"])
summary.lm(lm(prev$prev_obs[prev$Intervention == "PBO_ITN"]~prev$prev_pred_1[prev$Intervention == "PBO_ITN"]+0))
summary.lm(lm(prev$prev_obs[prev$Intervention == "PBO_ITN"]~prev$prev_pred_2[prev$Intervention == "PBO_ITN"]+0))
summary.lm(lm(prev$prev_obs[prev$Intervention == "PBO_ITN"]~prev$prev_pred_3[prev$Intervention == "PBO_ITN"]+0))
summary.lm(lm(prev$prev_obs[prev$Intervention == "PBO_ITN"]~prev$prev_pred_4[prev$Intervention == "PBO_ITN"]+0))
summary.lm(lm(prev$prev_obs[prev$Intervention == "PBO_ITN"]~prev$prev_pred_5[prev$Intervention == "PBO_ITN"]+0))
summary.lm(lm(prev$prev_obs[prev$Intervention == "PBO_ITN"]~prev$prev_pred_6[prev$Intervention == "PBO_ITN"]+0))

length(effi$Obs_mean_eff[effi$Intervention == "PBO_ITN"])
summary.lm(lm(effi$Obs_mean_eff[effi$Intervention == "PBO_ITN"]~effi$Pred_1_eff[effi$Intervention == "PBO_ITN"]+0))
summary.lm(lm(effi$Obs_mean_eff[effi$Intervention == "PBO_ITN"]~effi$Pred_2_eff[effi$Intervention == "PBO_ITN"]+0))
summary.lm(lm(effi$Obs_mean_eff[effi$Intervention == "PBO_ITN"]~effi$Pred_3_eff[effi$Intervention == "PBO_ITN"]+0))
summary.lm(lm(effi$Obs_mean_eff[effi$Intervention == "PBO_ITN"]~effi$Pred_4_eff[effi$Intervention == "PBO_ITN"]+0))
summary.lm(lm(effi$Obs_mean_eff[effi$Intervention == "PBO_ITN"]~effi$Pred_5_eff[effi$Intervention == "PBO_ITN"]+0))
summary.lm(lm(effi$Obs_mean_eff[effi$Intervention == "PBO_ITN"]~effi$Pred_6_eff[effi$Intervention == "PBO_ITN"]+0))

## tABLE 1
length(prev$prev_obs[prev$Intervention == "PBO_ITN" | prev$Intervention == "PBO_ITN_IRS"])
summary.lm(lm(prev$prev_obs[prev$Intervention == "PBO_ITN" | prev$Intervention == "PBO_ITN_IRS"]~
                prev$prev_pred_1[prev$Intervention == "PBO_ITN" | prev$Intervention == "PBO_ITN_IRS"]+0))
summary.lm(lm(prev$prev_obs[prev$Intervention == "PBO_ITN" | prev$Intervention == "PBO_ITN_IRS"]~
                prev$prev_pred_2[prev$Intervention == "PBO_ITN" | prev$Intervention == "PBO_ITN_IRS"]+0))
summary.lm(lm(prev$prev_obs[prev$Intervention == "PBO_ITN" | prev$Intervention == "PBO_ITN_IRS"]~
                prev$prev_pred_3[prev$Intervention == "PBO_ITN" | prev$Intervention == "PBO_ITN_IRS"]+0))
summary.lm(lm(prev$prev_obs[prev$Intervention == "PBO_ITN" | prev$Intervention == "PBO_ITN_IRS"]~
                prev$prev_pred_4[prev$Intervention == "PBO_ITN" | prev$Intervention == "PBO_ITN_IRS"]+0))
summary.lm(lm(prev$prev_obs[prev$Intervention == "PBO_ITN" | prev$Intervention == "PBO_ITN_IRS"]~
                prev$prev_pred_5[prev$Intervention == "PBO_ITN" | prev$Intervention == "PBO_ITN_IRS"]+0))
summary.lm(lm(prev$prev_obs[prev$Intervention == "PBO_ITN" | prev$Intervention == "PBO_ITN_IRS"]~
                prev$prev_pred_6[prev$Intervention == "PBO_ITN" | prev$Intervention == "PBO_ITN_IRS"]+0))



#############################################
##
##
## which studies are we poor at predicting?

## (at this point (=0) all data would be perfectly predicted)
effi$dist_from_perfect = abs(effi$Pred_2_eff - effi$Obs_mean_eff)
prev$dist_from_perfect = abs(prev$prev_pred_2 - prev$prev_obs)

effi[c(which(effi$dist_from_perfect > 0.2)),2:4]
prev[c(which(prev$dist_from_perfect > 0.2)),c(2,3,29)]


