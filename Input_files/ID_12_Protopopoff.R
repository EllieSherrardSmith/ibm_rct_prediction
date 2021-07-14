###################################
##
## input file for study ID 13 Protopopoff


## Protopopoff et al 2018

## 1. Pyrethroid resistance
## Estimated from the susceptibility bioassay test
## and the ratio of mosquito species at the trial site
## (each arm independently where data are available)

## Species ratio
an_gam = 0.916
an_fun = 0.038
an_ara = 0.046

## Species bioassay result (from % mortality)

an_gam_bio = 0.078
an_fun_bio = 0.545
an_ara_bio = 0.078 ## without specific data on arabiensis we use gamb sl

## Additional data from Natacha email: Wed 6/3/2020 9:53 AM
## Species bioassay result (from % mortality)

an_gam_bio = c(0.28+0.04)/2 ##lambdacyhalothrin and permethrin mortality
an_fun_bio = c(0.43+0.31)/2 ##lambdacyhalothrin and permethrin mortality
an_ara_bio = c(0.28+0.04)/2 ## without specific data on arabiensis we use gamb sl

## we do not have sufficient data to parameterise the model distinctly for species
## instead, we use a global weighted estimate for pyrethroid resistance

pyrethroid_resistance = (an_gam * (1 - an_gam_bio)) + 
  (an_fun * (1- an_fun_bio)) + 
  (an_ara * (1 - an_ara_bio))

pyrethroid_resistance

## res_istance params LLIN 
Act = read.csv("data/intervention_data/Actellic_uncertainty.csv",header=TRUE)


## Add in uncertainty for other params estimated from statistical work elsewhere

## itn_leave_dur
# dim(params_X1)

## Waning net usage 
library(rstan)
library(adegenet)

## DATA FROM MARK ROWLAND AND DAVID BATH (Unpublished)
## Different to VCAG report and not available for Actellic arms
# Survey A, 4 months	69.9 (61.1-77.5), 1080	74.6 (67.5-80.7), 1084
# Survey B, 9 months	81.0 (75.2-85.7), 1048	76.9 (72.0-81.2), 983
# Survey C,  16 months	66.2 (56.4-74.7), 1034	55.4 (47.5-63.0), 984
# SurveyD, 21 months	56.3 (51.3-61.3), 1044	44.3 (38.8-49.8), 958
# SurveyE, 28 months	54.6 (48.9-60.1), 1006	45.2 (39.1-51.5), 891
# SurveyF, 33 months	31.3 (28.5-34.4), 876	31.9 (25.3-39.3), 900

time_obs = c(4/12,9/12,16/12,21/12,28/12,33/12)
standard_net_usage = c(0.699,0.81,0.662,0.563,0.546,0.313)
standard_net_usage_min = c(0.611,0.752,0.564,0.513,0.489,0.285)
standard_net_usage_max = c(0.775,0.857,0.747,0.613,0.601,0.344)

pbo_net_usage = c(0.746,0.769,0.554,0.443,0.452,0.319)
pbo_net_usage_min = c(0.675,0.72,0.475,0.388,0.391,0.253)
pbo_net_usage_max = c(0.807,0.812,0.63,0.498,0.515,0.393)


y_standard_net_usage = log(standard_net_usage)

y_pbo_net_usage = log(pbo_net_usage)

time_m = seq(0,3,0.01)
# 

stanmodelcode <- "
data {
int<lower=0> N;
int<lower=0> N2; //the size of the new_X matrix
vector[N] y;
vector[N] x;
vector[N2] New_x;

}
parameters {
real beta0;
real beta1;
real sigma;
}
transformed parameters {
vector[N] m;
m = beta0 + beta1 * x;
} 
model {
// priors
beta0 ~ cauchy(0, 10); 
beta1 ~ cauchy(0, 10); 

// likelihood
y ~ normal(m, sigma);   
}
generated quantities {
vector[N2] y_pred;
y_pred = beta0 + beta1 * New_x; //the y values predicted by the model
}
"
stanDso <- stan_model(model_code = stanmodelcode) 

dat_standard <- list(N = length(time_obs), 
                     N2 = length(seq(0,3,0.01)),
                     y = y_standard_net_usage, 
                     x = time_obs,
                     New_x = seq(0,3,0.01)); 
fit <- sampling(stanDso, data = dat_standard, iter = 10000, warmup=2000) 
fit

#plotting the posterior distribution for the parameters
post_beta<-As.mcmc.list(fit,pars="beta0")
plot(post_beta)

## gradient is fit to the data for alpha
## standard_net_usage ~ exp(-alpha*time_obs)
b0 <- extract(fit, 'beta0')
b0<- unlist(b0, use.names=FALSE)
b1 <- extract(fit, 'beta1')
b1<- unlist(b1, use.names=FALSE)

y_predicted_stn_exp = exp(mean(b0)) * exp(mean(b1)*time_m)
y_predicted_stn_exp_min = exp(quantile(b0,0.45)) * exp(quantile(b1,0.25)*time_m)
y_predicted_stn_exp_max = exp(quantile(b0,0.55)) * exp(quantile(b1,0.75)*time_m)


dat_pbo <- list(N = length(time_obs), 
                N2 = length(seq(0,3,0.01)),
                y = y_pbo_net_usage, 
                x = time_obs,
                New_x = seq(0,3,0.01)); 

fit2 <- sampling(stanDso, data = dat_pbo, iter = 10000, warmup=2000) 
fit2
## Inference for Stan model: stanmodelcode.
# 4 chains, each with iter=10000; warmup=2000; thin=1; 
# post-warmup draws per chain=8000, total post-warmup draws=32000.
# 
# mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
# beta0       -0.11    0.00 0.17 -0.43 -0.19 -0.11 -0.03  0.23  5895    1
# beta1       -0.35    0.00 0.10 -0.54 -0.40 -0.35 -0.31 -0.17  6031    1
# sigma        0.17    0.00 0.11  0.07  0.11  0.14  0.20  0.46  3265    1
# 
# 
# Samples were drawn using NUTS(diag_e) at Wed Jan 23 11:17:41 2019.
# For each parameter, n_eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor on split chains (at convergence, Rhat=1).


#plotting the posterior distribution for the parameters
# post_beta<-As.mcmc.list(fit,pars="beta0")
# plot(post_beta)

## gradient is fit to the data for alpha
## standard_net_usage ~ exp(-alpha*time_obs)
b0_2 <- extract(fit2, 'beta0')
b0_2<- unlist(b0_2, use.names=FALSE)
b1_2 <- extract(fit2, 'beta1')
b1_2<- unlist(b1_2, use.names=FALSE)

## Translate back to exponential scale
y_predicted_pbo_exp = exp(mean(b0_2)) * exp(mean(b1_2)*time_m)

y_predicted_pbo_exp_min = exp(quantile(b0_2,0.45)) * exp(quantile(b1_2,0.25)*time_m)
y_predicted_pbo_exp_max = exp(quantile(b0_2,0.55)) * exp(quantile(b1_2,0.75)*time_m)

## Final output
par(mfrow = c(1,1))
plot(standard_net_usage ~ time_obs,ylab="Proportion of children using net (%)",
     xlab="Time in years",yaxt="n",ylim=c(0,1),cex.lab=1.6,cex.axis=1.4,xlim=c(0,3))
axis(2,las=2, at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.6,cex.axis=1.4)

# lines(net_usage_s ~ time_obs)
for(i in 1:length(standard_net_usage)){
  segments(x0=time_obs[i],x1 = time_obs[i], y0=standard_net_usage_min[i],y1=standard_net_usage_max[i])  
}
offset_time_obs = time_obs+0.04
points(pbo_net_usage ~ offset_time_obs,col="blue",pch=20)
for(i in 1:length(standard_net_usage)){
  segments(x0=offset_time_obs[i],x1 = offset_time_obs[i], y0=pbo_net_usage_min[i],y1=pbo_net_usage_max[i],col="blue",lty=2)  
}
lines(y_predicted_stn_exp ~ time_m,col="black",lty=2,lwd=2)
lines(y_predicted_pbo_exp ~ time_m,col="blue",lty=2,lwd=2)
polygon(c(time_m,rev(time_m)),c(y_predicted_stn_exp_min,rev(y_predicted_stn_exp_max)),border=NA,col=transp("grey","0.5"))
polygon(c(time_m,rev(time_m)),c(y_predicted_pbo_exp_min,rev(y_predicted_pbo_exp_max)),border=NA,col=transp("blue","0.5"))

parms_usage = data.frame(itn_leave_dur_standardLLIN = sample(b1[b1 >= quantile(b1,0.25) & b1 <= quantile(b1,0.75)],1000,replace=FALSE),
                         itn_leave_dur_standardLLIN = sample(b1_2[b1_2 >= quantile(b1_2,0.25) & b1_2 <= quantile(b1_2,0.75)],1000,replace=FALSE)) ##

parms_usage$itn_leave_dur_standardLLIN_bt = -1/parms_usage$itn_leave_dur_standardLLIN
parms_usage$itn_leave_dur_PBOLLIN_bt = -1/parms_usage$itn_leave_dur_standardLLIN.1

-1/colMeans(parms_usage)
-1/min(parms_usage[,1]);-1/max(parms_usage[,1])
-1/min(parms_usage[,2]);-1/max(parms_usage[,2])

parms_usage1 = sample(x = parms_usage[,3],size = 1000, replace = FALSE)
parms_usage2 = sample(x = parms_usage[,4],size = 1000, replace = FALSE)




## Protopopoff [X1] Standard LLIN
## create a parameter list of any params we can estimate uncertainty for:
params = c("base_prev","usage_t0","resistance")
param_best_est = c(0.68,0.391,0.904) ## Original
param_best_est = c(0.68,0.391,0.83202) ## Updated from the email from Natacha Wed 6/3/2020 (see above)
# param_min_est = c(0.61,0.358,0.80) ## We can discuss setting these!
# param_max_est = c(0.74,0.423,0.98) ## 

params_X1 = array(dim=c(1000,length(params)+16))
for(i in 1:length(params)){
  params_X1[,i] = uncertainty_fn(param_best_est[i])  
}

## Setting up for standard Pyrethroid LLINs
params_X1[,4] = parms_usage1

## res_istance params LLIN 
## Have options to estimate this here
## i think we should go back to the bayesian fits for impact and then track the uncertainty all the way through

## logistic all data
uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 3.571,param_b = 0.703,
                                                               param_c = 2.444,param_d = 0.422,param_e = 0.356,
                                                               param_f = 4.657,param_g = 0.040) ## and a range of resistance estimates
## log-logistic all data
uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 0.765,param_b = 0.468,
                                                               param_c = 2.444,param_d = 0.422,param_e = 0.356,
                                                               param_f = 4.657,param_g = 0.04) ## and a range of resistance estimates
## log-logistic west
uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 0.668,param_b = 0.870,
                                                               param_c = 3.045,param_d = 0.573,param_e = 0.39,
                                                               param_f = 5.129,param_g = 0.036) ## and a range of resistance estimates
## log-logistic east
uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 0.668,param_b = 0.870,
                                                               param_c = 1.546,param_d = 0.043,param_e = 0.454,
                                                               param_f = 3.43,param_g = 0.011)## and a range of resistance estimates

##west logistic
uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 3.784,param_b = 0.727,
                                                               param_c = 3.045,param_d = 0.573,param_e = 0.39,
                                                               param_f = 5.129,param_g = 0.036)
## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 3.28,param_b = 0.648,
                                                               param_c = 1.546,param_d = 0.043,param_e = 0.454,
                                                               param_f = 3.43,param_g = 0.011) ## and a range of resistance estimates


params_X1[,5] = uncertainty_resistance_params_nets[,2] ## itn_repel_fun
params_X1[,6] = uncertainty_resistance_params_nets[,2] ## itn_repel_gam
params_X1[,7] = uncertainty_resistance_params_nets[,2] ## itn_repel_ara

params_X1[,8] = uncertainty_resistance_params_nets[,1] ## itn_kill_fun
params_X1[,9] = uncertainty_resistance_params_nets[,1] ## itn_kill_gam
params_X1[,10] = uncertainty_resistance_params_nets[,1] ## itn_kill_ara

params_X1[,11] = uncertainty_resistance_params_nets[,3] ## itn_half_life

#Actellic base only
params_X1[,12] = Act$irs_decay_mort1
params_X1[,13] = Act$irs_decay_mort2
params_X1[,14] = Act$irs_decay_succ1
params_X1[,15] = Act$irs_decay_succ2
params_X1[,16] = Act$irs_decay_det1
params_X1[,17] = Act$irs_decay_det2 ## we use the same decay for mortality to fit deterrence

params_X1[,18] = rep(0,1000) ## IRS cover (none in X1 Standard LLIN only)
params_X1[,19] = 1
head(params_X1)


colnames(params_X1) = c("base_prev","usage_t0","resistance",
                        "itn_leave_dur",
                        "itn_repel_fun","itn_repel_gam","itn_repel_ara",
                        "itn_kill_fun","itn_kill_gam","itn_kill_ara",
                        "itn_half_life",
                        "irs_decay_mort1","irs_decay_mort2",
                        "irs_decay_succ1","irs_decay_succ2",
                        "irs_decay_det1","irs_decay_det2","irs_cover","run")


## Protopopoff [X2]PBO LLIN
## create a parameter list of any params we can estimate uncertainty for:
params = c("base_prev","usage_t0","resistance")
param_best_est = c(0.61,0.318,0.83202)
param_min_est = c(0.55,0.289,0.80) ## We can discuss setting these!
param_max_est = c(0.68,0.348,0.98) ## 

params_X2 = array(dim=c(1000,length(params)+16))
for(i in 1:length(params)){
  params_X2[,i] = uncertainty_fn(param_best_est[i])  
}

## Add in uncertainty for other params estimated from statistical work elsewhere
uncertainty_resistance_params_nets = resistance_ITN_params_f(is.pbo = 1, ## test for standard nets
                                                             species = 1, ## generic mosquitos
                                                             resistance_range = uncertainty_fn(pyrethroid_resistance)) ## and a range of resistance estimates


uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 1, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 3.57,param_b = 0.70,
                                                               param_c = 2.57,param_d = 0.49,param_e = 0.36,
                                                               param_f = 4.66,param_g = 0.04) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 1, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 0.89,param_b = 0.47,
                                                               param_c = 2.57,param_d = 0.49,param_e = 0.36,
                                                               param_f = 4.66,param_g = 0.04) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 1, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 0.72,param_b = 0.88,
                                                               param_c = 3.13,param_d = 0.66,param_e = 0.39,
                                                               param_f = 5.13,param_g = 0.04) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 1, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 1.05,param_b = 0.36,
                                                               param_c = 1.76,param_d = 0.07,param_e = 0.45,
                                                               param_f = 3.43,param_g = 0.01)## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 1, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 3.78,param_b = 0.73,
                                                               param_c = 3.13,param_d = 0.66,param_e = 0.39,
                                                               param_f = 5.13,param_g = 0.04) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 1, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 3.29,param_b = 0.66,
                                                               param_c = 1.76,param_d = 0.07,param_e = 0.45,
                                                               param_f = 3.43,param_g = 0.01) ## and a range of resistance estimates

## itn_leave_dur
params_X2[,4] = parms_usage2

params_X2[,5] = uncertainty_resistance_params_nets[,2] ## itn_repel_fun
params_X2[,6] = uncertainty_resistance_params_nets[,2] ## itn_repel_gam
params_X2[,7] = uncertainty_resistance_params_nets[,2] ## itn_repel_ara

params_X2[,8] = uncertainty_resistance_params_nets[,1] ## itn_kill_fun
params_X2[,9] = uncertainty_resistance_params_nets[,1] ## itn_kill_gam
params_X2[,10] = uncertainty_resistance_params_nets[,1] ## itn_kill_ara

params_X2[,11] = uncertainty_resistance_params_nets[,3] ## itn_half_life


#Actellic base only
params_X2[,12] = Act$irs_decay_mort1
params_X2[,13] = Act$irs_decay_mort2
params_X2[,14] = Act$irs_decay_succ1
params_X2[,15] = Act$irs_decay_succ2
params_X2[,16] = Act$irs_decay_det1
params_X2[,17] = Act$irs_decay_det2 ## we use the same decay for mortality to fit deterrence

params_X2[,18] = rep(0,1000) ## IRS cover (none in X1 Standard LLIN only)
params_X2[,19] = 2
head(params_X2)


colnames(params_X2) = c("base_prev","usage_t0","resistance",
                        "itn_leave_dur",
                        "itn_repel_fun","itn_repel_gam","itn_repel_ara",
                        "itn_kill_fun","itn_kill_gam","itn_kill_ara",
                        "itn_half_life",
                        "irs_decay_mort1","irs_decay_mort2",
                        "irs_decay_succ1","irs_decay_succ2",
                        "irs_decay_det1","irs_decay_det2","irs_cover","run")


## Protopopoff [X3] Standard LLIN + Actellic IRS
## create a parameter list of any params we can estimate uncertainty for:
params = c("base_prev","usage_t0","resistance")
param_best_est = c(0.67,0.309,0.83202)
param_min_est = c(0.61,0.281,0.80) ## We can discuss setting these!
param_max_est = c(0.74,0.338,0.98) ## 

params_X3 = array(dim=c(1000,length(params)+16))
for(i in 1:length(params)){
  params_X3[,i] = uncertainty_fn(param_best_est[i])  
}

## Add in uncertainty for other params estimated from statistical work elsewhere

## itn_leave_dur
dim(params_X3)
params_X3[,4] = parms_usage1
#sample(x = c(parms_usage[,1]),size = 1000,replace = FALSE)

uncertainty_resistance_params_nets = resistance_ITN_params_f(is.pbo = 0, ## test for standard nets
                                                             species = 1, ## generic mosquitos
                                                             resistance_range = uncertainty_fn(pyrethroid_resistance)) ## and a range of resistance estimates


uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 3.57,param_b = 0.70,
                                                               param_c = 2.57,param_d = 0.49,param_e = 0.36,
                                                               param_f = 4.66,param_g = 0.04) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 0.89,param_b = 0.47,
                                                               param_c = 2.57,param_d = 0.49,param_e = 0.36,
                                                               param_f = 4.66,param_g = 0.04) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 0.72,param_b = 0.88,
                                                               param_c = 3.13,param_d = 0.66,param_e = 0.39,
                                                               param_f = 5.13,param_g = 0.04) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 1.05,param_b = 0.36,
                                                               param_c = 1.76,param_d = 0.07,param_e = 0.45,
                                                               param_f = 3.43,param_g = 0.01) ## and a range of resistance estimates


params_X3[,5] = uncertainty_resistance_params_nets[,2] ## itn_repel_fun
params_X3[,6] = uncertainty_resistance_params_nets[,2] ## itn_repel_gam
params_X3[,7] = uncertainty_resistance_params_nets[,2] ## itn_repel_ara

params_X3[,8] = uncertainty_resistance_params_nets[,1] ## itn_kill_fun
params_X3[,9] = uncertainty_resistance_params_nets[,1] ## itn_kill_gam
params_X3[,10] = uncertainty_resistance_params_nets[,1] ## itn_kill_ara

params_X3[,11] = uncertainty_resistance_params_nets[,3] ## itn_half_life


#Actellic base only
params_X3[,12] = Act$irs_decay_mort1
params_X3[,13] = Act$irs_decay_mort2
params_X3[,14] = Act$irs_decay_succ1
params_X3[,15] = Act$irs_decay_succ2
params_X3[,16] = Act$irs_decay_det1
params_X3[,17] = Act$irs_decay_det2 ## we use the same decay for mortality to fit deterrence

params_X3[,18] = uncertainty_fn(0.94) ## IRS cover (none in X1 Standard LLIN only)
params_X3[,19] = 3
head(params_X3)


colnames(params_X3) = c("base_prev","usage_t0","resistance",
                        "itn_leave_dur",
                        "itn_repel_fun","itn_repel_gam","itn_repel_ara",
                        "itn_kill_fun","itn_kill_gam","itn_kill_ara",
                        "itn_half_life",
                        "irs_decay_mort1","irs_decay_mort2",
                        "irs_decay_succ1","irs_decay_succ2",
                        "irs_decay_det1","irs_decay_det2","irs_cover","run")


## Protopopoff [X4]PBO LLIN + IRS
## create a parameter list of any params we can estimate uncertainty for:
params = c("base_prev","usage_t0","resistance")
param_best_est = c(0.64,0.316,0.83202)
param_min_est = c(0.55,0.287,0.80) ## We can discuss setting these!
param_max_est = c(0.68,0.346,0.98) ## 

params_X4 = array(dim=c(1000,length(params)+16))
for(i in 1:length(params)){
  params_X4[,i] = uncertainty_fn(param_best_est[i])  
}

## Add in uncertainty for other params estimated from statistical work elsewhere

## itn_leave_dur
params_X4[,4] = parms_usage2

uncertainty_resistance_params_nets = resistance_ITN_params_f(is.pbo = 1, ## test for standard nets
                                                             species = 1, ## generic mosquitos
                                                             resistance_range = uncertainty_fn(pyrethroid_resistance)) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 1, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 3.57,param_b = 0.70,
                                                               param_c = 2.57,param_d = 0.49,param_e = 0.36,
                                                               param_f = 4.66,param_g = 0.04) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 1, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 0.89,param_b = 0.47,
                                                               param_c = 2.57,param_d = 0.49,param_e = 0.36,
                                                               param_f = 4.66,param_g = 0.04) ## and a range of resistance estimates


uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 1, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 0.72,param_b = 0.88,
                                                               param_c = 3.13,param_d = 0.66,param_e = 0.39,
                                                               param_f = 5.13,param_g = 0.04) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 1, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 1.05,param_b = 0.36,
                                                               param_c = 1.76,param_d = 0.07,param_e = 0.45,
                                                               param_f = 3.43,param_g = 0.01) ## and a range of resistance estimates


params_X4[,5] = uncertainty_resistance_params_nets[,2] ## itn_repel_fun
params_X4[,6] = uncertainty_resistance_params_nets[,2] ## itn_repel_gam
params_X4[,7] = uncertainty_resistance_params_nets[,2] ## itn_repel_ara

params_X4[,8] = uncertainty_resistance_params_nets[,1] ## itn_kill_fun
params_X4[,9] = uncertainty_resistance_params_nets[,1] ## itn_kill_gam
params_X4[,10] = uncertainty_resistance_params_nets[,1] ## itn_kill_ara

params_X4[,11] = uncertainty_resistance_params_nets[,3] ## itn_half_life

#Actellic base only
params_X4[,12] = Act$irs_decay_mort1
params_X4[,13] = Act$irs_decay_mort2
params_X4[,14] = Act$irs_decay_succ1
params_X4[,15] = Act$irs_decay_succ2
params_X4[,16] = Act$irs_decay_det1
params_X4[,17] = Act$irs_decay_det2 ## we use the same decay for mortality to fit deterrence

params_X4[,18] = uncertainty_fn(0.94) ## IRS cover (none in X1 Standard LLIN only)
params_X4[,19] = 4
head(params_X4)


colnames(params_X4) = c("base_prev","usage_t0","resistance",
                        "itn_leave_dur",
                        "itn_repel_fun","itn_repel_gam","itn_repel_ara",
                        "itn_kill_fun","itn_kill_gam","itn_kill_ara",
                        "itn_half_life",
                        "irs_decay_mort1","irs_decay_mort2",
                        "irs_decay_succ1","irs_decay_succ2",
                        "irs_decay_det1","irs_decay_det2","irs_cover","run")


params_all = rbind(params_X1,params_X2,params_X3,params_X4)

params_all = rbind(params_X1,params_X2)
# write.csv(params_all,"Input_files/data/input_files_13_update.csv")
# write.csv(params_all,"Input_files/data/input_files_13_update_net_parameters1_all_hutdata.csv")
# write.csv(params_all,"Input_files/data/input_files_13_update_net_parameters2_all_hutdata.csv")
# write.csv(params_all,"Input_files/data/input_files_13_update_net_parameters3_all_hutdata.csv")
# write.csv(params_all,"Input_files/data/input_files_13_update_net_parameters4_all_hutdata.csv")

