###################################
##
## input file for study ID 15 West


## West et al 2014

## 1. Pyrethroid resistance
## Estimated from the susceptibility bioassay test
## and the ratio of mosquito species at the trial site
## (each arm independently where data are available)

## Species ratio
an_gam = 0.8
an_fun = 0
an_ara = 0.2

## Species bioassay result (from % mortality)

an_gam_bio = 0.11
an_fun_bio = 0
an_ara_bio = 1 ## without specific data on arabiensis we use gamb sl

## we do not have sufficient data to parameterise the model distinctly for species
## instead, we use a global weighted estimate for pyrethroid resistance

pyrethroid_resistance = (an_gam * (1 - an_gam_bio)) + 
  (an_fun * (1- an_fun_bio)) + 
  (an_ara * (1 - an_ara_bio))

pyrethroid_resistance


## waning net use
time_obs = c(1/12,10/12,14/12,18/12,1/12,10/12,14/12,18/12,1/12,10/12,14/12,18/12)
standard_net_usage = c(0.533,0.466,0.407,0.360,
                       0.482,0.417,0.347,0.298,
                       0.583,0.516,0.47,0.426)
# standard_net_usage_min = c(0.482,0.417,0.347,0.298)
# standard_net_usage_max = c(0.583,0.516,0.47,0.426)

standard_net_usage_irs = c(0.582,0.530,0.441,0.361,
                           0.538,0.475,0.392,0.31,
                           0.625,0.583,0.492,0.415)
# standard_net_usage_irs_min = c(0.538,0.475,0.392,0.31)
# standard_net_usage_irs_max = c(0.625,0.583,0.492,0.415)


y_standard_net_usage = log(standard_net_usage)
# lm(y_standard_net_usage ~ time_obs)
# y_pred = exp(-0.07724) * exp(-0.30908*time_m)
# lines(y_pred ~ time_m)
## gradient is fit to the data for alpha
## standard_net_usage ~ exp(-alpha*time_obs)

y_net_irs_usage = log(standard_net_usage_irs)
# lm(y2~ time_obs)
## gradient is fit to the data for alpha
## standard_net_usage ~ exp(-alpha*time_obs)

time_m = seq(0,2,0.01)
# y_pred = exp(-0.1090) * exp(-0.3506*time_m)
# lines(y_pred ~ time_m)
## 

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
                     N2 = length(seq(0,2,0.01)),
                     y = y_standard_net_usage, 
                     x = time_obs,
                     New_x = seq(0,2,0.01)); 
fit <- sampling(stanDso, data = dat_standard, iter = 6000, warmup=2000) 
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

## Translate back to exponential scale
y_predicted_stn_exp = exp(mean(b0)) * exp(mean(b1)*time_m)

y_predicted_stn_exp_min = exp(quantile(b0,0.45)) * exp(quantile(b1,0.25)*time_m)
y_predicted_stn_exp_max = exp(quantile(b0,0.55)) * exp(quantile(b1,0.75)*time_m)


dat_pbo <- list(N = length(time_obs), 
                N2 = length(seq(0,2,0.01)),
                y = y_net_irs_usage, 
                x = time_obs,
                New_x = seq(0,2,0.01)); 

fit2 <- sampling(stanDso, data = dat_pbo, iter = 6000, warmup=2000) 
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

y_predicted_pbo_exp_min = exp(quantile(b0_2,0.2)) * exp(quantile(b1_2,0.2)*time_m)
y_predicted_pbo_exp_max = exp(quantile(b0_2,0.8)) * exp(quantile(b1_2,0.8)*time_m)

## Final output
par(mfrow = c(1,1))
plot(standard_net_usage ~ time_obs,ylab="Proportion of children using net (%)",
     xlab="Time in years",yaxt="n",ylim=c(0,1),cex.lab=1.6,cex.axis=1.4,xlim=c(0,2))
axis(2,las=2, at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.6,cex.axis=1.4)

# lines(net_usage_s ~ time_obs)
for(i in 1:length(standard_net_usage)){
  segments(x0=time_obs[i],x1 = time_obs[i], y0=standard_net_usage_min[i],y1=standard_net_usage_max[i])  
}
offset_time_obs = time_obs+0.04
points(standard_net_usage_irs ~ offset_time_obs,col="blue",pch=20)
# for(i in 1:length(standard_net_usage_irs)){
#   segments(x0=offset_time_obs[i],x1 = offset_time_obs[i], y0=standard_net_usage_irs_min[i],y1=standard_net_usage_irs_max[i],col="blue",lty=2)  
# }
lines(y_predicted_stn_exp ~ time_m,col="black",lty=2,lwd=2)
lines(y_predicted_pbo_exp ~ time_m,col="blue",lty=2,lwd=2)
polygon(c(time_m,rev(time_m)),c(y_predicted_stn_exp_min,rev(y_predicted_stn_exp_max)),border=NA,col=transp("grey","0.5"))
polygon(c(time_m,rev(time_m)),c(y_predicted_pbo_exp_min,rev(y_predicted_pbo_exp_max)),border=NA,col=transp("blue","0.5"))

parms_usage = data.frame(itn_leave_dur_standardLLIN = b1[b1 >= quantile(b1,0.2) & b1 <= quantile(b1,0.8)],
                         itn_leave_dur_standardLLIN = b1_2[b1_2 >= quantile(b1_2,0.2) & b1_2 <= quantile(b1_2,0.8)]) ##

colMeans(parms_usage)
min(parms_usage[,1]);max(parms_usage[,1])
min(parms_usage[,2]);max(parms_usage[,2])

parms_usage$itn_leave_dur_standardLLIN_bt = -1/parms_usage$itn_leave_dur_standardLLIN
parms_usage$itn_leave_dur_LLIN.IRS_bt = -1/parms_usage$itn_leave_dur_standardLLIN.1


parms_usage1_stn = sample(x = parms_usage[,3],size = 1000, replace = FALSE)
parms_usage2_irs = sample(x = parms_usage[,4],size = 1000, replace = FALSE)


Bn = read.csv("data\\intervention_data\\Bendiocarb_uncertainty.csv",header=TRUE)

## Adding uncertainty estimates to the model predictions
#######**************************UP TO HERE!
## West [X1] Standard LLIN
## create a parameter list of any params we can estimate uncertainty for:
params = c("base_prev","usage_t0","resistance")
param_best_est = c(0.103,0.533,0.712)
# param_min_est = c(0.052,0.482,0.662) ## Set min!
# param_max_est = c(0.193,0.583,0.758) ## Set max!

params_X1 = array(dim=c(1000,length(params)+17))
for(i in 1:length(params)){
  params_X1[,i] = uncertainty_fn(param_best_est[i])  
}

## Add in uncertainty for other params estimated from statistical work elsewhere

## itn_leave_dur
params_X1[,4] = parms_usage1_stn

## res_istance params LLIN 
uncertainty_resistance_params_nets = resistance_ITN_params_f(is.pbo = 0, ## test for standard nets
                                                             species = 1, ## generic mosquitos
                                                             resistance_range = params_X1[,3]) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X1[,3],
                                                               param_a = 3.57,param_b = 0.70,
                                                               param_c = 2.57,param_d = 0.49,param_e = 0.36,
                                                               param_f = 4.66,param_g = 0.04) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X1[,3],
                                                               param_a = 0.89,param_b = 0.47,
                                                               param_c = 2.57,param_d = 0.49,param_e = 0.36,
                                                               param_f = 4.66,param_g = 0.04) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X1[,3],
                                                               param_a = 0.72,param_b = 0.88,
                                                               param_c = 3.13,param_d = 0.66,param_e = 0.39,
                                                               param_f = 5.13,param_g = 0.04) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X1[,3],
                                                               param_a = 1.05,param_b = 0.36,
                                                               param_c = 1.76,param_d = 0.07,param_e = 0.45,
                                                               param_f = 3.43,param_g = 0.01) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X1[,3],
                                                               param_a = 3.78,param_b = 0.73,
                                                               param_c = 3.13,param_d = 0.66,param_e = 0.39,
                                                               param_f = 5.13,param_g = 0.04)## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X1[,3],
                                                               param_a = 3.29,param_b = 0.66,
                                                               param_c = 1.76,param_d = 0.07,param_e = 0.45,
                                                               param_f = 3.43,param_g = 0.01) ## and a range of resistance estimates

## 
params_X1[,5] = uncertainty_resistance_params_nets[,2] ## itn_repel_fun
params_X1[,6] = uncertainty_resistance_params_nets[,2] ## itn_repel_gam
params_X1[,7] = uncertainty_resistance_params_nets[,2] ## itn_repel_ara

params_X1[,8] = uncertainty_resistance_params_nets[,1] ## itn_kill_fun
params_X1[,9] = uncertainty_resistance_params_nets[,1] ## itn_kill_gam
params_X1[,10] = uncertainty_resistance_params_nets[,1] ## itn_kill_ara

params_X1[,11] = uncertainty_resistance_params_nets[,3] ## itn_half_life

head(params_X1)

#Bendiocarb### base only
params_X1[,12] = Bn$irs_decay_mort1
params_X1[,13] = Bn$irs_decay_mort2
params_X1[,14] = Bn$irs_decay_succ1
params_X1[,15] = Bn$irs_decay_succ2
params_X1[,16] = Bn$irs_decay_det1
params_X1[,17] = Bn$irs_decay_det2 ## we use the same decay for mortality to fit deterrence

params_X1[,18] = uncertainty_fn(0.033) ## IRS cover 1st spraying (none in X1 Standard LLIN only)
params_X1[,19] = uncertainty_fn(0.052) ## IRS cover 2nd spraying (none in X1 Standard LLIN only)

params_X1[,20] = 1

head(params_X1)


colnames(params_X1) = c("base_prev","usage_t0","resistance",
                        "itn_leave_dur",
                        "itn_repel_fun","itn_repel_gam","itn_repel_ara",
                        "itn_kill_fun","itn_kill_gam","itn_kill_ara",
                        "itn_half_life",
                        "irs_decay_mort1","irs_decay_mort2",
                        "irs_decay_succ1","irs_decay_succ2",
                        "irs_decay_det1","irs_decay_det2",
                        "irs_cover1", "irs_cover2",
                        "arm")


## West [X2]IRS + LLIN
## create a parameter list of any params we can estimate uncertainty for:
params = c("base_prev","usage_t0","resistance")
param_best_est = c(0.084,0.582,0.712)
# param_min_est = c(0.045,0.538,0.662) ## We can discuss setting these!
# param_max_est = c(0.153,0.625,0.758) ## 

params_X2 = array(dim=c(1000,length(params)+17))
for(i in 1:length(params)){
  params_X2[,i] = uncertainty_fn(param_best_est[i])  
}

## Add in uncertainty for other params estimated from statistical work elsewhere

## itn_leave_dur
dim(params_X2)
params_X2[,4] = parms_usage2_irs

## res_istance params LLIN 
uncertainty_resistance_params_nets = resistance_ITN_params_f(is.pbo = 0, ## test for standard nets
                                                             species = 1, ## generic mosquitos
                                                             resistance_range = params_X2[,3]) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X2[,3],
                                                               param_a = 3.57,param_b = 0.70,
                                                               param_c = 2.57,param_d = 0.49,param_e = 0.36,
                                                               param_f = 4.66,param_g = 0.04) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X2[,3],
                                                               param_a = 0.89,param_b = 0.47,
                                                               param_c = 2.57,param_d = 0.49,param_e = 0.36,
                                                               param_f = 4.66,param_g = 0.04) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X2[,3],
                                                               param_a = 0.72,param_b = 0.88,
                                                               param_c = 3.13,param_d = 0.66,param_e = 0.39,
                                                               param_f = 5.13,param_g = 0.04) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X2[,3],
                                                               param_a = 1.05,param_b = 0.36,
                                                               param_c = 1.76,param_d = 0.07,param_e = 0.45,
                                                               param_f = 3.43,param_g = 0.01)## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X2[,3],
                                                               param_a = 3.78,param_b = 0.73,
                                                               param_c = 3.13,param_d = 0.66,param_e = 0.39,
                                                               param_f = 5.13,param_g = 0.04)## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X2[,3],
                                                               param_a = 3.29,param_b = 0.66,
                                                               param_c = 1.76,param_d = 0.07,param_e = 0.45,
                                                               param_f = 3.43,param_g = 0.01) ## and a range of resistance estimates


## res_istance params LLIN 
## Have options to estimate this here
## i think we should go back to the bayesian fits for impact and then track the uncertainty all the way through
## files saved in H:\Ellie\Vector interventions Model validation\
## to do this (except for uncertainty on stand LLIN which i need to think about)
## From ERG Parameterisation Final we could also estimate uncertainty then use the beta params

## sWITCH UP TO LLIN+IRS
## e.g
params_X2[,5] = uncertainty_resistance_params_nets[,2] ## itn_repel_fun
params_X2[,6] = uncertainty_resistance_params_nets[,2] ## itn_repel_gam
params_X2[,7] = uncertainty_resistance_params_nets[,2] ## itn_repel_ara

params_X2[,8] = uncertainty_resistance_params_nets[,1] ## itn_kill_fun
params_X2[,9] = uncertainty_resistance_params_nets[,1] ## itn_kill_gam
params_X2[,10] = uncertainty_resistance_params_nets[,1] ## itn_kill_ara

params_X2[,11] = uncertainty_resistance_params_nets[,3] ## itn_half_life

head(params_X2)

#Bendiocarb base only
params_X2[,12] = Bn$irs_decay_mort1
params_X2[,13] = Bn$irs_decay_mort2
params_X2[,14] = Bn$irs_decay_succ1
params_X2[,15] = Bn$irs_decay_succ2
params_X2[,16] = Bn$irs_decay_det1
params_X2[,17] = Bn$irs_decay_det2 ## we use the same decay for mortality to fit deterrence

params_X2[,18] = uncertainty_fn(0.921) ## IRS cover 1st spraying
params_X2[,19] = uncertainty_fn(0.895) ## IRS cover 2nd spraying

params_X2[,20] = 2
  
head(params_X2)


colnames(params_X2) = c("base_prev","usage_t0","resistance",
                        "itn_leave_dur",
                        "itn_repel_fun","itn_repel_gam","itn_repel_ara",
                        "itn_kill_fun","itn_kill_gam","itn_kill_ara",
                        "itn_half_life",
                        "irs_decay_mort1","irs_decay_mort2",
                        "irs_decay_succ1","irs_decay_succ2",
                        "irs_decay_det1","irs_decay_det2","irs_cover1", "irs_cover2",
                        "arm")

params_all = rbind(params_X1,params_X2)

# write.csv(params_all,"Input_files/data/input_files_15.csv")
# write.csv(params_all,"Input_files/data/input_files_15_net_parameters1_all_hutdata.csv")
# write.csv(params_all,"Input_files/data/input_files_15_net_parameters2_all_hutdata.csv")
# write.csv(params_all,"Input_files/data/input_files_15_net_parameters3_all_hutdata.csv")
# write.csv(params_all,"Input_files/data/input_files_15_net_parameters4_all_hutdata.csv")

