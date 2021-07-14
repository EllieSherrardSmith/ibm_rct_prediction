###################################
##
## input file for study ID 8 Kafy

## Specific data needed to estimate parameters 

## 1. Pyrethroid resistance
## Estimated from the susceptibility bioassay test
## and the ratio of mosquito species at the trial site
## (each arm independently where data are available)

## Species ratio
an_gam = 0
an_fun = 0
an_ara = 1

## Species bioassay result (from % mortality)

an_gam_bio = 0
an_fun_bio = 0
an_ara_bio = 0.6 ## without specific data on arabiensis we use gamb sl

## we do not have sufficient data to parameterise the model distinctly for species
## instead, we use a global weighted estimate for pyrethroid resistance

pyrethroid_resistance = (an_gam * (1 - an_gam_bio)) + 
  (an_fun * (1- an_fun_bio)) + 
  (an_ara * (1 - an_ara_bio))

pyrethroid_resistance


params_estimated = read.csv("data/intervention_data/pyrethroid_IRS_resistance.csv",header=TRUE)[,2:7]
## Level of mortality in the LLIN only are 2012 is 60 and LLIN+IRS arm is 65
## much uncertainty (Figure 1: 30% to 80% for 75%iles of the two trial arms)
## this corersponds to 20% to 70% resistance
params_estimated[80,]#
params_estimated[30,]#

params_kafy = expand.grid(irs_mort_1 = c(params_estimated[30:80,1]))
params_kafy$irs_mort_2 = c(params_estimated[30:80,2])
params_kafy$irs_succ_1 = c(params_estimated[30:80,3])
params_kafy$irs_succ_2 = c(params_estimated[30:80,4])
params_kafy$irs_det_1 = c(params_estimated[30:80,5])
params_kafy$irs_det_2 = c(params_estimated[30:80,6])

params_irs_delta = dplyr::sample_n(params_kafy, 1000,replace=TRUE)
bendio = read.csv("data/intervention_data/Bendiocarb_uncertainty.csv",header=TRUE)

## Adding uncertainty estimates to the model predictions
## res_istance params LLIN 
                                                               ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               ...) ## and a range of resistance estimates
## Kafy [X1] Standard LLIN
## create a parameter list of any params we can estimate uncertainty for:
params = c("base_prev","usage_t0","resistance")
param_best_est = c(0.07,0.82,0.40)

params_X1 = expand.grid(runs = 1:1000)

##working out resistance effects for range 20 to 70 using Final_ERG_parameterisation.r
head(uncertainty_resistance_params_nets)
# uncertainty_resistance_params_nets2 <- uncertainty_resistance_params_nets[sample(nrow(uncertainty_resistance_params_nets)),]

params_X1[,1] = uncertainty_resistance_params_nets[,2]
params_X1[,2] = uncertainty_resistance_params_nets[,2]
params_X1[,3] = uncertainty_resistance_params_nets[,2]
params_X1[,4] = uncertainty_resistance_params_nets[,1]
params_X1[,5] = uncertainty_resistance_params_nets[,1]
params_X1[,6] = uncertainty_resistance_params_nets[,1]
params_X1[,7] = uncertainty_resistance_params_nets[,3]
#Deltamethrin only
params_X1[,8] = params_irs_delta$irs_mort_1
params_X1[,9] = params_irs_delta$irs_mort_2
params_X1[,10] = params_irs_delta$irs_succ_1
params_X1[,11] = params_irs_delta$irs_succ_2
params_X1[,12] = params_irs_delta$irs_det_1
params_X1[,13] = params_irs_delta$irs_mort_2 ## we use the same decay for mortality to fit deterrence

#bendiocarb
params_X1[,14] = bendio$irs_decay_mort1
params_X1[,15] = bendio$irs_decay_mort2
params_X1[,16] = bendio$irs_decay_succ1
params_X1[,17] = bendio$irs_decay_succ2
params_X1[,18] = bendio$irs_decay_det1
params_X1[,19] = bendio$irs_decay_det2 ## we use the same decay for mortality to fit deterrence



params_X1[,20] = uncertainty_fn(0.82)

params_X1[,21] = uncertainty_fn(0.09)
params_X1[,22] = uncertainty_fn(0.01)
params_X1[,23] = uncertainty_fn(0.04)

head(params_X1)


colnames(params_X1) = c("itn_repel_fun","itn_repel_gam","itn_repel_ara",
                        "itn_kill_fun","itn_kill_gam","itn_kill_ara",
                        "itn_half_life",
                        "irs_decay_mort1","irs_decay_mort2",
                        "irs_decay_succ1","irs_decay_succ2",
                        "irs_decay_det1","irs_decay_det2",
                        "irs_decay_mort1_2","irs_decay_mort2_2",
                        "irs_decay_succ1_2","irs_decay_succ2_2",
                        "irs_decay_det1_2","irs_decay_det2_2",
                        "itn_1","IRS_1","IRS_2","IRS_3")

params_X2 = params_X1
params_X2[,21] = uncertainty_fn(0.99)
params_X2[,22] = uncertainty_fn(0.82)
params_X2[,23] = uncertainty_fn(0.83)

head(params_X2)

params_all = rbind(params_X1,params_X2)
params_all$arm = rep(c(1,2),each=1000)
params_all$prevalence_baseline = c(uncertainty_fn(0.07),uncertainty_fn(0.10))

# 2. ITN leave duration
# The drop out rate of people using ITNs

library(rstan)
library(adegenet)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

time_obs = c(6/12,18/12,30/12)
standard_net_usage = c(0.82,0.79,0.74)

y_standard_net_usage = log(standard_net_usage)

time_m = seq(0,3,0.01)

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


y_predicted_stn_exp = exp(median(b0)) * exp(median(b1)*time_m)

y_predicted_stn_exp_min = exp(quantile(b0,0.45)) * exp(quantile(b1,0.45)*time_m)
y_predicted_stn_exp_max = exp(quantile(b0,0.55)) * exp(quantile(b1,0.55)*time_m)


## Final output
par(mfrow = c(1,1))
plot(standard_net_usage ~ time_obs,ylab="Proportion of children using net (%)",
     xlab="Time in years",yaxt="n",ylim=c(0,1),cex.lab=1.6,cex.axis=1.4,xlim=c(0,3))
axis(2,las=2, at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.6,cex.axis=1.4)

lines(y_predicted_stn_exp ~ time_m,col="black",lty=2,lwd=2)
polygon(c(time_m,rev(time_m)),c(y_predicted_stn_exp_min,rev(y_predicted_stn_exp_max)),border=FALSE,col=transp("grey",0.4))

parms_usage = data.frame(itn_leave_dur_standardLLIN = b1[b1 >= quantile(b1,0.45) & b1 <= quantile(b1,0.55)]) ##

parms_usage$itn_leave_dur_standardLLIN_bt = -1/parms_usage$itn_leave_dur_standardLLIN

parms_usage1 = sample(x = parms_usage[,2],size = 1000, replace = FALSE)

params_all$parms_usage1 = parms_usage1  




# write.csv(params_all,"Input_files/data/input_files_8.csv")
# write.csv(params_all,"Input_files/data/input_files_8_net_parameter1_all_hutdata.csv")
# write.csv(params_all,"Input_files/data/input_files_8_net_parameter2_all_hutdata.csv")
# write.csv(params_all,"Input_files/data/input_files_8_net_parameter3_all_hutdata.csv")
# write.csv(params_all,"Input_files/data/input_files_8_net_parameter4_all_hutdata.csv")
