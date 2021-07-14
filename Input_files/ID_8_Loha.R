###################################
##
## input file for study ID 9 Loha

# Re: An. pharoensis
# 
# From: Neil Lobo (17.33 18/07/2019)
# Pharoensis is outdoor biting and usually bites animals though we find a lot bite humans. A small number of a large population that bite humans can cause a lot of transmission.
# Also, a lot of the normal traps dont catch representative numbers (HLCs, CDC-LTs)…so if you catch a small number there are usually a LOT out there.
# 
# Zeimanni is part of the coustani complex and we have found 3-4 species in addition to ziemanni and tenebrous - all of which usually get lumped into coustani In Ethiopia (data being analyzed) we have 4(molecular)  species in the coustani group all with specific behaviors. 3 of these are found all over africa and all have had sporozoites.
# It is hard to parameterize them as they have regional differences in behavior… but overall, none are hugely indoor biting… they lean to outdoor biting but some will enter more than others I have data on these from Ethiopia, Kenya, Tanzania, DRC, Zim, Zambia, West Africa…..
# 
# ## Parameterising these as very outdoor biting (30%) i.e. phiI = 30%, phiB = 20%
# 
# 
# #########################
# ##
# ##   Gari paper
# ##
# ############################
# An arabiensis 71.1%: Q0 = 0.59, phi_I = 0.86, phi_B = 0.8 ## average age 14 days, reduce mosquito mortality 0.071 (?)
# An. pharoensis 21.1%: Q0 = 0.89, phi_I = 0.3, phi_B = 0.28 ## average age 14 days, reduce mosquito mortality 0.071 (?)
# An. zeimanni/An. funestus sl 5.1%.: Q0 = 0.59, phi_I = 0.86, phi_B = 0.8 ## average age 14 days, reduce mosquito mortality 0.071 (?)
# 
# 33/349 p. falciparum positive at baseline = 0.095 all age prevalence
# 23.4 cases per 10,000 person-weeks = 0.12168 cases per person per year
# 
# ## Pre trial net use = 11%

## Species ratio
an_pha = 0.211
an_ara = 0.711
an_oth = 1 - an_pha - an_ara

## Species bioassay result (from % mortality)

an_pha_bio = 1
an_ara_bio = 0.05
an_oth_bio = 1 ## without specific data on arabiensis we use gamb sl

## we do not have sufficient data to parameterise the model distinctly for species
## instead, we use a global weighted estimate for pyrethroid resistance

pyrethroid_resistance = (an_pha * (1 - an_pha_bio)) + 
  (an_oth * (1- an_oth_bio)) + 
  (an_ara * (1 - an_ara_bio))

pyrethroid_resistance




## Specific data needed to estimate parameters 
## res_istance params LLIN 
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

## create a parameter list of any params we can estimate uncertainty for:
params = c("base_prev","usage_t0","resistance")
param_best_est = c(33/349,0.49,0.67545)


params_X1 = expand.grid(runs = 1:1000)


params_X1[,1] = uncertainty_resistance_params_nets[,2]
params_X1[,2] = uncertainty_resistance_params_nets[,2]
params_X1[,3] = uncertainty_resistance_params_nets[,2]
params_X1[,4] = uncertainty_resistance_params_nets[,1]
params_X1[,5] = uncertainty_resistance_params_nets[,1]
params_X1[,6] = uncertainty_resistance_params_nets[,1]
params_X1[,7] = uncertainty_resistance_params_nets[,3]

params_X1[,8] = uncertainty_fn(param_best_est[1])## estimated base_prev
params_X1[,9] = uncertainty_fn(param_best_est[2])## usage_t0
params_X1[,10] = uncertainty_fn(param_best_est[3])## resistance

head(params_X1)


colnames(params_X1) = c("itn_repel_fun","itn_repel_gam","itn_repel_ara",
                        "itn_kill_fun","itn_kill_gam","itn_kill_ara",
                        "itn_half_life",
                        "estimated_baseline_prev",
                        "ITN_use_arm_2",
                        "resistance")


## 2. Waning net usage 
library(rstan)
library(adegenet)


options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

time_obs = c(13/52,39/52,66/52, 92/52,13/52,39/52,66/52, 92/52)
standard_net_usage = c(0.47,0.26,0.08,0.01, ##net use 1
                       0.49,0.27,0.06,0.01) ##net use 2

y_standard_net_usage = log(standard_net_usage)

time_m = seq(0,3,0.01)
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
                     N2 = length(seq(0,3,0.01)),
                     y = y_standard_net_usage, 
                     x = time_obs,
                     New_x = seq(0,3,0.01)); 
fit <- sampling(stanDso, data = dat_standard, iter = 3000, warmup=2000) 
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

y_predicted_stn_exp1 = exp(median(b0)) * exp(median(b1)*time_m)

y_predicted_stn_exp_min = exp(quantile(b0,0.20)) * exp(quantile(b1,0.20)*time_m)
y_predicted_stn_exp_max = exp(quantile(b0,0.80)) * exp(quantile(b1,0.80)*time_m)

## Final output
par(mfrow = c(1,1))
plot(standard_net_usage ~ time_obs,ylab="Proportion of all people using net (%)",
     xlab="Time in years",yaxt="n",ylim=c(0,1),cex.lab=1.6,cex.axis=1.4,xlim=c(0,3))
axis(2,las=2, at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.6,cex.axis=1.4)

lines(y_predicted_stn_exp1 ~ time_m,col="black",lty=2,lwd=2)

polygon(c(time_m,rev(time_m)),c(y_predicted_stn_exp_min,rev(y_predicted_stn_exp_max)),border=FALSE,col=transp("grey",0.4))

itn_leave_dur_standardLLIN1 = b1[b1 >= quantile(b1,0.35) & b1 <= quantile(b1,0.65)]

LLIN_IRSbt = -1/itn_leave_dur_standardLLIN1
parms_usage1 = sample(x = LLIN_IRSbt,size = 1000, replace = FALSE)

params_X1$itn_leave_dur = parms_usage1


write.csv(params_X1,"Input_files/data/input_files_9.csv")
write.csv(params_X1,"Input_files/data/input_files_9_net_parameters1_all_hutdata.csv")
write.csv(params_X1,"Input_files/data/input_files_9_net_parameters2_all_hutdata.csv")
