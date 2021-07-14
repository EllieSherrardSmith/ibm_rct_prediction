###################################
##
## input file for study ID 11 Nevill



## Nevill 1996

## 1. Pyrethroid resistance assume none (1993)
## Estimated from the susceptibility bioassay test
## and the ratio of mosquito species at the trial site
## (each arm independently where data are available)

## Species ratio
an_gam = 0.916
an_fun = 0.038
an_ara = 0.046

## Species bioassay result (from % mortality)

an_gam_bio = 1
an_fun_bio = 1
an_ara_bio = 1 ## without specific data on arabiensis we use gamb sl

## we do not have sufficient data to parameterise the model distinctly for species
## instead, we use a global weighted estimate for pyrethroid resistance

pyrethroid_resistance = (an_gam * (1 - an_gam_bio)) + 
  (an_fun * (1- an_fun_bio)) + 
  (an_ara * (1 - an_ara_bio))

pyrethroid_resistance


params = c("base_prev","resistance")
param_best_est = c(0.251,0)
# param_min_est = c(0.61,0.358,0.80) ## We can discuss setting these!
# param_max_est = c(0.74,0.423,0.98) ## 

params_X1 = params_X2 = expand.grid(runs = 1:1000)

## res_istance params LLIN 
uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = rep(0,1000),
                                                               param_a = 3.571,param_b = 0.703,
                                                               param_c = 2.444,param_d = 0.422,param_e = 0.356,
                                                               param_f = 4.657,param_g = 0.040)# and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = rep(0,1000),
                                                               param_a = 0.765,param_b = 0.468,
                                                               param_c = 2.444,param_d = 0.422,param_e = 0.356,
                                                               param_f = 4.657,param_g = 0.04)## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = rep(0,1000),
                                                               param_a = 0.72,param_b = 0.88,
                                                               param_c = 3.13,param_d = 0.66,param_e = 0.39,
                                                               param_f = 5.13,param_g = 0.04) ## and a range of resistance estimates

## .. 
params_X1[,1] = uncertainty_resistance_params_nets[,2]
params_X1[,2] = uncertainty_resistance_params_nets[,2]
params_X1[,3] = uncertainty_resistance_params_nets[,2]
params_X1[,4] = uncertainty_resistance_params_nets[,1]
params_X1[,5] = uncertainty_resistance_params_nets[,1]
params_X1[,6] = uncertainty_resistance_params_nets[,1]
params_X1[,7] = uncertainty_resistance_params_nets[,3]

params_X2[,1] = uncertainty_resistance_params_nets[,2]
params_X2[,2] = uncertainty_resistance_params_nets[,2]
params_X2[,3] = uncertainty_resistance_params_nets[,2]
params_X2[,4] = uncertainty_resistance_params_nets[,1]
params_X2[,5] = uncertainty_resistance_params_nets[,1]
params_X2[,6] = uncertainty_resistance_params_nets[,1]
params_X2[,7] = uncertainty_resistance_params_nets[,3]


## baseline prevalence
params_X1[,8] = uncertainty_fn(0.251)
params_X2[,8] = params_X1[,8]

## ITN use
## net use in U6s is noted, 25% of popn represented by pregnant and U6.
## so population net use for each trial arm is: 
## 1 0.38 * 0.25
## 2 0.65 * 0.25
## 3 0.45 * 0.25
## 4 0.7 * 0.25



## IRS coverage
params_X1[,9] = uncertainty_fn(0.06)
params_X2[,9] = uncertainty_fn(0.71)


colnames(params_X1) = c("itn_repel_fun","itn_repel_gam","itn_repel_ara",
                        "itn_kill_fun","itn_kill_gam","itn_kill_ara",
                        "itn_half_life",
                        "base_prev",
                        "itn_1")
colnames(params_X2) = c("itn_repel_fun","itn_repel_gam","itn_repel_ara",
                        "itn_kill_fun","itn_kill_gam","itn_kill_ara",
                        "itn_half_life",
                        "base_prev",
                        "itn_1")

head(params_X1)
head(params_X2)



## 2. WANING NET USE

time_obs = c(1/12,6/12,18/12,24/12,30/12)
standard_net_usage = c(0.71,0.6248,0.6461,0.6035,0.568)

y_standard_net_usage = log(standard_net_usage)
# lm(y_standard_net_usage ~ time_obs)
# y_pred = exp(-0.07724) * exp(-0.30908*time_m)
# lines(y_pred ~ time_m)
## gradient is fit to the data for alpha
## standard_net_usage ~ exp(-alpha*time_obs)

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


## Final output
par(mfrow = c(1,1))
plot(standard_net_usage ~ time_obs,ylab="Proportion of children using net (%)",
     xlab="Time in years",yaxt="n",ylim=c(0,1),cex.lab=1.6,cex.axis=1.4,xlim=c(0,3))
axis(2,las=2, at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.6,cex.axis=1.4)

lines(y_predicted_stn_exp ~ time_m,col="black",lty=2,lwd=2)
polygon(c(time_m,rev(time_m)),c(y_predicted_stn_exp_min,rev(y_predicted_stn_exp_max)),border=NA,col=transp("grey","0.5"))

parms_usage = data.frame(itn_leave_dur_standardLLIN = b1[b1 >= quantile(b1,0.25) & b1 <= quantile(b1,0.75)]) ##

parms_usage$itn_leave_dur_standardLLIN_bt = -1/parms_usage$itn_leave_dur_standardLLIN
parms_usage1 = sample(x = parms_usage[,2],size = 2000, replace = FALSE)


params_all = rbind(params_X1,params_X2)
params_all$itn_leave_dur = parms_usage1
dim(params_all)
params_all$arm = rep(c(1,2),each=1000)

# write.csv(params_all,"Input_files/data/input_files_11.csv")
# write.csv(params_all,"Input_files/data/input_files_11_net_parameters1_all_hutdata.csv")
# write.csv(params_all,"Input_files/data/input_files_11_net_parameters2_all_hutdata.csv")
# write.csv(params_all,"Input_files/data/input_files_11_net_parameters3_all_hutdata.csv")
write.csv(params_all,"Input_files/data/input_files_11_net_parameters4_all_hutdata.csv")
