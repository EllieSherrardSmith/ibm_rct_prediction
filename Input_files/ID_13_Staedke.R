##########################
##
## Staedke ID 14

## 1. Pyrethroid resistance
## Estimated from the susceptibility bioassay test
## and the ratio of mosquito species at the trial site
## (each arm independently where data are available)

## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5889576/#!po=31.8182
# EASTERN UGANDA 2015
# Deltamethrin: An. funestus = 82.9% (111), (Soroti)
# Deltamethrin: An. gambiae s.l. = 87% (95), (Soroti)
# Permethrin: An. gambiae s.l. = 61.3% (106), (Soroti)

# Permethrin: An. funestus = 20.8% (101), (Tororo)
# Permethrin: An. gambiae s.l. = 67.3% (147), (Tororo)
# Deltamethrin: An. gambiae s.l. = 81.6% (136), (Tororo)

# NORTHERN UGANDA 2015
# Permethrin: An. gambiae s.l. = 14% (100), (Apac)
# Deltamethrin: An. gambiae s.l. = 31.8% (101), (Apac)
# Deltamethrin: An. gambiae s.l. = 9.7% (82), (Apac)

# Permethrin: An. gambiae s.l. = 17.5% (137), (Lira)
# Deltamethrin: An. gambiae s.l. = 59.3% (91), (Lira)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3543752/ (figures not tables)
# EASTERN UGANDA 2011
# Permethrin: An. gambiae s.s. = 24.7% (166)
# Deltamethrin: An. gambiae s.s. = 33.0% (93)

# Permethrin: An. arabiensis = 56.9% (325)
# Deltamethrin: An. arabiensis = 84.6% (208)

## Species ratio
##ARMS 1, 2, 3, 4
#Proportion An. gambiae s.s.-like	1368 were An. gambiae (s.s.)	
#Proportion An. arabiensis-like	88 
#Proportion An. funestus s.s.-like	470

an_gam = 1368/(1368 + 88 + 470)
an_fun = 470/(1368 + 88 + 470)
an_ara = 88/(1368 + 88 + 470)

## Species bioassay result (from % mortality)

an_gam_bio = (0.87 * 95 + 
                0.613 * 106 +
                0.673 * 147 + 
                0.816 * 136 +
                0.14 * 100 + 
                0.318 * 101 +
                0.097 * 82 + 
                0.175 * 137 +
                0.593 * 91+ 
                0.247 * 166 +
                0.33 * 93)/
  (95+106+147+136+100+101+82+137+91 + 166 + 93)

an_fun_bio = (0.829 * 111 + 0.208 * 101) / (111 + 101) 

an_ara_bio = (0.569 * 325 + 0.846 * 208) / (325 + 208)

# ## OR arm1a, go with lowest resistance
# an_gam_bio = 0.097
# an_fun_bio = 0.208
# an_ara_bio = 0.569
# ## we do not have sufficient data to parameterise the model distinctly for species
# ## instead, we use a global weighted estimate for pyrethroid resistance

pyrethroid_resistance = (an_gam * (1 - an_gam_bio)) + 
  (an_fun * (1- an_fun_bio)) + 
  (an_ara * (1 - an_ara_bio))

pyrethroid_resistance


## Waning use of nets

time_obs = c(6/12,12/12,18/12,6/12,12/12,18/12,6/12,12/12,18/12)
standard_net_usage = c(0.854,0.786,0.731)
# Coverage 6-months: 71%
# Coverage 12-months: 63%
# Coverage 18-months: 51.1%
standard_net_usage = c(0.71,0.63,0.49,0.62,0.58,0.44,0.72,0.68,0.54)
y_standard_net_usage = log(standard_net_usage)

pbo_net_usage = c(0.73,0.63,0.511,0.64,0.58,0.461,0.74,0.68,0.561)
y_pbo_net_usage = log(pbo_net_usage)

time_m = seq(0,3,0.01)
# y_pred = exp(-0.1090) * exp(-0.3506*time_m)
# y_pred1 = exp(-0.50) * exp(-0.6*time_m)
# y_pred2 = exp(-0.01) * exp(-0.3506*time_m)
# lines(y_pred ~ time_m)
# lines(y_pred1 ~ time_m,lty=2)
# lines(y_pred2 ~ time_m,lty=2)
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
fit <- sampling(stanDso, data = dat_standard, iter = 4000, warmup=2000) 
fit

dat_pbo <- list(N = length(time_obs), 
                N2 = length(seq(0,3,0.01)),
                y = y_pbo_net_usage, 
                x = time_obs,
                New_x = seq(0,3,0.01)); 
fit2 <- sampling(stanDso, data = dat_pbo, iter = 4000, warmup=2000) 
fit2
# y_pred = extract(fit, 'y_pred')
#plotting the posterior distribution for the parameters
post_beta<-As.mcmc.list(fit,pars="beta0")
plot(post_beta)

## gradient is fit to the data for alpha
## standard_net_usage ~ exp(-alpha*time_obs)
b0 <- extract(fit, 'beta0')
b0<- unlist(b0, use.names=FALSE)
b1 <- extract(fit, 'beta1')
b1<- unlist(b1, use.names=FALSE)

b0p <- extract(fit2, 'beta0')
b0p<- unlist(b0p, use.names=FALSE)
b1p <- extract(fit2, 'beta1')
b1p<- unlist(b1p, use.names=FALSE)

plot(y_standard_net_usage ~ time_obs,ylab="Proportion of children using net (%)",
     xlab="Time in years",yaxt="n",cex.lab=1.6,cex.axis=1.4,xlim=c(0,3),pch=19)
axis(2,las=2, at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.6,cex.axis=1.4)

y_predicted = mean(b0) + mean(b1)*time_m; 
lines(y_predicted ~ time_m,col="black",lty=2,lwd=2)
## Translate back to exponential scale
y_predicted_stn_exp = exp(mean(b0)) * exp(mean(b1)*time_m)

y_predicted_stn_exp_min = exp(quantile(b0,0.25)) * exp(quantile(b1,0.25)*time_m)
y_predicted_stn_exp_max = exp(quantile(b0,0.75)) * exp(quantile(b1,0.75)*time_m)

y2_predicted = mean(b0p) + mean(b1p)*time_m; 
## Translate back to exponential scale
y2_predicted_stn_exp = exp(mean(b0p)) * exp(mean(b1p)*time_m)

y2_predicted_stn_exp_min = exp(quantile(b0p,0.25)) * exp(quantile(b1p,0.25)*time_m)
y2_predicted_stn_exp_max = exp(quantile(b0p,0.75)) * exp(quantile(b1p,0.75)*time_m)


## Final output
par(mfrow = c(1,1))
plot(standard_net_usage[1:3] ~ time_obs[1:3],ylab="Households with at least one net per two occupants (%)",
     xlab="Time in years",yaxt="n",ylim=c(0,1),cex.lab=1.4,cex.axis=1.4,xlim=c(0,3))
axis(2,las=2, at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.4,cex.axis=1.4)

time_obs_offset = time_obs + 0.02

polygon(c(time_m,rev(time_m)),c(y_predicted_stn_exp_min,rev(y_predicted_stn_exp_max)),border=NA,col=transp("grey","0.5"))
polygon(c(time_m,rev(time_m)),c(y2_predicted_stn_exp_min,rev(y2_predicted_stn_exp_max)),border=NA,col=transp("blue","0.5"))
lines(y_predicted_stn_exp ~ time_m,col="black",lty=2,lwd=1)
lines(y2_predicted_stn_exp ~ time_m,col="blue",lty=2,lwd=1)

for(i in 4:6){
  segments(x0=time_obs[i],x1=time_obs[i],
           y0=standard_net_usage[i],y1=standard_net_usage[i+3],lty=1)
  segments(x0=time_obs_offset[i],x1=time_obs_offset[i],
           y0=pbo_net_usage[i],y1=pbo_net_usage[i+3],lty=1,col="blue")
}
points(pbo_net_usage[1:3] ~ time_obs_offset[1:3],col="blue",pch=19)
points(standard_net_usage[1:3] ~ time_obs[1:3])

parms_usage = data.frame(itn_leave_dur_standardLLIN = sample(b1[b1 >= quantile(b1,0.25) & b1 <= quantile(b1,0.75)],1000,replace=FALSE),
                         itn_leave_dur_pbo_net = sample(b1p[b1p >= quantile(b1p,0.25) & b1p <= quantile(b1p,0.75)],1000,replace=FALSE)) ##

parms_usage$itn_leave_dur_standardLLIN_bt = -1/parms_usage$itn_leave_dur_standardLLIN
parms_usage$itn_leave_dur_pbo_bt = -1/parms_usage$itn_leave_dur_pbo_net

PARMS_USAGE = data.frame(parms_usage$itn_leave_dur_standardLLIN_bt,
                         parms_usage$itn_leave_dur_pbo_bt)



## Adding uncertainty estimates to the model predictions
# PREV IN 2 - 10 YEARS
# Arm 1 23.7% (12.8 - 42)
# Arm 2 25.1% (7.2 - 41.7)
# Arm 3 15.9% (1.5 - 37.5)
# Arm 4 8.2% (2.8 - 34.1)


## Staedke 2019/2020
## create a parameter list of any params we can estimate uncertainty for:
params = c("base_prev_arm1","base_prev_arm2","base_prev_arm3","base_prev_arm4","resistance","ITN_0")
param_best_est = c(0.237,0.251,0.159,0.082,0.5210686,0.854)
# param_best_est = c(0.237,0.251,0.159,0.082,0.8543468,0.854)

params_X1 = array(dim=c(1000,length(params)+16))
for(i in 1:length(params)){
  params_X1[,i] = uncertainty_fn(param_best_est[i])  
}

## Add in uncertainty for other params estimated from statistical work elsewhere

## itn_leave_dur
# dim(params_X1)
## PARMS_USAGE from uncertainty params
params_X1[,7] = sample(x = c(PARMS_USAGE[,1]),size = 1000,replace = FALSE)
params_X1[,8] = sample(x = c(PARMS_USAGE[,1]),size = 1000,replace = FALSE)


## res_istance params LLIN 
## res_istance params LLIN 
uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X1[,5],
                                                               param_a = 3.571,param_b = 0.703,
                                                               param_c = 2.444,param_d = 0.422,param_e = 0.356,
                                                               param_f = 4.657,param_g = 0.040) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X1[,5],
                                                               param_a = 0.765,param_b = 0.468,
                                                               param_c = 2.444,param_d = 0.422,param_e = 0.356,
                                                               param_f = 4.657,param_g = 0.04) ## and a range of resistance estimates
## These are arm 1 PermaNet 2.0	and  Arm 3 Olyset	
params_X1[,9] = uncertainty_resistance_params_nets[,2] ## itn_repel_fun
params_X1[,10] = uncertainty_resistance_params_nets[,2] ## itn_repel_fun
params_X1[,11] = uncertainty_resistance_params_nets[,2] ## itn_repel_fun
params_X1[,12] = uncertainty_resistance_params_nets[,1] ## itn_kill_fun
params_X1[,13] = uncertainty_resistance_params_nets[,1] ## itn_kill_fun
params_X1[,14] = uncertainty_resistance_params_nets[,1] ## itn_kill_fun
params_X1[,15] = uncertainty_resistance_params_nets[,3] ## itn_kill_fun

## ** re run ERG each time
uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 1, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X1[,5],
                                                               param_a = 3.571,param_b = 0.703,
                                                               param_c = 2.444,param_d = 0.422,param_e = 0.356,
                                                               param_f = 4.657,param_g = 0.040) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 1, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X1[,5],
                                                               param_a = 0.765,param_b = 0.468,
                                                               param_c = 2.444,param_d = 0.422,param_e = 0.356,
                                                               param_f = 4.657,param_g = 0.04)## and a range of resistance estimates
## These are arm 2 PermaNet 3.0	
params_X1[,16] = uncertainty_resistance_params_nets[,2] ## itn_repel_fun
params_X1[,17] = uncertainty_resistance_params_nets[,2] ## itn_repel_fun
params_X1[,18] = uncertainty_resistance_params_nets[,2] ## itn_repel_fun
params_X1[,19] = uncertainty_resistance_params_nets[,1] ## itn_kill_fun
params_X1[,20] = uncertainty_resistance_params_nets[,1] ## itn_kill_fun
params_X1[,21] = uncertainty_resistance_params_nets[,1] ## itn_kill_fun
params_X1[,22] = uncertainty_resistance_params_nets[,3] ## itn_kill_fun

head(params_X1)


colnames(params_X1) = c("base_prev_arm1","base_prev_arm2","base_prev_arm3","base_prev_arm4","resistance","ITN_0",
                        "itn_leave_dur_pyr","itn_leave_dur_pbo",
                        "itn_repel_fun_1","itn_repel_gam_1","itn_repel_ara_1",
                        "itn_kill_fun_1","itn_kill_gam_1","itn_kill_ara_1",
                        "itn_half_life_1",
                        "itn_repel_fun_pbo","itn_repel_gam_pbo","itn_repel_ara_pbo",
                        "itn_kill_fun_pbo","itn_kill_gam_pbo","itn_kill_ara_pbo",
                        "itn_half_life_pbo")

# write.csv(params_X1,"Input_files/data/input_files_14.csv")
# write.csv(params_X1,"Input_files/data/input_files_14_net_parameters1_all_hutdata.csv")
# write.csv(params_X1,"Input_files/data/input_files_14_net_parameters2_all_hutdata.csv")
# write.csv(params_X1,"Input_files/data/input_files_14_net_parameters3_all_hutdata.csv")
# write.csv(params_X1,"Input_files/data/input_files_14_net_parameters4_all_hutdata.csv")
