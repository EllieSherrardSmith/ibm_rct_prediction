###################################
##
## input file for study ID 4 Corbel

## Specific data needed to estimate parameters 

## Corbel 2012 Lancet

## 1. Pyrethroid resistance
## Estimated from the susceptibility bioassay test
## and the ratio of mosquito species at the trial site
## (each arm independently where data are available)

## Species ratio
## counts 
2508/(2508+1822)
an_gam = 0.579
an_fun = 1 - 0.579
an_ara = 0

## Species bioassay result (from % mortality)

an_gam_bio = 0.4
an_fun_bio = 0
an_ara_bio = 0 ## without specific data on arabiensis we use gamb sl

## we do not have sufficient data to parameterise the model distinctly for species
## instead, we use a global weighted estimate for pyrethroid resistance

pyrethroid_resistance = (an_gam * (1 - an_gam_bio)) + 
  (an_fun * (1- an_fun_bio)) + 
  (an_ara * (1 - an_ara_bio))

pyrethroid_resistance

bendio = read.csv("data/intervention_data/Bendiocarb_uncertainty.csv",header=TRUE)

## Adding uncertainty estimates to the model predictions

## Corbel [X1] Standard LLIN only
## create a parameter list of any params we can estimate uncertainty for:
params = c("base_prev","usage_t0","resistance")
param_best_est = c(0.253,0.25,0.7684)
# param_min_est = c(0.61,0.358,0.80) ## We can discuss setting these!
# param_max_est = c(0.74,0.423,0.98) ## 

resistance_range = uncertainty_fn(0.7684)


params_X1 = params_X2 = params_X3 = expand.grid(runs = 1:1000)

## res_istance params LLIN 
# logistic all data
uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = resistance_range,
                                                               param_a = 3.571,param_b = 0.703,
                                                               param_c = 2.444,param_d = 0.422,param_e = 0.356,
                                                               param_f = 4.657,param_g = 0.040) ## and a range of resistance estimates
# loglogistic all data
uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = resistance_range,
                                                               param_a = 0.765,param_b = 0.468,
                                                               param_c = 2.444,param_d = 0.422,param_e = 0.356,
                                                               param_f = 4.657,param_g = 0.04)## and a range of resistance estimates

# loglogistic west
uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = resistance_range,
                                                               param_a = 0.668,param_b = 0.870,
                                                               param_c = 3.045,param_d = 0.573,param_e = 0.39,
                                                               param_f = 5.129,param_g = 0.036) ## and a range of resistance estimates
# loglogistic east
uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = resistance_range,
                                                               param_a = 0.639,param_b = 0.355,
                                                               param_c = 1.76,param_d = 0.07,param_e = 0.45,
                                                               param_f = 3.43,param_g = 0.01)
## and a range of resistance estimates
## west logistic
uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = resistance_range,
                                                               param_a = 3.784,param_b = 0.727,
                                                               param_c = 3.045,param_d = 0.573,param_e = 0.39,
                                                               param_f = 5.129,param_g = 0.036) ## and a range of resistance estimates

## east logistic
uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range =resistance_range,
                                                               param_a = 3.784,param_b = 0.727,
                                                               param_c = 3.045,param_d = 0.573,param_e = 0.39,
                                                               param_f = 5.129,param_g = 0.036)## and a range of resistance estimates

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

params_X3[,1] = uncertainty_resistance_params_nets[,2]
params_X3[,2] = uncertainty_resistance_params_nets[,2]
params_X3[,3] = uncertainty_resistance_params_nets[,2]
params_X3[,4] = uncertainty_resistance_params_nets[,1]
params_X3[,5] = uncertainty_resistance_params_nets[,1]
params_X3[,6] = uncertainty_resistance_params_nets[,1]
params_X3[,7] = uncertainty_resistance_params_nets[,3]


#bendiocarb
params_X1[,8] = bendio$irs_decay_mort1
params_X1[,9] = bendio$irs_decay_mort2
params_X1[,10] = bendio$irs_decay_succ1
params_X1[,11] = bendio$irs_decay_succ2
params_X1[,12] = bendio$irs_decay_det1
params_X1[,13] = bendio$irs_decay_det2 ## we use the same decay for mortality to fit deterrence

params_X2[,8] = bendio$irs_decay_mort1
params_X2[,9] = bendio$irs_decay_mort2
params_X2[,10] = bendio$irs_decay_succ1
params_X2[,11] = bendio$irs_decay_succ2
params_X2[,12] = bendio$irs_decay_det1
params_X2[,13] = bendio$irs_decay_det2 ## we use the same decay for mortality to fit deterrence

params_X3[,8] = bendio$irs_decay_mort1
params_X3[,9] = bendio$irs_decay_mort2
params_X3[,10] = bendio$irs_decay_succ1
params_X3[,11] = bendio$irs_decay_succ2
params_X3[,12] = bendio$irs_decay_det1
params_X3[,13] = bendio$irs_decay_det2 ## we use the same decay for mortality to fit deterrence



## baseline prevalence
params_X1[,14] = uncertainty_fn(0.253)
params_X2[,14] = uncertainty_fn(0.301)
params_X3[,14] = uncertainty_fn(0.373)

## ITN use
## net use in U6s is noted, 25% of popn represented by pregnant and U6.
## so population net use for each trial arm is: 
## 1 0.38 * 0.25
## 2 0.65 * 0.25
## 3 0.45 * 0.25
## 4 0.7 * 0.25

params_X1[,15] = uncertainty_fn(0.43)
params_X2[,15] = uncertainty_fn(0.58)
params_X3[,15] = uncertainty_fn(0.43)

## IRS coverage
params_X1[,16] = 0
params_X2[,16] = 0
params_X3[,16] = uncertainty_fn(0.8)


colnames(params_X1) = c("itn_kill_fun","itn_kill_gam","itn_kill_ara",
                        "itn_repel_fun","itn_repel_gam","itn_repel_ara",
                        "itn_half_life",
                        "irs_decay_mort1","irs_decay_mort2",
                        "irs_decay_succ1","irs_decay_succ2",
                        "irs_decay_det1","irs_decay_det2",
                        "base_prev",
                        "itn_1","IRS_1")
colnames(params_X2) = c("itn_kill_fun","itn_kill_gam","itn_kill_ara",
                        "itn_repel_fun","itn_repel_gam","itn_repel_ara",
                        "itn_half_life",
                        "irs_decay_mort1","irs_decay_mort2",
                        "irs_decay_succ1","irs_decay_succ2",
                        "irs_decay_det1","irs_decay_det2",
                        "base_prev",
                        "itn_1","IRS_1")
colnames(params_X3) = c("itn_kill_fun","itn_kill_gam","itn_kill_ara",
                        "itn_repel_fun","itn_repel_gam","itn_repel_ara",
                        "itn_half_life",
                        "irs_decay_mort1","irs_decay_mort2",
                        "irs_decay_succ1","irs_decay_succ2",
                        "irs_decay_det1","irs_decay_det2",
                        "base_prev",
                        "itn_1","IRS_1")


params_all = rbind(params_X1,params_X2,params_X3)

params_all$arm = rep(c(1,2,3),each=1000)

# write.csv(params_all,"Input_files/data/input_files_4.csv")
# write.csv(params_all,"Input_files/data/input_files_4_net_parameters1_all_hutdata.csv")
# write.csv(params_all,"Input_files/data/input_files_4_net_parameters2_all_hutdata.csv")
# write.csv(params_all,"Input_files/data/input_files_4_net_parameters3_all_hutdata.csv")
# write.csv(params_all,"Input_files/data/input_files_4_net_parameters4_all_hutdata.csv")


