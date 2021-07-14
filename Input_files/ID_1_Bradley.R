###################################
##
## input file for study ID 2 Bradley

## Specific data needed to estimate parameters 

## 1. Pyrethroid resistance
## Estimated from the susceptibility bioassay test
## and the ratio of mosquito species at the trial site
## (each arm independently where data are available)

## Species ratio
an_gam = 1
an_fun = 0
an_ara = 0

## Species bioassay result (from % mortality)

an_gam_bio = 0.796
an_fun_bio = 1
an_ara_bio = 1 ## without specific data on arabiensis we use gamb sl

## we do not have sufficient data to parameterise the model distinctly for species
## instead, we use a global weighted estimate for pyrethroid resistance

pyrethroid_resistance = (an_gam * (1 - an_gam_bio)) + 
  (an_fun * (1- an_fun_bio)) + 
  (an_ara * (1 - an_ara_bio))

pyrethroid_resistance


## No other data

## 3. IRS product efficacies
params_estimated = read.csv("data/intervention_data/pyrethroid_IRS_resistance.csv",header=TRUE)[,2:7]

Bn = read.csv("data/intervention_data/Bendiocarb_uncertainty.csv",header=TRUE)

params_estimated[21,]##79.6% mortality at bioassay

params_bradley = expand.grid(irs_mort_1 = c(params_estimated[12:28,1]))
params_bradley$irs_mort_2 = c(params_estimated[12:28,2])
params_bradley$irs_succ_1 = c(params_estimated[12:28,3])
params_bradley$irs_succ_2 = c(params_estimated[12:28,4])
params_bradley$irs_det_1 = c(params_estimated[12:28,5])
params_bradley$irs_det_2 = c(params_estimated[12:28,6])

params_irs_delta = dplyr::sample_n(params_bradley, 1000,replace=TRUE)
##

## Bradley [X1] Standard LLIN
## create a parameter list of any params we can estimate uncertainty for:
params = c("base_prev_X1","base_prev_X2","usage_t0_x1","usage_t0_x2","resistance")
param_best_est = c(0.175,0.184,0.301,0.317,0.204)

params_X1 = array(dim=c(1000,length(params)+21))
for(i in 1:length(params)){
  params_X1[,i] = uncertainty_fn(param_best_est[i])  
}

## Add in uncertainty for other params estimated from statistical work elsewhere

## itn_leave_dur
## Use default

## res_istance params LLIN 
# logistic all data
uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X1[,5],
                                                               param_a = 3.571,param_b = 0.703,
                                                               param_c = 2.444,param_d = 0.422,param_e = 0.356,
                                                               param_f = 4.657,param_g = 0.040) ## and a range of resistance estimates
# loglogistic all data
uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X1[,5],
                                                               param_a = 0.765,param_b = 0.468,
                                                               param_c = 2.444,param_d = 0.422,param_e = 0.356,
                                                               param_f = 4.657,param_g = 0.04)## and a range of resistance estimates

# loglogistic west
uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X1[,5],
                                                               param_a = 0.668,param_b = 0.870,
                                                               param_c = 3.045,param_d = 0.573,param_e = 0.39,
                                                               param_f = 5.129,param_g = 0.036) ## and a range of resistance estimates
# loglogistic east
uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X1[,5],
                                                               param_a = 0.639,param_b = 0.355,
                                                               param_c = 1.76,param_d = 0.07,param_e = 0.45,
                                                               param_f = 3.43,param_g = 0.01)
## and a range of resistance estimates
## west logistic
uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X1[,5],
                                                               param_a = 3.784,param_b = 0.727,
                                                               param_c = 3.045,param_d = 0.573,param_e = 0.39,
                                                               param_f = 5.129,param_g = 0.036) ## and a range of resistance estimates

## east logistic
uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X1[,5],
                                                               param_a = 3.784,param_b = 0.727,
                                                               param_c = 3.045,param_d = 0.573,param_e = 0.39,
                                                               param_f = 5.129,param_g = 0.036)## and a range of resistance estimates


params_X1[,6] = uncertainty_resistance_params_nets[,2] ## itn_repel_fun
params_X1[,7] = uncertainty_resistance_params_nets[,2] ## itn_repel_gam
params_X1[,8] = uncertainty_resistance_params_nets[,2] ## itn_repel_ara

params_X1[,9] = uncertainty_resistance_params_nets[,1] ## itn_kill_fun
params_X1[,10] = uncertainty_resistance_params_nets[,1] ## itn_kill_gam
params_X1[,11] = uncertainty_resistance_params_nets[,1] ## itn_kill_ara

params_X1[,12] = uncertainty_resistance_params_nets[,3] ## itn_half_life

head(params_X1)

#Bendiocarb### base only
params_X1[,13] = Bn$irs_decay_mort1
params_X1[,14] = Bn$irs_decay_mort2
params_X1[,15] = Bn$irs_decay_succ1
params_X1[,16] = Bn$irs_decay_succ2
params_X1[,17] = Bn$irs_decay_det1
params_X1[,18] = Bn$irs_decay_det2 ## we use the same decay for mortality to fit deterrence

#Pyrethroid### base only
params_X1[,19] = params_irs_delta$irs_mort_1
params_X1[,20] = params_irs_delta$irs_mort_2
params_X1[,21] = params_irs_delta$irs_succ_1
params_X1[,22] = params_irs_delta$irs_succ_2
params_X1[,23] = params_irs_delta$irs_det_1
params_X1[,24] = params_irs_delta$irs_det_2 ## we use the same decay for mortality to fit deterrence

params_X1[,25] = uncertainty_fn(0.737) ## IRS cover 1st spraying (none in X1 Standard LLIN only)
params_X1[,26] = uncertainty_fn(0.765) ## IRS cover 2nd spraying (none in X1 Standard LLIN only)

head(params_X1)


colnames(params_X1) = c("base_prev_X1","base_prev_X2","usage_t0_x1","usage_t0_x2","resistance",
                        "itn_repel_fun","itn_repel_gam","itn_repel_ara",
                        "itn_kill_fun","itn_kill_gam","itn_kill_ara",
                        "itn_half_life",
                        "irs_decay_mort1_bn","irs_decay_mort2_bn",
                        "irs_decay_succ1_bn","irs_decay_succ2_bn",
                        "irs_decay_det1_bn","irs_decay_det2_bn",
                        "irs_decay_mort1_py","irs_decay_mort2_py",
                        "irs_decay_succ1_py","irs_decay_succ2_py",
                        "irs_decay_det1_py","irs_decay_det2_py",
                        "irs_cover1", "irs_cover2")

# write.csv(params_X1,"Input_files/data/input_files_2_net_parameters1_allhutdata.csv")
# write.csv(params_X1,"Input_files/data/input_files_2_net_parameters2_allhutdata.csv")
# write.csv(params_X1,"Input_files/data/input_files_2_net_parameters3_allhutdata.csv")
# write.csv(params_X1,"Input_files/data/input_files_2_net_parameters4_allhutdata.csv")


