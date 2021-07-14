###################################
##
## input file for study ID 3 Chaccour

## Specific data needed to estimate parameters 

## Chaccour - entomological data from Joe Wagman 

## 1. Pyrethroid resistance
## Estimated from the susceptibility bioassay test
## and the ratio of mosquito species at the trial site
## (each arm independently where data are available)

## Species ratio
## control arm
an_gam = 0.028
an_ara = 0.006
an_fun = 0.966

## Species bioassay result (from % mortality)

an_gam_bio = 1
an_fun_bio = 0.865
an_ara_bio = 1 ## without specific data on arabiensis we use gamb sl

## we do not have sufficient data to parameterise the model distinctly for species
## instead, we use a global weighted estimate for pyrethroid resistance

pyrethroid_resistance_1 = (an_gam * (1 - an_gam_bio)) + 
  (an_fun * (1- an_fun_bio)) + 
  (an_ara * (1 - an_ara_bio))

pyrethroid_resistance_1

## Actellic trial arm
an_gam = 0.053
an_ara = 0.03
an_fun = 0.917

## Species bioassay result (from % mortality)

an_gam_bio = 1
an_fun_bio = 0.865
an_ara_bio = 1 ## without specific data on arabiensis we use gamb sl

## we do not have sufficient data to parameterise the model distinctly for species
## instead, we use a global weighted estimate for pyrethroid resistance

pyrethroid_resistance_2 = (an_gam * (1 - an_gam_bio)) + 
  (an_fun * (1- an_fun_bio)) + 
  (an_ara * (1 - an_ara_bio))

pyrethroid_resistance_2


Act = read.csv("data/intervention_data/Actellic_uncertainty.csv",header=TRUE)

## Adding uncertainty estimates to the model predictions

## Chaccour [X1] Standard LLIN
## create a parameter list of any params we can estimate uncertainty for:
params = c("base_prev","ITN_1","ITN_2","resistance")
param_best_est = c(0.62,0.561,0.943,0.13041)


params_X1 = array(dim=c(1000,length(params)+15))
for(i in 1:length(params)){
  params_X1[,i] = uncertainty_fn(param_best_est[i])  
}

## Chaccour [X2] Standard LLIN + Actellic
## create a parameter list of any params we can estimate uncertainty for:
params = c("base_prev","ITN_1","ITN_2","resistance")
param_best_est2 = c(0.647,0.514,0.952,0.123795)


params_X2 = array(dim=c(1000,length(params)+15))
for(i in 1:length(params)){
  params_X2[,i] = uncertainty_fn(param_best_est2[i])  
}


## res_istance params LLIN 

## res_istance params LLIN 
# logistic all data
uncertainty_resistance_params_nets_x1 = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X1[,5],
                                                               param_a = 3.571,param_b = 0.703,
                                                               param_c = 2.444,param_d = 0.422,param_e = 0.356,
                                                               param_f = 4.657,param_g = 0.040) ## and a range of resistance estimates
# loglogistic all data
uncertainty_resistance_params_nets_x1 = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = params_X1[,5],
                                                               param_a = 0.765,param_b = 0.468,
                                                               param_c = 2.444,param_d = 0.422,param_e = 0.356,
                                                               param_f = 4.657,param_g = 0.04)## and a range of resistance estimates

## change as needed for east and west
params_X1[,5] =uncertainty_resistance_params_nets_x1[,2]# uncertainty_fn(ERG_r_ITN0[1]) ## itn_repel_fun
params_X1[,6] =uncertainty_resistance_params_nets_x1[,2]# uncertainty_fn(ERG_r_ITN0[1]) ## itn_repel_gam
params_X1[,7] =uncertainty_resistance_params_nets_x1[,2]# uncertainty_fn(ERG_r_ITN0[1]) ## itn_repel_ara

params_X1[,8] =uncertainty_resistance_params_nets_x1[,1]# uncertainty_fn(ERG_d_ITN0[1]) ## itn_kill_fun
params_X1[,9] =uncertainty_resistance_params_nets_x1[,1]# uncertainty_fn(ERG_d_ITN0[1]) ## itn_kill_gam
params_X1[,10] =uncertainty_resistance_params_nets_x1[,1]# uncertainty_fn(ERG_d_ITN0[1]) ## itn_kill_ara

params_X1[,11] = uncertainty_resistance_params_nets_x1[,3]#rnorm(itn_half_life[1],n = 1000,sd = 0.02) ## itn_half_life

#Actellic base only
params_X1[,12] = Act$irs_decay_mort1
params_X1[,13] = Act$irs_decay_mort2
params_X1[,14] = Act$irs_decay_succ1
params_X1[,15] = Act$irs_decay_succ2
params_X1[,16] = Act$irs_decay_det1
params_X1[,17] = Act$irs_decay_det2 ## we use the same decay for mortality to fit deterrence

params_X1[,18] = rep(0,1000) ## IRS cover (none in X1 Standard LLIN only)

params_X1[,19] = rep(1,1000) ## IRS cover (none in X1 Standard LLIN only)

head(params_X1)


colnames(params_X1) = c("base_prev","ITN_1","ITN_2","resistance",
                        "itn_repel_fun","itn_repel_gam","itn_repel_ara",
                        "itn_kill_fun","itn_kill_gam","itn_kill_ara",
                        "itn_half_life",
                        "irs_decay_mort1","irs_decay_mort2",
                        "irs_decay_succ1","irs_decay_succ2",
                        "irs_decay_det1","irs_decay_det2","irs_cover","arm")

params_X2[,5] = uncertainty_resistance_params_nets_x2[,2]#uncertainty_fn(ERG_r_ITN0[1]) ## itn_repel_fun
params_X2[,6] = uncertainty_resistance_params_nets_x2[,2]#uncertainty_fn(ERG_r_ITN0[1]) ## itn_repel_gam
params_X2[,7] = uncertainty_resistance_params_nets_x2[,2]#uncertainty_fn(ERG_r_ITN0[1]) ## itn_repel_ara

params_X2[,8] = uncertainty_resistance_params_nets_x2[,1]#uncertainty_fn(ERG_d_ITN0[1]) ## itn_kill_fun
params_X2[,9] = uncertainty_resistance_params_nets_x2[,1]#uncertainty_fn(ERG_d_ITN0[1]) ## itn_kill_gam
params_X2[,10] = uncertainty_resistance_params_nets_x2[,1]#uncertainty_fn(ERG_d_ITN0[1]) ## itn_kill_ara

params_X2[,11] = uncertainty_resistance_params_nets_x2[,3]#rnorm(itn_half_life[1],n = 1000,sd = 0.02) ## itn_half_life

#Actellic baSe only
params_X2[,12] = Act$irs_decay_mort1
params_X2[,13] = Act$irs_decay_mort2
params_X2[,14] = Act$irs_decay_succ1
params_X2[,15] = Act$irs_decay_succ2
params_X2[,16] = Act$irs_decay_det1
params_X2[,17] = Act$irs_decay_det2 ## we use the same decay for mortality to fit deterrence

params_X2[,18] = uncertainty_fn(0.83) ## IRS cover (none in X1 Standard LLIN only)

params_X2[,19] = rep(2,1000) ## IRS cover (none in X1 Standard LLIN only)

colnames(params_X2) = c("base_prev","ITN_1","ITN_2","resistance",
                        "itn_repel_fun","itn_repel_gam","itn_repel_ara",
                        "itn_kill_fun","itn_kill_gam","itn_kill_ara",
                        "itn_half_life",
                        "irs_decay_mort1","irs_decay_mort2",
                        "irs_decay_succ1","irs_decay_succ2",
                        "irs_decay_det1","irs_decay_det2","irs_cover","arm")

params = rbind(params_X1, params_X2)


# write.csv(params,"Input_files/data/input_files_[CHANGE AS NEEDED]_net_parameters1_all_hutdata.csv")
