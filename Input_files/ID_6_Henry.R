###################################
##
## input file for study ID 7 Henry

## Specific data needed to estimate parameters 

## 1. Pyrethroid resistance
## Estimated from the susceptibility bioassay test
## and the ratio of mosquito species at the trial site
## (each arm independently where data are available)

## Species ratio
## counts 
## uNSPECIFIed but roughly
an_gam = 0.66
an_fun = 0.34
an_ara = 0

## Species bioassay result (from % mortality)

an_gam_bio = 0.608
an_fun_bio = 1
an_ara_bio = 1 ## without specific data on arabiensis we use gamb sl

## we do not have sufficient data to parameterise the model distinctly for species
## instead, we use a global weighted estimate for pyrethroid resistance

pyrethroid_resistance = (an_gam * (1 - an_gam_bio)) + 
  (an_fun * (1- an_fun_bio)) + 
  (an_ara * (1 - an_ara_bio))

pyrethroid_resistance


## Adding uncertainty estimates to the model predictions
## untreated ITNs - assumed the equivalent impact from 100% resistance
uncertainty_resistance_params_untreated_nets = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = rep(1,1000),
                                                               param_a = 3.571,param_b = 0.703,
                                                               param_c = 2.444,param_d = 0.422,param_e = 0.356,
                                                               param_f = 4.657,param_g = 0.040) ## and a range of resistance estimates

uncertainty_resistance_params_untreated_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = rep(1,1000),
                                                               param_a = 0.765,param_b = 0.468,
                                                               param_c = 2.444,param_d = 0.422,param_e = 0.356,
                                                               param_f = 4.657,param_g = 0.04) ## and a range of resistance estimates
## resistance params LLIN 
uncertainty_resistance_params_nets = resistance_ITN_params_1_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 3.571,param_b = 0.703,
                                                               param_c = 2.444,param_d = 0.422,param_e = 0.356,
                                                               param_f = 4.657,param_g = 0.040) ## and a range of resistance estimates

uncertainty_resistance_params_nets = resistance_ITN_params_2_f(is.pbo = 0, ## test for standard nets
                                                               species = 1, ## generic mosquitos
                                                               resistance_range = uncertainty_fn(pyrethroid_resistance),
                                                               param_a = 0.765,param_b = 0.468,
                                                               param_c = 2.444,param_d = 0.422,param_e = 0.356,
                                                               param_f = 4.657,param_g = 0.04) ## and a range of resistance estimates

params = c("ITN_1")
param_best_est = c(0.8)
# param_min_est = c(0.61,0.358,0.80) ## We can discuss setting these!
# param_max_est = c(0.74,0.423,0.98) ## 

## How do we put in a bed net without insecticide

params_X1 = params_X2 = expand.grid(runs = 1:1000)


params_X1[,1] = uncertainty_resistance_params_untreated_nets[,2]
params_X1[,2] = uncertainty_resistance_params_untreated_nets[,2]
params_X1[,3] = uncertainty_resistance_params_untreated_nets[,2]
params_X1[,4] = uncertainty_resistance_params_untreated_nets[,1]
params_X1[,5] = uncertainty_resistance_params_untreated_nets[,1]
params_X1[,6] = uncertainty_resistance_params_untreated_nets[,1]
params_X1[,7] = uncertainty_resistance_params_untreated_nets[,3]

params_X2[,1] = uncertainty_resistance_params_nets[,2]
params_X2[,2] = uncertainty_resistance_params_nets[,2]
params_X2[,3] = uncertainty_resistance_params_nets[,2]
params_X2[,4] = uncertainty_resistance_params_nets[,1]
params_X2[,5] = uncertainty_resistance_params_nets[,1]
params_X2[,6] = uncertainty_resistance_params_nets[,1]
params_X2[,7] = uncertainty_resistance_params_nets[,3]




## ITN use for treated nets sites

params_X1[,8] = 0

params_X2[,8] = uncertainty_fn(0.8)

colnames(params_X1) = c("itn_repel_fun","itn_repel_gam","itn_repel_ara",
                        "itn_kill_fun","itn_kill_gam","itn_kill_ara",
                        "itn_half_life",
                        "ITN_1")
colnames(params_X2) = c("itn_repel_fun","itn_repel_gam","itn_repel_ara",
                        "itn_kill_fun","itn_kill_gam","itn_kill_ara",
                        "itn_half_life",
                        "ITN_1")

head(params_X1)
head(params_X2)


params_all = rbind(params_X1,params_X2)
params_all$arm = rep(c(1,2),each=1000)


# write.csv(params_all,"Input_files/data/input_files_7.csv")
# write.csv(params_all,"Input_files/data/input_files_7_net_parameters1_all_hutdata.csv")
# write.csv(params_all,"Input_files/data/input_files_7_net_parameters2_all_hutdata.csv")
# write.csv(params_all,"Input_files/data/input_files_7_net_parameters3_all_hutdata.csv")
# write.csv(params_all,"Input_files/data/input_files_7_net_parameters4_all_hutdata.csv")

