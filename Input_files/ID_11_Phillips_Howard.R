###################################
##
## input file for study ID 12 Phillips-Howard



## Phillips-Howard 2003

## 1. Pyrethroid resistance assume none (2003)
## 2. Net use information
ITN_1 = expand.grid(itns_asembo_control1 = uncertainty_fn(0.05)) ## prior to trial
ITN_1$itns_asembo_control2 = uncertainty_fn(0.774) ## during trial in Asembo
ITN_1$itns_asembo_t1 = uncertainty_fn(0.659) ## during trial in Asembo with the earlier ITN distribution
ITN_1$itns_asembo_t2 = uncertainty_fn(0.825) ## later distribution Asembo
ITN_1$itns_gem_t1 = uncertainty_fn(0.659) ## Gem distribution 1
ITN_1$itns_gem_t2 = uncertainty_fn(0.825) ## Gem distribution 2

ITN_1$total_MM  = uncertainty_fn(0.9)*150

write.csv(params_all,"Input_files/data/input_files_12.csv")
