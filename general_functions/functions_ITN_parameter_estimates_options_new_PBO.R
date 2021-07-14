## Function to translate some level of survival at bioassay into the transmission model parameters

## Nash et al updated the analysis of standard nets (Churcher et al 2016)
## in this analysis there are different functions that could
## potentially explain the association between 
## mortality in a bioassay and mortality in a hut trial

## Nash et al updated the association between net types
## so that we can rewrite those. We keep the assumption that the 
## association between blood feeding, mortality and deterrence is 
## consistent between net types (standard, PBO)

## Here we present alternative functions to estimate the efficacy parameters of the 
## respective net types.


## These are indexed as:


## 1 = logistic fit
## 2 = log-logistic fit


##***********************************************************
## Function 2 Index 1: logistic fit

resistance_ITN_params_1_f = function(is.pbo, species, resistance_range,
                                     param_a,param_b,
                                     param_c,param_d,param_e,
                                     param_f,param_g){
  
  metric = 1 
  
  ## The parameters included here are from the statistical analysis 
  ## determined by Rebecca Nash C:\Users\esherrar\Documents\LLINparameters\Parameter definitions and fitted values.doc
  ## Emails 03/08/2020 with Ben Lambert and Tom Churcher
  
  #Assay to hut mortality conversion		
  param_a = param_a ## All 3.57 ## West 3.784  ## East 3.28
  param_b = param_b ## All 0.70 ## West 0.727  ## East 0.648
  
  #Deterency from mortality		
  param_c = param_c ## All 2.444 ## West 3.045 ## East 1.546
  param_d = param_d ## All 0.422 ## West 0.573 ## East 0.043
  param_e = param_e ## All 0.356 ## West 0.390 ## East 0.454
  
  #Success from mortality		
  param_f = param_f ## All 4.657 ## West 5.129 ## East 3.430
  param_g = param_g ## All 0.040 ## West 0.036 ## East 0.011
  
  
  #Decay in insecticide non-PBO net		
  mup=	-2.429 ## ... gam.medians[9]
  rhop=	-3.007 ## ... gam.medians[10]
  
  
  ## The parameters here are from Griffin et al 2010
  ## originally. These will be updated following the durability studies 
  ## currently underway in Africa through New Nets Project
  
  ## The maximum successful feeding probability per feeding attempt 
  ## (feeding and not dying) in the absence of interventions 
  kp0=0.699 ## derived from Lines et al 1987 and Curtis et al 1990 
  
  ## The half-life of the net relative to it's capacity to kill mosquitoes
  ## with the insecticide active ingredient (a pyrethroid) when there is 
  ## no resistance in mosquitoes. 
  net_halflife=2.64
  
  
  ## The proportion of mosquitoes surviving in the susceptibility bioassay test
  surv_bioassay=resistance_range# a measure of resistance 
  # 0 = no resistance 
  # 1 = 100% survival in discriminating dose bioassay
  
  
  ## defining the associations for probable outcomes of feeding attempts
  #relationship mortality in bioassay -> hut trial, logit scale
  mort_assay = 1-surv_bioassay
  mort_huta   = 1/(1+exp(param_a*param_b-param_a*mort_assay))

  
  ## Next, we define the additional benefit from the next generation ITNs
  ## These are worked out in the ["prioritisations"] manuscript
  
  ## pyrethroid-PBO ITNs
  PBO_benefit = 1 / (1 + exp(-6.39 * mort_huta -  -1.82))

  
  ## Now we work through the probability steps to determine the key input parameters for the model
  ## These probability relationships are determined in Griffin et al 2010 and Churcher et al. 2016
  
  #specify whichever net is used in the RCT
  # mort_assay = if(is.pbo==0) 1-surv_bioassay else if(is.pbo==1) PBO_benefit else if(is.pbo==2) G2_benefit 
  
  mort_hut = if(is.pbo==0) mort_huta else if(is.pbo==1) PBO_benefit else if(is.pbo==2) G2_benefit 
  
  #relationship hut trial mortality -> deterrence
  det_hut = param_e * exp(param_d * (1 - exp(param_c * (1-mort_hut)))/param_c)
  
  #relationship hut trial mortality -> success
  suc_hut = 1-exp(param_g * (1 - exp(param_f * (1-mort_hut)))/param_f)
  
  rep_hut   = 1-suc_hut-mort_hut
  
  ## Combine to estimate the 3 key probable outcomes of feeding attempts
  ## Here we adjust for those mosquitoes not entering treated huts (determined by deterrence)
  n1n0 = 1-det_hut
  kp1  = n1n0*suc_hut
  jp1  = n1n0*rep_hut+(1-n1n0)
  lp1  = n1n0*mort_hut
  
  kp1 = ifelse(kp1 > kp0,kp0,kp1)
  ## (time = 0 time steps after net implementation)
  r_ITN0  = (1-kp1/kp0)*(jp1/(lp1+jp1))	#probability of dying with an encounter with ITN
  d_ITN0  = (1-kp1/kp0)*(lp1/(lp1+jp1))	#probability of repeating behaviour
  s_ITN0  = 1-d_ITN0-r_ITN0             #probability of successfully feeding (surviving and feeding)
  
  
  ## Repeat these to determin the maximum and minimum effects which combine to help determine ITN half life
  mort_maxA   = 1/(1+exp(param_a*param_b-param_a*(1)))
  
  mort_minA   = 1/(1+exp(param_a*param_b-param_a*(0)))
  
  ## pyrethroid-PBO ITNs
  PBO_benefitA = 1 / (1 + exp(-6.387 * mort_maxA -  -1.82))

  ## pyrethroid-PBO ITNs
  PBO_benefitB = 1 / (1 + exp(-6.387 * mort_minA -  -1.82))

  mort_max = if(is.pbo==0) mort_maxA else if(is.pbo==1) PBO_benefitA else if(is.pbo==2) G2_benefitA
  mort_mIN = if(is.pbo==0) mort_minA else if(is.pbo==1) PBO_benefitB else if(is.pbo==2) G2_benefitB
  #{halflife}
  my_max_washes_a = mup+rhop*(mort_max-0.5)		
  my_max_washes   = log(2)/(exp(my_max_washes_a)/(1+exp(my_max_washes_a)))
  
  wash_decay_rate_a = mup +rhop*(mort_hut-0.5)
  wash_decay_rate   = log(2)/(exp(wash_decay_rate_a)/(1+exp(wash_decay_rate_a)))
  itn_half_life     = wash_decay_rate/my_max_washes*net_halflife
  
  ##No need to re-adjust these anymore
  ##Final Parameter estimates for the transmission model
  ERG_d_ITN0 <- d_ITN0
  ERG_s_ITN0 <- s_ITN0
  ERG_r_ITN0 <- 1-ERG_d_ITN0-ERG_s_ITN0
  
  ## Print out these estimates to a data.frame as the function output
  uncertainty_resistance_params_nets = data.frame(ERG_d_ITN0,ERG_r_ITN0,itn_half_life)
  
  return(uncertainty_resistance_params_nets)
}

#Assay to hut mortality conversion		
param_a = param_a ## All 3.57 ## West 3.784  ## East 3.28
param_b = param_b ## All 0.70 ## West 0.727  ## East 0.648

#Deterency from mortality		
param_c = param_c ## All 2.444 ## West 3.045 ## East 1.546
param_d = param_d ## All 0.422 ## West 0.573 ## East 0.043
param_e = param_e ## All 0.356 ## West 0.390 ## East 0.454

#Success from mortality		
param_f = param_f ## All 4.657 ## West 5.129 ## East 3.430
param_g = param_g ## All 0.040 ## West 0.036 ## East 0.011


##All
pyr_itn = resistance_ITN_params_1_f(is.pbo=0, species=1, resistance_range=seq(0,1,0.1),
                          param_a = 3.571,param_b = 0.703,
                          param_c = 2.444,param_d = 0.422,param_e = 0.356,
                          param_f = 4.657,param_g = 0.040)
pyr_pbo_itn = resistance_ITN_params_1_f(is.pbo=1, species=1, resistance_range=seq(0,1,0.1),
                                        param_a = 3.571,param_b = 0.703,
                                        param_c = 2.444,param_d = 0.422,param_e = 0.356,
                                        param_f = 4.657,param_g = 0.040)


##West
pyr_itn = resistance_ITN_params_1_f(is.pbo=0, species=1, resistance_range=seq(0,1,0.01),
                                    param_a = 3.784,param_b = 0.727,
                                    param_c = 3.045,param_d = 0.573,param_e = 0.39,
                                    param_f = 5.129,param_g = 0.036)

pyr_pbo_itn  = resistance_ITN_params_1_f(is.pbo=1, species=1, resistance_range=seq(0,1,0.01),
                                         param_a = 3.784,param_b = 0.727,
                                         param_c = 3.045,param_d = 0.573,param_e = 0.39,
                                         param_f = 5.129,param_g = 0.036)

##East

pyr_itn = resistance_ITN_params_1_f(is.pbo=0, species=1, resistance_range=seq(0,1,0.01),
                                    param_a = 3.28,param_b = 0.648,
                                    param_c = 1.546,param_d = 0.043,param_e = 0.454,
                          param_f = 3.43,param_g = 0.011)

pyr_pbo_itn  = resistance_ITN_params_1_f(is.pbo=1, species=1, resistance_range=seq(0,1,0.01),
                                         param_a = 3.28,param_b = 0.648,
                                         param_c = 1.546,param_d = 0.043,param_e = 0.454,
                                         param_f = 3.43,param_g = 0.011)

##***********************************************************
## Function 2 Index 2: All data, log-logistic fit

resistance_ITN_params_2_f = function(is.pbo, species, resistance_range,
                                     param_a,param_b,
                                     param_c,param_d,param_e,
                                     param_f,param_g){
  
  metric = 1 
  
  ## The parameters included here are from the statistical analysis 
  ## determined by Rebecca Nash C:\Users\esherrar\Documents\LLINparameters\Parameter definitions and fitted values.doc
  ## Emails 03/08/2020 with Ben Lambert and Tom Churcher
  
  #Assay to hut mortality conversion		
  param_a = param_a ## All 0.765 ## West 0.668 ## East 0.639
  param_b = param_b ## All 0.468 ## West 0.870 ## East 0.355
  
  #Deterency from mortality		
  param_c = param_c ## All 2.444 ## West 3.045 ## East 1.546
  param_d = param_d ## All 0.422 ## West 0.573 ## East 0.043
  param_e = param_e ## All 0.356 ## West 0.390 ## East 0.454
  
  #Success from mortality		
  param_f = param_f ## All 4.657 ## West 5.129 ## East 3.430
  param_g = param_g ## All 0.040 ## West 0.036 ## East 0.011
  
  
  
  #Decay in insecticide non-PBO net		
  mup=	-2.429 ## ... gam.medians[9]
  rhop=	-3.007 ## ... gam.medians[10]
  
  
  ## The parameters here are from Griffin et al 2010
  ## originally. These will be updated following the durability studies 
  ## currently underway in Africa through New Nets Project
  
  ## The maximum successful feeding probability per feeding attempt 
  ## (feeding and not dying) in the absence of interventions 
  kp0=0.699 ## derived from Lines et al 1987 and Curtis et al 1990 
  
  ## The half-life of the net relative to it's capacity to kill mosquitoes
  ## with the insecticide active ingredient (a pyrethroid) when there is 
  ## no resistance in mosquitoes. 
  net_halflife=2.64
  
  
  ## The proportion of mosquitoes surviving in the susceptibility bioassay test
  surv_bioassay=resistance_range# a measure of resistance 
  # 0 = no resistance 
  # 1 = 100% survival in discriminating dose bioassay
  
  
  ## defining the associations for probable outcomes of feeding attempts
  #relationship mortality in bioassay -> hut trial, logit scale
  mort_assay = 1-surv_bioassay
  mort_huta   = 1/(1+(mort_assay/param_a)^-param_b)
  
  ## Next, we define the additional benefit from the next generation ITNs
  ## These are worked out in the ["prioritisations"] manuscript
  
  
  ## pyrethroid-PBO ITNs
  PBO_benefit = 1 / (1 + exp(-6.387 * mort_huta -  -1.82))
  PBO_benefit_log_upp = 1 / (1 + exp(-6.16 * mort_huta - -1.88))
  PBO_benefit_log_low = 1 / (1 + exp(-6.61 * mort_huta - -1.76))
  
  # ## pyrethroid-chlorfenapyr ITNs
  # G2_benefit = 1 / (1 + exp(-4.6637 * mort_huta -  -0.5537 ))
  # G2_benefit_log_upp = 1 / (1 + exp(-5.4243 * mort_huta - 0.3408))
  # G2_benefit_log_low = 1 / (1 + exp(-3.9817 * mort_huta - -0.7598))
  
  
  ## Now we work through the probability steps to determine the key input parameters for the model
  ## These probability relationships are determined in Griffin et al 2010 and Churcher et al. 2016
  
  #specify whichever net is used in the RCT
  # mort_assay = if(is.pbo==0) 1-surv_bioassay else if(is.pbo==1) PBO_benefit else if(is.pbo==2) G2_benefit 
  
  mort_hut = if(is.pbo==0) mort_huta else if(is.pbo==1) PBO_benefit #else if(is.pbo==2) G2_benefit 
  
  
  #relationship hut trial mortality -> deterrence
  det_hut = param_e * exp(param_d * (1 - exp(param_c * (1-mort_hut)))/param_c)
  
  #relationship hut trial mortality -> success
  suc_hut = 1-exp(param_g * (1 - exp(param_f * (1-mort_hut)))/param_f)
  
  rep_hut   = 1-suc_hut-mort_hut
  
  ## Combine to estimate the 3 key probable outcomes of feeding attempts
  ## Here we adjust for those mosquitoes not entering treated huts (determined by deterrence)
  n1n0 = 1-det_hut
  kp1  = n1n0*suc_hut
  jp1  = n1n0*rep_hut+(1-n1n0)
  lp1  = n1n0*mort_hut
  
  kp1 = ifelse(kp1 > kp0,kp0,kp1)
  ## (time = 0 time steps after net implementation)
  r_ITN0  = (1-kp1/kp0)*(jp1/(lp1+jp1))	#probability of dying with an encounter with ITN
  d_ITN0  = (1-kp1/kp0)*(lp1/(lp1+jp1))	#probability of repeating behaviour
  s_ITN0  = 1-d_ITN0-r_ITN0             #probability of successfully feeding (surviving and feeding)
  
  
  ## Repeat these to determin the maximum and minimum effects which combine to help determine ITN half life
  mort_maxA   = 1/(1+((1)/param_a)^-param_b)
   
  mort_minA   = 1/(1+((0)/param_a)^-param_b)

    ## pyrethroid-PBO ITNs
  PBO_benefitA = 1 / (1 + exp(-6.387 * mort_maxA -  -1.82))
  
  ## pyrethroid-PBO ITNs
  PBO_benefitB = 1 / (1 + exp(-6.387 * mort_minA -  -1.82))
  
  mort_max = if(is.pbo==0) mort_maxA else if(is.pbo==1) PBO_benefitA else if(is.pbo==2) G2_benefitA
  mort_mIN = if(is.pbo==0) mort_minA else if(is.pbo==1) PBO_benefitB else if(is.pbo==2) G2_benefitB
  
  #{halflife}
  my_max_washes_a = mup +rhop*(mort_max-0.5)		
  my_max_washes   = log(2)/(exp(my_max_washes_a)/(1+exp(my_max_washes_a)))
  
  wash_decay_rate_a = mup +rhop*(mort_hut-0.5)
  wash_decay_rate   = log(2)/(exp(wash_decay_rate_a)/(1+exp(wash_decay_rate_a)))
  itn_half_life     = wash_decay_rate/my_max_washes*net_halflife
  
  ##No need to re-adjust these anymore
  ##Final Parameter estimates for the transmission model
  ERG_d_ITN0 <- d_ITN0
  ERG_s_ITN0 <- s_ITN0
  ERG_r_ITN0 <- 1-ERG_d_ITN0-ERG_s_ITN0
  
  ## Print out these estimates to a data.frame as the function output
  uncertainty_resistance_params_nets = data.frame(ERG_d_ITN0,ERG_r_ITN0,itn_half_life)
  
  return(uncertainty_resistance_params_nets)
}

#Assay to hut mortality conversion		
param_a = param_a ## All 0.765 ## West 0.668 ## East 0.639
param_b = param_b ## All 0.468 ## West 0.870 ## East 0.355

#Deterency from mortality		
param_c = param_c ## All 2.444 ## West 3.045 ## East 1.546
param_d = param_d ## All 0.422 ## West 0.573 ## East 0.043
param_e = param_e ## All 0.356 ## West 0.390 ## East 0.454

#Success from mortality		
param_f = param_f ## All 4.657 ## West 5.129 ## East 3.430
param_g = param_g ## All 0.040 ## West 0.036 ## East 0.011
##all
pyr_itn = resistance_ITN_params_2_f(is.pbo=0, species=1, resistance_range=seq(0,1,0.1),
                          param_a = 0.765,param_b = 0.468,
                          param_c = 2.444,param_d = 0.422,param_e = 0.356,
                          param_f = 4.657,param_g = 0.04)
pyr_pbo_itn = resistance_ITN_params_2_f(is.pbo=1, species=1, resistance_range=seq(0,1,0.1),
                                        param_a = 0.765,param_b = 0.468,
                                        param_c = 2.444,param_d = 0.422,param_e = 0.356,
                                        param_f = 4.657,param_g = 0.04)


################
## West parameters

pyr_itn = resistance_ITN_params_2_f(is.pbo=0, species=1, resistance_range=seq(0,1,0.01),
                                    param_a = 0.668,param_b = 0.870,
                                    param_c = 3.045,param_d = 0.573,param_e = 0.39,
                                    param_f = 5.129,param_g = 0.036)

pyr_pbo_itn = resistance_ITN_params_2_f(is.pbo=1, species=1, resistance_range=seq(0,1,0.01),
                                        param_a = 0.668,param_b = 0.870,
                                        param_c = 3.045,param_d = 0.573,param_e = 0.39,
                                        param_f = 5.129,param_g = 0.036)

#Assay to hut mortality conversion		
param_a = param_a ## All 0.765 ## West 0.668  ## East 0.639
param_b = param_b ## All 0.468 ## West 0.870  ## East 0.355

#Deterency from mortality		
param_c = param_c ## All 2.444 ## West 3.045 ## East 1.546
param_d = param_d ## All 0.422 ## West 0.573 ## East 0.043
param_e = param_e ## All 0.356 ## West 0.390 ## East 0.454

#Success from mortality		
param_f = param_f ## All 4.657 ## West 5.129 ## East 3.430
param_g = param_g ## All 0.040 ## West 0.036 ## East 0.011


################
## East parameters

pyr_itn =resistance_ITN_params_2_f(is.pbo=0, species=1, resistance_range=seq(0,1,0.01),
                          param_a = 0.639,param_b = 0.355,
                          param_c = 1.76,param_d = 0.07,param_e = 0.45,
                          param_f = 3.43,param_g = 0.01)


pyr_pbo_itn =resistance_ITN_params_2_f(is.pbo=1, species=1, resistance_range=seq(0,1,0.01),
                                       param_a = 0.639,param_b = 0.355,
                                       param_c = 1.76,param_d = 0.07,param_e = 0.45,
                                       param_f = 3.43,param_g = 0.01)
