#############################################################
##
## Final script to determine parameters for net efficacy
##
##############################################################



####################################################
##
## Fits look ok so now extract the uncertainti using the function 

## Now we need to add in the net parameter estimates 
## working from the function the estimate crude estimates and
## using the parameters determined in Table 2 main manuscipt

## Critical is to keep, for each fit
## the same row of parameters for that simualtion

library(rstan)

## keep rows associated from the bayes posterior draws
data_picker = sample(1:4000,size = 1000,replace=TRUE)
data_picker2 =  sample(1:1000,size = 100,replace=FALSE)

## b is when we use just pyrethroid and pyrethroid-PBO nets
## h is just permaNet 3
## j is just Olyset Plus

## draw from the posterior distribution with the respective inputs
## 
## THIS IS THE SAME FOR ALL COMBINATIONS
## ASSOCIATION 1 MORTALITY TO BIOASSAY
setwd("C:/Users/esherrar/Documents/Rprojects/ibm_rct_prediction")
getwd() ## confirm

# Global set
ll_1 = readRDS("data/Entomological model fits/rds files/Combined/Bioassay/fit_ew_comb_log_logistic.rds")
# ll_1 = readRDS("data/Entomological model fits/rds files/Huts separate/Bioassay/fit_west_log_logistic.rds")
# ll_1 = readRDS("data/Entomological model fits/rds files/Huts separate/Bioassay/fit_east_log_logistic.rds")

# ll_1 = readRDS("data/Entomological model fits/rds files/Combined/Bioassay/fit_ew_comb_logistic.rds")
# ll_1 = readRDS("data/Entomological model fits/rds files/Huts separate/Bioassay/fit_west_logistic.rds")
# ll_1 = readRDS("data/Entomological model fits/rds files/Huts separate/Bioassay/fit_east_logistic.rds")

LL_fit <- extract(ll_1, permuted = TRUE)
# L_fit <- extract(ll_1, permuted = TRUE)

#function shape:
f_LOG_logistic <- function(x, a, c){
  mort = 1/(1+(x/a)^(-c))
  surv = (1-mort)
  return(surv)
}
f_logistic <- function(x, a, c){
  mort = 1/(1+exp(-(x-c)*a))
  surv = (1-mort)
  return(surv)
}

LL_a_full <- LL_fit$a[data_picker]
LL_c_full <- LL_fit$c[data_picker]

quantile(LL_fit$a,c(0.05,0.5,0.950))
quantile(LL_fit$c,c(0.05,0.5,0.950))

# quantile(L_fit$a,c(0.05,0.5,0.950))
# quantile(L_fit$c,c(0.05,0.5,0.950))

## Just confirm nothing spurious here
hut_surv_LL <- f_LOG_logistic(x = seq(0,1,0.01), LL_a_full[1], LL_c_full[1])
plot(hut_surv_LL ~ seq(1,0,-0.01),ylim=c(0,1),
)
for(i in 1:10400){
  hut_surv_LL <- f_LOG_logistic(x = seq(0,1,0.01), LL_a_full[i], LL_c_full[i])
  lines(hut_surv_LL ~ seq(1,0,-0.01))
  
}

## 
## ASSOCIATION 2 BENEFIT OF PBO
## THIS REMAINS THE SAME FOR ALL COMBINATIONS
benefitall <- readRDS("data/Entomological model fits/rds files/ento_pbo_benefit.rds")
pbo_bene <- extract(benefitall, permuted = TRUE)


fitbene_1 <- pbo_bene$alpha1[data_picker]
fitbene_2 <- pbo_bene$alpha2[data_picker]

quantile(pbo_bene$alpha1,c(0.05,0.5,0.950))
quantile(pbo_bene$alpha2,c(0.05,0.5,0.950))

# PBO_benefit = 1 / (1 + exp(-fitbene_2[1] * seq(0,1,0.01) -  fitbene_1[1]))
# plot(PBO_benefit,ylim=c(0,1))
# 
# PBO_benefit = array(dim=c(101,100))
# for(i in 1:100){
#   PBO_benefit[,i] = 1 / (1 + exp(-fitbene_2[i] * seq(0,1,0.01) -  fitbene_1[i]))
#   
# }
# for(i in 1:100) {
#   lines(PBO_benefit[,i])
# }

## 
## ASSOCIATION 3 MORTALITY TO DETERRENCE
## THIS CAN BE EITHER (GENERIC TO 'NETS')
##                    (SPECIFIC TO 'NET CLASS': Pyr-ITN or Pyr-PBO-ITN)
##                    (PRODUCT SPECIFIC)

# load fit
fit1_a <- readRDS("data/Entomological model fits/rds files/Combined/Feeding attempt/deterrence_ew_fit.rds")
# fit1_a <- readRDS("data/Entomological model fits/rds files/Huts separate/Feeding attempt/deterrence_west_fit_corrected.rds")
# fit1_a <- readRDS("data/Entomological model fits/rds files/Huts separate/Feeding attempt/deterrence_east_fit_corrected.rds")

fit1_a_fit <- extract(fit1_a, permuted = TRUE)

## Just pyrethroid and pyrethroid-PBO 
fit1_a_c <- fit1_a_fit$c[c(data_picker)]
fit1_a_d <- fit1_a_fit$d[c(data_picker)]
fit1_a_e <- fit1_a_fit$e[c(data_picker)]

quantile(fit1_a_fit$c,0.5)
quantile(fit1_a_fit$d,0.5)
quantile(fit1_a_fit$e,0.5)

deterrence = fit1_a_e[1] * exp(fit1_a_d[1] * (1 - exp(fit1_a_c[1] * hut_surv_LL)) / fit1_a_c[1])
plot(deterrence ~ hut_surv_LL,ylim=c(0,1),xlim=c(0,1))
for(i in 1:100){
  hut_surv_LL <- 1 - f_LOG_logistic(x = seq(0,1,0.01), LL_a_full[i], LL_c_full[i])
  hut_det_LL <- fit1_a_e[i] * exp(fit1_a_d[i] * (1 - exp(fit1_a_c[i] * hut_surv_LL)) / fit1_a_c[i])
  lines(hut_det_LL ~ hut_surv_LL)
  
}

## 
## ASSOCIATION 4 MORTALITY TO feed
## THIS CAN BE EITHER (GENERIC TO 'NETS')
##                    (SPECIFIC TO 'NET CLASS': Pyr-ITN or Pyr-PBO-ITN)
##                    (PRODUCT SPECIFIC)
#full_model<- stan_model("ento analysis/Nash et al/Model_succfed_WithREs.stan")

fit3_a <- readRDS("data/Entomological model fits/rds files/Combined/Feeding attempt/Succfed_ew_fit.rds")
# fit3_a <- readRDS("data/Entomological model fits/rds files/Huts separate/Feeding attempt/Succfed_west_fit_corrected.rds")
# fit3_a <- readRDS("data/Entomological model fits/rds files/Huts separate/Feeding attempt/Succfed_east_fit_corrected.rds")

fit3_a_fit <- extract(fit3_a, permuted = TRUE)

fit3_a_f <- fit3_a_fit$a[c(data_picker)]
fit3_a_g <- fit3_a_fit$b[c(data_picker)]

quantile(fit3_a_fit$a,0.5)
quantile(fit3_a_fit$b,0.5)

feeding = (1 - (exp(fit3_a_g[1] * (1 - exp(fit3_a_f[1] * hut_surv_LL))))/fit3_a_f[1])
plot(feeding ~ hut_surv_LL,ylim=c(0,1),xlim=c(0,1))
for(i in 1:1000){
  hut_surv_LL <- 1 - f_LOG_logistic(x = seq(0,1,0.01), LL_a_full[i], LL_c_full[i])
  hut_fed_LL <- (1 - (exp(fit3_a_g[i] * (1 - exp(fit3_a_f[i] * hut_surv_LL))))/fit3_a_f[i])
  lines(hut_fed_LL ~ hut_surv_LL)
  
}


##
## ASSOCIATION 5 MORTALITY TO halflife
## THIS CAN BE EITHER (GENERIC TO 'NETS')
##                    (SPECIFIC TO 'NET CLASS': Pyr-ITN or Pyr-PBO-ITN)
##                    (PRODUCT SPECIFIC)
# load fit
#full_model<- stan_model("ento analysis/Model_halflife_0.stan") # flex params
#saveRDS(fit_full, paste0("ento analysis/statistical fits/ento_fit1_half_life_",dataset,".rds"))
# setwd("C:/Users/esherrar/Documents/Rprojects/LLINEUP_Uganda/ento analysis/statistical fits/")
# fit_a <- readRDS("ento_fit1_half_life_a.rds") ## All data
# fit_b <- readRDS("ento_fit1_half_life_b.rds")

# testa = rstan::extract(fit_a)
# testb = rstan::extract(fit_b)
# for half life we are continuing initially with Tom's eLife estimates'

# names(test)
# sim_X = seq(0,1,0.01)
# P_hl1 <- 1/(1 + exp(-(mean(test$a) + mean(test$b) * sim_X)))
# P_hl1_low <- 1/(1 + exp(-(quantile(test$a,0.025) + quantile(test$b,0.025) * sim_X)))
# P_hl1_upp <- 1/(1 + exp(-(quantile(test$a,0.975) + quantile(test$b,0.975) * sim_X)))
# plot(P_hl1 ~ sim_X)


resistance_ITN_default_params_2_f = function(shape, product, resistance_range,
                                             ## inputs
                                             LL_a_full,LL_c_full,
                                             fitbene_1,fitbene_2,
                                             fit1_c,fit1_d,fit1_e,
                                             fit3_f,fit3_g){
  ## The parameters included here are from the statistical analysis 
  ## following on from Nash et al 2021 but focusing on PermaNet and Olyset Nets
  ## These include all West and East Africa hut design experimental huts
  ## The associations between mortality for the different net types
  ## DEFAULT PARAMETERS
  
  
  #Assay to hut mortality conversion - median estimates	
  param_a = LL_a_full ## 0.89 ## All 0.89 ## West 0.72 ## East 1.05
  param_b = LL_c_full ## 0.47 ## All 0.47 ## West 0.88 ## East 0.36
  
  #Deterency from mortality		
  param_c = fit1_c  ## 2.57 ## All 2.57 ## West 3.04 ## East 1.76
  param_d = fit1_d  ## 0.49 ## All 0.49 ## West 0.62 ## East 0.10
  param_e = fit1_e  ## 0.36 ## All 0.36 ## West 0.39 ## East 0.44
  
  #Success from mortality		
  param_f = fit3_f ## 4.66 ## All 4.66 ## West 5.18 ## East 3.42
  param_g = fit3_g ## 0.04 ## All 0.04 ## West 0.03 ## East 0.01
  
  ## AND THESE ARE THE ENTIRE DATABASE FOR NASH UNCERTAINTY
  ## I AM USING ONLY THOSE FOR THE NET TYPES IN THE UGANDA TRIAL
  #a = 0.89 (0.27 - 2.23)
  #b = 0.47 (0.38 - 0.56)
  
  #c = 2.57 (0.27 - 5.65)
  #d = 0.49 (0.06 - 1.28)
  #e = 0.36 (0.26 - 0.46)
  
  #f = 4.66 (4.48 - 4.84)
  #g = 0.04 (0.03 - 0.05)
  
  
  #Decay in insecticide non-PBO net		
  mup =	-2.429 #mu_p ## sample(test$a,size=1) ##array(c(rep(-2.429,3),rep(-2.984025,3),rep(-1.866,3)),c(3,3)) ## ... gam.medians[9]
  rhop = -3.007# rho_p ##  sample(test$b,size=1)	##array(c(rep(-3.007,3),rep(-3.74,3),rep(-2.295,3)),c(3,3)) ## ... gam.medians[10]
  
  
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
  mort_huta   = if(shape =="log-log") 1/(1+(mort_assay/param_a)^-param_b) else if(shape=="logistic") 1/(1+exp(param_a*param_b-param_a*mort_assay))
  
  ## Next, we define the additional benefit from the next generation ITNs
  ## These are worked out in the ["prioritisations"] manuscript
  
  
  ## pyrethroid-PBO ITNs
  # PBO_benefit = 1 / (1 + exp(-5.603 * mort_huta -  -1.431))
  # PBO_benefit_log_upp = 1 / (1 + exp(-5.252 * mort_huta - -1.538))
  # PBO_benefit_log_low = 1 / (1 + exp(-5.924 * mort_huta -1.320))
  # 
  PBO_benefit = 1 / (1 + exp(-fitbene_2 * mort_huta -  fitbene_1))
  # PBO_benefit_log_upp = 1 / (1 + exp(-quantile(fitbene_2,0.975) * mort_huta - quantile(fitbene_1,0.975)))
  # PBO_benefit_log_low = 1 / (1 + exp(-quantile(fitbene_2,0.025) * mort_huta - quantile(fitbene_1,0.025)))
  # 
  ## pyrethroid-chlorfenapyr ITNs
  # G2_benefit = 1 / (1 + exp(-4.6637 * mort_huta -  -0.5537 ))
  # G2_benefit_log_upp = 1 / (1 + exp(-5.4243 * mort_huta - 0.3408))
  # G2_benefit_log_low = 1 / (1 + exp(-3.9817 * mort_huta - -0.7598))
  # 
  
  ## Now we work through the probability steps to determine the key input parameters for the model
  ## These probability relationships are determined by Rebecca Nash, Ben and Tom see email notes above
  
  #specify whichever net is used in the RCT
  mort_hut = if(product==0) mort_huta else if(product==1) PBO_benefit else if(product==2) G2_benefit 
  
  
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
  
  kp1 = ifelse(kp1 > kp0,kp0,kp1) ## Capping impact so max feeding is no bigger than assumed 
  ## max feeding for no interventions (kp0 = 0.699, Griffin et al 2010)
  ## (time = 0 time steps after net implementation)
  r_ITN0  = (1-kp1/kp0)*(jp1/(lp1+jp1))	#probability of repeating behaviour
  d_ITN0  = (1-kp1/kp0)*(lp1/(lp1+jp1))	#probability of dying with an encounter with ITN
  s_ITN0  = 1-d_ITN0-r_ITN0             #probability of successfully feeding (surviving and feeding)
  
  
  ## Repeat these to determin the maximum and minimum effects which combine to help determine ITN half life
  mort_maxA   = if (shape=="log-log") 1/(1+((1)/param_a)^-param_b) else if (shape=="logistic") 1/(1+exp(param_a*param_b-param_a*(1)))
  
  mort_minA   = if (shape=="log-log") 1/(1+((0)/param_a)^-param_b) else if (shape=="logistic") 1/(1+exp(param_a*param_b-param_a*(0)))
  
  
  
  ## pyrethroid-PBO ITNs  ********************** update to not the mean
  PBO_benefitA = 1 / (1 + exp(-fitbene_2 * mort_maxA -  fitbene_1))
  
  ## pyrethroid-PBO ITNs  ********************** update to not the mean
  PBO_benefitB = 1 / (1 + exp(-fitbene_2 * mort_minA -  fitbene_1))
  
  mort_max = if(product==0) mort_maxA else if(product==1) PBO_benefitA else if(product==2) G2_benefitA
  mort_min = if(product==0) mort_minA else if(product==1) PBO_benefitB else if(product==2) G2_benefitB
  
  
  #{halflife}
  my_max_washes_a = mup +rhop*(mort_max-0.5)		
  my_max_washes   = log(2)/(exp(my_max_washes_a)/(1+exp(my_max_washes_a)))
  
  
  ## Uncertainty
  net_half_life_min = 2
  net_half_life_max = 3
  
  wash_decay_rate_a = mup +rhop*(mort_hut-0.5)
  wash_decay_rate   = log(2)/(exp(wash_decay_rate_a)/(1+exp(wash_decay_rate_a)))
  itn_half_life     = wash_decay_rate/my_max_washes*net_halflife
  itn_half_life_max     = wash_decay_rate/my_max_washes*net_half_life_max  
  itn_half_life_min     = wash_decay_rate/my_max_washes*net_half_life_min
  ##No need to re-adjust these anymore
  ##Final Parameter estimates for the transmission model
  ERG_d_ITN0 <- d_ITN0
  ERG_s_ITN0 <- s_ITN0
  ERG_r_ITN0 <- 1-ERG_d_ITN0-ERG_s_ITN0
  
  ## Print out these estimates to a data.frame as the function output
  uncertainty_resistance_params_nets = data.frame(ERG_d_ITN0,ERG_r_ITN0,itn_half_life,itn_half_life_min,itn_half_life_max)
  
  return(uncertainty_resistance_params_nets)
}


## To include some uncertainty around estimates
se = function(data){
  x = sd(data) / sqrt(1000) 
  return(x)
}

uncertainty_fn = function(param_best){
  se_dat = se(rnorm(n = n_uncertainty, mean = param_best, sd = 0.5))
  N_dat = (param_best*(1-param_best)/se_dat^2)-1
  
  alpha = param_best * N_dat
  beta = N_dat - alpha
  
  reps_out = rbeta(n = n_uncertainty,alpha,beta)
  
  return(reps_out)
}



data_picker2 =  sample(1:1000,size = 100,replace=FALSE)
extract_uncertainty_nets_f = function(res,net_product){
  
  nets_ALL = array(dim=c(1000,5))
  for(i in 1:1000){
    nets_ALL[i,] = as.numeric(resistance_ITN_default_params_2_f(shape = "logistic",
                                                                product = net_product, ## Always pyrethroid LLINs
                                                                resistance_range = res,
                                                                LL_a_full=LL_a_full[i],
                                                                LL_c_full=LL_c_full[i],
                                                                fitbene_1=fitbene_1[i],
                                                                fitbene_2=fitbene_2[i],
                                                                fit1_c=fit1_a_c[i],
                                                                fit1_d=fit1_a_d[i],
                                                                fit1_e=fit1_a_e[i],
                                                                fit3_f=fit3_a_f[i],
                                                                fit3_g=fit3_a_g[i]))
    
  }
  
  nets_ALL = nets_ALL[order(nets_ALL[,1]), ]
  
  ## Returns the median, 95% Credible intervals and mean estimates
  return(list(data.frame(net1_dITN = c(nets_ALL[c(500,50,950),1],mean(nets_ALL[,1])),
                         net1_rITN = c(nets_ALL[c(500,50,950),2],mean(nets_ALL[,2])),
                         net1_hlf = c(nets_ALL[500,3],nets_ALL[50,4],nets_ALL[950,5],2.64)),
              B_all_nets = data.frame(nets_ALL[data_picker2,])
  )
  
  
  )
  
}
RESIST0 = extract_uncertainty_nets_f(res = 0,net_product = 0) ## 0 for standard nets, 1 for pbo
RESIST0.2 = extract_uncertainty_nets_f(res = 0.2,net_product = 0)
RESIST0.4 = extract_uncertainty_nets_f(res = 0.4,net_product = 0)
RESIST0.6 = extract_uncertainty_nets_f(res = 0.6,net_product = 0)
RESIST0.8 = extract_uncertainty_nets_f(res = 0.8,net_product = 0)
RESIST1 = extract_uncertainty_nets_f(res = 1,net_product = 0)





res = seq(0,1,0.01)
nets_ALL = array(dim=c(1000,5,length(res)))
for(j in 1:length(res)){
  for(i in 1:1000){
    nets_ALL[i,,j] = as.numeric(resistance_ITN_default_params_2_f(shape = "log-log",
                                                                  product = 0, ## Always pyrethroid LLINs
                                                                  resistance_range = res[j],
                                                                  LL_a_full=LL_a_full[i],
                                                                  LL_c_full=LL_c_full[i],
                                                                  fitbene_1=fitbene_1[i],
                                                                  fitbene_2=fitbene_2[i],
                                                                  fit1_c=fit1_a_c[i],
                                                                  fit1_d=fit1_a_d[i],
                                                                  fit1_e=fit1_a_e[i],
                                                                  fit3_f=fit3_a_f[i],
                                                                  fit3_g=fit3_a_g[i]))
    
  }
  
}

for(j in 1:length(res)){
  nets_ALL[,,j] = nets_ALL[order(nets_ALL[,1,j]),,j]
}

to_ret_halflife = to_ret_dITN = to_ret_rITN = array(dim=c(length(res),4))
for(j in 1:length(res)){
  to_ret_dITN[j,] = c(nets_ALL[c(500,50,950),1,j],mean(nets_ALL[,1,j]))
  to_ret_rITN[j,] = c(nets_ALL[c(500,50,950),2,j],mean(nets_ALL[,2,j]))
  to_ret_halflife[j,] = c(nets_ALL[500,3,j],nets_ALL[50,4,j],nets_ALL[950,5,j],NA)
}

params = cbind(to_ret_dITN,to_ret_rITN,to_ret_halflife)
colnames(params) = c("median_dITN","min_dITN","max_dITN","mean_dITN",
                     "median_rITN","min_rITN","max_rITN","mean_rITN",
                     "median_hl","min_hl","max_hl","mean_hl")
write.csv(params,"data/uncertainty_parameters_pyr_nets.csv")




res = seq(0,1,0.01)
nets_ALL = array(dim=c(1000,5,length(res)))
for(j in 1:length(res)){
  for(i in 1:1000){
    nets_ALL[i,,j] = as.numeric(resistance_ITN_default_params_2_f(shape = "log-log",
                                                                  product = 1, ## Always pyrethroid LLINs
                                                                  resistance_range = res[j],
                                                                  LL_a_full=LL_a_full[i],
                                                                  LL_c_full=LL_c_full[i],
                                                                  fitbene_1=fitbene_1[i],
                                                                  fitbene_2=fitbene_2[i],
                                                                  fit1_c=fit1_a_c[i],
                                                                  fit1_d=fit1_a_d[i],
                                                                  fit1_e=fit1_a_e[i],
                                                                  fit3_f=fit3_a_f[i],
                                                                  fit3_g=fit3_a_g[i]))
    
  }
  
}

for(j in 1:length(res)){
  nets_ALL[,,j] = nets_ALL[order(nets_ALL[,1,j]),,j]
}

to_ret_halflife = to_ret_dITN = to_ret_rITN = array(dim=c(length(res),4))
for(j in 1:length(res)){
  to_ret_dITN[j,] = c(nets_ALL[c(500,50,950),1,j],mean(nets_ALL[,1,j]))
  to_ret_rITN[j,] = c(nets_ALL[c(500,50,950),2,j],mean(nets_ALL[,2,j]))
  to_ret_halflife[j,] = c(nets_ALL[500,3,j],nets_ALL[50,4,j],nets_ALL[950,5,j],NA)
}

params = cbind(to_ret_dITN,to_ret_rITN,to_ret_halflife)
colnames(params) = c("median_dITN","min_dITN","max_dITN","mean_dITN",
                     "median_rITN","min_rITN","max_rITN","mean_rITN",
                     "median_hl","min_hl","max_hl","mean_hl")
write.csv(params,"data/uncertainty_parameters_pbo_nets.csv")
