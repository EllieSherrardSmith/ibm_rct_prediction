## Functions and script to estimate the IRS impact

## Taking data and analysis from Sherrard-Smith et al 2018

## These are provided in the data/intervention_data folder. 

resistance_IRS_product_f = function(product, resistance_range){


  ## IRS parameters and pyrethroid IRS effect with resistance
  
  ## Bioassay mortality relationships
  ##Using the relationships defined as other chemistries and not accounting for Controls
  alpha1 = -1.059923 ##THIS IS IF DOING DATA REARRANGING AS FUNCTION imp_list_cleaner_f -1.20088# -1.942755
  alpha2 = 0.02387883##0.01878267# 0.03361892
  beta1 =  0.7674771 ##0.8682404# 1.849397
  beta2 =  -0.03166185 ##-0.01735203# -0.04921294
  theta1 = -1.673563 ##-1.672587# -2.071873
  theta2 = -0.0007511497 ##-0.0007890287# 0.02004906
  
  # alpha1qu = -1.1387719 ## -1.285526 # -2.112937
  # alpha2qu =0.02226354 ##  0.01705237#0.03114543 
  # beta1qu = 0.6885264  ##0.7894646  # 1.706689 
  # beta2qu = -0.03346408 ##-0.01889573  # -0.05236427 
  # theta1qu =-1.794143 ##-1.777039 # -2.225711
  # theta2qu = -0.003069259  ##-0.002947599   #0.01760429
  # 
  # alpha1ql = -0.9785583##  -1.098071 #-1.793964 
  # alpha2ql =0.02556298##   0.02030704#0.03624486
  # beta1ql = 0.8468829 ## 0.9470872 # 1.996350
  # beta2ql =  -0.02982156 ##  -0.01598732 #-0.04645465
  # theta1ql = -1.566623 ##  -1.566337# -1.908551 
  # theta2ql = 0.001473022 ##0.001528303 # 0.02252894
  
  # ##If AandB              UPDATE THESE...***** ONCE PYRETHROID INDIVIDUAL STUDIES HAVE RUN
  int1 = -2.587528
  int2 =  -0.002989876
  int3 =  -2.95455
  int4 =   0.01200367
  int5 =  -3.917539
  int6 = 0.002780567 ##same as mortality decay...
  
  grad1 =5.777369
  grad2 = -0.01362388
  grad3 = 5.231271
  grad4 =-0.004870386
  grad5 =8.542869
  grad6 =  0.005316776
  
  ## Estimating parameters where there is pyrethroid resistance
  ## For any given level of mortality in a bioassay (x) where x=1 indicates all mosquitoes die,
  ## these are linked by the relationship:
  ##  e.g.
  
  
  time = 1:365
  
  x=1 - resistance_range ##when x is 1 all mosquitoes die
  
  # death = feed = rep = deter = array(dim=c(365,length(x)))
  
  irs_params = expand.grid(res_range = resistance_range)
  
  for(i in 1:length(x)){
    temp_mort_1 = grad1 * (1 / (1 + exp(-alpha1 - alpha2 * 100 * x[i]))) + int1
    temp_mort_2 = grad2 * (1 / (1 + exp(-alpha1 - alpha2 * 100 * x[i]))) + int2
    
    # death[,i] = 1 / (1 + exp(-(temp_mort_1+ temp_mort_2 * time)))
    
    temp_succ_1 = grad3 * (1 / (1 + exp(-beta1 - beta2 * 100 * x[i]))) + int3
    temp_succ_2 = grad4 * (1 / (1 + exp(-beta1 - beta2 * 100 * x[i]))) + int4
    
    # feed[,i] = 1 / (1 + exp(-(temp_succ_1 + temp_succ_2 * time)))
    
    temp_det_1 = grad5 * (1 / (1 + exp(-theta1 - theta2 * 100 * x[i]))) + int5
    temp_det_2 = grad6 * (1 / (1 + exp(-theta1 - theta2 * 100 * x[i]))) + int6
    
    # deter[,i] = 1 / (1 + exp(-(temp_det_1 + temp_mort_2 * time)))
    
    irs_params[i,2] = temp_mort_1
    irs_params[i,3] = temp_mort_2
    irs_params[i,4] = temp_succ_1
    irs_params[i,5] = temp_succ_2
    irs_params[i,6] = temp_det_1
    irs_params[i,7] = temp_mort_2
  }
  
  colnames(irs_params) = c("res_range",
                           "irs_decay_mort1",
                           "irs_decay_mort2",
                           "irs_decay_succ1",
                           "irs_decay_succ2",
                           "irs_decay_det1",
                           "irs_decay_det2")
  
                              ##pyrethroid         ## Actellic                             ## Sumi                               ##Bendio                                ## DDT (assumed to be as good as pyrethroid with no res in absence of data)          
  mort1 = data.frame(PYR=irs_params[,2], ACT=rep(2.024412039,nrow(irs_params)),  SUM=rep(1.735248,nrow(irs_params)),   BEN=rep(1.094161249,nrow(irs_params)),  DDT=rep(irs_params[1,2],nrow(irs_params)))
  mort2 = data.frame(PYR=irs_params[,3], ACT=rep(-0.008887407,nrow(irs_params)), SUM=rep(-0.01318361,nrow(irs_params)),BEN=rep(-0.024488819,nrow(irs_params)), DDT=rep(irs_params[1,3],nrow(irs_params)))
  succ1 = data.frame(PYR=irs_params[,4], ACT=rep(-2.222632837,nrow(irs_params)), SUM=rep(-2.374664,nrow(irs_params)),  BEN=rep(-1.267712211,nrow(irs_params)), DDT=rep(irs_params[1,4],nrow(irs_params)))
  succ2 = data.frame(PYR=irs_params[,5], ACT=rep(0.008469928,nrow(irs_params)),  SUM=rep(0.01146379,nrow(irs_params)), BEN=rep(0.019556452,nrow(irs_params)),  DDT=rep(irs_params[1,5],nrow(irs_params)))
  det1 =  data.frame(PYR=irs_params[,6], ACT=rep(1.231733302,nrow(irs_params)),  SUM=rep(-2.643243,nrow(irs_params)),  BEN=rep(-1.700676448,nrow(irs_params)), DDT=rep(irs_params[1,6],nrow(irs_params)))
  det2 =  data.frame(PYR=irs_params[,7], ACT=rep(-0.008887407,nrow(irs_params)), SUM=rep(-0.01318361,nrow(irs_params)),BEN=rep(-0.024488819,nrow(irs_params)), DDT=rep(irs_params[1,7],nrow(irs_params)))
  

  product = product ## the column number matches the respective produc
  IRS_params = data.frame(irs_decay_mort1 = mort1[,product],irs_decay_mort2 = mort2[,product],
                          irs_decay_succ1 = succ1[,product],irs_decay_succ2 = succ2[,product],
                          irs_decay_det1 = det1[,product],irs_decay_det2 = det2[,product])
  
  return(IRS_params)
  
  
}

resistance_IRS_product_f(product = 1,            ## 1 is pyrethroid products, 2 is Actellic, 2 is Sumishield, 4 is bendiocarb, 5 is DDT 
                         resistance_range = seq(0,1,0.2)) ##


resistance_IRS_product_f(product = 2,            ## 1 is pyrethroid products, 2 is Actellic, 2 is Sumishield, 4 is bendiocarb, 5 is DDT 
                         resistance_range = seq(0,1,0.2)) ##


resistance_IRS_product_f(product = 3,            ## 1 is pyrethroid products, 2 is Actellic, 2 is Sumishield, 4 is bendiocarb, 5 is DDT 
                         resistance_range = seq(0,1,0.2)) ##


resistance_IRS_product_f(product = 4,            ## 1 is pyrethroid products, 2 is Actellic, 2 is Sumishield, 4 is bendiocarb, 5 is DDT 
                         resistance_range = seq(0,1,0.2)) ##


resistance_IRS_product_f(product = 5,            ## 1 is pyrethroid products, 2 is Actellic, 2 is Sumishield, 4 is bendiocarb, 5 is DDT 
                         resistance_range = seq(0,1,0.2)) ##
