#


## Bradley ID 1
## U-Bradley

# ## CHECK THESE
# interventions_file	interventions.txt
# sim_file	sim_parms_U.txt
# mosq_file	mosq_species_parms.txt
# corr_file	corr_file_0.txt
# itn_file	itn_parms_2010.txt
# irs_file	irs_parms_ellie_pyrethroid.txt
# mda_file	mda_parms.txt
# vacc_file	vacc_parms.txt
# output_file	output_vars_U.txt
# model_file	model_parms.txt

functions_model_runs_bradley<-function(DIDE_CODE,
                                       arm,run,
                                       
                                       ITN_1, 
                                       IRS_1, 
                                       
                                       itn_repel_fun,itn_repel_gam,itn_repel_ara,
                                       itn_kill_fun,itn_kill_gam,itn_kill_ara,
                                       itn_half_life,
                                       
                                       irs_decay_mort1,irs_decay_mort2,
                                       irs_decay_succ1,irs_decay_succ2,
                                       irs_decay_det1,irs_decay_det2,
                                       
                                       total_MM,
                                       
                                       sites){
  
  arm <- arm
  run <- run
  
  pop_size<- 30000 #Sim_pop_size(site_file[site_file[,1]=='prev',2])
  total_M = total_MM
  
  Int_set_up<-paste('num_people', pop_size, ##This is arbitrary and simply speeds up or smooths out runs (start with small popn size to make sure it runs fine then up to about 30,000 - 80,000)
                    'total_M', total_M,
                    'itn_irs_corr',	1, ##same people receive nets and spray unless trial specifies otherwise
                    'output_type 0 recalculate 0 itn_usage 1 add itn 0 itn_start 0 continue_itn 0 itn_flexible 1', ## for calibration we need this to be output_type 0
                    
                    'num_runs 1 itn_cov_0_0', ITN_1, 'itn_cov_1_0 0', ## year 0 is march 2013
                    'itn_cov_2_0 0 itn_cov_3_0', ITN_1, 'itn_cov_4_0 0 itn_cov_5_0 0 itn_cov_6_0 0 itn_cov_7_0 0 itn_cov_8_0 0 itn_cov_9_0 0 itn_cov_10_0 0',
                    'itn_cov_11_0 0 itn_cov_12_0 0 itn_cov_13_0 0 itn_cov_14_0 0 itn_cov_15_0 0 itn_cov_16_0 0 itn_cov_17_0 0 itn_cov_18_0 0 itn_cov_19_0 0',
                    
                    ## Assuming resistance is unchanging over time *and was present pre the trials at whatever rate determined
                    'itn_repel_fun_1', itn_repel_fun, 'itn_repel_gamb_ss_1', itn_repel_gam, 'itn_repel_arab_1', itn_repel_ara,
                    'itn_kill_fun_1', itn_kill_fun, 'itn_kill_gamb_ss_1', itn_kill_gam, 'itn_kill_arab_1', itn_kill_ara,
                    'itn_half_life_1',itn_half_life,
                    
                    'irs 1 irs_coverage', IRS_1, 'irs_start 1 irs_max_rounds 2 irs_frequency 1 irs_offset_absolute 1 irs_offset 0.25',
                    'change_irs 1 change_irs_time 1',
                    'irs_decay_mort1', irs_decay_mort1, 'irs_decay_mort2', irs_decay_mort2,
                    'irs_decay_succ1', irs_decay_succ1,'irs_decay_succ2', irs_decay_succ2,
                    'irs_decay_det1', irs_decay_det1,'irs_decay_det2', irs_decay_det2,
                    
                    'irs_ellie 1')
  
  Options<-paste(Int_set_up)
  
  sites <- sites
  draw<-0
  
  ## Run the simulation
  Model_launcher(OutputName = paste("VALIDATION_trial_arm",arm,"net_1",draw,run, sep = "_"),
                 OutputRoot = paste0("P:/validating_modelling/output_files/FINAL_netparm2/bradley/arm", arm),
                 Options=Options,
                 Exe = "P:/validating_modelling/bin/irs_restructure.exe",
                 Root="P:/validating_modelling/input_files",
                 Demog = "demog_GADM18_101013.txt",  ## Match up for  Bioko Island, North, Equatorial Guinea (using Eq Guinea... no site file for bioko)
                 Pop = "pop_GADM18_101013.txt",        ## Match up for  Bioko Island, North, Equatorial Guinea
                 Site=sites,
                 Parameter_draw = draw,
                 Return_output = FALSE)
}


#############################
##
## Run the model 2


## CHECK THE CORRECT INPUT FILES ARE ADDED
# sim_file	sim_parms_P.txt               ##Check year of trial, check age group of prevalence cohort etc unless on command line
# mosq_file	mosq_species_parms.txt        ##Check Q0 and phi unless in site file or command lines
# corr_file	corr_file_0.txt               ##check correlated
# itn_file	itn_parms_2010_T.txt          ##check on command line and input folder
# irs_file	irs_parms_ellie_pyrethroid.txt  ##use default pyrethroid and then switch up to actellic as pyrethroid in 2014
# mda_file	mda_parms.txt                 ##prob wont use - interventions.txt drug 0 
# vacc_file	vacc_parms.txt                ##wont use - interventions.txt epi_[] 0
# output_file	output_vars_P.txt           ##specify for the trial (Incidence 0 to 5 year olds; prevalence all age)
# model_file	model_parms.txt             ##will not change this

## Function to run the model across parameter draws 
## Site data, Zambesia Province 101979
## Match site data to that estimated from the study for mosquito outdoor and human biting.


chac_trial_f<-function(DIDE_CODE, ## Zambezia list 1:2000
                       
                       arm,
                       run, ## either 1 for ITN only, or 2 for ITN+IRS
                       
                       # ITN_1, ## this is the cover prior to the net distribution campaign
                       ITN_2, ## this is the universal net distribution in June July 2017
                       
                       IRS_1, ## this is the Oct/Nov 2017
                       
                       itn_repel_fun,itn_repel_gam,itn_repel_ara, ## figure out the pyrethroid resistance estimates using funestus 2018 data
                       itn_kill_fun,itn_kill_gam,itn_kill_ara,
                       itn_half_life,
                       
                       irs_decay_mort1,irs_decay_mort2, ## actellic so assuming no resistance
                       irs_decay_succ1,irs_decay_succ2,
                       irs_decay_det1,irs_decay_det2,
                       
                       total_MM,
                       # IRS_2, ## this is the Oct/Nov 2017 assuming the same cover of IRS in second spray campaign
                       sites){
  
  arm <- arm
  run <- run
  
  pop_size<- 30000 #Sim_pop_size(site_file[site_file[,1]=='prev',2])
  
  Int_set_up<-paste('num_people', pop_size, 'total_M', total_MM,##This is arbitrary and simply speeds up or smooths out runs (start with small popn size to make sure it runs fine then up to about 30,000 - 80,000)
                    'itn_irs_corr',	1, ##same people receive nets and spray unless trial specifies otherwise
                    
                    # 'output_type 0 recalculate 0 itn_usage 1 add itn 0 itn_start 1.5 continue_itn 0 itn_flexible 1',
                    'itn_leave_dur 3',
                    
                    'output_type 0 recalculate 0 itn_usage 1 add itn 0 itn_start 0 continue_itn 0 itn_flexible 1', ## for calibration we need this to be output_type 0
                    
                    'num_runs 1 itn_cov_0_0 0.05 itn_cov_1_0', ITN_2, ## year 0 is june 2016
                    'itn_cov_2_0 0 itn_cov_3_0 0 itn_cov_4_0 0 itn_cov_5_0', ITN_2, 'itn_cov_6_0 0 itn_cov_7_0 0 itn_cov_8_0 0 itn_cov_9_0 0 itn_cov_10_0 0',
                    'itn_cov_11_0 0 itn_cov_12_0 0 itn_cov_13_0 0 itn_cov_14_0 0 itn_cov_15_0 0 itn_cov_16_0 0 itn_cov_17_0 0 itn_cov_18_0 0 itn_cov_19_0 0',
                    
                    'itn_repel_fun_1', itn_repel_fun, 'itn_repel_gamb_ss_1', itn_repel_gam, 'itn_repel_arab_1', itn_repel_ara,
                    'itn_kill_fun_1', itn_kill_fun, 'itn_kill_gamb_ss_1', itn_kill_gam, 'itn_kill_arab_1', itn_kill_ara,
                    'itn_half_life_1',itn_half_life, ## This is constant given estimate of resistance
                    
                    # 'num_runs 1 itn_cov_0_0 0.05 itn_cov_1_0', ITN_2,
                    # 'itn_cov_2_0 0.05 itn_cov_3_0 0.05 itn_cov_4_0 0.05 itn_cov_5_0 0.05 itn_cov_6_0 0.05 itn_cov_7_0 0.05 itn_cov_8_0 0.05 itn_cov_9_0 0 itn_cov_10_0 0',
                    # 'itn_cov_11_0 0 itn_cov_12_0 0 itn_cov_13_0 0 itn_cov_14_0 0 itn_cov_15_0 0 itn_cov_16_0 0 itn_cov_17_0 0 itn_cov_18_0 0 itn_cov_19_0 0',
                    # 
                    'irs 1 irs_coverage', IRS_1, 'irs_start 0 irs_max_rounds 2 irs_offset_absolute 1 irs_offset 0.416667', ##This is for ITS to start in November 2016
                    'irs_2 1 irs_coverage_2', IRS_1, 'irs_start_2 0 irs_max_rounds_2 2', 
                    'change_irs 1 change_irs_time 0',
                    'irs_decay_mort1', irs_decay_mort1, 'irs_decay_mort2', irs_decay_mort2,
                    'irs_decay_succ1', irs_decay_succ1,'irs_decay_succ2', irs_decay_succ2,
                    'irs_decay_det1', irs_decay_det1,'irs_decay_det2', irs_decay_det2,
                    
                    'irs_ellie 1')
  
  Options<-paste(Int_set_up)
  
  sites <- sites
  draw<-0
  
  ## Run the simulation
  Model_launcher(OutputName = paste("VALIDATION_trial_arm",arm,"net_1",draw,run, sep = "_"),
                 OutputRoot = paste0("P:/validating_modelling/output_files/FINAL_netparm2/chaccour/arm",arm),
                 Options=Options,
                 Exe = "P:/validating_modelling/bin/irs_restructure.exe",
                 Root="P:/validating_modelling/input_files",
                 Demog = "demog_GADM18_101979.txt",  ## Match up for Zambezia
                 Pop = "pop_GADM18_101979.txt",        ## Match up for Zambezia
                 Site=sites,
                 Parameter_draw = draw,
                 Return_output = FALSE)
}


#############################
##
## Run the model

######################################
##
## Modelling Corbel 2012
# 3
## 2008 LLIN distribution permanet 2
## 2008 bendiocarb IRS sprayed every 8 months
## 2008 bendiocarb CTPS every 4 months

## under6 YEAR OLDS prevalence to hit 

## DELTA mortality 76.84% LLIN only at year 2008


#############################
##
## Run the model


## CHECK THE CORRECT INPUT FILES ARE ADDED
# sim_file	sim_parms_4.txt                 ##Check year of trial, check age group of prevalence cohort etc unless on command line
# mosq_file	mosq_species_parms.txt        ##Using matched to study G 
# corr_file	corr_file_0.txt               ##check correlated
# itn_file	itn_parms_2010.txt          ##check on command line and input folder
# irs_file	irs_parms_ellie_pyrethroid.txt  ##use study default and check no IRS on if not implemented
# mda_file	mda_parms.txt                 ##prob wont use - interventions.txt drug 0 
# vacc_file	vacc_parms.txt                ##wont use - interventions.txt epi_[] 0
# output_file	output_vars_4.txt           ##specify for the trial (prev 0 to 14 year olds)
# model_file	model_parms.txt             ##will not change this


#############################
##
## Calibrating to prevalence estimate (already done ... need to do for txt site files then translate into model)



## Function to run the model across parameter draws 
## Site data, Kagera Region 

corb_trial_f1<-function(DIDE_CODE,
                        run,
                        
                        itn_repel_fun,itn_repel_gam,itn_repel_ara,
                        itn_kill_fun,itn_kill_gam,itn_kill_ara,
                        itn_half_life,
                        
                        irs_decay_mort1,irs_decay_mort2,
                        irs_decay_succ1,irs_decay_succ2,
                        irs_decay_det1,irs_decay_det2,
                        ITN_1,
                        IRS_1, 
                        irs_freq,
                        arm,
                        total_MM,
                        sites){
  
  arm <- arm
  run <- run
  
  # Load the site_file file
  # site_file<-read.table(paste0('P:/Ellies_cool_model_folder2/model_files/sites/Africa_Sites_Ellie_Copy/Africa_sites_0/Tanz_Kagera_677_', site, '.txt'))
  ## THIS IS Tanz_Kagera_677_dbHGH_1.txt for the high baseline coverage
  ## THIS IS Tanz_Kagera_677_dbMID_1.txt for the median baseline coverage
  ## THIS IS Tanz_Kagera_677_dbLOW_1.txt for the low baseline coverage
  
  pop_size<- 20000 #Sim_pop_size(site_file[site_file[,1]=='prev',2])
  # total_M = total_MM
  # 
  Int_set_up<-paste('num_people', pop_size, ##This is arbitrary and simply speeds up or smooths out runs (start with small popn size to make sure it runs fine then up to about 30,000 - 80,000)
                    'total_M', total_MM,
                    'itn_irs_corr',	1, ##same people receive nets and spray unless trial specifies otherwise
                    
                    ## Assuming resistance is unchanging over time *and was present pre the trials at whatever rate determined
                    'itn_repel_fun', itn_repel_fun, 'itn_repel_gamb_ss', itn_repel_gam, 'itn_repel_arab', itn_repel_ara,
                    'itn_kill_fun', itn_kill_fun, 'itn_kill_gamb_ss', itn_kill_gam, 'itn_kill_arab', itn_kill_ara,
                    'itn_half_life',itn_half_life,
                    
                    'output_type 0 recalculate 0 itn_usage 1 add itn 0 itn_start 0 continue_itn 0 itn_flexible 1',## for calibration we need this to be output_type 0
                    
                    'num_runs 1 itn_cov_0_0', ITN_1, 'itn_cov_1_0 0', ## year 0 is 2008 and aiming for under 6 year old prev of c25-30%
                    'itn_cov_2_0 0 itn_cov_3_0', ITN_1, 'itn_cov_4_0 0 itn_cov_5_0 0 itn_cov_6_0 0 itn_cov_7_0 0 itn_cov_8_0 0 itn_cov_9_0 0 itn_cov_10_0 0',
                    'itn_cov_11_0 0 itn_cov_12_0 0 itn_cov_13_0 0 itn_cov_14_0 0 itn_cov_15_0 0 itn_cov_16_0 0 itn_cov_17_0 0 itn_cov_18_0 0 itn_cov_19_0 0',
                    
                    'irs 1 irs_coverage', IRS_1, 'irs_start 0 irs_max_rounds 20 irs_frequency', irs_freq,'irs_offset_absolute 1 irs_offset 0',
                    'change_irs 1 change_irs_time 0',
                    'irs_decay_mort1', irs_decay_mort1, 'irs_decay_mort2', irs_decay_mort2,
                    'irs_decay_succ1', irs_decay_succ1,'irs_decay_succ2', irs_decay_succ2,
                    'irs_decay_det1', irs_decay_det1,'irs_decay_det2', irs_decay_det2,
                    
                    'irs_ellie 1')
  
  Options<-paste(Int_set_up)
  
  sites <- sites
  draw<-0
  
  ## Run the simulation
  Model_launcher(OutputName = paste("VALIDATION_trial_arm",arm,"net_1",draw,run, sep = "_"),
                 OutputRoot = paste0("P:/validating_modelling/output_files/FINAL_netparm2/corbel/arm",arm),
                 Options=Options,
                 Exe = "P:/validating_modelling/bin/irs_restructure.exe",
                 Root="P:/validating_modelling/input_files",
                 Demog = "demog_GADM18_100198.txt",  ## Match up for  Atlantique, Benin
                 Pop = "pop_GADM18_100198.txt",        ## Match up for  Atlantique, Benin
                 Site=sites,
                 Parameter_draw = draw,
                 Return_output = FALSE)
}

## CHECK THE CORRECT INPUT FILES ARE ADDED
# sim_file	sim_parms_B.txt                 ##Check year of trial, check age group of prevalence cohort etc unless on command line
# mosq_file	mosq_species_parms.txt        ##Using matched to study G 
# corr_file	corr_file_0.txt               ##check correlated
# itn_file	itn_parms_2010.txt          ##check on command line and input folder
# irs_file	irs_parms_ellie_pyrethroid.txt  ##use study default and check no IRS on if not implemented
# mda_file	mda_parms.txt                 ##prob wont use - interventions.txt drug 0 
# vacc_file	vacc_parms.txt                ##wont use - interventions.txt epi_[] 0
# output_file	output_vars_B.txt           ##specify for the trial (prev 0 to 14 year olds)
# model_file	model_parms.txt             ##will not change this


#############################
##
## Calibrating to prevalence estimate (already done ... need to do for txt site files then translate into model)

## 4

## Function to run the model across parameter draws 
## Site data, Kagera Region 

curt_trial_f1<-function(DIDE_CODE,
                        run,
                        
                        itn_repel_fun,itn_repel_gam,itn_repel_ara,
                        itn_kill_fun,itn_kill_gam,itn_kill_ara,
                        itn_half_life,
                        
                        irs_decay_mort1,irs_decay_mort2,
                        irs_decay_succ1,irs_decay_succ2,
                        irs_decay_det1,irs_decay_det2,
                        
                        ITN_1,
                        IRS_1, 
                        
                        arm,
                        total_MM,
                        sites){
  
  arm <- arm
  run <- run
  
  # Load the site_file file
  # site_file<-read.table(paste0('P:/Ellies_cool_model_folder2/model_files/sites/Africa_Sites_Ellie_Copy/Africa_sites_0/Tanz_Kagera_677_', site, '.txt'))
  ## THIS IS Tanz_Kagera_677_dbHGH_1.txt for the high baseline coverage
  ## THIS IS Tanz_Kagera_677_dbMID_1.txt for the median baseline coverage
  ## THIS IS Tanz_Kagera_677_dbLOW_1.txt for the low baseline coverage
  
  pop_size<- 10000 #Sim_pop_size(site_file[site_file[,1]=='prev',2])
  # total_M = total_MM
  # 
  Int_set_up<-paste('num_people', pop_size, ##This is arbitrary and simply speeds up or smooths out runs (start with small popn size to make sure it runs fine then up to about 30,000 - 80,000)
                    'total_M', total_MM,
                    'itn_irs_corr',	1, ##same people receive nets and spray unless trial specifies otherwise
                    'output_type 0 recalculate 0 itn_leave 1 add itn 1 itn_start 0 continue_itn 0 itn_flexible 1',## for calibration we need this to be output_type 0
                    'itn_leave_dur 150 itn_frequency 3 itn_coverage', ITN_1,
                    
                    # 'num_runs 1 itn_cov_0_0',ITN_1,'itn_cov_1_0 0',## year 0 is 1995 and aiming for under 6 year old prev of c25-30%
                    # 'itn_cov_2_0 0 itn_cov_3_0',ITN_1,'itn_cov_4_0 0 itn_cov_5_0 0 itn_cov_6_0 0 itn_cov_7_0 0 itn_cov_8_0 0 itn_cov_9_0 0 itn_cov_10_0 0',
                    # 'itn_cov_11_0 0 itn_cov_12_0 0 itn_cov_13_0 0 itn_cov_14_0 0 itn_cov_15_0 0 itn_cov_16_0 0 itn_cov_17_0 0 itn_cov_18_0 0 itn_cov_19_0 0',
                    # 
                    ## Assuming resistance is unchanging over time *and was present pre the trials at whatever rate determined
                    'itn_repel_fun_1', itn_repel_fun, 'itn_repel_gamb_ss_1', itn_repel_gam, 'itn_repel_arab_1', itn_repel_ara,
                    'itn_kill_fun_1', itn_kill_fun, 'itn_kill_gamb_ss_1', itn_kill_gam, 'itn_kill_arab_1', itn_kill_ara,
                    'itn_half_life_1',itn_half_life,
                    
                    'irs 1 irs_coverage', IRS_1, 'irs_start 0 irs_max_rounds 2 irs_frequency 0.67 irs_offset_absolute 1 irs_offset 0',
                    'change_irs 1 change_irs_time 0',
                    'irs_decay_mort1', irs_decay_mort1, 'irs_decay_mort2', irs_decay_mort2,
                    'irs_decay_succ1', irs_decay_succ1,'irs_decay_succ2', irs_decay_succ2,
                    'irs_decay_det1', irs_decay_det1,'irs_decay_det2', irs_decay_det2,
                    
                    'irs_2 1 irs_coverage_2 0 irs_start_2 2 irs_max_rounds_2 20',
                    
                    'irs_ellie 1')
  
  Options<-paste(Int_set_up)
  
  sites <- sites
  draw<-0
  
  ## Run the simulation
  Model_launcher(OutputName = paste("VALIDATION_trial_arm",arm,"net_1",draw,run, sep = "_"),
                 OutputRoot = paste0("P:/validating_modelling/output_files/FINAL_netparm2/curtis/arm",arm),
                 Options=Options,
                 Exe = "P:/validating_modelling/bin/irs_restructure.exe",
                 Root="P:/validating_modelling/input_files",
                 Demog = "demog_GADM18_103243.txt",  ## Match up for  
                 Pop = "pop_GADM18_103243.txt",        ## Match up for  
                 Site=sites,
                 Parameter_draw = draw,
                 Return_output = FALSE)
}


######################################
## 
## Modelling D'Alessandro ID 5

## No nets distribution permanet 2
## ITN dip nets 
## under6 YEAR OLDS prevalence to hit 


## CHECK THE CORRECT INPUT FILES ARE ADDED
# sim_file	sim_parms_C.txt                 ##Check year of trial, check age group of prevalence cohort etc unless on command line
# mosq_file	mosq_species_parms.txt        ##Using matched to study G 
# corr_file	corr_file_0.txt               ##check correlated
# itn_file	itn_parms_2010.txt          ##check on command line and input folder
# irs_file	irs_parms_ellie_pyrethroid.txt  ##use study default and check no IRS on if not implemented
# mda_file	mda_parms.txt                 ##prob wont use - interventions.txt drug 0 
# vacc_file	vacc_parms.txt                ##wont use - interventions.txt epi_[] 0
# output_file	output_vars_T.txt           ##specify for the trial (prev 0 to 14 year olds)
# model_file	model_parms.txt             ##will not change this


#############################
##
## Calibrating to prevalence estimate (already done ... need to do for txt site files then translate into model)



## Function to run the model across parameter draws 
## Site data, Kagera Region 

DALESS_trial_f1<-function(DIDE_CODE,
                          run,
                          
                          itn_repel_fun,itn_repel_gam,itn_repel_ara,
                          itn_kill_fun,itn_kill_gam,itn_kill_ara,
                          itn_half_life,
                          
                          itn_repel_fun_1,itn_repel_gam_1,itn_repel_ara_1,
                          itn_kill_fun_1,itn_kill_gam_1,itn_kill_ara_1,
                          itn_half_life_1,
                          
                          ITN_1, 
                          
                          sitnum,arm,dide_code_match,
                          
                          total_MM,
                          sites){
  sitnum <- sitnum
  arm <- arm
  run <- run
  
  pop_size<- 10000 
  # 
  Int_set_up<-paste('num_people', pop_size, ##This is arbitrary and simply speeds up or smooths out runs (start with small popn size to make sure it runs fine then up to about 30,000 - 80,000)
                    'total_M', total_MM,
                    'itn_irs_corr',	1, ##same people receive nets and spray unless trial specifies otherwise
                    'output_type 0 recalculate 0 itn_usage 0 add itn 1 itn_start 0 continue_itn 0 itn_flexible 1',## for calibration we need this to be output_type 0
                    
                    'itn_coverage', ITN_1,
                    # 'num_runs 1 itn_cov_0_0 0 itn_cov_1_0', ITN_1, ## year 0 is JULY 1995 and aiming for under 6 year old prev of c25-30%
                    # 'itn_cov_2_0 0 itn_cov_3_0', ITN_1, 'itn_cov_4_0 0 itn_cov_5_0 0 itn_cov_6_0 0 itn_cov_7_0 0 itn_cov_8_0 0 itn_cov_9_0 0 itn_cov_10_0 0',
                    # 'itn_cov_11_0 0 itn_cov_12_0 0 itn_cov_13_0 0 itn_cov_14_0 0 itn_cov_15_0 0 itn_cov_16_0 0 itn_cov_17_0 0 itn_cov_18_0 0 itn_cov_19_0 0',
                    
                    ## Assuming resistance is unchanging over time *and was present pre the trials at whatever rate determined
                    'itn_repel_fun', itn_repel_fun, 'itn_repel_gamb_ss', itn_repel_gam, 'itn_repel_arab', itn_repel_ara,
                    'itn_kill_fun', itn_kill_fun, 'itn_kill_gamb_ss', itn_kill_gam, 'itn_kill_arab', itn_kill_ara,
                    'itn_half_life',itn_half_life,
                    
                    'change_itn 1 change_itn_time 0',
                    'itn_repel_fun_1', itn_repel_fun_1, 'itn_repel_gamb_ss_1', itn_repel_gam_1, 'itn_repel_arab_1', itn_repel_ara_1,
                    'itn_kill_fun_1', itn_kill_fun_1, 'itn_kill_gamb_ss_1', itn_kill_gam_1, 'itn_kill_arab_1', itn_kill_ara_1,
                    'itn_half_life_1',itn_half_life_1)
  
  Options<-paste(Int_set_up)
  
  draw<-0
  dide_code_match <- dide_code_match
  sites <- sites
  
  ## Run the simulation
  Model_launcher(OutputName = paste("VALIDATION_trial_arm",arm,sitnum,"net_1",draw,run, sep = "_"),
                 OutputRoot = paste0("P:/validating_modelling/output_files/FINAL_netparm2/dalessandro/arm",arm),
                 Options=Options,
                 Exe = "P:/validating_modelling/bin/irs_restructure.exe", ## The model executable
                 Root="P:/validating_modelling/input_files",
                 Demog = paste0("demog_GADM18_",dide_code_match,".txt"),  ## Match up for  The Gambia
                 Pop = paste0("pop_GADM18_",dide_code_match,".txt"),        ## Match up for The Gambia
                 Site=sites,
                 Parameter_draw = draw,
                 Return_output = FALSE)
}


################
##
## Henry ID 6

#############################
##
## Run the model


## CHECK THE CORRECT INPUT FILES ARE ADDED
# sim_file	sim_parms_E.txt                 ##Check year of trial, check age group of prevalence cohort etc unless on command line
# mosq_file	mosq_species_parms.txt        ##Using matched to study G 
# corr_file	corr_file_0.txt               ##check correlated
# itn_file	itn_parms_2010.txt          ##check on command line and input folder
# irs_file	irs_parms_ellie_pyrethroid.txt  ##use study default and check no IRS on if not implemented
# mda_file	mda_parms.txt                 ##prob wont use - interventions.txt drug 0 
# vacc_file	vacc_parms.txt                ##wont use - interventions.txt epi_[] 0
# output_file	output_vars_E.txt           ##specify for the trial (prev 0 to 14 year olds)
# model_file	model_parms.txt             ##will not change this


#############################
##
## Calibrating to prevalence estimate (already done ... need to do for txt site files then translate into model)



## Function to run the model across parameter draws 
## Site data, Kagera Region 

HENRY_trial_f1<-function(DIDE_CODE,
                         run,
                         
                         itn_repel_fun,itn_repel_gam,itn_repel_ara,
                         itn_kill_fun,itn_kill_gam,itn_kill_ara,
                         itn_half_life,
                         
                         ITN_1,
                         
                         arm,
                         total_MM,
                         sites){
  
  arm <- arm
  run <- run
  
  # Load the site_file file
  # site_file<-read.table(paste0('P:/Ellies_cool_model_folder2/model_files/sites/Africa_Sites_Ellie_Copy/Africa_sites_0/Tanz_Kagera_677_', site, '.txt'))
  ## THIS IS Tanz_Kagera_677_dbHGH_1.txt for the high baseline coverage
  ## THIS IS Tanz_Kagera_677_dbMID_1.txt for the median baseline coverage
  ## THIS IS Tanz_Kagera_677_dbLOW_1.txt for the low baseline coverage
  
  pop_size<- 20000 #Sim_pop_size(site_file[site_file[,1]=='prev',2])
  # total_M = total_MM
  # 
  Int_set_up<-paste('num_people', pop_size, ##This is arbitrary and simply speeds up or smooths out runs (start with small popn size to make sure it runs fine then up to about 30,000 - 80,000)
                    'total_M', total_MM,
                    'itn_irs_corr',	1, ##same people receive nets and spray unless trial specifies otherwise
                    'output_type 0 recalculate 0 itn_usage 1 add itn 0 itn_start 0 continue_itn 0 itn_flexible 1',## for calibration we need this to be output_type 0
                    
                    'num_runs 1 itn_cov_0_0', ITN_1, 'itn_cov_1_0', ITN_1, ## year 0 is 1995 and aiming for under 6 year old prev of c25-30%
                    'itn_cov_2_0 0 itn_cov_3_0 0 itn_cov_4_0 0 itn_cov_5_0 0 itn_cov_6_0 0 itn_cov_7_0 0 itn_cov_8_0 0 itn_cov_9_0 0 itn_cov_10_0 0',
                    'itn_cov_11_0 0 itn_cov_12_0 0 itn_cov_13_0 0 itn_cov_14_0 0 itn_cov_15_0 0 itn_cov_16_0 0 itn_cov_17_0 0 itn_cov_18_0 0 itn_cov_19_0 0',
                    
                    ## Assuming resistance is unchanging over time *and was present pre the trials at whatever rate determined
                    'itn_repel_fun_1', itn_repel_fun, 'itn_repel_gamb_ss_1', itn_repel_gam, 'itn_repel_arab_1', itn_repel_ara,
                    'itn_kill_fun_1', itn_kill_fun, 'itn_kill_gamb_ss_1', itn_kill_gam, 'itn_kill_arab_1', itn_kill_ara,
                    'itn_half_life_1',itn_half_life)
  
  Options<-paste(Int_set_up)
  
  sites <- sites
  draw<-0
  
  ## Run the simulation
  Model_launcher(OutputName = paste("VALIDATION_trial_arm",arm,"net_1",draw,run, sep = "_"),
                 OutputRoot = paste0("P:/validating_modelling/output_files/FINAL_netparm2/henry/arm",arm),
                 Options=Options,
                 Exe = "P:/validating_modelling/bin/irs_restructure.exe", ## The model executable
                 Root="P:/validating_modelling/input_files",
                 Demog = "demog_GADM18_100521.txt",  ## Match up for Korhogo, Savanes, Cote D'Ivoire
                 Pop = "pop_GADM18_100521.txt",        ## Match up for  Korhogo, Savanes, Cote D'Ivoire
                 Site=sites,
                 Parameter_draw = draw,
                 Return_output = FALSE)
}


#####################
##
## Kafy 7

## F-Kafy
functions_model_runs_kafy<-function(DIDE_CODE,
                                    arm,run,
                                    
                                    ITN_1, ITN_2, 
                                    IRS_1, IRS_2, IRS_3, 
                                    
                                    itn_leave_dur,
                                    itn_repel_fun,itn_repel_gam,itn_repel_ara,
                                    itn_kill_fun,itn_kill_gam,itn_kill_ara,
                                    itn_half_life,
                                    
                                    irs_decay_mort1,irs_decay_mort2,
                                    irs_decay_succ1,irs_decay_succ2,
                                    irs_decay_det1,irs_decay_det2,
                                    
                                    irs_decay_mort1_2,irs_decay_mort2_2,
                                    irs_decay_succ1_2,irs_decay_succ2_2,
                                    irs_decay_det1_2,irs_decay_det2_2,
                                    total_MM,
                                    
                                    sites){
  
  arm <- arm
  run <- run
  
  pop_size<- 30000 #Sim_pop_size(site_file[site_file[,1]=='prev',2])
  total_M = total_MM
  
  Int_set_up<-paste('num_people', pop_size, ##This is arbitrary and simply speeds up or smooths out runs (start with small popn size to make sure it runs fine then up to about 30,000 - 80,000)
                    'total_M', total_M,
                    'itn_irs_corr',	1, ##same people receive nets and spray unless trial specifies otherwise
                    'output_type 0 recalculate 0 itn_usage 1 add itn 0 itn_start 0 continue_itn 0 itn_flexible 1 itn_leave_dur',	itn_leave_dur, ## for calibration we need this to be output_type 0
                    
                    'num_runs 1 itn_cov_0_0', ITN_1, 'itn_cov_1_0 0', ## year 0 is 2011 and aiming for 2012 4-6 year old prev of 7 or 10%
                    'itn_cov_2_0 0 itn_cov_3_0', ITN_2, 'itn_cov_4_0 0 itn_cov_5_0 0 itn_cov_6_0 0 itn_cov_7_0 0 itn_cov_8_0 0 itn_cov_9_0 0 itn_cov_10_0 0',
                    'itn_cov_11_0 0 itn_cov_12_0 0 itn_cov_13_0 0 itn_cov_14_0 0 itn_cov_15_0 0 itn_cov_16_0 0 itn_cov_17_0 0 itn_cov_18_0 0 itn_cov_19_0 0',
                    
                    ## Assuming resistance is unchanging over time *and was present pre the trials at whatever rate determined
                    'itn_repel_fun_1', itn_repel_fun, 'itn_repel_gamb_ss_1', itn_repel_gam, 'itn_repel_arab_1', itn_repel_ara,
                    'itn_kill_fun_1', itn_kill_fun, 'itn_kill_gamb_ss_1', itn_kill_gam, 'itn_kill_arab_1', itn_kill_ara,
                    'itn_half_life_1',itn_half_life,
                    
                    'irs 1 irs_coverage', IRS_1, 'irs_start 0 irs_max_rounds 2 irs_frequency 0.5 irs_offset_absolute 1 irs_offset 0.45',
                    'change_irs 1 change_irs_time 0',
                    'irs_decay_mort1', irs_decay_mort1, 'irs_decay_mort2', irs_decay_mort2,
                    'irs_decay_succ1', irs_decay_succ1,'irs_decay_succ2', irs_decay_succ2,
                    'irs_decay_det1', irs_decay_det1,'irs_decay_det2', irs_decay_det2,
                    
                    'irs_2 1 irs_coverage_2', IRS_2, 'irs_start_2 1 irs_max_rounds_2 4 irs_frequency_2 0.5 irs_offset_absolute_2 1 irs_offset_2 0.45',
                    'change_irs_2 1 change_irs_time_2 1',
                    'irs_decay_mort1_2', irs_decay_mort1_2, 'irs_decay_mort2_2', irs_decay_mort2_2,
                    'irs_decay_succ1_2', irs_decay_succ1_2,'irs_decay_succ2_2', irs_decay_succ2_2,
                    'irs_decay_det1_2', irs_decay_det1_2,'irs_decay_det2_2', irs_decay_det2_2,
                    
                    'irs_3 1 irs_coverage_3', IRS_3, 'irs_start_3 2 irs_max_rounds_3 4 irs_frequency_3 0.5 irs_offset_absolute_3 1 irs_offset_3 0.45',
                    'change_irs_3 1 change_irs_time_3 2',
                    'irs_decay_mort1_2', irs_decay_mort1_2, 'irs_decay_mort2_2', irs_decay_mort2_2,
                    'irs_decay_succ1_2', irs_decay_succ1_2,'irs_decay_succ2_2', irs_decay_succ2_2,
                    'irs_decay_det1_2', irs_decay_det1_2,'irs_decay_det2_2', irs_decay_det2_2,
                    
                    'irs_4 1 irs_coverage_4 0 irs_start_4 3 irs_max_rounds_4 20',
                    
                    'irs_ellie 1')
  
  Options<-paste(Int_set_up)
  
  sites <- sites
  draw<-0
  
  ## Run the simulation
  Model_launcher(OutputName = paste("VALIDATION_trial_arm",arm,"net_1",draw,run, sep = "_"),
                 OutputRoot = paste0("P:/validating_modelling/output_files/FINAL_netparm2/kafy/arm", arm),
                 Options=Options,
                 Exe = "P:/validating_modelling/bin/irs_restructure.exe",
                 Root="P:/validating_modelling/input_files",
                 Demog = "demog_GADM18_102708.txt",  ## Match up for  Galabat, Al Qadarif Sudan
                 Pop = "pop_GADM18_102708.txt",        ## Match up for Galabat, Al Qadarif Sudan
                 Site=sites,
                 Parameter_draw = draw,
                 Return_output = FALSE)
}


#########################
## Marbiah ID 8

############################
##
## Run the model


## CHECK THE CORRECT INPUT FILES ARE ADDED
# sim_file	sim_parms_H.txt               ##Check year of trial, check age group of prevalence cohort etc unless on command line
# mosq_file	mosq_species_parms.txt        ##Check Q0 and phi unless in site file or command lines
# corr_file	corr_file_0.txt               ##check correlated
# itn_file	itn_parms_2010_H2.txt          ##check on command line and input folder
# irs_file	irs_parms_ellie_pyrethroid.txt  ##use default pyrethroid assuming no resistance
# mda_file	mda_parms.txt                 ##prob wont use - interventions.txt drug 0 
# vacc_file	vacc_parms.txt                ##wont use - interventions.txt epi_[] 0
# output_file	output_vars_H.txt           ##specify for the trial (Incidence 0 to 6 year olds; prevalence under 7s and under 6)
# model_file	model_parms.txt             ##will not change this

## Function to run the model across parameter draws 
## Site data, Zambesia Province 101979
## Match site data to that estimated from the study for mosquito outdoor and human biting.


marb_trial_f<-function(DIDE_CODE, ## Southern Province Sierre Leone 102761
                       
                       arm, ## 1 or 2
                       run, ## either 1 for control or 2 for ITNs
                       
                       ITN_1, ## this is 'all people have nets in the intervention...'
                       itn_repel_fun,itn_repel_gam,itn_repel_ara,
                       itn_kill_fun,itn_kill_gam,itn_kill_ara,
                       itn_half_life,
                       
                       sites){
  
  arm <- arm
  run <- run
  
  pop_size<- 20000 #Sim_pop_size(site_file[site_file[,1]=='prev',2])
  
  Int_set_up<-paste('num_people', pop_size, ##This is arbitrary and simply speeds up or smooths out runs (start with small popn size to make sure it runs fine then up to about 30,000 - 80,000)
                    'itn_irs_corr',	1, ##same people receive nets and spray unless trial specifies otherwise
                    
                    'itn_repel_fun', itn_repel_fun, 'itn_repel_gamb_ss', itn_repel_gam, 'itn_repel_arab', itn_repel_ara,
                    'itn_kill_fun', itn_kill_fun, 'itn_kill_gamb_ss', itn_kill_gam, 'itn_kill_arab', itn_kill_ara,
                    'itn_half_life',itn_half_life,
                    
                    # 'output_type 0 recalculate 0 itn_usage 1 add itn 0 itn_start 1.5 continue_itn 0 itn_flexible 1',
                    'output_type 0 recalculate 0 add itn 1 itn_start 0 itn_coverage', ITN_1,
                    'itn_max_rounds 1 itn_leave_dur 3',
                    
                    'irs 0',
                    
                    'irs_ellie 1')
  
  Options<-paste(Int_set_up)
  
  sites <- sites
  draw<-0
  
  ## Run the simulation
  Model_launcher(OutputName = paste("VALIDATION_trial_arm",arm,"net_1",draw,run, sep = "_"),
                 OutputRoot = paste0("P:/validating_modelling/output_files/FINAL_netparm2/marbiah/arm", arm),
                 Options=Options,
                 Exe = "P:/validating_modelling/bin/irs_restructure.exe",
                 Root="P:/validating_modelling/input_files",
                 Demog = "demog_GADM18_102761.txt",  ## Match up for SOUTHERN Province Sierre Leone
                 Pop = "pop_GADM18_102761.txt",        ## Match up for Southern Province Sierre Leone
                 Site=sites,
                 Parameter_draw = draw,
                 Return_output = FALSE)
}

#########################
## Nevill ID 10

############################
##
## Run the model


## CHECK THE CORRECT INPUT FILES ARE ADDED
# sim_file	sim_parms_S.txt               ##Check year of trial, check age group of prevalence cohort etc unless on command line
# mosq_file	mosq_species_parms.txt        ##Check Q0 and phi unless in site file or command lines
# corr_file	corr_file_0.txt               ##check correlated
# itn_file	itn_parms_2010.txt          ##check on command line and input folder
# irs_file	irs_parms_ellie_pyrethroid.txt  ##use default pyrethroid assuming no resistance
# mda_file	mda_parms.txt                 ##prob wont use - interventions.txt drug 0 
# vacc_file	vacc_parms.txt                ##wont use - interventions.txt epi_[] 0
# output_file	output_vars_S.txt           ##specify for the trial (Incidence 0 to 6 year olds; prevalence under 7s and under 6)
# model_file	model_parms.txt             ##will not change this

## Function to run the model across parameter draws 
## Site data, Kilifi Kenya
## Match site data to that estimated from the study for mosquito outdoor and human biting.


nevi_trial_f<-function(DIDE_CODE, ## Southern Province Sierre Leone 102761
                       
                       itn_repel_fun,itn_repel_gam,itn_repel_ara,
                       itn_kill_fun,itn_kill_gam,itn_kill_ara,
                       itn_half_life,
                       
                       arm, ## 1 or 2
                       run, ## either 1 for control or 2 for ITNs
                       
                       ITN_1, ## this is 'all people have nets in the intervention...'
                       itn_leave_dur,
                       
                       sites){
  
  arm <- arm
  run <- run
  
  pop_size<- 20000 #Sim_pop_size(site_file[site_file[,1]=='prev',2])
  
  Int_set_up<-paste('num_people', pop_size, ##This is arbitrary and simply speeds up or smooths out runs (start with small popn size to make sure it runs fine then up to about 30,000 - 80,000)
                    'itn_irs_corr',	1, ##same people receive nets and spray unless trial specifies otherwise
                    
                    'itn_repel_fun', itn_repel_fun, 'itn_repel_gamb_ss', itn_repel_gam, 'itn_repel_arab', itn_repel_ara,
                    'itn_kill_fun', itn_kill_fun, 'itn_kill_gamb_ss', itn_kill_gam, 'itn_kill_arab', itn_kill_ara,
                    'itn_half_life',itn_half_life,
                    # 'output_type 0 recalculate 0 itn_usage 1 add itn 0 itn_start 1.5 continue_itn 0 itn_flexible 1',
                    'output_type 0 recalculate 0 add itn 1 itn_start 0 itn_coverage', ITN_1,
                    'itn_max_rounds 8 itn_leave_dur', itn_leave_dur, ## year 0 is jan 1993
                    'itn_2 1 itn_start_2 0.5 itn_coverage_2', ITN_1, 'itn_leave_dur_2', itn_leave_dur,
                    'itn_3 1 itn_start_3 1 itn_coverage_3', ITN_1, 'itn_leave_dur_3', itn_leave_dur,
                    'itn_4 1 itn_start_4 1.5 itn_coverage_4', ITN_1, 'itn_leave_dur_4', itn_leave_dur,
                    'irs 0',
                    
                    'irs_ellie 1')
  
  Options<-paste(Int_set_up)
  
  sites <- sites
  draw<-0
  
  ## Run the simulation
  Model_launcher(OutputName = paste("VALIDATION_trial_arm",arm,"net_1",draw,run, sep = "_"),
                 # OutputRoot = paste0("P:/Ellies_output_folder/validations/nevill/arm",arm),
                 OutputRoot = paste0("P:/validating_modelling/output_files/FINAL_netparm2/nevill/arm", arm),
                 Options=Options,
                 Exe = "P:/validating_modelling/bin/irs_restructure.exe",
                 Root="P:/validating_modelling/input_files",
                 Demog = "demog_GADM18_101474.txt",  ## Match up for Kilifi Kenya
                 Pop = "pop_GADM18_101474.txt",        ## Match up for Kilifi Kenya
                 Site=sites,
                 Parameter_draw = draw,
                 Return_output = FALSE)
}


##################################
##
## Philips-Howard 11


## CHECK THE CORRECT INPUT FILES ARE ADDED
# sim_file	sim_parms_12.txt               ##Check year of trial, check age group of prevalence cohort etc unless on command line
# mosq_file	mosq_species_parms.txt        ##Check Q0 and phi unless in site file or command lines
# corr_file	corr_file_0.txt               ##check correlated
# itn_file	itn_parms_2010.txt          ##check on command line and input folder
# irs_file	irs_parms_ellie_pyrethroid.txt  ##use default pyrethroid assuming no resistance
# mda_file	mda_parms.txt                 ##prob wont use - interventions.txt drug 0 
# vacc_file	vacc_parms.txt                ##wont use - interventions.txt epi_[] 0
# output_file	output_vars_12.txt           ##specify for the trial (Incidence 0 to 6 year olds; prevalence under 7s and under 6)
# model_file	model_parms.txt             ##will not change this


PHILHO_trial_f1<-function(DIDE_CODE,
                          run,
                          
                          itn_repel_fun,itn_repel_gam,itn_repel_ara,
                          itn_kill_fun,itn_kill_gam,itn_kill_ara,
                          itn_half_life,
                          
                          # itn_repel_fun_1,itn_repel_gam_1,itn_repel_ara_1,
                          # itn_kill_fun_1,itn_kill_gam_1,itn_kill_ara_1,
                          # itn_half_life_1,
                          # 
                          ITN_1,ITN_2,ITN_3,ITN_4,
                          
                          sitnum,arm,dide_code_match,sites){
  sitnum <- sitnum
  arm <- arm
  run <- run
  
  # Load the site_file file
  # site_file<-western kenya
  
  pop_size<- 20000 #Sim_pop_size(site_file[site_file[,1]=='prev',2])
  # total_M = total_MM
  # 
  Int_set_up<-paste('num_people', pop_size, ##This is arbitrary and simply speeds up or smooths out runs (start with small popn size to make sure it runs fine then up to about 30,000 - 80,000)
                    # 'total_M', total_M,
                    'itn_irs_corr',	1, ##same people receive nets and spray unless trial specifies otherwise
                    'output_type 0 recalculate 0 itn_usage 1 add itn 0 itn_start 0.05 continue_itn 0 itn_flexible 1',## for calibration we need this to be output_type 0
                    
                    'num_runs 1 itn_cov_0_0', ITN_1, 'itn_cov_1_0', ITN_2, ## year 0 is JAN 1997 and aiming for under 5 year old prev of c75%
                    'itn_cov_2_0', ITN_3, 'itn_cov_3_0', ITN_4, ## year 2 is March 1999 and control arms now get ITNs
                    'itn_cov_4_0 0 itn_cov_5_0 0', ## End of trial
                    'itn_cov_6_0 0 itn_cov_7_0 0 itn_cov_8_0 0 itn_cov_9_0 0',
                    'itn_cov_10_0 0 itn_cov_11_0 0 itn_cov_12_0 0 itn_cov_13_0 0',
                    'itn_cov_14_0 0 itn_cov_15_0 0 itn_cov_16_0 0 itn_cov_17_0 0',
                    'itn_cov_18_0 0 itn_cov_19_0 0',
                    
                    ## Assuming resistance is unchanging over time *and was present pre the trials at whatever rate determined
                    'itn_repel_fun', itn_repel_fun, 'itn_repel_gamb_ss', itn_repel_gam, 'itn_repel_arab', itn_repel_ara,
                    'itn_kill_fun', itn_kill_fun, 'itn_kill_gamb_ss', itn_kill_gam, 'itn_kill_arab', itn_kill_ara,
                    'itn_half_life',itn_half_life) ## dipped nets are continually dipped throughout so should retain use well
  
  
  Options<-paste(Int_set_up)
  
  draw<-0
  dide_code_match <- dide_code_match
  sites <- sites
  
  ## Run the simulation
  Model_launcher(OutputName = paste("VALIDATION_trial_arm",arm,sitnum,"net_1",draw,run, sep = "_"),
                 OutputRoot = paste0("P:/validating_modelling/output_files/FINAL_netparm2/philipshoward/arm",arm),
                 Options=Options,
                 Exe = "P:/validating_modelling/bin/irs_restructure.exe", ## The model executable
                 Root="P:/validating_modelling/input_files",
                 Demog = paste0("demog_GADM18_",dide_code_match,".txt"),  ## Match up for Western Kenya, Asembo and Gem
                 Pop = paste0("pop_GADM18_",dide_code_match,".txt"),      ## Match up for Western Kenya, Asembo and Gem
                 Site=sites,
                 Parameter_draw = draw,
                 Return_output = FALSE)
}



##################################
##
## Protopopoff ID 12

#############################
##
## Run the model


## CHECK THE CORRECT INPUT FILES ARE ADDED
# sim_file	sim_parms_L.txt                 ##Check year of trial, check age group of prevalence cohort etc unless on command line
# mosq_file	mosq_species_parms.txt        ##Check Q0 and phi unless in site file or command lines
# corr_file	corr_file_0.txt               ##check correlated
# itn_file	itn_parms_2010.txt          ##check on command line and input folder
# irs_file	irs_parms_ellie_actellic.txt  ##use study default and check no IRS on if not implemented
# mda_file	mda_parms.txt                 ##prob wont use - interventions.txt drug 0 
# vacc_file	vacc_parms.txt                ##wont use - interventions.txt epi_[] 0
# output_file	output_vars_L.txt           ##specify for the trial (prev 0 to 14 year olds)
# model_file	model_parms.txt             ##will not change this


pbo_trial_f2<-function(DIDE_CODE,
                       ITN_1, IRS_1, itn_leave_dur,
                       itn_repel_fun,itn_repel_gam,itn_repel_ara,
                       itn_kill_fun,itn_kill_gam,itn_kill_ara,
                       itn_half_life,
                       irs_decay_mort1,irs_decay_mort2,
                       irs_decay_succ1,irs_decay_succ2,
                       irs_decay_det1,irs_decay_det2,
                       total_M,
                       arm,run,net_type,
                       sites){
  
  
  # Load the site_file file
  # site_file<-read.table(paste0('P:/Ellies_cool_model_folder2/model_files/sites/Africa_Sites_Ellie_Copy/Africa_sites_0/Tanz_Kagera_677_', site, '.txt'))
  ## THIS IS Tanz_Kagera_677_dbHGH_1.txt for the high baseline coverage
  ## THIS IS Tanz_Kagera_677_dbMID_1.txt for the median baseline coverage
  ## THIS IS Tanz_Kagera_677_dbLOW_1.txt for the low baseline coverage
  
  pop_size<- 30000 #Sim_pop_size(site_file[site_file[,1]=='prev',2])
  
  Int_set_up<-paste('num_people', pop_size, ##This is arbitrary and simply speeds up or smooths out runs (start with small popn size to make sure it runs fine then up to about 30,000 - 80,000)
                    'total_M', total_M,
                    'itn_irs_corr',	1, ##same people receive nets and spray unless trial specifies otherwise
                    'output_type 0 recalculate 0 itn_usage 1 add itn 0 itn_start 0 continue_itn 0 itn_flexible 1 itn_leave_dur',	itn_leave_dur, ## for calibration we need this to be output_type 0
                    
                    'num_runs 1 itn_cov_0_0', ITN_1, 'itn_cov_1_0 0',
                    'itn_cov_2_0 0 itn_cov_3_0 0 itn_cov_4_0 0 itn_cov_5_0 0 itn_cov_6_0 0 itn_cov_7_0 0 itn_cov_8_0 0 itn_cov_9_0 0 itn_cov_10_0 0',
                    'itn_cov_11_0 0 itn_cov_12_0 0 itn_cov_13_0 0 itn_cov_14_0 0 itn_cov_15_0 0 itn_cov_16_0 0 itn_cov_17_0 0 itn_cov_18_0 0 itn_cov_19_0 0',
                    
                    'change_itn 1 change_itn_time 0',
                    'itn_repel_fun_1', itn_repel_fun, 'itn_repel_gamb_ss_1', itn_repel_gam, 'itn_repel_arab_1', itn_repel_ara,
                    'itn_kill_fun_1', itn_kill_fun, 'itn_kill_gamb_ss_1', itn_kill_gam, 'itn_kill_arab_1', itn_kill_ara,
                    'itn_half_life_1',itn_half_life,
                    
                    'irs 1 irs_coverage', IRS_1, 'irs_start 0 irs_max_rounds 1 irs_offset_absolute 1 irs_offset 0.083 irs_2 1 irs_coverage_2 0 irs_start_2 0 irs_max_rounds_2 20', 
                    'change_irs 1 change_irs_time 0',
                    'irs_decay_mort1', irs_decay_mort1, 'irs_decay_mort2', irs_decay_mort2,
                    'irs_decay_succ1', irs_decay_succ1,'irs_decay_succ2', irs_decay_succ2,
                    'irs_decay_det1', irs_decay_det1,'irs_decay_det2', irs_decay_det2,
                    
                    'irs_ellie 1')
  
  Options<-paste(Int_set_up)
  
  sites <- sites
  draw<-0
  
  ## Run the simulation
  Model_launcher(OutputName = paste("VALIDATION_trial_arm",arm,"net",net_type,draw,run, sep = "_"),
                 OutputRoot = paste0("P:/validating_modelling/output_files/FINAL_netparm2/protopopoff/arm", arm),
                 Options=Options,
                 Exe = "P:/validating_modelling/bin/irs_restructure.exe",
                 Root="P:/validating_modelling/input_files",
                 Demog = "demog_GADM18_103222.txt",  ## Match up for Kagera
                 Pop = "pop_GADM18_103222.txt",        ## Match up for Kagera
                 Site=sites,
                 Parameter_draw = draw,
                 Return_output = FALSE)
}

##########################
##
## Staedke ID 13

## Modelling Staedke

############################
##
## Run the model


## CHECK THE CORRECT INPUT FILES ARE ADDED
# sim_file	sim_parms_Q.txt               ##Check year of trial, check age group of prevalence cohort etc unless on command line
# mosq_file	mosq_species_parms.txt        ##Check Q0 and phi unless in site file or command lines
# corr_file	corr_file_0.txt               ##check correlated
# itn_file	itn_parms_2010.txt          ##check on command line and input folder
# irs_file	irs_parms_ellie_pyrethroid.txt  ##use default pyrethroid assuming no resistance
# mda_file	mda_parms.txt                 ##prob wont use - interventions.txt drug 0 
# vacc_file	vacc_parms.txt                ##wont use - interventions.txt epi_[] 0
# output_file	output_vars_Q.txt           ##specify for the trial (Incidence 0 to 6 year olds; prevalence under 7s and under 6)
# model_file	model_parms.txt             ##will not change this

## Function to run the model across parameter draws 
## Site data, Kibale Uganda
## Match site data to that estimated from the study for mosquito outdoor and human biting.


STAEDKE_trial_f<-function(DIDE_CODE, ## Uganda 103270
                          
                          arm, ## 1 or 2
                          run, ## either 1 for control or 2 for ITNs
                          net_type,
                          
                          ITN_1, ## figure 3 in Staedke et al 2020; historic ITN use was 39.1% (Gonahasa et al 2018)
                          itn_leave_dur, ## can work this out from figures in Staedke et al 2020
                          itn_repel_fun,itn_repel_gam,itn_repel_ara,
                          itn_kill_fun,itn_kill_gam,itn_kill_ara,
                          itn_half_life,
                          total_M,
                          
                          sites){
  
  arm <- arm
  run <- run
  net_type <- net_type
  
  pop_size<- 20000 #Sim_pop_size(site_file[site_file[,1]=='prev',2])
  
  Int_set_up<-paste('num_people', pop_size, ##This is arbitrary and simply speeds up or smooths out runs (start with small popn size to make sure it runs fine then up to about 30,000 - 80,000)
                    'total_M',total_M,'itn_irs_corr',	1, ##same people receive nets and spray unless trial specifies otherwise
                    ## JULY - SEPTEMBER 2017 NET DISTRIBUTION
                    'output_type 0 recalculate 0 itn_usage 0 add itn 1 itn_start 0.67 continue_itn 0 itn_flexible 1 itn_leave_dur',	itn_leave_dur, ## for calibration we need this to be output_type 0
                    
                    'num_runs 1 itn_cov_0_0', ITN_1, 'itn_cov_1_0 0',
                    'itn_cov_2_0 0 itn_cov_3_0 0 itn_cov_4_0 0 itn_cov_5_0 0 itn_cov_6_0 0 itn_cov_7_0 0 itn_cov_8_0 0 itn_cov_9_0 0 itn_cov_10_0 0',
                    'itn_cov_11_0 0 itn_cov_12_0 0 itn_cov_13_0 0 itn_cov_14_0 0 itn_cov_15_0 0 itn_cov_16_0 0 itn_cov_17_0 0 itn_cov_18_0 0 itn_cov_19_0 0',
                    
                    'change_itn 1 change_itn_time 0.67',
                    'itn_repel_fun_1', itn_repel_fun, 'itn_repel_gamb_ss_1', itn_repel_gam, 'itn_repel_arab_1', itn_repel_ara,
                    'itn_kill_fun_1', itn_kill_fun, 'itn_kill_gamb_ss_1', itn_kill_gam, 'itn_kill_arab_1', itn_kill_ara,
                    'itn_half_life_1',itn_half_life,
                    
                    'irs 0', ## no sites had historic IRS (Staedke et al 2020)
                    
                    'irs_ellie 1')
  
  Options<-paste(Int_set_up)
  
  sites <- sites
  draw<-0
  
  ## Run the simulation
  Model_launcher(OutputName = paste("VALIDATION_trial_arm",arm,"net",net_type,draw,run, sep = "_"),
                 OutputRoot = paste0("P:/validating_modelling/output_files/FINAL_netparm2/staedkev2/arm",arm),
                 Options=Options,
                 Exe = "P:/validating_modelling/bin/irs_restructure.exe", ## The model executable
                 Root="P:/validating_modelling/input_files",
                 Demog = "demog_GADM18_103270.txt",  ## Match up for Kibale Uganda
                 Pop = "pop_GADM18_103270.txt",        ## Match up for Kibale Uganda
                 Site=sites,
                 Parameter_draw = draw,
                 Return_output = FALSE)
}

##########################
##
## West ID 14


#############################
##
## Run the model


## CHECK THE CORRECT INPUT FILES ARE ADDED
# sim_file	  sim_parms_O.txt                   ##Check year of trial, check age group of prevalence cohort etc unless on command line
# mosq_file 	mosq_species_parms_O.txt        ##Check Q0 and phi unless in site file or command lines
# corr_file	  corr_file_0.txt                 ##check correlated
# itn_file  	itn_parms_2010_O.txt            ##check on command line and input folder
# irs_file	  irs_parms_bendiocarb.txt  ##use study default and check no IRS on if not implemented
# mda_file	  mda_parms.txt                   ##prob wont use - interventions.txt drug 0 
# vacc_file	  vacc_parms.txt                  ##wont use - interventions.txt epi_[] 0
# output_file	output_vars_O.txt               ##specify for the trial (prev 0 to 14 year olds)
# model_file	model_parms.txt                 ##will not change this


west_trial_f1<-function(arm,run,
                        DIDE_CODE,
                        ITN_1, IRS_1, itn_leave_dur,
                        itn_repel_fun,itn_repel_gam,itn_repel_ara,
                        itn_kill_fun,itn_kill_gam,itn_kill_ara,
                        itn_half_life,
                        irs_decay_mort1,irs_decay_mort2,
                        irs_decay_succ1,irs_decay_succ2,
                        irs_decay_det1,irs_decay_det2,
                        IRS_2,
                        net_type,
                        total_M,
                        sites){
  
  arm <- arm
  run <- run
  net_type <- 0
  total_M <- total_M
  # Load the site_file file
  # site_file<-read.table(paste0('P:/Ellies_cool_model_folder2/model_files/sites/Africa_Sites_Ellie_Copy/Africa_sites_0/Tanz_Kagera_677_', site, '.txt'))
  ## THIS IS Tanz_Kagera_677_dbHGH_1.txt for the high baseline coverage
  ## THIS IS Tanz_Kagera_677_dbMID_1.txt for the median baseline coverage
  ## THIS IS Tanz_Kagera_677_dbLOW_1.txt for the low baseline coverage
  
  pop_size<- 30000 #Sim_pop_size(site_file[site_file[,1]=='prev',2])
  
  Int_set_up<-paste('num_people', pop_size, 'itn_irs_corr',	1, 'total_M',	total_M, 
                    'add output_type 0 itn 1 itn_start 0 continue_itn 0 itn_flexible 1 itn_leave_dur',	itn_leave_dur, 
                    'num_runs 1 itn_cov_0_0', ITN_1, 'itn_cov_1_0 0',
                    'itn_cov_2_0 0 itn_cov_3_0 0 itn_cov_4_0 0 itn_cov_5_0 0 itn_cov_6_0 0 itn_cov_7_0 0 itn_cov_8_0 0 itn_cov_9_0 0 itn_cov_10_0 0',
                    'itn_cov_11_0 0 itn_cov_12_0 0 itn_cov_13_0 0 itn_cov_14_0 0 itn_cov_15_0 0 itn_cov_16_0 0 itn_cov_17_0 0 itn_cov_18_0 0 itn_cov_19_0 0',
                    'change_itn 1 change_itn_time 0',
                    'itn_repel_fun_1', itn_repel_fun, 'itn_repel_gamb_ss_1', itn_repel_gam, 'itn_repel_arab_1', itn_repel_ara,
                    'itn_kill_fun_1', itn_kill_fun, 'itn_kill_gamb_ss_1', itn_kill_gam, 'itn_kill_arab_1', itn_kill_ara,
                    'itn_half_life_1',itn_half_life,
                    
                    'irs 1 irs_coverage', IRS_1, 'irs_start 0 irs_max_rounds 2 irs_frequency 0.33 irs_offset_absolute 1 irs_offset 0.92', 
                    'change_irs 1 change_irs_time 0',
                    'irs_decay_mort1', irs_decay_mort1, 'irs_decay_mort2', irs_decay_mort2,
                    'irs_decay_succ1', irs_decay_succ1,'irs_decay_succ2', irs_decay_succ2,
                    'irs_decay_det1', irs_decay_det1,'irs_decay_det2', irs_decay_det2,
                    'irs_ellie 1')
  
  Options<-paste(Int_set_up)
  
  sites <- sites
  draw<-0
  
  ## Run the simulation
  Model_launcher(OutputName = paste("VALIDATION_trial_arm",arm,"net",net_type,draw,run, sep = "_"), ## filepath
                 OutputRoot = paste0("P:/validating_modelling/output_files/FINAL_netparm2/west/arm",arm),
                 Options=Options,
                 Exe = "P:/validating_modelling/bin/irs_restructure.exe", ## The model executable
                 Root="P:/validating_modelling/input_files",
                 Demog = "demog_GADM18_103222.txt",    ## Match up for Kagera
                 Pop = "pop_GADM18_103222.txt",        ## Match up for Kagera
                 Site=sites,
                 Parameter_draw = draw,
                 Return_output = FALSE)
}

