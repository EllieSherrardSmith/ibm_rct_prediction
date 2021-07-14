## vALIDATIONS SET UP#
# source("https://dide-tools.github.io/didehpc/install")
share <- didehpc::path_mapping("malaria", "P:", "//projects.dide.ic.ac.uk/malaria/Ellie/Rprojects/Malaria", "P:")
# share <- didehpc::path_mapping("malaria", "K:", "//fi--didenas1/malaria/Ellie", "K:")


#file.exists() vcheck the file is correctly named
#dir() lists the files in a directory
#options(didehpc.cluster = "fi--didemrchnb")#   
#config <- didehpc::didehpc_config(shares = share)#small

#config <- didehpc::didehpc_config(shares = share, use_rrq = TRUE, cluster = "fi--didemrchnb")#  fi--dideclusthn #This is if you want to block out some cores and circle through jobs on those
config <- didehpc::didehpc_config(shares = share, cluster = "fi--didemrchnb")#  fi--dideclusthn # fi--didemrchnb
context::context_log_start()
root <- "contextsA"
src = provisionr::package_sources(repos = "file:///h:/drat")

ctx <- context::context_save(root, 
                             packages = c("MalariaLaunchR"), ##Any packages I need
                             sources = "functions_validations_netparms_3.r", ##Any functions I need
                             package_sources = src) ##Any location of resources I need (not on CRAN stuff)

obj = didehpc::queue_didehpc(ctx,config = config)
##Mdog29d

obj$cluster_load() ## check for space
didehpc::didehpc_config()

##############################
## Bradley ID 1
## Prepare data appropriately

library("dplyr")
## Data parameter input files
arm1 = read.csv("Input_files/data/Model_params/input_params_1_net_parameters3_all_hutdata.csv",header=TRUE)
arm1$total_MM = rnorm(mean=18.75829*2/5,sd = 0.8,n = 2000)
head(arm1)
## Site data, 
site_data <- read.csv("Site_files/site_file_1.csv", stringsAsFactors = FALSE)

# Create list columns
sites <- as_tibble(select(site_data, DIDE_CODE)) %>%
  mutate(site = apply(select(site_data, - DIDE_CODE), 1,
                      function(x){ as.list(x)}))

arm1_a = arm1[c(1:10),]#[1001:1003,]

v1 <- left_join(arm1, sites, by = "DIDE_CODE")


# purrr::pmap(v1, functions_model_runs_bradley)
Bdd = obj$enqueue_bulk(v1, functions_model_runs_bradley)

###############################
##
## Chaccour ID 2
##

## Data parameter input files
arm1 = read.csv("Input_files/data/Model_params/input_params_2_net_parameters3_all_hutdata.csv",header=TRUE)
arm1$total_MM = rnorm(mean=110,sd = 30,n = 2000)
head(arm1)
## Site data, N
site_data <- read.csv("Site_files/site_file_2.csv", stringsAsFactors = FALSE)

# Create list columns
sites <- as_tibble(select(site_data, DIDE_CODE)) %>%
  mutate(site = apply(select(site_data, - DIDE_CODE), 1,
                      function(x){ as.list(x)}))

arm1_a = arm1[c(1:3),]#[1001:1003,]

v1 <- left_join(arm1, sites, by = "DIDE_CODE")


# purrr::pmap(v1, chac_trial_f)
Bdd = obj$enqueue_bulk(v1, chac_trial_f)


###############################
##
## Corbel ID 3
##

## Data parameter input files
arm1 = read.csv("Input_files/data/Model_params/input_params_3_net_parameters3_all_hutdata.csv",header=TRUE)
arm1$total_MM = c(rnorm(mean=4.8,sd = 0.3,n = 1000),
                  rnorm(mean=6,sd = 0.3,n = 1000),
                  rnorm(mean=8.5,sd = 0.3,n = 1000),
                  rnorm(mean=7.5,sd = 0.3,n = 1000))
head(arm1)
## Site data, N
site_data <- read.csv("Site_files/site_file_3.csv", stringsAsFactors = FALSE)

# Create list columns
sites <- as_tibble(select(site_data, DIDE_CODE)) %>%
  mutate(site = apply(select(site_data, - DIDE_CODE), 1,
                      function(x){ as.list(x)}))

arm1_a = arm1[c(1:3000),]#[1001:1003,]

v1 <- left_join(arm1_a, sites, by = "DIDE_CODE")


# purrr::pmap(v1[1,], corb_trial_f1)
Bdd = obj$enqueue_bulk(v1, corb_trial_f1)



###############################
##
## Curtis ID 4
##

## Data parameter input files
arm1 = read.csv("Input_files/data/Model_params/input_params_4_net_parameters3_all_hutdata.csv",header=TRUE)
arm1$total_MM = rnorm(mean=500,sd = 30,n = 3000)
head(arm1)
## Site data,
site_data <- read.csv("Site_files/site_file_4.csv", stringsAsFactors = FALSE)

# Create list columns
sites <- as_tibble(select(site_data, DIDE_CODE)) %>%
  mutate(site = apply(select(site_data, - DIDE_CODE), 1,
                      function(x){ as.list(x)}))

arm1_a = arm1[c(1001:1003,2001:2003),]#[1001:1003,]

v1 <- left_join(arm1, sites, by = "DIDE_CODE")


# purrr::pmap(v1, curt_trial_f1)
Bdd = obj$enqueue_bulk(v1, curt_trial_f1)


################################
##
## D'Alessandro ID 5
##

## Data parameter input files
arm1 = read.csv("Input_files/data/Model_params/input_params_5_net_parameters3_all_hutdata.csv",header=TRUE)
arm1$total_MM = rep(c(rnorm(mean=17.06387,sd = 0.1,n = 1000),
                      rnorm(mean=10.90881,sd = 0.1,n = 1000),
                      rnorm(mean=17.06387,sd = 0.1,n = 1000),
                      rnorm(mean=69.9967,sd = 0.1,n = 1000),
                      rnorm(mean=229.7857,sd = 0.1,n = 1000)),2)
head(arm1)
## Site data,
site_data <- read.csv("Site_files/site_file_5.csv", stringsAsFactors = FALSE)

# Create list columns
sites <- as_tibble(select(site_data, DIDE_CODE)) %>%
  mutate(site = apply(select(site_data, - DIDE_CODE), 1,
                      function(x){ as.list(x)}))

arm1_a = arm1[c(1:2000),]#[1001:1003,]
arm1_b = arm1[c(2001:4000),]#[1001:1003,]
arm1_c = arm1[c(4001:8000),]#[1001:1003,]
arm1_d = arm1[c(8001:10000),]#[1001:1003,]

v1 <- left_join(arm1_a, sites, by = "DIDE_CODE")
v2 <- left_join(arm1_b, sites, by = "DIDE_CODE")
v3 <- left_join(arm1_c, sites, by = "DIDE_CODE")
v4 <- left_join(arm1_d, sites, by = "DIDE_CODE")


# purrr::pmap(v1, curt_trial_f1)
Bdd = obj$enqueue_bulk(v1, DALESS_trial_f1)
BEd = obj$enqueue_bulk(v2, DALESS_trial_f1)
BFd = obj$enqueue_bulk(v3, DALESS_trial_f1)
BGd = obj$enqueue_bulk(v4, DALESS_trial_f1)


################################
##
## Henry ID 6
##

## Data parameter input files
arm1 = read.csv("Input_files/data/Model_params/input_params_6_net_parameters3_all_hutdata.csv",header=TRUE)
arm1$total_MM = c(rnorm(mean=234,sd=1.2,n=1000),
                  rnorm(mean=284,sd=1.2,n=1000))
## Site data, 
site_data <- read.csv("Site_files/site_file_6.csv", stringsAsFactors = FALSE)

# Create list columns
sites <- as_tibble(select(site_data, DIDE_CODE)) %>%
  mutate(site = apply(select(site_data, - DIDE_CODE), 1,
                      function(x){ as.list(x)}))

arm1_a = arm1[c(1:5,1001:1005),]#[1001:1003,]

v1 <- left_join(arm1, sites, by = "DIDE_CODE")


# purrr::pmap(v1, HENRY_trial_f1)
Bdd = obj$enqueue_bulk(v1, HENRY_trial_f1)


################################
##
## Kafy ID 7
##

## Data parameter input files
arm1 = read.csv("Input_files/data/Model_params/input_params_7_net_parameters3_all_hutdata.csv",header=TRUE)
arm1$total_MM = c(rnorm(mean=5.6,sd=0.2,n=1000),
                  rnorm(mean=7,sd=0.2,n=1000))
## Site data, 
site_data <- read.csv("Site_files/site_file_7.csv", stringsAsFactors = FALSE)

# Create list columns
sites <- as_tibble(select(site_data, DIDE_CODE)) %>%
  mutate(site = apply(select(site_data, - DIDE_CODE), 1,
                      function(x){ as.list(x)}))

arm1_a = arm1[c(1:5,1001:1005),]#[1001:1003,]

v1 <- left_join(arm1, sites, by = "DIDE_CODE")


# purrr::pmap(v1, functions_model_runs_kafy)
Bdd = obj$enqueue_bulk(v1, functions_model_runs_kafy)



##############################
## Marbiah ID 9
## Prepare data appropriately

library("dplyr")
## Data parameter input files
arm1 = read.csv("Input_files/data/Model_params/input_params_9_net_parameters3_all_hutdata.csv",header=TRUE)
## Site data,  
site_data <- read.csv("Site_files/site_file_9.csv", stringsAsFactors = FALSE)

# Create list columns
sites <- as_tibble(select(site_data, DIDE_CODE)) %>%
  mutate(site = apply(select(site_data, - DIDE_CODE), 1,
                      function(x){ as.list(x)}))

arm1_a = arm1[c(1:10),]#[1001:1003,]

v1 <- left_join(arm1, sites, by = "DIDE_CODE")


# purrr::pmap(v1, marb_trial_f)
Bdd = obj$enqueue_bulk(v1, marb_trial_f)


##############################
## Nevill ID 10
## Prepare data appropriately

library("dplyr")
## Data parameter input files
arm1 = read.csv("Input_files/data/Model_params/input_params_10_net_parameters3_all_hutdata.csv",header=TRUE)
# arm1$total_M = rnorm(mean=81,sd=4,n=2000)
## Site data,  
site_data <- read.csv("Site_files/site_file_10.csv", stringsAsFactors = FALSE)

# Create list columns
sites <- as_tibble(select(site_data, DIDE_CODE)) %>%
  mutate(site = apply(select(site_data, - DIDE_CODE), 1,
                      function(x){ as.list(x)}))

arm1_a = arm1[c(1:10),]#[1001:1003,]

v1 <- left_join(arm1, sites, by = "DIDE_CODE")


# purrr::pmap(v1, nevi_trial_f)
Bdd = obj$enqueue_bulk(v1, nevi_trial_f)


##############################
## PhilipsHoward ID 11
## Prepare data appropriately

library("dplyr")
## Data parameter input files
arm1 = read.csv("Input_files/data/Model_params/input_params_11_net_parameters3_all_hutdata.csv",header=TRUE)
# arm1$total_M = rnorm(mean=81,sd=4,n=2000)
## Site data,  
site_data <- read.csv("Site_files/site_file_11.csv", stringsAsFactors = FALSE)

# Create list columns
sites <- as_tibble(select(site_data, DIDE_CODE)) %>%
  mutate(site = apply(select(site_data, - DIDE_CODE), 1,
                      function(x){ as.list(x)}))

arm1_a = arm1[c(1:10),]#[1001:1003,]

v1 <- left_join(arm1, sites, by = "DIDE_CODE")


# purrr::pmap(v1, PHILHO_trial_f1)
Bdd = obj$enqueue_bulk(v1, PHILHO_trial_f1)


##############################
## Protopopoff ID 12
## Prepare data appropriately

library("dplyr")
## Data parameter input files
arm1 = read.csv("Q:/RProjects/Model_validation_ITN_IRS/Input_files/data/Model_params/input_params_12_update_net_parameters3_all_hutdataB.csv",header=TRUE)

## Site data,  
site_data <- read.csv("Site_files/site_file_12.csv", stringsAsFactors = FALSE)

# Create list columns
sites <- as_tibble(select(site_data, DIDE_CODE)) %>%
  mutate(site = apply(select(site_data, - DIDE_CODE), 1,
                      function(x){ as.list(x)}))

arm1_a = arm1[c(1001:2000,3001:4000),]#[1001:1003,]

v1 <- left_join(arm1_a, sites, by = "DIDE_CODE")


# purrr::pmap(v1, pbo_trial_f2)
Bdd = obj$enqueue_bulk(v1, pbo_trial_f2)


##############################
## Staedke ID 13
## Prepare data appropriately

library("dplyr")
## Data parameter input files
arm1 = read.csv("Q:/RProjects/Model_validation_ITN_IRS/Input_files/data/Model_params/input_params_13_net_parameters3_all_hutdata_B.csv",header=TRUE)

## Site data,  
site_data <- read.csv("Site_files/site_file_13.csv", stringsAsFactors = FALSE)

# Create list columns
sites <- as_tibble(select(site_data, DIDE_CODE)) %>%
  mutate(site = apply(select(site_data, - DIDE_CODE), 1,
                      function(x){ as.list(x)}))

arm1_a = arm1[c(1001:2000),]#[1001:1003,]

v1 <- left_join(arm1_a, sites, by = "DIDE_CODE")


# purrr::pmap(v1, STAEDKE_trial_f)
Bdd = obj$enqueue_bulk(v1, STAEDKE_trial_f)


##############################
## West ID 14
## Prepare data appropriately

library("dplyr")
## Data parameter input files
arm1 = read.csv("Input_files/data/Model_params/input_params_14_net_parameters3_all_hutdata.csv",header=TRUE)
arm1$total_M = c(rnorm(n=1000,mean=12.5,sd=0.7),
                 rnorm(n=1000,mean=11,sd=0.7))
## Site data,  
site_data <- read.csv("Site_files/site_file_14.csv", stringsAsFactors = FALSE)

# Create list columns
sites <- as_tibble(select(site_data, DIDE_CODE)) %>%
  mutate(site = apply(select(site_data, - DIDE_CODE), 1,
                      function(x){ as.list(x)}))

arm1_a = arm1[c(1:2),]#[1001:1003,]

v1 <- left_join(arm1, sites, by = "DIDE_CODE")


# purrr::pmap(v1, west_trial_f1)
Bdd = obj$enqueue_bulk(v1, west_trial_f1)
