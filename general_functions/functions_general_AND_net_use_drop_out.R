write.csv(params_all,"temp.csv")

## Functions to estimate standard error and uncertainty around data points

se = function(data){
  x = sd(data) / sqrt(1000) 
  return(x)
}

uncertainty_fn = function(param_best){
  se_dat = se(rnorm(n = 1000, mean = param_best, sd = 0.5))
  N_dat = (param_best*(1-param_best)/se_dat^2)-1
  
  alpha = param_best * N_dat
  beta = N_dat - alpha
  
  reps_out = rbeta(n = 1000,alpha,beta)
  
  return(reps_out)
}


## Function and statistical model written in stan to estimate the drop out rate of people using nets

library(rstan)
library(adegenet)

## Using Staedke et al as an example

## Specify the times when data were observed
time_obs = c(6/12,12/12,18/12,
             6/12,12/12,18/12,
             6/12,12/12,18/12)

# Coverage 6-months: 71% ## From supplementary figures Staedke et al 2020
# Coverage 12-months: 63%
# Coverage 18-months: 49%
## Specify the data observed (including the uncertainties as this provides more data points for a better fit)
## where more specific data for different clusters are available this can be improved.
standard_net_usage = c(0.71,0.63,0.49,
                       0.62,0.58,0.44, ##these are the estimated uncertainty ranges
                       0.72,0.68,0.54) ##these are the estimated uncertainty ranges

## Transform data for the model fitting
y_standard_net_usage = log(standard_net_usage)

## Specify a time series to project outcomes across
time_m = seq(0,3,0.01)

## Specify the exponential decay model in Rstan
stanmodelcode <- "
data {
int<lower=0> N;
int<lower=0> N2; //the size of the new_X matrix
vector[N] y;
vector[N] x;
vector[N2] New_x;

}
parameters {
real beta0;
real beta1;
real sigma;
}
transformed parameters {
vector[N] m;
m = beta0 + beta1 * x;
} 
model {
// priors
beta0 ~ cauchy(0, 10); 
beta1 ~ cauchy(0, 10); 

// likelihood
y ~ normal(m, sigma);   
}
generated quantities {
vector[N2] y_pred;
y_pred = beta0 + beta1 * New_x; //the y values predicted by the model
}
"
stanDso <- stan_model(model_code = stanmodelcode) 

## Put data into a list
dat_standard <- list(N = length(time_obs), 
                     N2 = length(seq(0,3,0.01)),
                     y = y_standard_net_usage, 
                     x = time_obs,
                     New_x = seq(0,3,0.01))

## Run the statistical model
fit <- sampling(stanDso, 
                data = dat_standard, ## specify data
                iter = 4000,         ## speify iterations
                warmup=2000)         ## specify warm up (default is 4 chains here)


## plotting the posterior distribution for the parameters
post_beta<-As.mcmc.list(fit,pars="beta0")
plot(post_beta)

## gradient is fit to the data for alpha
## standard_net_usage ~ exp(-alpha*time_obs)
b0 <- extract(fit, 'beta0')
b0<- unlist(b0, use.names=FALSE)
b1 <- extract(fit, 'beta1')
b1<- unlist(b1, use.names=FALSE)

## Back translate output estimates
y_predicted = mean(b0) + mean(b1)*time_m; 
y_predicted_stn_exp = exp(mean(b0)) * exp(mean(b1)*time_m)

## and select some uncertainty for the Bayesian estimates
y_predicted_stn_exp_min = exp(quantile(b0,0.25)) * exp(quantile(b1,0.25)*time_m)
y_predicted_stn_exp_max = exp(quantile(b0,0.75)) * exp(quantile(b1,0.75)*time_m)


## Final output plotted to confirm the fit
par(mfrow = c(1,1))
plot(standard_net_usage[1:3] ~ ## plotting the mean estimates first
       time_obs[1:3],
     ylab="Households with at least one net per two occupants (%)", ## specify data used
     xlab="Time in years",yaxt="n",ylim=c(0,1),cex.lab=1.4,cex.axis=1.4,xlim=c(0,3))
axis(2,las=2, at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.4,cex.axis=1.4)

polygon(c(time_m,rev(time_m)),c(y_predicted_stn_exp_min,rev(y_predicted_stn_exp_max)),border=NA,col=transp("grey","0.5"))
lines(y_predicted_stn_exp ~ time_m,col="black",lty=2,lwd=1)

## add the range in the uncertainty from the observed data for the trial if available
for(i in 4:6){
  segments(x0=time_obs[i],x1=time_obs[i],
           y0=standard_net_usage[i],y1=standard_net_usage[i+3],lty=1)
  segments(x0=time_obs_offset[i],x1=time_obs_offset[i],
           y0=pbo_net_usage[i],y1=pbo_net_usage[i+3],lty=1,col="blue")
}

points(pbo_net_usage[1:3] ~ time_obs_offset[1:3],col="blue",pch=19)
points(standard_net_usage[1:3] ~ time_obs[1:3])

## Pull out the parameter required to specify drop out from using ITNs in the transmission model 
parms_usage = data.frame(itn_leave_dur_standardLLIN = b1[b1 >= quantile(b1,0.25) & b1 <= quantile(b1,0.75)]) ##
parms_usage$itn_leave_dur_standardLLIN_bt = -1/parms_usage$itn_leave_dur_standardLLIN

PARMS_USAGE = data.frame(parms_usage$itn_leave_dur_standardLLIN_bt)

# Confirm this parameter is specifying what we expect
D = median(PARMS_USAGE[,1])
D_LOW = quantile(PARMS_USAGE[,1],0.01)
D_UPP = quantile(PARMS_USAGE[,1],0.99)

cover_at_start = 0.85
cover_at_start_UPP = 0.88
cover_at_start_LOW = 0.8

## The below replicate the way this parameter is included in the transmission model
aa = cover_at_start * exp((-1/D)*time_m)
aL = cover_at_start_LOW * exp((-1/D_LOW)*time_m)
aU = cover_at_start_UPP * exp((-1/D_UPP)*time_m)

lines(aa ~ time_m,col="blue")
lines(aL ~ time_m,col="blue",lty=3)
lines(aU ~ time_m,col="blue",lty=3)
