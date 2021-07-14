## Function to adjust the seasonal pattern so that ITNs can be distributed in the appropriate month
## The transmission model default is for January distribution. 
## To switch this we can adjust the Fourrier function so that the first month of the 
## time series post implementation corresponds to the distribution month

# Obtain fourier coefficients

## This is a file with the original site parameters for the setting 
## following Garske et al 2013 and Walker et al 2016

data = read.csv("data\\seasonal_data\\season_original.csv", stringsAsFactors = FALSE)
data$start_month[4:5] = c("June","June") ## ITN distribution in June 2017 in Mopeia
data$itn_start[4:5] = c(0.5,0.5)         ## ITN distribution in June 2017 in Mopeia
seasonal_fourier_coefs = array(dim=c(nrow(data),7))
offset = data$itn_start

for(i in 1:nrow(data)){
  
  ##Re-ordering to match the site file we already have
  seasonal_fourier_coefs[i,1] = data$seasonal_a0[i]
  
  seasonal_fourier_coefs[i,2] = data$seasonal_a1[i] * cos(1*2*pi*offset[i]) - data$seasonal_b1[i] * sin(1*2*pi*offset[i])
  seasonal_fourier_coefs[i,3] = data$seasonal_a1[i] * sin(1*2*pi*offset[i]) + data$seasonal_b1[i] * cos(1*2*pi*offset[i])
  seasonal_fourier_coefs[i,4] = data$seasonal_a2[i] * cos(2*2*pi*offset[i]) - data$seasonal_b2[i] * sin(2*2*pi*offset[i])
  seasonal_fourier_coefs[i,5] = data$seasonal_a2[i] * sin(2*2*pi*offset[i]) + data$seasonal_b2[i] * cos(2*2*pi*offset[i])
  seasonal_fourier_coefs[i,6] = data$seasonal_a3[i] * cos(3*2*pi*offset[i]) - data$seasonal_b3[i] * sin(3*2*pi*offset[i])
  seasonal_fourier_coefs[i,7] = data$seasonal_a3[i] * sin(3*2*pi*offset[i]) + data$seasonal_b3[i] * cos(3*2*pi*offset[i])
  
}

colnames(seasonal_fourier_coefs) = names(data[,5:11])
head(seasonal_fourier_coefs)

## We can double check that the original data align with these new data

ssa0 = 	data$seasonal_a0[1]
ssa1 = data$seasonal_a1[1]
ssb1 = data$seasonal_b1[1]
ssa2 = data$seasonal_a2[1]
ssb2 =  data$seasonal_b2[1]
ssa3 = 	data$seasonal_a3[1]
ssb3 =	data$seasonal_b3[1]

#((ssa0+ssa1*cos(2*pi*TIME/365)+ssa2*cos(2*2*pi*TIME/365)+ssa3*cos(3*2*pi*TIME/365)+ssb1*sin(2*pi*TIME/365)+ssb2*sin(2*2*pi*TIME/365)+ ssb3*sin(3*2*pi*TIME/365) ) /theta_c,0.001)

###divide by theta_c to make mean 1
TIME = 1:365
data = (ssa0+ssa1*cos(2*pi*TIME/365)+ssa2*cos(2*2*pi*TIME/365)+
          ssa3*cos(3*2*pi*TIME/365)+ssb1*sin(2*pi*TIME/365)+
          ssb2*sin(2*2*pi*TIME/365)+ ssb3*sin(3*2*pi*TIME/365) )

data[which(data < 0)] = 0.001

theta_c =  mean(data)
mal1 = data/theta_c
par(mfrow=c(2,1))
plot(TIME,mal1)

ssa0 = 	as.numeric(seasonal_fourier_coefs[1,1])
ssa1 = as.numeric(seasonal_fourier_coefs[1,2])
ssb1 = as.numeric(seasonal_fourier_coefs[1,3])
ssa2 = as.numeric(seasonal_fourier_coefs[1,4])
ssb2 =  as.numeric(seasonal_fourier_coefs[1,5])
ssa3 = 	as.numeric(seasonal_fourier_coefs[1,6])
ssb3 =	as.numeric(seasonal_fourier_coefs[1,7])

#((ssa0+ssa1*cos(2*pi*TIME/365)+ssa2*cos(2*2*pi*TIME/365)+ssa3*cos(3*2*pi*TIME/365)+ssb1*sin(2*pi*TIME/365)+ssb2*sin(2*2*pi*TIME/365)+ ssb3*sin(3*2*pi*TIME/365) ) /theta_c,0.001)

###divide by theta_c to make mean 1
TIME = 1:365
data = (ssa0+ssa1*cos(2*pi*TIME/365)+ssa2*cos(2*2*pi*TIME/365)+
          ssa3*cos(3*2*pi*TIME/365)+ssb1*sin(2*pi*TIME/365)+
          ssb2*sin(2*2*pi*TIME/365)+ ssb3*sin(3*2*pi*TIME/365) )

data[which(data < 0)] = 0.001

theta_c =  mean(data)
mal2 = data/theta_c

lines(mal2)
offset_check=365*offset[1]
abline(v=offset_check,lty=2)

## As expected
## Now create site_parameters_offset.csv
write.csv(seasonal_fourier_coefs,"data\\site_files_offset_seasonality.csv")


