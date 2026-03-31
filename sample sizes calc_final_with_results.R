#sample size calc and associated results 
library (rethinking)

set.seed (100)
N <- 100

#doing computing for mean gradient values first 
mean <- 40 
sd <- 15
sdlog <- sqrt (log(1 + (sd^2 / mean^2)))
meanlog <- log(mean) - (sdlog^2 / 2)

sim_detect_effect_ulam_v2 <- function (N) {
  sim_data <- list (
    gender = rbinom (N, 1, 0.5), 
    age = rnorm (N, 70, 10), 
    mean_gradient = rlnorm (N, meanlog, sdlog)
  )
  sim_data$age_c <- sim_data$age - 70
  sim_data$mean_gradient_c <- sim_data$mean_gradient - 40 
  
  #simulate the "true" parameters 
  a <- 15
  b_gender <-  1.4
  b_age <- -0.2
  b_meanGradient <-  -0.15
  sigma <- 4
  
  #build mu 
  mu <- a + b_gender*sim_data$gender + b_age*sim_data$age_c + b_meanGradient*sim_data$mean_gradient_c
  sim_data$VO2 <- rnorm (N, mu, sigma)
  
  #fit the model with ulam MCMC
  m <- ulam (
    alist(
      VO2 ~ dnorm (mu, sigma), 
      mu <- a + b_gender*gender + b_age*age_c + b_meanGradient*mean_gradient_c,
      a ~ dnorm (15, 2), 
      b_gender ~ dnorm (1.4, 0.5), 
      b_age ~ dnorm (-0.2, 0.05), 
      b_meanGradient ~ dnorm (-0.15, 0.05), 
      sigma ~ dunif (1,5)
    ),
    data = sim_data,
    chains = 4,
    cores = 4, 
    iter = 1000, 
    log_lik = FALSE
  )
  
  post <- extract.samples (m)
  
  #clinically meaningful thresholds have been set  
  detected_age <- mean (post$b_age < -0.15) > 0.90
  detected_gender <- mean (post$b_gender > 1) > 0.90 
  detected_gradient <- mean (post$b_meanGradient < -0.1) > 0.90
  
  return (c(detected_age, detected_gender, detected_gradient)) #this has been changes to return power per parameter 
}
#testing a single run of the simulation function 
sim_detect_effect_ulam_v2(50)

sample_sizes <- c(30, 40, 60, 80, 100)
n_sims <- 100 #will take a while to do, because of the number of simulations. Last time took around 40 minutes 

#loop coding below will display power calculation for each parameter's effect size 
power_results <- matrix (NA, nrow = length(sample_sizes), ncol = 3)
colnames(power_results) <- c("age", "gender", "gradient")

for (i in seq_along(sample_sizes)) {
  detections <- replicate (n_sims, sim_detect_effect_ulam_v2(sample_sizes[i]))
  
  power_results[i, ] <- rowMeans(detections, na.rm = TRUE)
}


#table of results 
results_table <- data.frame (
  sample_size = sample_sizes,
  power_age = power_results[, "age"], 
  power_gender = power_results[, "gender"], 
  power_gradient = power_results[, "gradient"]
)

print (results_table)

#results table below, sample size 80 is identical to previous runs with this set of sample sizes (30, 40, 60, 80, 100). with the 0.85 power, but it is different to the simulation for the sample size set (30, 50, 70, 90, 100)

#sample_size power_age power_gender power_gradient
#1          30      0.53         0.13           0.53
#2          40      0.53         0.11           0.65
#3          60      0.62         0.15           0.70
#4          80      0.65         0.29           0.85
#5         100      0.67         0.24           0.81



#set 30, 50, 70, 80, 100
#now running sample size calculation for the set: 30, 50, 70, 80, 100
library (rethinking)

set.seed (100)
N <- 100

#doing computing for mean gradient values first 
mean <- 40 
sd <- 15
sdlog <- sqrt (log(1 + (sd^2 / mean^2)))
meanlog <- log(mean) - (sdlog^2 / 2)

sim_detect_effect_ulam_v2 <- function (N) {
  sim_data <- list (
    gender = rbinom (N, 1, 0.5), 
    age = rnorm (N, 70, 10), 
    mean_gradient = rlnorm (N, meanlog, sdlog)
  )
  sim_data$age_c <- sim_data$age - 70
  sim_data$mean_gradient_c <- sim_data$mean_gradient - 40 
  
  #simulate the "true" parameters 
  a <- 15
  b_gender <-  1.4
  b_age <- -0.2
  b_meanGradient <-  -0.15
  sigma <- 4
  
  #build mu 
  mu <- a + b_gender*sim_data$gender + b_age*sim_data$age_c + b_meanGradient*sim_data$mean_gradient_c
  sim_data$VO2 <- rnorm (N, mu, sigma)
  
  #fit the model with ulam MCMC
  m <- ulam (
    alist(
      VO2 ~ dnorm (mu, sigma), 
      mu <- a + b_gender*gender + b_age*age_c + b_meanGradient*mean_gradient_c,
      a ~ dnorm (15, 2), 
      b_gender ~ dnorm (1.4, 0.5), 
      b_age ~ dnorm (-0.2, 0.05), 
      b_meanGradient ~ dnorm (-0.15, 0.05), 
      sigma ~ dunif (1,5)
    ),
    data = sim_data,
    chains = 4,
    cores = 4, 
    iter = 1000, 
    log_lik = FALSE
  )
  
  post <- extract.samples (m)
  
  #clinically meaningful thresholds have been set  
  detected_age <- mean (post$b_age < -0.15) > 0.90
  detected_gender <- mean (post$b_gender > 1) > 0.90 
  detected_gradient <- mean (post$b_meanGradient < -0.1) > 0.90
  
  return (c(detected_age, detected_gender, detected_gradient)) #this has been changes to return power per parameter 
}
#testing a single run of the simulation function 
sim_detect_effect_ulam_v2(50)

sample_sizes <- c(30, 50, 70, 80, 100)
n_sims <- 100 #will take a while to do, because of the number of simulations. Last time took around 40 minutes 

#loop coding below will display power calculation for each parameter's effect size  
power_results <- matrix (NA, nrow = length(sample_sizes), ncol = 3)
colnames(power_results) <- c("age", "gender", "gradient")

for (i in seq_along(sample_sizes)) {
  detections <- replicate (n_sims, sim_detect_effect_ulam_v2(sample_sizes[i]))
  
  power_results[i, ] <- rowMeans(detections, na.rm = TRUE)
}


#table of results 
results_table <- data.frame (
  sample_size = sample_sizes,
  power_age = power_results[, "age"], 
  power_gender = power_results[, "gender"], 
  power_gradient = power_results[, "gradient"]
)

print (results_table)
#sample_size power_age power_gender power_gradient
#1          30      0.39         0.07           0.57
#2          50      0.58         0.16           0.69
#3          70      0.53         0.15           0.73
#4          80      0.58         0.22           0.76
#5         100      0.60         0.21           0.86

#could be simulation noise, but now 80 sample size has fallen below 0.80 power with simulation. On a separate run as a double check of this result, it was at 0.79 power for mean gradient detection


#set 30, 50, 70, 90, 100
library (rethinking)

set.seed (100)
N <- 100

#doing computing for mean gradient values first 
mean <- 40 
sd <- 15
sdlog <- sqrt (log(1 + (sd^2 / mean^2)))
meanlog <- log(mean) - (sdlog^2 / 2)

sim_detect_effect_ulam_v2 <- function (N) {
  sim_data <- list (
    gender = rbinom (N, 1, 0.5), 
    age = rnorm (N, 70, 10), 
    mean_gradient = rlnorm (N, meanlog, sdlog)
  )
  sim_data$age_c <- sim_data$age - 70
  sim_data$mean_gradient_c <- sim_data$mean_gradient - 40 
  
  #simulate the "true" parameters 
  a <- 15
  b_gender <-  1.4
  b_age <- -0.2
  b_meanGradient <-  -0.15
  sigma <- 4
  
  #build mu 
  mu <- a + b_gender*sim_data$gender + b_age*sim_data$age_c + b_meanGradient*sim_data$mean_gradient_c
  sim_data$VO2 <- rnorm (N, mu, sigma)
  
  #fit the model with ulam MCMC
  m <- ulam (
    alist(
      VO2 ~ dnorm (mu, sigma), 
      mu <- a + b_gender*gender + b_age*age_c + b_meanGradient*mean_gradient_c,
      a ~ dnorm (15, 2), 
      b_gender ~ dnorm (1.4, 0.5), 
      b_age ~ dnorm (-0.2, 0.05), 
      b_meanGradient ~ dnorm (-0.15, 0.05), 
      sigma ~ dunif (1,5)
    ),
    data = sim_data,
    chains = 4,
    cores = 4, 
    iter = 1000, 
    log_lik = FALSE
  )
  
  post <- extract.samples (m)
  
  #clinically meaningful thresholds have been set  
  detected_age <- mean (post$b_age < -0.15) > 0.90
  detected_gender <- mean (post$b_gender > 1) > 0.90 
  detected_gradient <- mean (post$b_meanGradient < -0.1) > 0.90
  
  return (c(detected_age, detected_gender, detected_gradient)) #this has been changes to return power per parameter 
}
#testing a single run of the simulation function 
sim_detect_effect_ulam_v2(50)

sample_sizes <- c(30, 50, 70, 90, 100)
n_sims <- 100 #will take a while to do, because of the number of simulations. Last time took around 40 minutes 

#loop coding below will display power calculation for each parameter's effect size  
power_results <- matrix (NA, nrow = length(sample_sizes), ncol = 3)
colnames(power_results) <- c("age", "gender", "gradient")

for (i in seq_along(sample_sizes)) {
  detections <- replicate (n_sims, sim_detect_effect_ulam_v2(sample_sizes[i]))
  
  power_results[i, ] <- rowMeans(detections, na.rm = TRUE)
}


#table of results 
results_table <- data.frame (
  sample_size = sample_sizes,
  power_age = power_results[, "age"], 
  power_gender = power_results[, "gender"], 
  power_gradient = power_results[, "gradient"]
)

print (results_table)
#sample_size power_age power_gender power_gradient
#1          30      0.39         0.07           0.57
#2          50      0.58         0.16           0.69
#3          70      0.53         0.15           0.73
#4          90      0.70         0.20           0.83
#5         100      0.65         0.26           0.80
#90 may be the sweet spot, one further simulation below with the in between numbers

#simulation with 30, 55, 75, 85, 95, 100
library (rethinking)

set.seed (100)
N <- 100

#doing computing for mean gradient values first 
mean <- 40 
sd <- 15
sdlog <- sqrt (log(1 + (sd^2 / mean^2)))
meanlog <- log(mean) - (sdlog^2 / 2)

sim_detect_effect_ulam_v2 <- function (N) {
  sim_data <- list (
    gender = rbinom (N, 1, 0.5), 
    age = rnorm (N, 70, 10), 
    mean_gradient = rlnorm (N, meanlog, sdlog)
  )
  sim_data$age_c <- sim_data$age - 70
  sim_data$mean_gradient_c <- sim_data$mean_gradient - 40 
  
  #simulate the "true" parameters 
  a <- 15
  b_gender <-  1.4
  b_age <- -0.2
  b_meanGradient <-  -0.15
  sigma <- 4
  
  #build mu 
  mu <- a + b_gender*sim_data$gender + b_age*sim_data$age_c + b_meanGradient*sim_data$mean_gradient_c
  sim_data$VO2 <- rnorm (N, mu, sigma)
  
  #fit the model with ulam MCMC
  m <- ulam (
    alist(
      VO2 ~ dnorm (mu, sigma), 
      mu <- a + b_gender*gender + b_age*age_c + b_meanGradient*mean_gradient_c,
      a ~ dnorm (15, 2), 
      b_gender ~ dnorm (1.4, 0.5), 
      b_age ~ dnorm (-0.2, 0.05), 
      b_meanGradient ~ dnorm (-0.15, 0.05), 
      sigma ~ dunif (1,5)
    ),
    data = sim_data,
    chains = 4,
    cores = 4, 
    iter = 1000, 
    log_lik = FALSE
  )
  
  post <- extract.samples (m)
  
  #clinically meaningful thresholds have been set  
  detected_age <- mean (post$b_age < -0.15) > 0.90
  detected_gender <- mean (post$b_gender > 1) > 0.90 
  detected_gradient <- mean (post$b_meanGradient < -0.1) > 0.90
  
  return (c(detected_age, detected_gender, detected_gradient)) #this has been changes to return power per parameter 
}
#testing a single run of the simulation function 
sim_detect_effect_ulam_v2(50)

sample_sizes <- c(30, 55, 75, 85, 95, 100)
n_sims <- 100 #will take a while to do, because of the number of simulations. Last time took around 40 minutes 

#loop coding below will display power calculation for each parameter's effect size  
power_results <- matrix (NA, nrow = length(sample_sizes), ncol = 3)
colnames(power_results) <- c("age", "gender", "gradient")

for (i in seq_along(sample_sizes)) {
  detections <- replicate (n_sims, sim_detect_effect_ulam_v2(sample_sizes[i]))
  
  power_results[i, ] <- rowMeans(detections, na.rm = TRUE)
}


#table of results 
results_table <- data.frame (
  sample_size = sample_sizes,
  power_age = power_results[, "age"], 
  power_gender = power_results[, "gender"], 
  power_gradient = power_results[, "gradient"]
)

print (results_table)

# sample_size power_age power_gender power_gradient
#1          30      0.39         0.07           0.57
#2          55      0.59         0.16           0.68
#3          75      0.52         0.17           0.78
#4          85      0.56         0.27           0.74
#5          95      0.70         0.21           0.85
#6         100      0.71         0.23           0.85

#95 and a 100 are almost identical , 80s range still looking unable to hold 0.80 threshold with simulation noise. 90 may be the safest option for a sample size