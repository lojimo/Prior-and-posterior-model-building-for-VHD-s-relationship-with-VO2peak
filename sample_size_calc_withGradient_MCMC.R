#sample size calculation and simulation using MCMC approximation for the posterior with the mean gradient added in as a predictor variable now 
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
  a <- rnorm (1, 15, 2)
  b_gender <- rnorm (1, 1.4, 0.5)
  b_age <- rnorm (1, -0.2, 0.05)
  b_meanGradient <- rnorm (1, -0.15, 0.05)
  sigma <- runif (1, 1, 5)
  
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
  
  #detection criteria , for now we will make this relative to 0, won't make it clinically meaningful 
  detected_age <- mean (post$b_age < 0) > 0.9
  detected_gender <- mean (post$b_gender > 0) > 0.9 
  detected_gradient <- mean (post$b_meanGradient < 0) > 0.9
  
  return (detected_age & detected_gender & detected_gradient)
}
#testing a single run of the simulation function 
sim_detect_effect_ulam_v2(50)

sample_sizes <- c(30, 40, 60, 80, 100)
n_sims <- 100 #will take a while to do, because of the number of simulations. Last time took around 40 minutes 

power_results <- numeric (length(sample_sizes))

for (i in seq_along(sample_sizes)) {
  detections <- replicate (n_sims, sim_detect_effect_ulam_v2(sample_sizes[i]))
  power_results[i] <- mean (detections, na.rm = TRUE)
}

#results table code 
results_table <- data.frame(
  sample_size = sample_sizes,
  power_to_detect = power_results
)
print(results_table)

table (detections, useNA = "always")

#next step here would be to stress test it against different effect sizes. Currently it is outputting that detecting the direction of change, whether the effect size is greater or less than 0 is easy. 
#this is the table it outputs for effect size detection of >0 or <0 90% of the time: 
#sample_size power_to_detect
#1          30            1.00
#2          40            0.99
#3          60            0.99
#4          80            1.00
#5         100            1.00

#so should try testing to see whether this changes at all when we change >0.90 to >0.95 , and what happens when we make it a meaningful effect size detection, rather than just less than or more than zero 
#the other thing we could consider here would be going back to adding more variation in the priors 



#trying different effect sizes to detect 
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
  a <- rnorm (1, 15, 2)
  b_gender <- rnorm (1, 1.4, 0.5)
  b_age <- rnorm (1, -0.2, 0.05)
  b_meanGradient <- rnorm (1, -0.15, 0.05)
  sigma <- runif (1, 1, 5)
  
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
  
  #detection criteria , for now we will make this relative to 0, won't make it clinically meaningful 
  detected_age <- mean (post$b_age < -0.1) > 0.9
  detected_gender <- mean (post$b_gender > 1.0) > 0.9 
  detected_gradient <- mean (post$b_meanGradient < -0.05) > 0.9
  
  return (detected_age & detected_gender & detected_gradient)
}
#testing a single run of the simulation function 
sim_detect_effect_ulam_v2(50)

sample_sizes <- c(30, 40, 60, 80, 100)
n_sims <- 100 #will take a while to do, because of the number of simulations. Last time took around 40 minutes 

power_results <- numeric (length(sample_sizes))

for (i in seq_along(sample_sizes)) {
  detections <- replicate (n_sims, sim_detect_effect_ulam_v2(sample_sizes[i]))
  power_results[i] <- mean (detections, na.rm = TRUE)
}

#results table code 
results_table <- data.frame(
  sample_size = sample_sizes,
  power_to_detect = power_results
)
print(results_table)
#slightly weird output when you do the "meaningful change" , where 80 sample size has higher power than the 100 sample size 
#sample_size power_to_detect
#1          30            0.22
#2          40            0.23
#3          60            0.26
#4          80            0.42
#5         100            0.36


#now going to go back to the greater than or less than 0 as our threshold for detection, but going to change from 0.90 to 0.95 , make it stricter 
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
  a <- rnorm (1, 15, 2)
  b_gender <- rnorm (1, 1.4, 0.5)
  b_age <- rnorm (1, -0.2, 0.05)
  b_meanGradient <- rnorm (1, -0.15, 0.05)
  sigma <- runif (1, 1, 5)
  
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
  
  #detection criteria , for now we will make this relative to 0, won't make it clinically meaningful 
  detected_age <- mean (post$b_age < 0) > 0.95
  detected_gender <- mean (post$b_gender > 0) > 0.95 
  detected_gradient <- mean (post$b_meanGradient < 0) > 0.95
  
  return (detected_age & detected_gender & detected_gradient)
}
#testing a single run of the simulation function 
sim_detect_effect_ulam_v2(50)

sample_sizes <- c(30, 40, 60, 80, 100)
n_sims <- 100 #will take a while to do, because of the number of simulations. Last time took around 40 minutes 

power_results <- numeric (length(sample_sizes))

for (i in seq_along(sample_sizes)) {
  detections <- replicate (n_sims, sim_detect_effect_ulam_v2(sample_sizes[i]))
  power_results[i] <- mean (detections, na.rm = TRUE)
}

#results table code 
results_table <- data.frame(
  sample_size = sample_sizes,
  power_to_detect = power_results
)
print(results_table)
#this is the output, again a bit weird in that it shows 30 as having sliightlyy more power than 40, but still showing almost all of them to have power all the time 
#sample_size power_to_detect
#1          30            0.99
#2          40            0.98
#3          60            0.99
#4          80            1.00
#5         100            1.00
#the next step is to run it with the 0.95 and the meaningful effect size detection


#meaningful effect size for age, if decreases by 0.2 units per year roughly , then we would want a meaningful effect size detection to be around -0.05 or -0.10 
#meaningful effect size for gender, would be more than one, since we are saying increases by 40% so anything more than 1 is showing a gender effect size 
#meaningful effect size for the gradient, if decreases by 0.15 units per year (for a 1.5 unit change per decade, being conservative) then we would want to say a detection of anything more than 0.05 change is meaningful. 

#now with the 0.95 and the meaningful effect sizes 
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
  a <- rnorm (1, 15, 2)
  b_gender <- rnorm (1, 1.4, 0.5)
  b_age <- rnorm (1, -0.2, 0.05)
  b_meanGradient <- rnorm (1, -0.15, 0.05)
  sigma <- runif (1, 1, 5)
  
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
  
  #detection criteria , for now we will make this relative to 0, won't make it clinically meaningful 
  detected_age <- mean (post$b_age < -0.05) > 0.95
  detected_gender <- mean (post$b_gender > 1) > 0.95 
  detected_gradient <- mean (post$b_meanGradient < -0.05) > 0.95
  
  return (detected_age & detected_gender & detected_gradient)
}
#testing a single run of the simulation function 
sim_detect_effect_ulam_v2(50)

sample_sizes <- c(30, 40, 60, 80, 100)
n_sims <- 100 #will take a while to do, because of the number of simulations. Last time took around 40 minutes 

power_results <- numeric (length(sample_sizes))

for (i in seq_along(sample_sizes)) {
  detections <- replicate (n_sims, sim_detect_effect_ulam_v2(sample_sizes[i]))
  power_results[i] <- mean (detections, na.rm = TRUE)
}

#results table code 
results_table <- data.frame(
  sample_size = sample_sizes,
  power_to_detect = power_results
)
print(results_table)
#new results look weak across all sample sizes :(
#sample_size power_to_detect
#1          30            0.06
#2          40            0.17
#3          60            0.17
#4          80            0.29
#5         100            0.24

#simulation with these newly defined meaingful effect sizes, but with detection >0.90 rather than >0.95 to see if there is any change at all to these results 

#meaningful effect size for age, if decreases by 0.2 units per year roughly , then we would want a meaningful effect size detection to be around -0.05 or -0.10 
#meaningful effect size for gender, would be more than one, since we are saying increases by 40% so anything more than 1 is showing a gender effect size 
#meaningful effect size for the gradient, if decreases by 0.15 units per year (for a 1.5 unit change per decade, being conservative) then we would want to say a detection of anything more than 0.05 change is meaningful. 

#now with meaningful effect sizes as before, but 0.90 instead of 0.95  
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
  a <- rnorm (1, 15, 2)
  b_gender <- rnorm (1, 1.4, 0.5)
  b_age <- rnorm (1, -0.2, 0.05)
  b_meanGradient <- rnorm (1, -0.15, 0.05)
  sigma <- runif (1, 1, 5)
  
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
  
  #detection criteria , for now we will make this relative to 0, won't make it clinically meaningful 
  detected_age <- mean (post$b_age < -0.05) > 0.90
  detected_gender <- mean (post$b_gender > 1) > 0.90 
  detected_gradient <- mean (post$b_meanGradient < -0.05) > 0.90
  
  return (detected_age & detected_gender & detected_gradient) 
}
#testing a single run of the simulation function 
sim_detect_effect_ulam_v2(50)

sample_sizes <- c(30, 40, 60, 80, 100)
n_sims <- 100 #will take a while to do, because of the number of simulations. Last time took around 40 minutes 

power_results <- numeric (length(sample_sizes))

for (i in seq_along(sample_sizes)) {
  detections <- replicate (n_sims, sim_detect_effect_ulam_v2(sample_sizes[i]))
  power_results[i] <- mean (detections, na.rm = TRUE)
}

mean (is.na(detections)) #this is to see how many failed simulations are there 
#results table code 
results_table <- data.frame(
  sample_size = sample_sizes,
  power_to_detect = power_results
)
print(results_table)
#again some weird results where 80 has more power than 100, but the reason for that could be the fact that our "true parameters" themselves are normally distributed, so there is variation each time the simulation runs, the way to fix this would be to assign them as "true" parameters from the start. 
#sample size      powe_to_detect 
#1 30               0.24 
#2 40               0.24 
#3 60               0.31 
#4 80               0.43 
#5 100              0.38 

#instead of doing that straight away, going to run one where it will return power per parameter to see if there is a specific parameter that is causing issue
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
  a <- rnorm (1, 15, 2)
  b_gender <- rnorm (1, 1.4, 0.5)
  b_age <- rnorm (1, -0.2, 0.05)
  b_meanGradient <- rnorm (1, -0.15, 0.05)
  sigma <- runif (1, 1, 5)
  
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
  
  #detection criteria , for now we will make this relative to 0, won't make it clinically meaningful 
  detected_age <- mean (post$b_age < -0.05) > 0.90
  detected_gender <- mean (post$b_gender > 1) > 0.90 
  detected_gradient <- mean (post$b_meanGradient < -0.05) > 0.90
  
  return (c(detected_age, detected_gender, detected_gradient)) #this has been changes to return power per parameter 
}
#testing a single run of the simulation function 
sim_detect_effect_ulam_v2(50)

sample_sizes <- c(30, 40, 60, 80, 100)
n_sims <- 100 #will take a while to do, because of the number of simulations. Last time took around 40 minutes 

#need to correct the loop below to display the power per sample size 
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

#meaningful effect size detection 90% of the time, power per parameter these are the results. So the main issue seems to be happening with gender effect size detection.  
#sample_size power_age power_gender power_gradient
#1          30         1         0.25           0.93
#2          40         1         0.26           0.94
#3          60         1         0.32           0.94
#4          80         1         0.44           0.99
#5         100         1         0.40           0.96

#an option here could be to change the parameter definition inside the simulation from the simulated rnorm parameters to "true" parameters. 
#currently each simulation has a slightly different simulated parameter which could mean simulation noise accounts for the 0.44 vs. 0.40 in the 80 vs. 100 we are seeing 

#going to run a sequential analysis with true parameters but also with higher thresholds for meaningful detections 
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
  
  detected_age <- mean (post$b_age < -0.15) > 0.90
  detected_gender <- mean (post$b_gender > 1) > 0.90 
  detected_gradient <- mean (post$b_meanGradient < -0.1) > 0.90
  
  return (c(detected_age, detected_gender, detected_gradient)) #this has been changes to return power per parameter 
}
#testing a single run of the simulation function 
sim_detect_effect_ulam_v2(50)

sample_sizes <- c(30, 40, 60, 80, 100)
n_sims <- 100 #will take a while to do, because of the number of simulations. Last time took around 40 minutes 

#need to correct the loop below to display the power per sample size 
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
#the results from this are solid, no nonsensical 1s , still doing the weird thing with 80 as the sample size being higher than 100 but better 
#sample_size power_age power_gender power_gradient
#1          30      0.53         0.13           0.53
#2          40      0.53         0.11           0.65
#3          60      0.62         0.15           0.70
#4          80      0.65         0.29           0.85
#5         100      0.67         0.24           0.81

#Doing this again, but increasing the threshold for detection to be very close to the proposed likeliest effect size 
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
  
  #detection criteria , for now we will make this relative to 0, won't make it clinically meaningful 
  detected_age <- mean (post$b_age < -0.19) > 0.90
  detected_gender <- mean (post$b_gender > 1) > 0.90 
  detected_gradient <- mean (post$b_meanGradient < -0.14) > 0.90
  
  return (c(detected_age, detected_gender, detected_gradient)) #this has been changes to return power per parameter 
}
#testing a single run of the simulation function 
sim_detect_effect_ulam_v2(50)

sample_sizes <- c(30, 40, 60, 80, 100)
n_sims <- 100 #will take a while to do, because of the number of simulations. Last time took around 40 minutes 

#need to correct the loop below to display the power per sample size 
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

#weird results here with smaller sample sizes having higher power sometimes, but similar to what you would expect, it is not detecting much when it comes to being higher than the proposed effect size 
#sample_size power_age power_gender power_gradient
#1          30      0.04         0.13           0.12
#2          40      0.03         0.11           0.07
#3          60      0.13         0.15           0.13
#4          80      0.08         0.29           0.14
#5         100      0.10         0.24           0.13

#now need to run the sample size calculations for the in between ones , 50, 70, 80, and 90 
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
  
  #detection criteria , for now we will make this relative to 0, won't make it clinically meaningful 
  detected_age <- mean (post$b_age < -0.15) > 0.90
  detected_gender <- mean (post$b_gender > 1) > 0.90 
  detected_gradient <- mean (post$b_meanGradient < -0.1) > 0.90
  
  return (c(detected_age, detected_gender, detected_gradient)) #this has been changes to return power per parameter 
}
#testing a single run of the simulation function 
sim_detect_effect_ulam_v2(50)

sample_sizes <- c(30, 50, 70, 80, 90)
n_sims <- 100 #will take a while to do, because of the number of simulations. Last time took around 40 minutes 

#need to correct the loop below to display the power per sample size 
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

#these were the results, some things here that don't quite make sense: 
#sample_size power_age power_gender power_gradient
#1          30      0.53         0.13           0.53
#2          50      0.58         0.16           0.68
#3          70      0.49         0.28           0.74
#4          80      0.63         0.22           0.76
#5          90      0.59         0.33           0.76
#80 consistently lower with the power across all three parameters now that sample sizes have changed. same sample size of 80 had a 0.85 power to detect the mean gradient before. 
#going to run the same simulation but with 80 as the only sample size to see what that on its own comes out to. 
detections <- replicate (n_sims, sim_detect_effect_ulam_v2(sample_sizes[4]))
power_results [4, ] <- rowMeans(detections, na.rm = TRUE)
results_table <- data.frame (
  sample_size = sample_sizes,
  power_age = power_results[, "age"], 
  power_gender = power_results[, "gender"], 
  power_gradient = power_results[, "gradient"]
)

print (results_table)
#these are the results with running it just for 80: 
#sample_size power_age power_gender power_gradient
#1          30        NA           NA             NA
#2          50        NA           NA             NA
#3          70        NA           NA             NA
#4          80      0.63         0.25           0.79
#5          90        NA           NA             NA
#so the mean gradient outcme is inching closer to 0.80 , still a way off from 0.85 that we saw , might be useful to run the original set again 
#running that on a quarto document, will print the outcomes below here. 

#same simulation but for the following sample size set: 30, 55, 75, 85, 100
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
  
  #clinically meaningful thresholds have been set!  
  detected_age <- mean (post$b_age < -0.15) > 0.90
  detected_gender <- mean (post$b_gender > 1) > 0.90 
  detected_gradient <- mean (post$b_meanGradient < -0.1) > 0.90
  
  return (c(detected_age, detected_gender, detected_gradient)) #this has been changes to return power per parameter 
}
#testing a single run of the simulation function 
sim_detect_effect_ulam_v2(50)

sample_sizes <- c(30, 55, 75, 85, 100)
n_sims <- 100 #will take a while to do, because of the number of simulations. Last time took around 40 minutes 

#need to correct the loop below to display the power per sample size 
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
#again results here are very weird in that they are substantially different from when you run it with the sample size set of 30, 40, 60, 80, and 100 
#sample_size power_age power_gender power_gradient
#1          30      0.39         0.07           0.57
#2          55      0.59         0.16           0.68
#3          75      0.52         0.17           0.78
#4          85      0.56         0.27           0.74
#5         100      0.67         0.25           0.84