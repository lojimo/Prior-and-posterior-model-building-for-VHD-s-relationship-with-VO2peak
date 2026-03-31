#now we need to do prior predictive simulation again to check what our prior (and subsequently posterior) distributions look like when we add in the AS effect size variable 
#this will be described in terms of amount of change in VO2 per mmHg change in mean gradient (mean gradient as our continuous measure of severity)
#previously decided of 1.5 unit change in VO2 per mmHg of mean gradient change 
#according to literature, MG mean should be around 40 mmHg, and most papers had ±15 , so rnorm (1, 40, 15)
library (rethinking)
#first step is to do predictive prior simulation to check that priors make sense 
set.seed (1e4)
a <- rnorm (1e4, 15, 2)
b_gender <- rnorm (1e4, 1.4, 0.5)
b_age <- rnorm (1e4, -0.2, 0.05)
sigma <- runif (1e4, 1, 5)
b_meanGradient <- rnorm (1e4, -1.5, 0.5)

#now simulating the predictors, the hypotehtical individuals 
gender <- rbinom (1e4, 1 , 0.5)
age <- rnorm (1e4, 70, 10)
age_c <- age - 70
#the difference here is that we need to do the mean gradient on a log-normal distribution because it is a right-skewed variable, and is not normally distributed 
#we can't just put 40 and 15 in as the mean and sd into rlnorm, need to convert them into log-scale parameters using the following formulas: 
#SD = log (1 + (SD^2/mean^2))
#mean = log(mean) - (SD^2/2)
mean <- 40 
sd <- 15
sdlog <- sqrt (log(1 + (sd^2 / mean^2)))
meanlog <- log(mean) - (sdlog^2 / 2)
mean_gradient <- rlnorm (1e4, meanlog, sdlog)

#building mu in the model, same as before except now we add the mean gradient variable 
mu <- a + b_gender*gender + b_age*age_c + b_meanGradient*mean_gradient

#simulate the outcomes 

VO2_prior <- rnorm (1e4, mu, sigma)
dens (VO2_prior)
precis (VO2_prior)

#the above didn't work because the b_meanGradient*mean_gradient was dominating the model, as roughly -60 contributing (-1.5 x 40)
#it was also messing with the intercept, as at mean gradient 0, the intercept was 15. This is not logically sound as mean gradient 0 is not meaningful
#what we need to do is center the mean gradient exactly as we did with age, so that at mean gradient 40, VO2 is at 15, and that is our intercept. 
#below is the above code but with mean gradient changed to be centred , helpful to run both for comparison 

library (rethinking)
#first step is to do predictive prior simulation to check that priors make sense 
set.seed (1e4)
a <- rnorm (1e4, 15, 2)
b_gender <- rnorm (1e4, 1.4, 0.5)
b_age <- rnorm (1e4, -0.2, 0.05)
sigma <- runif (1e4, 1, 5)
b_meanGradient <- rnorm (1e4, -0.15, 0.05)

#now simulating the predictors, the hypotehtical individuals 
gender <- rbinom (1e4, 1 , 0.5)
age <- rnorm (1e4, 70, 10)
age_c <- age - 70
#the difference here is that we need to do the mean gradient on a log-normal distribution because it is a right-skewed variable, and is not normally distributed 
#we can't just put 40 and 15 in as the mean and sd into rlnorm, need to convert them into log-scale parameters using the following formulas: 
#SD = log (1 + (SD^2/mean^2))
#mean = log(mean) - (SD^2/2)
mean <- 40 
sd <- 15
sdlog <- sqrt (log(1 + (sd^2 / mean^2)))
meanlog <- log(mean) - (sdlog^2 / 2)
mean_gradient <- rlnorm (1e4, meanlog, sdlog)
mean_gradient_c <- mean_gradient - 40 

#building mu in the model, same as before except now we add the mean gradient variable 
mu <- a + b_gender*gender + b_age*age_c + b_meanGradient*mean_gradient_c

#simulate the outcomes 

VO2_prior <- rnorm (1e4, mu, sigma)
dens (VO2_prior)
precis (VO2_prior)
#this is looking better than the previous one. Note also changed the effect size of mean gradient to be -0.15. Before it was computing -1.5 change ber mmHg, which is too much, and logically wrong 
#if we are predicting a conservative -1.5 change per 10 mmHg (say if someone was at 30mmHg in modrate then jumps to severe at 40 mmHg), then -0.15 change per mmHg would make more sense 

#can plot variables against each other to check scatter and relationship 
plot (mean_gradient, VO2_prior, 
      main = "Scatter plot of mean gradient vs. VO2max",
      xlab = "mean gradient in mmHg", 
      ylab = "VO2 max", 
      pch = 19, 
      col = rgb (0, 0, 1, alpha = 0.1)
      )
abline (lm(VO2_prior ~ mean_gradient), col = "red", lwd = 1)

plot (age , VO2_prior, 
      main = "Scatter plot of age vs. VO2max", 
      xlab = "age", 
      ylab = "VO2max",
      pch = 19, 
      col = rgb (0, 0, 1, alpha = 0.1)
      )
abline (lm(VO2_prior ~ age), col = "red", lwd =1)
#both of these plots show VO2 max concentrated values around 15-7 but could be erring on the side of higehr than 20 because of the standard deviation
#assume that is okay? healthy amount of scatter, not too concentrated 

#now onto building the posterior
set.seed (100)
N <- 100 

#doing the computing for the mean gradient values first 
mean <- 40 
sd <- 15
sdlog <- sqrt (log(1 + (sd^2 / mean^2)))
meanlog <- log(mean) - (sdlog^2 / 2)

#put into sim_data
sim_data <- list (
  gender = rbinom (N, 1, 0.5), 
  age = rnorm (N, 70, 10), 
  mean_gradient = rlnorm (N, meanlog, sdlog)
)

#do the centering 
#doing the computing for the mean gradient values first 
sim_data$mean_gradient_c <- sim_data$mean_gradient - 40 
sim_data$age_c <- sim_data$age - 70

#now for the assumed parameters, but still keeping uncertainty in them rnorm function rather than using true values 
a <- rnorm (1, 15, 2)
b_gender <- rnorm (1, 1.4, 0.5)
b_age <- rnorm (1, -0.2, 0.05)
b_meanGradient <- rnorm (1, -0.15, 0.05)
sigma <- runif (1, 1, 5)

#building mu in the model 
mu <- a + b_gender*sim_data$gender + b_age*sim_data$age_c + b_meanGradient*sim_data$mean_gradient_c

sim_data$VO2 <- rnorm (N, mu, sigma)

#MCMC approximation 
mv5 <- ulam (
  alist ( 
    VO2 ~ dnorm (mu , sigma), 
    mu <- a + b_gender*gender + b_age*age_c + b_meanGradient*mean_gradient_c,
    a ~ dnorm (15 , 2), 
    b_gender ~ dnorm (1.4, 0.5), 
    b_age ~ dnorm (-0.2, 0.05), 
    b_meanGradient ~ dnorm (-0.15, 0.05), 
    sigma ~ dunif (1, 5)
    ),
  data = sim_data, 
  chains = 4, 
  cores = 4, 
  iter = 1000, 
  log_lik = FALSE
)
precis (mv5, depth = 2)
pairs (mv5)

post <- extract.samples (mv5, 1000)
precis (post)
diag (vcov(mv5))
cov2cor (vcov(mv5))
plot (VO2 ~ gender + age + mean_gradient , data = sim_data, col = rangi2)
#all of the above are checks for the posterior to make sure it is making sense, and it is 
#time to move onto simulation :)
