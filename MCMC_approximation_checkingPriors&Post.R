#this is going to be posterior sampling from the MCMC model to check that it is sound 
library (rethinking)
#first step is to do preditive prior simulation to check that make sense (essentially the same code as before)
set.seed (1e4)
a <- rnorm (1e4, 15, 2)
b_gender <- rnorm (1e4, 1.4, 0.5)
b_age <- rnorm (1e4, -0.2, 0.05)
sigma <- runif (1e4, 1, 5)

#now simulate the predictors, the hypotehtical individuals 
gender <- rbinom (1e4, 1, 0.5)
age <- rnorm (1e4, 70, 10)
age_c <- age - 70

#build mu in the model 
mu <- a + b_gender*gender + b_age*age_c

#simulate the outcomes 
VO2_prior <- rnorm (1e4, mu, sigma)
dens(VO2_prior)
precis (VO2_prior)
#this is a prior that makes sense with a healthy amount of variation and a mean that makes sense 
#can still plot two variables against each other to see how their relationship is, for example plotting age against VO2 
plot (age, VO2_prior, 
      main = "Scatter plot of age vs. VO2max",
      xlab = "age",
      ylab = "VO2max",
      pch=19,
      col = rgb (0, 0, 1, alpha = 0.1)
      )
abline (lm(VO2_prior ~ age), col = "red", lwd = 1)
#this shows that the priors are having the correct relationship with each other. Age is concentrated at 70, there is a healthy amount of scatter, negative decline in V02 as age increases which is what we want to see 


#now onto building the posterior 
set.seed (100)
N <- 100 
sim_data <- list(
  gender = rbinom (N, 1, 0.5),
  age = rnorm (N, 70, 10)
)
sim_data$age_c <- sim_data$age - 70

#now for the assumed parameters but still keeping uncertainty in them with rnorm rather than using true values 
a <- rnorm (1, 15, 2)
b_gender <- rnorm (1, 1.4, 0.5)
b_age <- rnorm (1, -0.2, 0.05)
sigma <- runif (1, 1, 5)

#building mu into the model 
mu <- a + b_gender*sim_data$gender + b_age*sim_data$age_c

sim_data$VO2 <- rnorm (N, mu, sigma)

#MCMC approximation 
mv4 <- ulam (
  alist ( 
    VO2 ~ dnorm (mu, sigma),
    mu <- a + b_gender*gender + b_age*age_c,
    a ~ dnorm (15,2),
    b_gender ~ dnorm (1.4, 0.5),
    b_age ~ dnorm (-0.2, 0.05),
    sigma ~ dunif (1,5)
    ),
  data = sim_data,
  chains = 4, 
  cores = 4, 
  iter = 1000, 
  log_lik = FALSE
)

precis (mv4, depth = 2)

post <- extract.samples (mv4, 1000)
precis (post)

show (mv4)

pairs (mv4)


