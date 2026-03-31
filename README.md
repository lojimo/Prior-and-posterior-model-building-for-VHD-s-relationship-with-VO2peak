# Prior and posterior model building 
This repository contains the R scripts for the determination of priors and subsequent posterior distribution for the effect of age, gender, and mean gradient in aortic stenosis (AS) on VO2max as determined by cardiopulmonary exercise testing (CPET). The priors and posterior distribution were then used to inform sample size calculations for our study examining physical capacity as assessed by CPETs for AS patients. We were looking for the most appropriate sample size with enough power to detect the effect size of mean gradient on VO2max. 

We first started of with determining prior distribution values based off current literature. To check the behaviour of the priors within the model, we conducted prior predictive simulation. After checking the priors were sound, we moved onto posterior disribution approximation. We first started conducting this using quadratic approximation, however this later proved ineffective. Quadratic approximation of the posterior failed to deal with increasing errors rates when we moved onto testing the model with sample size simulations. 

As such, we shifted from quadratic approximation of the posterior to approximating it via Hamiltonian Markov Chain Monte Carlo (MCMC). The new posterior distribution model handled simulation noise well and was sound when tested via posterior sampling. 
The final sample size calculations indicate that a sample of around 90 is the safest option when wanting to detect the effect size of mean gradient on VO2max. 
MCMC approximation for the posterior was done the ulam function in the rethinking package on R. 

The priors included an effect size of -0.15 per mmHg change, with a hypothesised 1.5 unit decrease in VO2max with 10 mmHg change in mean gradient. A sample size of 90 would have 83% power to detect this effect size according to our simulations. 

