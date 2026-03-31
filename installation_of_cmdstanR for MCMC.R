#installing cmdstanr, needed for running MCMC ulam approximation 
remotes::install_github("stan-dev/cmdstanr")
cmdstanr::install_cmdstan(overwrite=TRUE)
cmdstanr::set_cmdstan_path()
cmdstanr::cmdstan_path()
