## ----libraries, echo=T, results='hide', message=F, warning=F------------------------------------------------------------
library(tidyverse)
library(broom)
library(rstan)
library(gdata)
library(bayesplot)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

knitr::purl(input = "nacaAnalysis.Rmd",output = "nacaAnalysis.R")


## ----priors-------------------------------------------------------------------------------------------------------------
prior_logkChemInt <- -3
priorSD_logkChemInt <- 0.5
prior_logkBleachInt <- -4
priorSD_logkBleachInt <- 0.5
prior_mChem <- -0.05
priorSD_mChem <- 0.5
prior_mBleach <- 0
priorSD_mBleach <- 0.5
prior_ext <- 0.02
priorSD_ext <- 0.01
prior_delay <- 2
priorSD_delay <- 0.01

priorSD_sigma <- 0.05

debug <- 0; #0 for no debug output, 1 for printing gen values, 2 for printing params

init_params <- list(list(logkChemInt=prior_logkChemInt,
                    logkBleachInt=prior_logkBleachInt,
                    mChem=prior_mChem,
                    mBleach=prior_mBleach,
                    ext=prior_ext,
                    delay=prior_delay,
                    sigma = 0.02))
init_params_4x <- rep(init_params,4)


## ----conditions---------------------------------------------------------------------------------------------------------
initCys <- 50
initDTP <- 100


## ----nullData-----------------------------------------------------------------------------------------------------------
N <- 2
time <- c(5,6)
gdn <- c(0,0)
absorbance <- c(0,0)
numTrials <- 1
trialStarts <- c(1,3)


## ----nData--------------------------------------------------------------------------------------------------------------
nDat_gen <- 100
deadTime <- 5


## ----formatStanSim------------------------------------------------------------------------------------------------------
stan_data <-list(prior_logkChemInt=prior_logkChemInt, priorSD_logkChemInt=priorSD_logkChemInt,
                 prior_logkBleachInt=prior_logkBleachInt, priorSD_logkBleachInt=priorSD_logkBleachInt,
                 prior_mChem=prior_mChem, priorSD_mChem=priorSD_mChem, prior_mBleach=prior_mBleach,
                 priorSD_mBleach=priorSD_mBleach, prior_ext=prior_ext,priorSD_ext=priorSD_ext,
                 prior_delay=prior_delay, priorSD_delay=priorSD_delay, initCys=initCys, initDTP=initDTP, N=N,
                 time=time, gdn=gdn, absorbance=absorbance,
                 numTrials=numTrials,trialStarts=trialStarts,nDat_gen=nDat_gen,deadTime=deadTime,
                 priorSD_sigma=priorSD_sigma,debug=debug)


## ----stanSimulate-------------------------------------------------------------------------------------------------------
stan_model<- "nacaModel.stan"
#stanc(file = stan_model, verbose = TRUE)
sim_data <- stan(file = stan_model,
             data = stan_data, init = init_params,
             chains = 1, iter = 1)


## ----plot-sim-----------------------------------------------------------------------------------------------------------
posterior <- rstan::extract(sim_data)

extractData <- function(fitObject,name,gdn){
  vals <- as.matrix(fitObject,pars = name)
  times <- as.matrix(fitObject, pars = "times_gen")
  vals <- as.vector(vals)
  times <- as.vector(times)
  pred_data <- data.frame("Time" = times, "Gdn" = gdn, "Abs" = vals)
  return(pred_data)
}
  
pred_data_0 <- extractData(sim_data,"pred_gdn_0",0)
pred_data_3 <- extractData(sim_data,"pred_gdn_3",3)
pred_data_6 <- extractData(sim_data,"pred_gdn_6",6)

pred_trial_starts = c(1,
                      nDat_gen+1,
                      nDat_gen*2+1,
                      nDat_gen*3+1)

numPredTrials = 3

pred_data_ALL <- rbind(pred_data_0,pred_data_3,pred_data_6)
ggplot(pred_data_ALL,aes(x = Time, y = Abs,color = factor(Gdn))) +
  geom_point() +
  labs(x = "Time (s)", y = "Absorbance", title = "Predicted Values") +
  scale_x_continuous(trans='log10', labels = scales::label_comma(), limits = c(1,1800))



## ----format-sim---------------------------------------------------------------------------------------------------------
pred_stan_data <- stan_data
# the same as before except for these changes
pred_stan_data$N <- 3*nDat_gen
pred_stan_data$time <- pull(pred_data_ALL,Time)
pred_stan_data$gdn <- pull(pred_data_ALL,Gdn)
pred_stan_data$absorbance <- pred_data_ALL %>%
  mutate(Abs = if_else(Abs<=0,0.0001,Abs)) %>%
  pull(Abs)
pred_stan_data$numTrials <- numPredTrials
pred_stan_data$trialStarts <- pred_trial_starts


## ----fit-same-priors----------------------------------------------------------------------------------------------------
sim_fit <- stan(file = stan_model,
             data = pred_stan_data, init = init_params_4x,
             chains = 4,iter = 2000)


## ----posterior-plots----------------------------------------------------------------------------------------------------
sim_fit_posterior <- rstan::extract(sim_fit)

par(mfrow = c(2,2))

plot(density(sim_fit_posterior$logkChemInt), main = "logkChemInt")
abline(v = prior_logkChemInt, col = 4, lty = 2)

plot(density(sim_fit_posterior$logkBleachInt), main = "logkBleachInt")
abline(v = prior_logkBleachInt, col = 4, lty = 2)

plot(density(sim_fit_posterior$mChem), main = "mChem")
abline(v = prior_mChem, col = 4, lty = 2)

plot(density(sim_fit_posterior$mBleach), main = "mBleach")
abline(v = prior_mBleach, col = 4, lty = 2)


## ----diagnostic-plots---------------------------------------------------------------------------------------------------
traceplot(sim_fit,pars=c("logkChemInt","logkBleachInt","mChem","mBleach"))


## ----save-targets-------------------------------------------------------------------------------------------------------
target_logkChemInt <- prior_logkChemInt
target_logkBleachInt <- prior_logkBleachInt
target_mChem <- prior_mChem
target_mBleach <- prior_mBleach
target_ext <- prior_ext
target_delay <- prior_delay


## ----change-priors------------------------------------------------------------------------------------------------------
prior_logkChemInt <- -2
priorSD_logkChemInt <- 2
prior_logkBleachInt <- -1
priorSD_logkBleachInt <- 2
prior_mChem <- 0.05
priorSD_mChem <- 2
prior_mBleach <- 0.1
priorSD_mBleach <- 2
prior_ext <- 0.01
priorSD_ext <- 0.02
prior_delay <- 5
priorSD_delay <- 1

priorSD_sigma <- 0.05

debug <- 1; #0 for no debug output, 1 for printing gen values, 2 for printing params

init_params <- list(list(logkChemInt=prior_logkChemInt,
                    logkBleachInt=prior_logkBleachInt,
                    mChem=prior_mChem,
                    mBleach=prior_mBleach,
                    ext=prior_ext,
                    delay=prior_delay,
                    sigma = 0.02))
init_params_4x <- rep(init_params,4)


## ----format-sim-2-------------------------------------------------------------------------------------------------------
pred_stan_data_2 <- pred_stan_data
# the same as the previous fit except for changing the priors

pred_stan_data_2$prior_logkChemInt <- prior_logkChemInt
pred_stan_data_2$priorSD_logkChemInt <- priorSD_logkChemInt
pred_stan_data_2$prior_logkBleachInt <- prior_logkBleachInt
pred_stan_data_2$priorSD_logkBleachInt <- priorSD_logkBleachInt
pred_stan_data_2$prior_mChem <- prior_mChem
pred_stan_data_2$priorSD_mChem <- priorSD_mChem
pred_stan_data_2$prior_mBleach <- prior_mBleach
pred_stan_data_2$priorSD_mBleach <- priorSD_mBleach
pred_stan_data_2$prior_ext <- prior_ext
pred_stan_data_2$priorSD_ext <- priorSD_ext
pred_stan_data_2$prior_delay <- prior_delay
pred_stan_data_2$priorSD_delay <- priorSD_delay
pred_stan_data_2$debug <- debug


## ----fit-diff-priors----------------------------------------------------------------------------------------------------
sim_fit_diff <- stan(file = stan_model,
             data = pred_stan_data_2, init = init_params_4x,
             chains = 4,iter = 2000)


## ----posterior-plots-2--------------------------------------------------------------------------------------------------
sim_fit_diff_posterior <- rstan::extract(sim_fit_diff)

par(mfrow = c(2,2))

plot(density(sim_fit_diff_posterior$logkChemInt), main = "logkChemInt",xlim = c(-3.1,-1.9))
abline(v = prior_logkChemInt, col = 4, lty = 2)
abline(v = target_logkChemInt, col = "red", lty = 2)

plot(density(sim_fit_diff_posterior$logkBleachInt), main = "logkBleachInt",xlim = c(-5,-0.9))
abline(v = prior_logkBleachInt, col = 4, lty = 2)
abline(v = target_logkBleachInt, col = "red", lty = 2)

plot(density(sim_fit_diff_posterior$mChem), main = "mChem",xlim = c(-0.08,0.06))
abline(v = prior_mChem, col = 4, lty = 2)
abline(v = target_mChem, col = "red", lty = 2)

plot(density(sim_fit_diff_posterior$mBleach), main = "mBleach",xlim = c(-5,0.12))
abline(v = prior_mBleach, col = 4, lty = 2)
abline(v = target_mBleach, col = "red", lty = 2)


## ----diagnostic-plots-2-------------------------------------------------------------------------------------------------
traceplot(sim_fit_diff,pars=c("logkChemInt","logkBleachInt","mChem","mBleach"))

