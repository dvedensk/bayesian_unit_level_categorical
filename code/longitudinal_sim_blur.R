library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)
library(Matrix)
library(extraDistr)
library(BayesLogit)
library(mvtnorm)
library(survey)
library(mase)


here::i_am("code/longitudinal_sim_blur.R")

devtools::load_all(here::here("blur"))
source(here::here("code", "helper_functions.R"))

##########################
### Simulation Settings ##
##########################
set.seed(555)

out_dir <- here::here("results", "GAD2", "longitudinal")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

data_dir <- here::here("data", "ordinal")
load(here::here(data_dir, "empirical_samples_GAD2.RData"))
load(here::here("data", "unscaled_basis_functions.RData"))
Nsim <- length(samples)

pop_covars <- select(HPS_pop_long, starts_with("COVAR")) %>% names %>% paste(collapse="+")
covars <-  gsub("\\+COVAR.NEW", "", pop_covars)

pop_grouped <- HPS_pop_long %>%
                      group_by(WEEK, AREA, COVAR.1, COVAR.2, COVAR.3, COVAR.NEW) %>%
                      summarize(popsize=n())

popPsi <-  basis_funcs[pop_grouped$AREA, ]
popX <- model.matrix(as.formula(paste("~ -1 + ", covars)), data=pop_grouped)[,-1]

population_counts_by_time <- HPS_pop_long %>%
    group_by(WEEK) %>%
    summarize(N=n()) %>%
    select(N) %>%
    unlist()

n_mcmc <- 2000
n_burn <- 500
nsamp <- n_mcmc - n_burn

#Parallelization parameters 
numCores <- 10
cl<-makeCluster(numCores, type="FORK", outfile="") 
registerDoParallel(cl)

lon_fit <- foreach(i=1:Nsim,
                      .combine="rbind",
                      .verbose=TRUE,
                      .packages=c('dplyr')) %dopar% {

  set.seed(i + 15) 
  print(i)

  HPS_sample_long <- samples[[i]]
  areas <- HPS_sample_long$AREA
  timepoints <- HPS_sample_long$WEEK
  n_areas <- nlevels(areas)
  n_weeks <- max(HPS_sample_long$WEEK)

  svy_dir_est <- get_svy_direst(HPS_pop_long, HPS_sample_long)
  svy_dir_est$model <- "svy_direst"
  svy_dir_est$runtime <- NA

  scale_weights <- HPS_sample_long$SCALE_WEIGHT
  Psi <- basis_funcs[areas, ]
  X <- model.matrix(as.formula(paste("~ -1 + ", covars)), data=HPS_sample_long)[,-1]
  Y <- model.matrix( ~ 0 + RESPONSE,data=HPS_sample_long)
  K <- ncol(Y)

  vb_out <-  ulm(Y=Y,
                 X=X,
                 Psi=Psi,
                 weights=scale_weights,
                 n_samples=nsamp,
                 epsilon=1e-3,
                 response_type="ordinal",
                 algorithm="VB",
                 timepoints=timepoints,
                 prev_covar=HPS_sample_long$COVAR.NEW,
                 longitudinal=TRUE)

  vb_pp <- agg_predict(vb_out,
                       predX=popX,
                       predPsi=popPsi,
                       predTimes=pop_grouped$WEEK,
                       counts=pop_grouped$popsize,
                       prev_covar=pop_grouped$COVAR.NEW,
                       alpha=.05,
                       grouping_vars=c("AREA", "WEEK"),
                       pop_df=pop_grouped,
                       K=K)

  mcmc_out <-  ulm(Y=Y,
                   X=X,
                   Psi=Psi,
                   weights=scale_weights,
                   n_iter=n_mcmc,
                   n_burn=n_burn,
                   response_type="ordinal",
                   algorithm="MCMC",
                   timepoints=timepoints,
                   prev_covar=HPS_sample_long$COVAR.NEW,
                   longitudinal=TRUE)

  mcmc_pp <- agg_predict(mcmc_out,
                         predX=popX,
                         predPsi=popPsi,
                         predTimes=pop_grouped$WEEK,
                         counts=pop_grouped$popsize,
                         prev_covar=pop_grouped$COVAR.NEW,
                         alpha=.05,
                         grouping_vars=c("AREA", "WEEK"),
                         pop_df=pop_grouped,
                         K=K)
    
  mcmc_pp$model <- "mcmc-lon"
  vb_pp$model <- "vb-lon"
  mcmc_pp$runtime <- attr(mcmc_out, "runtime")
  vb_pp$runtime <- attr(vb_out, "runtime")

  results_df <- rbind(mcmc_pp, vb_pp, svy_dir_est)
  results_df$sim_num <- i

  save(
    results_df,
    file = file.path(out_dir, paste0(i, ".RData"))
  )

  return(results_df)
}
stopCluster(cl)

save(lon_fit, file=here::here("results", "GAD2" ,"longitudinal_simulation_results.RData"))
