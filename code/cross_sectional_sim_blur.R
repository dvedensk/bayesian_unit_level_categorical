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

here::i_am("code/cross_sectional_sim_blur.R")

devtools::load_all(here::here("blur"))
source(here::here("code", "helper_functions.R"))

##########################
### Simulation Settings ##
##########################
set.seed(555)

out_dir <- here::here("results", "GAD2", "cross_sectional")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data_dir <- here::here("data", "ordinal")
load(file.path(data_dir, "empirical_samples_GAD2.RData"))
load(here::here("data", "unscaled_basis_functions.RData"))
Nsim <- length(samples)

pop_covars <- select(HPS_pop_long, starts_with("COVAR")) %>% names %>% paste(collapse="+")
week_covars <-  gsub("\\+COVAR.NEW", "", pop_covars)

pop_grouped_bulm <- HPS_pop_long %>%
                      group_by(WEEK, AREA, COVAR.1, COVAR.2, COVAR.3) %>%
                         summarize(popsize=n())

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

cs_fit <-foreach(i=1:Nsim, .combine="rbind", .verbose=TRUE, .packages=c('dplyr')) %dopar% {
############################
### Generate Sample Data ###
############################
  set.seed(i + 15) 
  print(i)

  ### Process sample
  HPS_sample_long <- samples[[i]]
  areas <- HPS_sample_long$AREA
  n_areas <- nlevels(areas)
  n_weeks <- max(HPS_sample_long$WEEK)

  ##Fit non-time-dependent
  results_df <- c()
  for(tt in 1:n_weeks){
    print(paste("Fitting week", tt))

    HPS_week <- filter(HPS_sample_long, WEEK==tt) %>% select(-COVAR.NEW)
    N_week <- nrow(HPS_week)
    scale_weights <- HPS_week$PWEIGHT * N_week/sum(HPS_week$PWEIGHT)
    Y_week <- model.matrix( ~ 0 + RESPONSE,data=HPS_week)
    Psi_week <- basis_funcs[HPS_week$AREA, ]
    X_week <- model.matrix(as.formula(paste("~", week_covars)), data=HPS_week)[,-1]

    pop_grouped_week <- filter(pop_grouped_bulm, WEEK==tt) %>%
        group_by(WEEK, AREA, across(starts_with("COVAR"))) 
    
    popX_week <- model.matrix(as.formula(paste("~", week_covars)),
                              data=pop_grouped_week)[,-1]
    popPsi_week <- basis_funcs[pop_grouped_week$AREA, ]

    vb_out <- ulm(Y=Y_week,
                  X=X_week,
                  Psi=Psi_week,
                  weights=scale_weights,
                  n_samples=nsamp,
                  epsilon=1e-3,
                  response_type="ordinal",
                  algorithm="VB")

    vb_pp <- agg_predict(vb_out,
                         predX=popX_week,
                         predPsi=popPsi_week,
                         counts=pop_grouped_week$popsize,
                         alpha=.05,
                         grouping_vars="AREA",
                         pop_df=pop_grouped_week)

    mcmc_out <- ulm(Y=Y_week,
                  X=X_week,
                  Psi=Psi_week,
                  weights=scale_weights,
                  n_samples=nsamp,
                  n_iter=n_mcmc,
                  n_burn=n_burn,
                  response_type="ordinal",
                  algorithm="MCMC")

    mcmc_pp <- agg_predict(mcmc_out,
                           predX=popX_week,
                           predPsi=popPsi_week,
                           counts=pop_grouped_week$popsize,
                           alpha=.05,
                           grouping_vars="AREA",
                           pop_df=pop_grouped_week)

    mcmc_pp$model <- "mcmc-cs"
    vb_pp$model <- "vb-cs"

    mcmc_pp$runtime <- attr(mcmc_out, "runtime")
    vb_pp$runtime <- attr(vb_out, "runtime")

    week_out <- rbind(mcmc_pp, vb_pp) %>%
        tibble::add_column(.before="AREA", WEEK=tt)

    results_df <- rbind(results_df, week_out)
  }
  results_df$sim_num <- i

  save(
    results_df,
    file = file.path(out_dir, paste0(i, ".RData"))
  )

  return(results_df)
}

stopCluster(cl)

save(cs_fit, file=here::here("results", "GAD2" ,"cross_sectional_simulation_results.RData"))
