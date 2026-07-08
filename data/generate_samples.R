library(dplyr)
library(sampling)
library(tidyr)
source("code/helper_functions.R")

response_type <- "ordinal"
save_dir <- paste0("data/", response_type, "/") 

set.seed(555)
Nsim <- 100

load(paste0(save_dir, "HPS_empirical_pop_df_GAD2.RData"))

HPS_pop_long <- HPS_df_long
HPS_pop_wide <- HPS_df_wide

state_exclude <- c("Alaska", "Hawaii")

HPS_pop_long <- filter(HPS_pop_long, !AREA %in% state_exclude) %>% 
                     mutate(AREA=droplevels(AREA)) %>% 
                     arrange(PREV_RESPONSE)

HPS_pop_wide <- filter(HPS_pop_wide, !AREA %in% state_exclude) %>% 
                     mutate(AREA=droplevels(AREA))

n_weeks <- max(HPS_pop_long$WEEK)
n_areas <- nlevels(HPS_pop_long$AREA)

pop_size <- nrow(HPS_pop_wide)
sample_size <- floor(.05 * pop_size)

samples <- list()
for(i in 1:Nsim){
    set.seed(i + 15)
    print(paste("Taking sample...", i))
    sample_out <- get_sample(HPS_pop_long=HPS_pop_long, HPS_pop_wide=HPS_pop_wide, 
                             sample_size=sample_size, n_weeks=n_weeks, n_areas=n_areas)
    HPS_sample_long <- sample_out$HPS_sample_long
    HPS_sample_long <- HPS_sample_long %>% group_by(WEEK) %>%
                                           mutate(SCALE_WEIGHT=
                                                  PWEIGHT*n()/sum(PWEIGHT)) %>%
                                           ungroup() %>%
                                           arrange(WEEK)
    samples[[i]] <- HPS_sample_long
}

save(samples, HPS_pop_long, file=paste0(save_dir, "empirical_samples_GAD2.RData"))
