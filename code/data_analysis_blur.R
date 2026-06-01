integer_breaks <- function(x)
  seq(floor(min(x)), ceiling(max(x)))

library(ggplot2)
library(ggh4x)
library(dplyr)
library(tidyr)
library(Matrix)
library(extraDistr)
library(BayesLogit)
library(survey)
library(ggrepel)
library(readr)
library(tibble)
library(ggpattern)
library(viridis)
source("code/helper_functions.R")
devtools::load_all("blur/")

output_dir <- file.path("output", "GAD2", "data_analysis")
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive=TRUE)
}

##########################
### Simulation Settings ##
##########################
set.seed(555)

load(file.path("data", "ordinal", "HPS_empirical_pop_df_GAD2.RData"))
load("data/unscaled_basis_functions.RData")

HPS_df_long <- filter(HPS_df_long, !AREA %in% c("Alaska", "Hawaii")) %>%
                 mutate(AREA=droplevels(AREA)) %>%
                 arrange(PREV_RESPONSE) %>%
                 group_by(WEEK) %>%
                 mutate(SCALE_WEIGHT = n()*ORIG_WEIGHT/sum(ORIG_WEIGHT)) %>%
                 ungroup()

covars <- select(HPS_df_long, matches("COVAR.\\d")) %>%
                 names %>%
                 paste(collapse="+")

n_areas <- nlevels(HPS_df_long$AREA)
n_weeks <- max(HPS_df_long$WEEK)

direst_file <- file.path(output_dir, "data_analysis_direst.RData")
if(file.exists(direst_file)) {
    load(direst_file)
} else {
  files <- sprintf(
    "data/HPS_PUFs/HPS_Week%02d_PUF_CSV/pulse2020_repwgt_puf_%02d.csv",
    1:n_weeks,
    1:n_weeks)
  repwgt_df <- do.call(rbind, lapply(files, read_csv))
  svy_df <- HPS_df_long %>% left_join(repwgt_df, by=c("SCRAM"="SCRAM", "WEEK"="WEEK"))
  hps_svy <- svrepdesign(data = svy_df,
                         weights = ~ORIG_WEIGHT,
                         type="Fay",
                         rho=0.5,
                         repweights = "PWEIGHT[1-9]+")

  direst_props_by_area <- svyby(~RESPONSE, by=~AREA + WEEK,# + COVAR.1 + COVAR.2 + COVAR.3,
                            design=hps_svy, FUN=svymean)
  rownames(direst_props_by_area)  <- c()

  direst_props_finegrained <- c()
  for(week in 1:n_weeks){
      svy_week <- svrepdesign(data = filter(svy_df, WEEK==week),
                              weights = ~ORIG_WEIGHT,
                              repweights = "PWEIGHT[1-9]+")

      direst_week <- svyby(~RESPONSE, by=~AREA + COVAR.1 + COVAR.2 + COVAR.3,
                           design=svy_week, FUN=svymean)

      direst_week <- direst_week %>%
          add_column(WEEK=week, .after="AREA")
      rownames(direst_week) <- c()
      direst_props_finegrained <- rbind(direst_props_finegrained, direst_week)
  }

  direst_finegrained_long <- direst_props_finegrained %>%
     pivot_longer(
       cols = c(starts_with("RESPONSE"), starts_with("se")),
       names_to = c(".value", "CATEGORY"),
       names_pattern = "(RESPONSE|se)(\\d+)"
     ) %>%
      rename(point_est = RESPONSE)

  save(direst_props_by_area, direst_props_finegrained, direst_finegrained_long,
     file=direst_file)
}

n_iter <- 1500
n_burn <- 500
nsamp <- n_iter - n_burn

X <- model.matrix(as.formula(paste("~ -1 + ", covars)), data=HPS_df_long)[,-1]
Psi <- basis_funcs[HPS_df_long$AREA, ]
Y <- model.matrix( ~ 0 + RESPONSE,data=HPS_df_long)
K <- ncol(Y)

census.raw <- readr::read_csv(file.path("data", "census_tables.csv"))

age_breaks <- c(17, seq(25,65,5),100)

census.raw <- census.raw %>%
    select(STATE, NAME, SEX, ORIGIN, RACE, AGE,
           POPESTIMATE2020) %>%
    filter(NAME %in% levels(HPS_df_long$AREA)) %>%
    filter(SEX!=0, AGE >= 18, ORIGIN!=0) %>%
    mutate(SEX=recode(as.factor(SEX),
                      "1" = "MALE", "2" = "FEMALE")) %>%
    mutate(RACE=recode(as.factor(RACE),
                       "1" = "White", "2" = "Black",
                       "4" = "Asian", "3" = "Other",
                       "5" = "Other", "6" = "Other")) %>%
    mutate(AGE_CAT=cut(AGE, breaks=age_breaks)) %>%
    group_by(NAME, SEX, AGE_CAT, RACE) %>%
    summarize(COUNT=sum(POPESTIMATE2020)) %>%
    mutate(AREA=as.factor(NAME)) %>%
    ungroup() %>%
    select(AREA, SEX, AGE_CAT, RACE, COUNT) 

census.raw <- census.raw %>% mutate(CELL.ID = row_number(), .before=AREA)
cell_pops <- census.raw$COUNT #keep a count of how many people are in each cell to start
n_cells <- nrow(census.raw)

prev_resp_fact <- c("PREV_RESP_NA", paste0("PREV_RESP_", 1:K))

census.rep <- c() #for no prev covar we need to repeat all weeks
for(week in 1:n_weeks) {
    census.rep <- rbind(census.rep,
                        cbind(census.raw, WEEK=week))
}

census.table <- expand_grid(census.raw,
                            COVAR.NEW=factor(prev_resp_fact, levels=prev_resp_fact)) %>%
                 arrange(AREA, SEX, AGE_CAT, RACE) %>%
                 arrange(COVAR.NEW)

popPsi <- Matrix(basis_funcs[census.table$AREA,])
names(popPsi) <- c()
n_bf <- ncol(popPsi)

popX <- Matrix(model.matrix(~SEX+AGE_CAT+RACE,data=census.table))[,-1]
p <- ncol(popX)
#need to swap RACEASIAN and RACEOTHER to match betas from model
popX <- popX[,c(1:(p-2), p, (p-1))]

vb_out <-  ulm(Y=Y,
               X=X,
               Psi=Psi,
               weights=HPS_df_long$SCALE_WEIGHT,
               n_samples=nsamp,
               epsilon=.001,
               response_type="ordinal",
               algorithm="VB",
               timepoints=HPS_df_long$WEEK,
               prev_covar=HPS_df_long$COVAR.NEW,
               longitudinal=TRUE)
save(vb_out, file=file.path(output_dir, "data_analysis_vb_output_newer_test.RData"))

vb_out_no_prev <-  ulm(Y=Y,
                       X=X,
                       Psi=Psi,
                       weights=HPS_df_long$SCALE_WEIGHT,
                       n_samples=nsamp,
                       epsilon=.001,
                       response_type="ordinal",
                       algorithm="VB",
                       timepoints=HPS_df_long$WEEK,
                       longitudinal=TRUE)
save(vb_out_no_prev, file=file.path(output_dir, "data_analysis_vb_output_noprev.RData"))

popPsi_rep <- Matrix(basis_funcs[census.rep$AREA,])
names(popPsi_rep) <- c()
popX_rep <- Matrix(model.matrix(~SEX+AGE_CAT+RACE,data=census.rep))[,-1]

vb_no_prev_pp <- agg_predict(vb_out_no_prev,
                             predX=popX_rep,
                             predPsi=popPsi_rep,
                             predTimes=census.rep$WEEK,
                             counts=census.rep$COUNT,                             
                             alpha=.05,
                             grouping_vars=c("WEEK", "AREA", "SEX", "AGE_CAT", "RACE"), 
                             pop_df=census.rep,
                             K=K)

save(vb_no_prev_pp, file=file.path(output_dir, "data_analysis_vb_noprev_pp.RData"))
#load(file.path(output_dir, "data_analysis_model_output_oldbf.RData"))

## mcmc_out <-  ulm(Y=Y,
##               X=X,
##               Psi=Psi,
##               weights=HPS_df_long$SCALE_WEIGHT,
##               n_samples=nsamp,
##               n_iter=n_iter,
##               n_burn=n_burn,
##               response_type="ordinal",
##               algorithm="MCMC",
##               timepoints=HPS_df_long$WEEK,
##               prev_covar=HPS_df_long$COVAR.NEW,
##               longitudinal=TRUE)
## save(mcmc_out, file=file.path(output_dir, "data_analysis_mcmc_output.RData"))

betas <- t(vb_out$posteriors$beta)
etas <- vb_out$posteriors$eta
gammas <- t(vb_out$posteriors$gamma)

#betas <- t(mcmc_out$posteriors$beta)
#etas <- mcmc_out$posteriors$eta
#gammas <- t(mcmc_out$posteriors$gamma)

##posterior predictions
n_Z <- n_weeks + K*(n_weeks-1) #width of Z before KR product
n_gammas <- (K-1)*n_Z
gamma_id2 <- c()
gamma_id2_tmp <- (2*(K+1)):(3*(K+1)-1) - 2*K
for(kk in 1:(K-1)){
  gamma_id2 <- c(gamma_id2, gamma_id2_tmp + (kk-1)*n_Z)
}

#post preds for time 1 

idx <- which(census.table$COVAR.NEW == "PREV_RESP_NA")
cell_pops <- rep(census.table[idx,]$COUNT, nsamp) 
gamma_tt_id <- seq(1, n_gammas, n_Z)
gamma_mat <- gammas[gamma_tt_id, ] 
probs <- matrix(NA, nrow=n_cells*nsamp, ncol=K-1) 
linpred <- popX[idx,] %*% betas + popPsi[idx,] %*% t(etas[,,1])
for(kk in 1:(K-1)){
  gamma_mat_kk <- matrix(rep(gamma_mat[kk, ], each=n_cells), nrow=n_cells)
  stopifnot(length(gamma_mat_kk) == nrow(probs))
  probs[, kk] <- as.vector(gamma_mat_kk - linpred) 
}
probs <- stick(plogis(probs)) 
preds <- rmnom(n=n_cells*nsamp,
               size=cell_pops,  
               prob=probs)

preds <- tibble(cbind(census.table[idx,],
                        MCMC.iter=rep(1:nsamp, each=n_cells),
                        WEEK=1,
                      PRED=preds)) %>%
         select(-CELL.ID, -COUNT, -COVAR.NEW)

#cell_pops ordering needs to follow the order of probs
cell_pop.df <- preds %>%
                 pivot_longer(cols=starts_with("PRED"),
                   names_to=c(".value", "COVAR.NEW"),
                   values_to="PREV_RESP_{.value}",
                   names_sep="\\.") %>%
    pivot_wider(names_from="MCMC.iter",
                values_from="PRED",
                names_prefix="MCMC.") %>%
    mutate(COVAR.NEW=paste0("PREV_RESP_",COVAR.NEW))
                             
for(tt in 2:n_weeks){
  print(tt)
  
  idx <- which(census.table$COVAR.NEW != "PREV_RESP_NA")
  gamma_tt_id <- gamma_id2 + (K+1)*(tt-2)

  gamma_mat <- gammas[gamma_tt_id,] 
  n_cens <- length(idx)  
  probs <- matrix(NA, nrow=n_cens*nsamp, ncol=K-1) 
  
  linpred <- popX[idx,] %*% betas + popPsi[idx,] %*% t(etas[,,tt])      

  #here we have K options per category depending on prev. response
  for(kk in 1:(K-1)){
    g_kk_id <- 1:(K+1) + (K+1)*(kk-1) #subset of columns to iterate over (each set of
                                      # K+1 columns corresponds to a category) so this is
                                      #just sliding a window of width K+1 along within each
                                      #category
    g_kk <- c()
    for(jj in g_kk_id[-1]){ #always skip 1 b/c PREV_RESP_NA never occurs in the population
                            #each iteration here is its own census.table
        g_jj <- matrix(rep(gamma_mat[jj, ], each=n_cells), nrow=n_cells) #
        g_kk <- rbind(g_kk, g_jj)
    }
    if(length(g_kk) != nrow(probs) || !all.equal(dim(g_kk), dim(linpred))){print("ERROR")}
    probs[,kk] <- as.vector(g_kk - linpred) 
  }
  probs <- stick(plogis(probs)) 

  cell_pops <- census.table[idx,] %>%
      right_join(cell_pop.df, by=c("AREA",
                                   "SEX",
                                   "AGE_CAT",
                                   "RACE",
                                   "COVAR.NEW")) %>%
      select(starts_with("MCMC")) %>%
      as.matrix %>%
      as.vector  

  stopifnot(n_cens*nsamp == length(cell_pops) && n_cens*nsamp == nrow(probs))
  
  preds_t <- rmnom(n=n_cens*nsamp,
                   size=cell_pops,  
                   prob=probs)

  pred.names <- tibble(data.frame(census.table[idx,],
                                  MCMC=as.matrix(g_kk-linpred))) %>%
                   pivot_longer(cols=starts_with("MCMC"), 
                   names_to=c(".value", "MCMC.iter"),
                   names_sep="\\.")  %>%
                mutate(MCMC.iter=as.integer(MCMC.iter)) %>%
                arrange(MCMC.iter, COVAR.NEW, CELL.ID) %>% 
                add_column(WEEK=tt, .after="MCMC.iter") %>%
                select(-MCMC)
  
  preds_t <- tibble(cbind(pred.names,
                          PRED=as.matrix(preds_t)))
 
  #need to summarize preds_t to group by category and marginalize out PREV_RESP's
  preds_t <- preds_t %>% group_by(AREA, SEX, AGE_CAT, RACE, MCMC.iter, WEEK) %>%
                                  summarize(across(starts_with("PRED"), sum)) %>%
                                  ungroup()

  #update cell_pops for next iter, which means flattening preds_t correctly
  cell_pop.df <- preds_t %>%
                 pivot_longer(cols=starts_with("PRED"),
                   names_to=c(".value", "COVAR.NEW"),
                   values_to="PREV_RESP_{.value}",
                   names_sep="\\.") %>%
    pivot_wider(names_from="MCMC.iter",
                values_from="PRED",
                names_prefix="MCMC.") %>%
    mutate(COVAR.NEW=paste0("PREV_RESP_",COVAR.NEW))

  preds <- rbind(preds, preds_t)
}

#get means of proportions
blur_preds <- preds %>%
    mutate(TOTAL=rowSums(select(. ,starts_with("PRED")))) %>%
    mutate(across(starts_with("PRED"), ~ .x/TOTAL)) %>% #proportions
    group_by(AREA, SEX, AGE_CAT, RACE, WEEK) %>%
    summarize(across(starts_with("PRED"), list(prop=mean, stderr=sd))) %>%
    pivot_longer(cols=starts_with("PRED"), 
                 names_to=c(".value", "CATEGORY"),
                 names_sep="\\.") %>%
    separate(CATEGORY, into = c("CATEGORY", "fun"), sep = "_") %>%
    pivot_wider(names_from=fun, values_from=PRED) 

#for GAD2, need first category to be zero
blur_preds <- blur_preds %>%
    mutate(CATEGORY=as.integer(CATEGORY) - 1) 

save(blur_preds, file=file.path(output_dir, "blur_preds.RData"))

##compare direst and model preds 
#by sample size 
tmp_df <- blur_preds %>%
    rename(blur_point_est = prop) %>%
    mutate(CATEGORY=as.character(CATEGORY)) %>%
    left_join(direst_finegrained_long,
              by=c("AREA", "WEEK", "CATEGORY",
                   "SEX"="COVAR.1", "AGE_CAT"="COVAR.2",
                   "RACE"="COVAR.3"))

sample_size_df <-
    HPS_df_long %>%
    group_by(WEEK, AREA, RESPONSE, COVAR.1, COVAR.2, COVAR.3) %>%
    summarize("sample_size"=n()) %>%
    left_join(tmp_df,
              by=c("AREA", "WEEK", "RESPONSE"="CATEGORY",
                   "COVAR.1"="SEX", "COVAR.2"="AGE_CAT",
                   "COVAR.3"="RACE"))

sample_size_df %>%
    mutate(RESPONSE=as.factor(RESPONSE)) %>%
    ggplot() +
    geom_point(aes(x=sample_size, y=blur_point_est/point_est,
                   alpha=.3)) +
    geom_abline(aes(slope=0, intercept=1))+# , color="red") +
    guides(alpha="none") + 
    labs(x="Sample size", y="Ratio of model est. to dir. est.") +
    ylim(c(0, 10)) +
    theme_bw(base_size=30)

ggsave(file.path(output_dir, "point_est_comparison.pdf"))#,

# Figure 3) plot direct estimates by state, faceted by week
plot_df <- direst_finegrained_long %>%
    complete(AREA, WEEK, COVAR.1,
             COVAR.2, COVAR.3, CATEGORY) %>%
    filter(COVAR.1=="MALE",
           COVAR.2=="(35,40]",
           COVAR.3=="Asian",
           CATEGORY=="1") %>%
    rename(prop=point_est)

make_areal_plot_DE(plot_df=plot_df,
                   filepath=file.path(output_dir, "DE_areal_estimates_over_time.pdf"))

# Figure 4) plot ratio of model standard errors to DE standard errors
direst_finegrained_long %>%
    rename(direst_se=se) %>%
    mutate(CATEGORY=as.numeric(CATEGORY)) %>%
    right_join(rename(vb_no_prev_pp, model_se=se),
               by=c("AREA", "WEEK", "CATEGORY",
                    "COVAR.1"="SEX",
                    "COVAR.2"="AGE_CAT",
                    "COVAR.3"="RACE")) %>%
    filter(direst_se > 0, !is.na(direst_se)) %>%
    ggplot() +
    geom_point(aes(x=model_se, y=direst_se), size=2, alpha=1) +
    geom_abline(aes(intercept=0, slope=1)) +
    theme_bw(base_size=35) +
    theme(aspect.ratio=.5) +
    theme(strip.text = element_text(size = 20)) +
    xlab("Model estimate standard error") +
    ylab("Direct estimate standard error")

ggsave(file.path(output_dir, "SE_compare.pdf"))

# Figure 5) same plot as 4 for VB estimates
plot_df <- filter(blur_preds, SEX=="MALE", AGE_CAT=="(35,40]", RACE=="Asian", CATEGORY==1)
make_areal_plot_model(plot_df=plot_df,
                      filepath=file.path(output_dir, "TEST_VB_areal_estimates_over_time.pdf"))

# Figure 6) plot same subdemographic, but trajectories for each category by census division
blur_preds$abbr <-  state.abb[match(blur_preds$AREA, state.name)]
blur_preds$DIVISION <- state.division[match(blur_preds$AREA, state.name)]
blur_preds$WEEK <- as.integer(blur_preds$WEEK)

blur_preds <- blur_preds %>% mutate(REGION = case_when(
                                        DIVISION %in% c("New England",
                                                        "Middle Atlantic") ~ "Northeast",
                                        DIVISION %in% c("East North Central",
                                                        "West North Central") ~ "Midwest",
                                      DIVISION %in% c("South Atlantic", "East South Central",
                                                        "West South Central") ~ "South",
                                        DIVISION %in% c("Mountain", "Pacific") ~ "West"
                                    ))
blur_preds$REGION <- as.factor(blur_preds$REGION)

plot_df <- filter(blur_preds,
                  AREA != "District of Columbia",
                  SEX=="MALE",
                  AGE_CAT=="(35,40]",
                  RACE=="Asian")
                  
plot_df  <- plot_df %>%
    group_by(DIVISION) %>%
    arrange(AREA) %>%
    mutate(color_dummy=dense_rank(AREA))

ggplot(data=plot_df, mapping = aes(x=WEEK, y=prop)) +
    geom_line(mapping=aes(group=AREA, color=factor(color_dummy)), linewidth=0.8) +
    geom_point(aes(color=factor(color_dummy)), size=.3) +
    geom_text_repel(data = subset(plot_df, WEEK==max(WEEK)),
                    mapping = aes(x = WEEK, y = prop,
                                  color=factor(color_dummy), label = abbr),
                    size = 9,
                    segment.color = "grey50",
                    nudge_x = 2,              
                    direction = "y",
                    hjust = 0,
                    min.segment.length = 0,
                    max.overlaps = Inf) +
    ggthemes::scale_color_colorblind() +
    coord_cartesian(
        xlim = c(min(plot_df$WEEK), max(plot_df$WEEK) + 3), 
        clip = "off"
    ) +
    facet_grid(CATEGORY~DIVISION, scales="free") +
    labs(x="Week", y="Proportion") +
    guides(color="none") +
    theme_bw(base_size=32) +
    theme(aspect.ratio=.5,
          plot.margin = margin(t = 10, r = 20, b = 10, l = 10))

ggsave(file.path(output_dir, "trajectories_by_census_division.pdf"),
       width=48,
       height=32,
       limitsize=FALSE)

# Figure 7) Compare Male vs. Female all four response categories for just New England
plot_df <- filter(blur_preds,
                  DIVISION == "New England",
                  AGE_CAT=="(35,40]",
                  RACE=="Asian")

plot_df$grouping_var <- paste0(plot_df$AREA, plot_df$SEX)

design <- "
  AABB
  CCDD
  EEFF
  #GG#
"

ggplot(data=plot_df, mapping = aes(x=WEEK, y=prop)) +
    geom_line(mapping=aes(group=grouping_var, color=AREA, linetype=SEX)) +
    geom_point(aes(color=AREA, shape=SEX), size=2.5) +
    scale_x_continuous(breaks=integer_breaks) +
    geom_text_repel(data = subset(plot_df, WEEK==max(WEEK)),
                    mapping = aes(x = WEEK, y = prop,
                                  color=AREA, label = abbr),
                    size = 9,
                    segment.color = "grey50",
                    nudge_x = 2,
                    direction = "y",
                    hjust = 0,
                    max.overlaps = Inf) +
    ggthemes::scale_color_colorblind() +
    coord_cartesian(xlim = c(min(plot_df$WEEK), max(plot_df$WEEK)), clip = "off") +
    labs(x="Week", y="Proportion", linetype="Sex", shape="Sex") +
    guides(color = "none", 
           shape = guide_legend(override.aes = list(size = 5)), 
           linetype = guide_legend()) +
    theme_bw(base_size=50) +
    theme(aspect.ratio=.5,
          text=element_text(face="bold"),
          legend.key.width = unit(2.5, "cm"),
          plot.margin = margin(t = 10, r = 200, b = 10, l = 10)) +
    facet_manual(~CATEGORY, design = design, scales="free_y")

ggsave(file.path(output_dir, "sex_comparison_new_england.pdf"),
       width=48,
       height=36,
       limitsize=FALSE)


#domain sample size summaries for Section 5.3
domains <- expand_grid(
  WEEK     = unique(HPS_df_long$WEEK),
  AREA     = unique(HPS_df_long$AREA),
  COVAR.1  = unique(HPS_df_long$COVAR.1),
  COVAR.2  = unique(HPS_df_long$COVAR.2),
  COVAR.3  = unique(HPS_df_long$COVAR.3),
  RESPONSE = unique(HPS_df_long$RESPONSE)
)

domain_counts <- HPS_df_long %>%
  count(WEEK, AREA, COVAR.1, COVAR.2, COVAR.3, RESPONSE, name = "n") %>%
  right_join(
    domains,
    by = c("WEEK", "AREA", "COVAR.1", "COVAR.2", "COVAR.3", "RESPONSE")
  ) %>%
  mutate(n = if_else(is.na(n), 0L, n))

#count domains where SE couldn't be produced
sum(is.na(direst_finegrained_long$se) |
    near(direst_finegrained_long$se, 0, tol=1e-8))/nrow(domain_counts)

se_compare_df <-  direst_finegrained_long %>%
    rename(direst_se=se) %>%
    mutate(CATEGORY=as.numeric(CATEGORY)) %>%
    right_join(rename(blur_preds, model_se=stderr),
               by=c("AREA", "WEEK", "CATEGORY",
                    "COVAR.1"="SEX",
                    "COVAR.2"="AGE_CAT",
                    "COVAR.3"="RACE"))  %>%
    filter(direst_se > 0, !is.na(direst_se)) 

#count percent of domains where model SE was not lower than DE SE
sum(se_compare_df$direst_se <= se_compare_df$model_se)/nrow(se_compare_df)

vb_no_prev_pp <- mutate(vb_no_prev_pp, CATEGORY=CATEGORY-1)

prev_no_prev_joined_df <- blur_preds %>%
    rename(prev_se=stderr,
           prev_est=prop) %>%
    right_join(rename(vb_no_prev_pp,
                      noprev_se=se,
                      noprev_est=point_est),
               by=c("AREA",
                    "WEEK",
                    "CATEGORY",
                    "SEX",
                    "AGE_CAT",
                    "RACE"))

# compare model std err with full model and noprev
prev_no_prev_joined_df %>%
    ggplot() +
    geom_point(aes(x=prev_se, y=noprev_se), size=2, alpha=.6) +
    geom_abline(aes(intercept=0, slope=1)) +
    theme_bw(base_size=35) +
    theme(aspect.ratio=.5) +
    theme(strip.text = element_text(size = 20)) +
    xlab("Full model estimate standard error") +
    ylab("No previous response model estimate standard error")+
    facet_wrap(~CATEGORY)

ggsave(file.path(output_dir, "prev_noprev_se_comparison.pdf"))

# compare model std err with full model and noprev
prev_no_prev_joined_df %>%
    ggplot() +
    geom_point(aes(x=prev_est, y=noprev_est),
               size=2, alpha=.6) +
    geom_abline(aes(intercept=0, slope=1)) +
    theme_bw(base_size=35) +
    theme(aspect.ratio=.5) +
    theme(strip.text = element_text(size = 20)) +
    xlab("Full model point estimate") +
    ylab("No previous response model point estimate") +
    facet_wrap(~CATEGORY)

ggsave(file.path(output_dir, "prev_noprev_point_est_comparison.pdf"))
