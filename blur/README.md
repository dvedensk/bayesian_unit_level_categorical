# Installation 

To install the package, run 

```
devtools::install_all("blur")
library(dplyr)
```

# Example with cross-sectional ordinal model

We illustrate the usage with an example from the 2018 Survey of Income and Program Participation (SIPP). 
First we process the provided `sipp` and `census_table_2018` datasets to make sure all the factor levels agree. 
We also bin age into groups for convenience.

```
levels(sipp$health_status) <- rev(levels(sipp$health_status))

 census_df <- census_table_2018 %>%
    filter(sex != "total",
           age > 15) %>%  #min age in SIPP data                                                                                                                                                           mutate(sex = droplevels(sex))

 age_groups <- c(15, 30, 40, 50, 65, Inf)
 census_df <- mutate(census_df, age=cut(age, age_groups))
 sipp <- sipp %>% mutate(age=cut(age, age_groups))
                                                                
 levels(census_df$sex) <- levels(sipp$sex)
 levels(census_df$race) <- levels(sipp$race)
 census_df <- census_df %>%
    mutate(state=as.factor(tolower(state)))
```

We can fit both a binary and an ordinal model, using race, sex, age as covariates and assuming an iid areal random effect

```
 Y_bin <- abs(as.integer(sipp$medicaid)-2)
 Y_ord <- model.matrix(~ -1 + sipp$health_status)
 X <- model.matrix(~ -1 + race + sex + age, data=sipp)
 Psi <- model.matrix(~ -1 + state, data=sipp)
```

After re-scaling the survey weights to sum to the sample size, we can run both models. (Note that the ordinal model requires dropping the first column of the X matrix for identifiability.)

```
 scale_weights <- sipp$weight/sum(sipp$weight) * nrow(sipp)

 bin_mod_mcmc <- ulm(Y=Y_bin,
                    X=X,
                    Psi=Psi,
                    response_type="binary",
                    algorithm="VB",
                    weights=scale_weights,
                    n_samples=1000,
                    epsilon=.001)

 ord_mod_vb <- ulm(Y=Y_ord,
               X=X[,-1],
               Psi=Psi,
               response_type="ordinal",
               algorithm="VB",
               weights=scale_weights,
               n_samples=1000,
               epsilon=.01)
```
Lastly, calculate posterior predictions  

```
 popX <- model.matrix(~ -1 + race + sex + age, data=census_df)
 popPsi <- model.matrix(~ -1 + state, data=census_df)

 ord_pp <- agg_predict(ord_mod_vb,
                      predX=popX[,-1],
                      predPsi=popPsi,
                      alpha=.05,
                      counts=census_df$count,
                      pop_df=census_df,
                      grouping_vars="state")
```
