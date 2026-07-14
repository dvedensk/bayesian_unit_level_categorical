library(dplyr)
library(knitr)
library(scales)
library(ggforce)

data_dir <- here::here("data", "ordinal")
output_dir <- file.path("output", "GAD2", "empirical_sim")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

load(here::here(data_dir, "empirical_samples_GAD2.RData"))

int_score <- function(Q, L, U, alpha){
  return(
    (U - L) + 2/alpha*(Q < L)*(L - Q) + 2/alpha*(Q > U)*(Q - U)
  )
}

rbind_results <- function(path = ".") {
  files <- list.files(
    path,
    pattern = "^\\d+\\.RData$",
    full.names = TRUE
  )

  res_list <- vector("list", length(files))

  for (i in seq_along(files)) {
    e <- new.env()
    load(files[i], envir = e)
    res_list[[i]] <- e$results_df
  }

  dplyr::bind_rows(res_list)
}

cs_results <- rbind_results("results/GAD2/cross_sectional")
lon_results <- rbind_results("results/GAD2/longitudinal")
all_results <- rbind(cs_results, lon_results)

#temporarily replace svy_direst until next re-run of CS file
#load("svy_direst.RData")
#all_results <- all_results %>% filter(model!="svy_direst")
#all_results <- rbind(all_results, results_df)

truth <- HPS_pop_long %>%
           group_by(AREA, WEEK, RESPONSE) %>%
           summarize(Q=n()) %>%
           group_by(AREA, WEEK) %>%
           mutate(P=Q/sum(Q)) %>%
    arrange(WEEK, AREA) %>%
    mutate(CATEGORY=as.integer(RESPONSE)) %>%
    select(-RESPONSE, -Q)

pred_avg <- all_results %>%
  group_by(WEEK, AREA, CATEGORY, model) %>%
  mutate(model=ifelse(model=="mcmc", "mcmc-cs", model)) %>%
    mutate(model=ifelse(model=="vb", "vb-cs", model)) #%>%
#    filter(!(point_est==0 & ci_lower==0 & ci_upper==0))


pred_avg <- ungroup(pred_avg) %>%
    filter(model %in% c("mcmc-cs", "vb-cs")) %>%
    select(model, sim_num, WEEK, runtime) %>%
    unique %>%
    group_by(model, sim_num) %>%
    summarize(cs_runtime = sum(runtime)) %>%
    right_join(pred_avg, by=c("model", "sim_num")) %>%
    mutate(runtime=ifelse(model %in% c("mcmc-cs", "vb-cs"), cs_runtime, runtime)) %>%
    select(-cs_runtime)
 
compare_df <- pred_avg %>%
  inner_join(truth,
             by = c("WEEK", "AREA", "CATEGORY"))

summary_df <- compare_df %>%
  group_by(model) %>% #, WEEK, AREA, CATEGORY) %>%
  summarise(
      mse = mean((point_est - P)^2),
      rel_mse = mean((point_est - P)^2/P^2), 
      coverage = mean(between(P, ci_lower, ci_upper)),
#      abs_bias = mean(abs(point_est - P)),
      rel_bias = mean( (point_est - P)/P ),
      IS = mean(int_score(P, ci_lower, ci_upper, .05)),
      med_runtime = median(runtime, na.rm=TRUE),
    .groups = "drop"
  )


summary_df <- summary_df %>%
  mutate(
    mse_rel_to_DE = mse / mse[model == "svy_direst"]
  ) %>%
    relocate(mse_rel_to_DE, .before = med_runtime)  %>%
    select(-mse)

compare_df %>%
#  filter(model != "svy_direst") %>%
  group_by(model, WEEK) %>%
    summarize(mse = mean((point_est - P)^2),
              rel_mse = mean((point_est - P)^2/P^2), 
             .groups = "drop") %>%
  mutate(
    cs_dummy = factor(ifelse(model %in% c("mcmc-cs", "vb-cs"), 1, 0)),
    model = recode(
      model,
      "mcmc-cs"  = "Gibbs-CS",
      "mcmc-lon" = "Gibbs-Lon",
      "svy_direst" = "Dir. Est.",
      "vb-cs"    = "VB-CS",
      "vb-lon"   = "VB-Lon"
    )
  ) %>%
  ggplot(aes(x = WEEK, y = rel_mse, color = model)) +
  geom_line(aes(linetype = cs_dummy), show.legend = FALSE) +
  geom_point() +
    ggthemes::scale_color_colorblind(name = "Estimator") +
    facet_zoom(ylim=c(0.04,.01)) +
  scale_x_continuous(
    name = "Week",
    breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)
  ) +
    ylab("Relative MSE") + 
#  scale_y_log10(name = "Relative MSE (log scale)") +
  theme_bw(base_size = 30)

ggsave(file.path(output_dir, "log_MSE_time_series_plot.pdf"),
       width=48,
       height=36,
       limitsize=FALSE)

plot_df <- compare_df %>%
  group_by(model, WEEK, AREA, CATEGORY) %>%
    summarise(
#        mse = mean((point_est - P)^2),
        rel_mse = mean((point_est - P)^2/P^2), 
        .groups = "drop") %>%
  tidyr::pivot_wider(names_from = model, values_from = rel_mse) %>%
  tidyr::pivot_longer(
    -c(WEEK, AREA, CATEGORY, svy_direst),
    names_to = "model",
    values_to = "rel_mse_model"
  )# %>%
#  filter(!is.na(svy_direst), !is.na(mse_model))

plot_df %>%
  mutate(
    model = recode(
      model,
      "mcmc-cs"  = "Gibbs-CS",
      "mcmc-lon" = "Gibbs-Lon",
      "svy_direst" = "Dir. Est.",
      "vb-cs"    = "VB-CS",
      "vb-lon"   = "VB-Lon"
    )
  ) %>%
  ggplot(aes(x = rel_mse_model, y = svy_direst)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ model, scales = "free") +
  labs(
    x = "Model Rel. MSE",
    y = "Direct estimate Rel. MSE"
  ) +
  theme_bw(base_size=30)

ggsave(file.path(output_dir, "MSE_ratios.pdf"),
       width=48,
       height=32,
       limitsize=FALSE)

table_df <- summary_df %>%
  mutate(
    Method = recode(
      model,
      "svy_direst" = "Direct",
      "vb-cs"      = "VB-CS",
      "vb-lon"     = "VB-Lon",
      "mcmc-cs"    = "Gibbs-CS",
      "mcmc-lon"   = "Gibbs-Lon"
    ),
#    `Rel. MSE` = scientific(rel_mse, digits = 1),
#    `Rel. Bias` = scientific(rel_bias, digits = 1),
    `Rel. MSE` = round(rel_mse, 3),
    `Rel. Bias` = round(rel_bias, 3),
    `Cov.` = percent(coverage, accuracy = 1),
    IS = formatC(IS, format = "f", digits = 3),
    `Runtime (s)` = ifelse(
        Method == "Direct",
        "-",
        format(round(med_runtime), big.mark = ",")
    )) %>%
    select(
        Method,
        `Rel. MSE`,
        `Rel. Bias`,
        `Cov.`,
        IS,
        `Runtime (s)`
  ) %>%
  arrange(factor(
    Method,
    levels = c("Direct", "VB-CS", "VB-Lon", "Gibbs-CS", "Gibbs-Lon")
  ))

kable(
  table_df,
  align = "lrrrrr",
  format="latex",
  caption = "Performance comparison across estimators"
)
