agg_predict_ordinal_cs <- function(object,
                                   predX,
                                   predPsi,
                                   counts,
                                   pop_df,
                                   alpha,
                                   grouping_vars) {
  #ERROR CHECK that all grouping_Vars exist in pop_df
  pop_size <- nrow(pop_df)
  nsamp <- nrow(object$posteriors$beta)
  K <- ncol(object$posteriors$gamma) + 1

  linpred <- predX %*% t(object$posteriors$beta) +
        predPsi %*% t(object$posteriors$eta)

  probs <-  plogis(object$posteriors$gamma[rep(1:nsamp, each=pop_size),] -
                     as.vector((linpred)))
      
  preds <- rmnom(n=pop_size * nsamp,
                 size=rep(counts, nsamp), 
                 prob=stick(probs))

  #reshape into n_pop by n_iter matrix
  preds <- mat_split(preds, pop_size)
  preds <- matrix(aperm(preds, c(3,2,1)), ncol=nsamp, byrow=F)

  #get total sums for each category in each area
  preds <- data.frame(group_ids = rep(pop_df$group_ids, each=K),
                      CATEGORY = rep(1:K, pop_size),
                      MCMC_iter = preds) %>% 
            group_by(group_ids, CATEGORY) %>%
            summarize_all(sum) 

  #rescale to props.
  preds <- preds %>% mutate(across(starts_with("MCMC_iter"), ~ .x/sum(.x)))

  preds <- cbind(preds[, 1:2], point_est=rowMeans(preds[, -c(1:2)]),
             ci_lower=apply(preds[, -c(1:2)], 1, \(x)(quantile(x, alpha/2))),
             ci_upper=apply(preds[, -c(1:2)], 1, \(x)(quantile(x, 1-alpha/2))))

  preds <- add_grouping_names(preds,
                              select(pop_df, all_of(c("group_ids", grouping_vars))), 
			      grouping_vars)
}
