agg_predict_ordinal_lon <- function(object,
                                    predX,
                                    predPsi,
                                    predTimes,
                                    prev_covar,
                                    counts,
                                    grouping_vars,
                                    alpha,
                                    pop_df,
                                    K) {
  pop_size <- nrow(pop_df)
  nsamp <- nrow(object$posteriors$beta)
  n_tp <- dim(object$posteriors$eta)[3]
  A <- model.matrix(~ 0 + as.factor(predTimes))
  B <- predPsi # basis functions for everyone in pop
  predPsiTime <- t(Matrix::KhatriRao(t(A), t(B)))

  popZ <- model.matrix(~ prev_covar:factor(predTimes) -1)
  popZ <- popZ[,-c(2:(K+1))]
  gamma_col_id <- max.col(popZ)
  t_gamma <- t(object$posteriors$gamma)

  linpred <- predX %*% t(object$posteriors$beta) +
      predPsiTime %*% matrix(object$posteriors$eta, ncol=nsamp, byrow=T)

  gamma_reduced <- c()  #stack gamma into a matrix that matches linpred
  for(k in 1:(K-1)){
    width <- (K+1)*n_tp-K 
    start_id <- (k-1)*width+1
    end_id <- start_id + width - 1
    gamma_col_k <- t_gamma[start_id:end_id,][gamma_col_id, 1:nsamp] 
    gamma_reduced <- cbind(gamma_reduced, as.vector(gamma_col_k))
  }
  rm(gamma_col_k)

  preds <- plogis(gamma_reduced - as.vector(linpred))
  rm(gamma_reduced, linpred)
  gc()
      
  preds <- rmnom(n=pop_size * nsamp,
                 size=rep(counts, nsamp), 
                 prob=stick(preds))

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
