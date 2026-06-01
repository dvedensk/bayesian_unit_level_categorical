agg_predict_ordinal_lon_noprev <- function(object,
                                    predX,
                                    predPsi,
                                    predTimes,
                                    counts,
                                    grouping_vars,
                                    alpha,
                                    pop_df,
                                    K) {

  pop_size <- nrow(pop_df)
  nsamp <- nrow(object$posteriors$beta)
  n_tp <- dim(object$posteriors$eta)[3]
  A <- model.matrix(~ 0 + as.factor(predTimes))
  B <- predPsi 
  predPsiTime <- t(Matrix::KhatriRao(t(A), t(B)))

  linpred <- predX %*% t(object$posteriors$beta) +
      predPsiTime %*% matrix(object$posteriors$eta, ncol=nsamp, byrow=T)

#Can verify it matches this loop
#  linpred <- c()
#  for(tt in 1:n_tp) {
#    idx_tt <- which(predTimes == tt)
#    eta_tt <- object$posteriors$eta[, , tt]
#    tmp <- predX[idx_tt, ] %*% t(object$posteriors$beta) + predPsi[idx_tt, ] %*% t(eta_tt)
#    linpred <- rbind(linpred, tmp)
#  }

  probs <-  plogis(as.matrix(object$posteriors$gamma[rep(1:nsamp, each=pop_size), ]) -
                     as.vector((linpred)))

  preds <- rmnom(n=pop_size * nsamp,
                 size=rep(counts, nsamp), 
                 prob=stick(probs))

  preds <- mat_split(preds, pop_size)
  preds <- matrix(aperm(preds, c(3,2,1)), ncol=nsamp, byrow=F)
  
  preds <- preds / rep(counts, each = K) #rescale to proportions

  preds <- data.frame(
    group_ids = rep(pop_df$group_ids, each = K),
    CATEGORY  = rep(1:K, pop_size),
    point_est = matrixStats::rowMeans2(preds),
    se        = matrixStats::rowSds(preds),
    ci_lower  = matrixStats::rowQuantiles(preds, probs = alpha/2),
    ci_upper  = matrixStats::rowQuantiles(preds, probs = 1 - alpha/2)
  )

  preds <- add_grouping_names(preds,
                              pop_df,
                              grouping_vars)
}



## linpred <- predX %*% t(object$posteriors$beta) +
##     predPsiTime %*% matrix(object$posteriors$eta, ncol=nsamp, byrow=T)

## #currently do
## tstmat <- matrix(object$posteriors$eta, ncol=nsamp, byrow=T)
## #this does all areas in time 1 first, then all areas in time 2
## tstmat_llm <- matrix(aperm(object$posteriors$eta, c(1, 3, 2)), ncol=nsamp, byrow=FALSE)
