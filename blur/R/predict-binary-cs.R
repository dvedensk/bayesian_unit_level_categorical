predict_binary_cs <- function(object, predX, predPsi, counts) {

  beta <- object$posteriors$beta
  eta <- object$posteriors$eta

  n_pred <- nrow(predX)
  n_mcmc <- nrow(beta)
  
  preds <- rbinom(n=(n_pred * n_mcmc),
                  size=rep(counts, n_mcmc), 
                  prob=plogis(as.numeric(predX %*% t(beta) + predPsi %*% t(eta))))
    
  preds <- matrix(preds, nrow=n_pred, ncol=n_mcmc)

}
