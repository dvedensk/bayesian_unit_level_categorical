predict_ordinal_cs <- function(object, predX, predPsi, counts) {
  betas <- object$posteriors$betas
  etas <- object$posteriors$etas
  gammas <- object$posteriors$gammas

  K <- ncol(gammas) + 1 
  n_pred <- nrow(predX)
  n_mcmc <- nrow(betas)

  linpred <- as.vector(cbind(predX, predPsi) %*% t(cbind(betas, etas)))
  linpred <- gammas[rep(1:n_mcmc, each=n_pred), ] - linpred #needs to be matrix w/ K columns

  preds <- rmnom(n=(n_pred * n_mcmc),
                 size=counts, #rep(counts, n_mcmc),
                 prob=stick(plogis(linpred)))
  
  preds <- array(preds, dim = c(n_pred, n_mcmc, K))
  
  preds
}
