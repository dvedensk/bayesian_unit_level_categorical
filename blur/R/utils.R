q025 <- function(x){quantile(x, .025)}
q975 <- function(x){quantile(x, .975)}

lse <- function(x){
  m <- max(x)
  m + log(sum(exp(x-m)))
}

tr <- function(x) {
  sum(diag(x))
}

softmax <- function(x){
  exp(x - lse(x))
}

sample_MVN <- function(mu, prec){
  U <- base::chol(prec)
  tmp_norm <- rnorm(nrow(prec))
  ret_mat <- backsolve(U, backsolve(U, mu, transpose=TRUE) + tmp_norm)
  return(ret_mat)
}

validate_list_arg <- function(defaults, args=NULL) {
  default_names_str <- paste(names(defaults), collapse=", ")
    
  if(is.null(args)) {
     args <- defaults
  }
  
  if(!all( names(args) %in% names(defaults) )) {
    stop(paste("Priors for this model are ", default_names_str))
  }
  
  for (name in names(defaults)) {
    if (is.null(args[[name]])) {
      args[[name]] <- defaults[[name]]
    }
  }

  return(args)
}

#takes an N x K matrix of probabilities and returns an N x (K+1) matrix
stick <- function(Theta){
    K <- ncol(Theta) + 1
    pred_probs <- matrix(NA, nrow=nrow(Theta), ncol=K)
    pred_probs[, 1] <- Theta[, 1]
    for(kk in 1:(K-2)){
      pred_probs[ , kk + 1] <- Theta[ , kk + 1]/Theta[ , kk] * (1 - Theta[ , kk]) * pred_probs[, kk]
    }
    pred_probs[ , K] <- (1 - Theta[ , (K-1)])/Theta[,(K-1)] * pred_probs[, (K-1)]
    return(pred_probs)
}

#function for reconstructing stick-breaking probabilities output by nominal model
#input is an array of size n_pred x (K-1) x n_mcmc
#output flattens this array into a matrix of dimension  (n_pred*n_mcmc) x K
array_stick <- function(pred_probs, n_pred, n_mcmc, K){
    Theta <- aperm(pred_probs[,,], c(1,3,2))
    dim(Theta) <- c(n_pred*n_mcmc, K - 1)
    pred_probs <- matrix(NA, nrow=n_pred*n_mcmc, ncol=K)
    pred_probs[,1] <- Theta[,1]
    for(kk in 1:(K-2)){
      pred_probs[ , kk + 1] <- Theta[ , kk + 1]/Theta[ , kk] * (1 - Theta[ , kk]) * pred_probs[,kk]
    }
    pred_probs[ , K] <- (1 - Theta[ , (K-1)])/Theta[,(K-1)] * pred_probs[, (K-1)]
    return(pred_probs)
}

agg_df <- function(group_ids, mcmc_mat, alpha, pop_df, grouping_vars) {
    #turn mcmc_mat into proportions
    tmp <- data.frame(group_ids=group_ids, MCMC=mcmc_mat) |>
        group_by(group_ids) |>
        summarize_all(sum) |>
        mutate(across(starts_with("MCMC"), ~ .x/sum(.x))) 

     preds <- cbind(preds[, 1], point_est=rowMeans(preds[, -1]),
             ci_lower=apply(preds[, -1], 1, \(x)(quantile(x, alpha/2))),
             ci_upper=apply(preds[, -1], 1, \(x)(quantile(x, 1-alpha/2))))

     preds <- add_grouping_names(preds,
                                select(pop_df, all_of(c("group_ids", grouping_vars))))

}

# Draw from MVN by specifying a covariance matrix
rmvnorm = function(n, mean, covar)
{
  k <- length(mean)
  stopifnot(k == nrow(covar) && k == ncol(covar))
  Z <- matrix(rnorm(n*k), k, n)
  A <- t(chol(covar))
  out <- A %*% Z + mean

  if(n==1){ return(as.vector(out)) } else { return(out) }
}

# Draw from MVN by specifying a precision matrix
rmvnorm_prec = function(n, mean, prec)
{
  k <- length(mean)
  stopifnot(k == nrow(prec) && k == ncol(prec))
  Z <- matrix(rnorm(n*k), k, n)

  # Note that Ainv %*% t(Ainv) is the Cholesky decomposition of the covariance
  # matrix solve(Omega)
  A <- chol(prec)
  Ainv <- backsolve(A, diag(1,k,k))
  out <- Ainv %*% Z + mean

  if(n==1){ return(as.vector(out)) } else { return(out) }
}


add_grouping_names <- function(pp_df, grouping_df, grouping_vars) {
  unique(select(grouping_df, all_of(grouping_vars), group_ids)) %>%
      right_join(by=c("group_ids"), pp_df) %>%
      select(-group_ids)
}


#takes in an n x m matrix
#returns an array of n/block_size matrices each of dim. block_size x m
mat_split <- function(input.matrix, block.size){
  n.rows <- nrow(input.matrix)
  n.cols <- ncol(input.matrix)
  if(n.rows %% block.size != 0){stop("Not a valid block size")}

  n.blocks <- n.rows/block.size
  output.array <- array(NA, dim=c(n.blocks, block.size, n.cols))
  for(i in 1:n.blocks){
    start_ind <- (i-1)*block.size +1
    end_ind <- i*block.size
    output.array[i,,]  <- input.matrix[start_ind:end_ind,]
  }
  return(output.array)
}

