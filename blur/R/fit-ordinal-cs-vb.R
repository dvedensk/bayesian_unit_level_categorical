#' This function fits the ordinal cross-sectional model via variational Bayes
#' @import Matrix 
#' @param Y is a vector of binary survey responses
#' @param X is the input covariate matrix, one row per sample unit
#' @param Psi is the matrix of spatial basis functions, one row per sample unit
#' @param weights is the vector of survey weights
#' @param epsilon is the convergence tolerance
#' @param n_samples number of samples to draw from variational posterior
#' @param init_vals Initial values
#' @param hyperparams Fixed hyperparameters for fixed effects variance and inverse gamma priors

fit_ordinal_cs_vb <- function(Y,
                             X,
                             Psi,
                             epsilon=0.001,
                             n_samples,
                             weights=NULL,
                             init_vals=NULL,
                             hyperparams=NULL) {

  n <- nrow(Y)
  K <- ncol(Y)
  p <- ncol(X)
  r <- ncol(Psi)
  
  if(is.null(weights)){weights <- rep(1, n)}

  categories <- apply(Y, 1, which.max)

  default_init_vals <- list(beta = rep(0, p),
                            eta = rep(0, r),
                            sigma2_eta = 1)
    
  init_vals <- validate_list_arg(default_init_vals, init_vals)
    
  default_hyperparams <- list(sigma2_beta = 3,
                              sigma2_gamma = 3,
                              sigma2_eta = 1,
                              a=.1,
                              b=.1)
  
  hyperparams <- validate_list_arg(default_hyperparams, hyperparams)
  
  keep <- 1:n
  G_ind <- matrix(0, nrow=n, ncol=(K-1))
  G_ind[,1] <- 1
  kappa <- as.matrix(weights * (Y[,-K] - 0.5)) #column for each cat
  kappa_block <- kappa[ , 1]
  for(k in 2:(K-1)){
    #indices for stacked matrices
    id_k <- which(categories >= k)
    keep <- c(keep, id_k) 

    #matrix for picking out the indices of gamma
    tmp <- matrix(0, nrow=length(id_k), ncol=(K-1))
    tmp[,k] <- 1
    G_ind <- rbind(G_ind, tmp)
 
    kappa_block <- c(kappa_block, kappa[id_k, k])    
  }
  kappa <- kappa_block
  X <- X[keep, ]
  Psi <- Psi[keep, ]
  weights <- weights[keep]
  n <- length(weights)

  D <- cbind(-X, -Psi, G_ind)
  Bsigma2_u <- 1
  n_param <- p + r + K - 1
  id_omit <- c(1:p, (p+r+1):n_param)
  muZeta <- rep(1, n_param)
  sigmaZeta <- Diagonal(n_param)
  checkOld <- Inf
  iter <- 1
  repeat{
    ## Latent Variables
    xi <- as.numeric(sqrt(colSums(t(D)*(sigmaZeta%*%t(D))) + (D%*%muZeta)^2))
    
    ## Regression coefficients
    Omega <- Diagonal(,(weights*0.5/xi)*tanh(0.5*xi))
    precZeta <- bdiag((1/hyperparams$sigma2_beta)*Diagonal(p), 
                         as.numeric((hyperparams$a+r/2)/(Bsigma2_u))*Diagonal(r),
                         (1/hyperparams$sigma2_gamma)*Diagonal(K-1)) +
                       t(D)%*%Omega%*%D
    sigmaZeta <- solve(precZeta)
    muZeta <- sigmaZeta %*% t(D) %*% kappa_block
    
    ## RE Variance
    Bsigma2_u <- hyperparams$b + 0.5*(t(muZeta[-id_omit])%*%(muZeta[-id_omit]) +
                              sum(diag(sigmaZeta[-id_omit,-id_omit])))
    
    PPrec <- bdiag((1/hyperparams$sigma2_beta)*Diagonal(p), 
                   as.numeric((hyperparams$a+r/2)/(Bsigma2_u))*Diagonal(r),
                   (1/hyperparams$sigma2_gamma)*Diagonal(K-1))
    
    ## Check for convergence
    checkNew  <- 0.5*(n_param) +
         0.5*determinant(sigmaZeta, logarithm = T)$modulus
        
    if (as.logical(abs(checkOld - checkNew) < epsilon)) break
    checkOld <- checkNew
    iter<-iter+1
  }

  variational_samples <- rmvnorm_prec(n=n_samples,
                                      mean=muZeta,
                                      prec=precZeta)

  betas <- t(variational_samples[1:p, ])
  etas <- t(variational_samples[c((p+1):(p+r)), ])
  gammas <- t(variational_samples[c((p+r+1):n_param), ])   
  sigma2_etas <- rinvgamma(n=n_samples, hyperparams$a+r/2, Bsigma2_u)

  return(list(posteriors=list(beta=as.matrix(betas),
                              eta=as.matrix(etas),
                              gamma=as.matrix(gammas),
                              sigma2_eta=sigma2_etas),
              hyperparams=hyperparams))
}
