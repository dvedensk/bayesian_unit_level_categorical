#' This function fits the binary cross-sectional model via variational Bayes
#' @import Matrix 
#' @import mvtnorm
#' @param Y is a vector of binary survey responses
#' @param X is the input covariate matrix, one row per sample unit
#' @param Psi is the matrix of spatial basis functions, one row per sample unit
#' @param weights is the vector of survey weights
#' @param epsilon is the convergence tolerance
#' @param n_samples number of samples to draw from variational posterior
#' @param init_vals Initial values
#' @param hyperparams Fixed hyperparameters for fixed effects variance and inverse gamma priors

fit_binary_cs_vb <- function(Y,
                             X,
                             Psi,
                             epsilon=0.001,
                             n_samples,
                             weights=NULL,
                             init_vals,
                             hyperparams){

  n <- length(Y)
  p <- ncol(X)
  r <- ncol(Psi)
  C <- cbind(X, Psi)
  I_p <- Diagonal(p)
  I_r <- Diagonal(r)

  default_init_vals <- list(beta = rep(0, p),
                            eta = rep(0, r),
                            sigma2_eta = 1)
    
  init_vals <- validate_list_arg(default_init_vals, init_vals)
    
  default_hyperparams <- list(sigma2_beta = 3,
                              sigma2_eta = 1,
                              a=.1,
                              b=.1)
  
  hyperparams <- validate_list_arg(default_hyperparams, hyperparams)
  
  Bsigma2_u <- 1
  muZeta <- c(init_vals$beta, init_vals$eta) #Zeta is beta and eta blocked together
  sigmaZeta <- Diagonal(p+r)
  Binv <- (1/hyperparams$sigma2_beta) * I_p

  checkOld <- Inf
  checkVec <- c()
  iter <- 1

  repeat{
    ## Latent Variables
    xi <- as.numeric(sqrt(colSums(t(C)*(sigmaZeta%*%t(C))) + (C%*%muZeta)^2))
    
    ## Regression coefficients
    Zbar <- Diagonal(n, (weights*0.5/xi) * tanh(0.5*xi))
    Rinv <- as.numeric((hyperparams$a+r/2)/(Bsigma2_u))*I_r
    precZeta <- bdiag(Binv, Rinv) + t(C)%*%Zbar%*%C
    sigmaZeta <- solve(precZeta)
    muZeta <-  sigmaZeta %*% t(C) %*% (weights * (Y - 0.5))
           
    ## RE Variance
    Bsigma2_u <- hyperparams$b + 0.5*(t(muZeta[-c(1:p)])%*%(muZeta[-c(1:p)]) +
                                     sum(diag(sigmaZeta[-c(1:p),-c(1:p)])))
    
    PPrec <- bdiag(Binv,
                   as.numeric((hyperparams$a+r/2)/(Bsigma2_u))*I_r)
   
    ## Check for convergence
    checkNew  <- 0.5*(p+r) + 0.5*determinant(sigmaZeta, logarithm = T)$modulus  - 
     0.5*t(muZeta)%*%PPrec%*%(muZeta) + sum(weights*(Y-0.5)*as.numeric(C%*%muZeta) + weights*log(stats::plogis(xi)) - 0.5*weights*xi) - 
      0.5*sum(diag(PPrec %*% sigmaZeta)) - log(Bsigma2_u)
      
    checkVec <- c(checkVec, as.numeric(checkNew))
        
    if (as.logical(abs(checkOld - checkNew) < epsilon)) break
    iter <- iter + 1
    checkOld <- checkNew
  }

  #variational_samples <- rmvnorm_prec(n=n_samples, mean=muZeta, prec=precZeta)
  variational_samples <- rmvnorm_prec(n=n_samples,
                                      mean=muZeta,
                                      prec=precZeta) 
  beta <- t(variational_samples[1:p, ])
  eta <- t(variational_samples[-c(1:p), ])
  
  return(list(posteriors = list(beta=beta, eta=eta),
              hyperparams=hyperparams))
}
