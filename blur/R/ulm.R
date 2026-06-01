#' Constructor for unit level model
#'
#' @description Fits a Bayesian pseudolikelihood unit-level model for the given response type using the specified sampling algorithm
#' @param response_type data model, one of "binary" or "ordinal"
#' @param algorithm sampling algorithm to employ, either "MCMC" or "VB"
#' @param Y vector or matrix of response values
#' @param X n x p design matrix of covariates
#' @param Psi n x r matrix of basis functions
#' @param weights n-dimensional vector of weights
#' @param n_iter number of iterations
#' @param n_burn number of iterations to discard for burn-in
#' @param init_vals list of initial values
#' @param hyperparams list of hyperparams
#' @param longitudinal logical value for whether to fit longitudinal model
#' @param epsilon tolerance for VB (VB only)
#' @param n_samples number of samples (VB only)
#' @param prev_covar response values from the previous timepoint (long format)
#' @param timepoints list of timepoints (integer valued with first time point equal to 1)
#' @param fix_sigma optional numeric value to fix sigma to in AL sampler (can be null)
#' @returns a list of matrices of MCMC chains for each parameter
#' @export

ulm <- function(Y,
                X,
                Psi,
                weights = NULL,
                timepoints = NULL,
                n_iter=2000,
                n_burn=500,
                init_vals=NULL,
                hyperparams=NULL,
                longitudinal=FALSE,
                epsilon=NULL,
                n_samples=NULL,
                fix_sigma=NULL,
                prev_covar=NULL,
                response_type="ordinal",
                algorithm=c("MCMC", "VB")) {

    if(!is.null(hyperparams) & !is.list(hyperparams)) {
        stop("hyperparams must be a list argument")
    }
    if(!is.null(init_vals) & !is.list(init_vals)) {
        stop("init_vals must be a list argument")
    }
                     
    response_type <- match.arg(response_type,
                               c( "binary", "ordinal"))                                 
    
    algorithm <- match.arg(algorithm, c("MCMC", "VB"))

    if(algorithm == "MCMC") {
      args <- list(Y=Y, X=X, Psi=Psi, weights=weights, n_iter=n_iter, 
                   n_burn=n_burn, init_vals=init_vals, hyperparams=hyperparams)
    }
    
    if(algorithm == "VB") {
      args <- list(Y=Y, X=X, Psi=Psi, weights=weights, epsilon=epsilon, 
                   n_samples=n_samples, init_vals=init_vals, hyperparams=hyperparams)
    }

    if(algorithm == "MCMC" & (n_iter <= n_burn)) {
        stop("n_iter must be greater than n_burn")        
    }

    if(algorithm == "VB" & (is.null(epsilon) | is.null(n_samples))){
        stop("Running VB requires specifying a tolerance `epsilon`
              and a number of samples to draw, `n_samples`")
    }

    if(longitudinal) {
      args$prev_covar <- prev_covar
      if(is.null(timepoints)) {
        stop("timepoints must be provided for longitudinal model")
      } else {
        args$timepoints <- timepoints
      }
    }
    
    function_name <- paste("fit",
                            response_type,
                            ifelse(longitudinal, "lon", "cs"),
                            tolower(algorithm),
                            sep="_")

    new_ulm(args, function_name, longitudinal, response_type, algorithm)
}

new_ulm <- function(args, function_name, longitudinal, response_type, algorithm){

    starttime <- Sys.time()
    out <- do.call(function_name, args)

    #set attributes
    class(out) <- "ulm"
    attr(out, "response_type") <- response_type
    attr(out, "algorithm") <- algorithm
    attr(out, "longitudinal") <- longitudinal
    attr(out, "coeff_names") <- c(colnames(args$X), colnames(args$Psi))
    if(is.null(args$prev_coovar)){
      attr(out, "has_prev_covar") <- TRUE
    } else { attr(out, "has_prev_covar") <- FALSE }

    attr(out, "runtime") <- Sys.time() - starttime

    out
}

validate_ulm <- function(object) {
  if(attr(object, "algorithm") == "VB"){} #need to make sure arguments include epsilon and n_samples
  if(attr(object, "algorithm") == "MCMC"){} #need to make sure arguments include n_burn and n_iter
}

#' Summarize unit-level model fit
#' @param include_reff indicate whether random effects should be included in summary output
#' @param object ulm model output object
#' @param ... other params
#' @method summary ulm
#' @export
#' 
summary.ulm <- function(object, ..., include_reff = FALSE) {
  out <- cbind(mean=round(coef(object, include_reff=include_reff), 3),
               round(confint(object, include_reff=include_reff), 3))

  out
}

#' Print a ulm object
#' @param x ulm model output object
#' @param ... other params
#' @method print ulm
#' @export
print.ulm <- function(x, ...) {
  class <- "ulm"
  summary.ulm(x)
}


#' @method confint ulm
#' @export
#' 
confint.ulm <- function(object, ..., include_reff=FALSE) {
  longitudinal <- attr(object, "longitudinal")
  p <- ncol(object$posteriors$beta)

  if(longitudinal){ #only print fixed effects for longitudinal for now
    beta_eta <- object$posteriors$beta 
  } else {
    beta_eta <- cbind(object$posteriors$beta, object$posteriors$eta)
#    print(beta_eta)
  }
  
  lower_CI <- apply(beta_eta, 2, q025)
  upper_CI <- apply(beta_eta, 2, q975)
  out <- cbind(`2.5%`=lower_CI, `97.5%`=upper_CI)

  if(longitudinal){
    rownames(out) <- attr(object, "coeff_names")[1:p]
  } else {
    rownames(out) <- attr(object, "coeff_names")
  }
  
  if(!include_reff){ out <- out[1:p, ] }
  
  out
}

#' Get fixed effect coefficients (posterior means)
#'
#' @param object ulm model output object
#' @param ... other params
#' @param include_reff logical value indicating whether to print random effect coefficients
#' @method coef ulm
#' @export

coef.ulm <- function(object, ..., include_reff=FALSE) {
  longitudinal <- attr(object, "longitudinal")
  p <- ncol(object$posteriors$beta)

  out <- c(colMeans(object$posteriors$beta), colMeans(object$posteriors$eta))
  names(out) <- attr(object, "coeff_names")
  if(!include_reff){ out <- out[1:p] }

  out
}

## fitted.ulm <- function(object, ...) {}

#' Make model predictions
#' 
#' @description Generates predictions for the population given a unit-level model
#' @param object ulm model output object
#' @param ... other params
#' @param predX X matrix for the population
#' @param predPsi X matrix for the population
#' @param predTimes list of timepoints
#' @param counts number of units in each cell
#' @param K number of categories for categorical models
#' @method predict ulm
#' @export

predict.ulm <- function(object, ..., predX, predPsi, predTimes=NULL, counts=NULL, K=NULL) {
    if(nrow(predX) != nrow(predPsi)) {
      stop("predX and predPsi must have the same number of rows")
    }
    if (!is.null(counts) & (nrow(predX) != length(counts))) {
      stop("length(counts) must match nrow(predX")
    }
    if(is.null(counts) & attr(object, "response_type") %in% c("binary",                                
                                                              "ordinal")) {
      counts <- rep(1, nrow(predX))
        
    }
    if(is.null(predTimes) & attr(object, "longitudinal")) {
      stop("Predictions for longitudinal model require predTimes argument")
    }    
    function_name <- paste("predict",
                            attr(object, "response_type"),
                            ifelse(attr(object, "longitudinal"), "lon", "cs"),
                            sep="_")
   
    args <- list(object=object,
                 predX=predX,
                 predPsi=predPsi,
                 counts=counts)
    if(attr(object, "longitudinal")) { args$predTimes <- predTimes }

    out <- do.call(function_name, args)
    
    out
}

#' Make model predictions aggregated by cell
#' 
#' @description Generates predictions for the population given a unit-level model
#' @param object ulm model output object
#' @param ... other params
#' @param predX X matrix for the population
#' @param predPsi X matrix for the population
#' @param predTimes list of timepoints
#' @param prev_covar vector of previous responses, if known
#' @param counts number of units in each cell
#' @param alpha confidence level for CIs
#' @param group_ids numeric ID to group and aggregate by 
#' @param K number of categories for categorical models
#' @export

agg_predict <- function(object,
                        ...,
                        predX,
                        predPsi,
                        predTimes=NULL,
                        prev_covar=NULL,
                        counts=NULL,
                        grouping_vars,
                        alpha=.05,
                        pop_df,
                        K=NULL) {
    
    if(!all.equal(nrow(predX), nrow(predPsi), nrow(pop_df))) {
        stop("predX, predPsi, and group_ids dimensions must match")
    }

    pop_df <- pop_df %>%
        group_by(across(all_of(grouping_vars))) %>%
        mutate(group_ids=cur_group_id())

    group_ids <- pop_df$group_ids

    if(attr(object, "response_type") == "ordinal") {
        if(attr(object, "longitudinal")) {
		if(is.null(prev_covar)) {
		   return(agg_predict_ordinal_lon_noprev(object,
                                           grouping_vars=grouping_vars,
                                           predX=predX,
                                           predPsi=predPsi,
                                           predTimes=predTimes,
                                           alpha=alpha,
                                           counts=counts,
                                           pop_df=pop_df,
                                           K=K)) 
		} else {
	           return(agg_predict_ordinal_lon(object,
                                           grouping_vars=grouping_vars,
                                           predX=predX,
                                           predPsi=predPsi,
                                           predTimes=predTimes,
                                           prev_covar=prev_covar,
                                           alpha=alpha,
                                           counts=counts,
                                           pop_df=pop_df,
                                           K=K)) 
			}
        } else {
            return(agg_predict_ordinal_cs(object,
                                          predX=predX,
                                          grouping_vars=grouping_vars,
                                          predPsi=predPsi, alpha=alpha,
                                          counts=counts, pop_df=pop_df))
        }
    }

    raw_preds <- predict.ulm(object, predX=predX, predPsi=predPsi,
                             predTimes=predTimes, counts=counts)

    if(attr(object, "response_type") == "binary") {
        agg_pred <- agg_df(group_ids=group_ids, mcmc_mat=raw_preds, alpha=alpha,
                           pop_df=pop_df, grouping_vars=grouping_vars)
    }

    return(agg_pred)
}
