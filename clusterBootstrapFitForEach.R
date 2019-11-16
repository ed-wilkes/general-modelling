#' @name clusterBootstrapFitParallel
#' @author Ed Wilkes
#' 
#' @description fits new rlmer models to bootstrap resampled data in parallel
#' 
#' @param model lmer or rlmer object containing fitted model
#' @param data Data frame containing original data to be sampled
#' @param x_var String denoting the name of the x variable to model
#' @param id_var Vector of ID variables corresponding to nesting levels
#' @param seed Random seed number
#' @param max_boot Numeric denoting the maximum number of iterations to perform
#' @param type String denoting model type "random" or "mixed"
#' @param data_new Data from which to draw response predictions
#' @param n_cores Numeric denoting number of cores (virtual) to use for parallel
#' processing. Note you may quickly run out of RAM if you use too many cores!
#'
#' @return List containing the following:
#'           (i) Data frame of resampled effect estimates
#'           (ii) List of predicted responses for each bootstrapped model
#'
clusterBootstrapFitForEach <- function(data
                                       ,x_var = NULL
                                       ,id_var
                                       ,seed
                                       ,response
                                       ,max_boot = 999
                                       ,type = "mixed"
                                       ,data_new
                                       ,n_cores) {
  
  ## Packages ----
  require(doParallel)
  require(doRNG)
  require(dplyr)
  require(foreach)
  require(iterators)
  require(lme4)
  require(robustlmm)
  
  ## Formula definition ----
  if (type == "mixed") {
    
    form <- as.formula(paste0(response, " ~1 +", x_var, "+ (1|",id_var,")"))
    
  } else if (type == "random") {
    
    form <- as.formula(paste0(response, "~1 + (1|", id_var, ")"))
    
  }
  
  ## Generate bootstrap samples in list ----
  set.seed(seed)
  list_samples <- list()
  
  for (i in 1:max_boot) {
    list_samples[[i]] <- clusterBootSample(data = data, id_var = id_var)
  }
  
  ## Set up clusters ----
  cores_max <- detectCores()
  if (n_cores > cores_max) {
    stop("Number of cores detected is lower than n_cores!")
  }
  cluster <- makeCluster(n_cores)
  registerDoParallel(cluster)
  
  ## Foreach parallel loop over list_samples ----
  boot_results <- foreach(i = iter(list_samples)
                          ,.combine = rbind
                          ,.packages = c("lme4", "dplyr")
                          ,.options.RNG = seed) %dorng% {
    
    # Fit model to bootstrap sample
    boot_fit <- robustlmm::rlmerRcpp(form, data = i)
    
    # Return coefficient estimates
    if (type == "mixed") {
      
      df_coef <- data.frame(u0 = attr(lme4::VarCorr(boot_fit)[[1]], "stddev")
                            ,resid = sigma(boot_fit)
                            ,intercept = lme4::fixef(boot_fit)[1]
                            ,slope = lme4::fixef(boot_fit)[2])
      
    } else if (type == "random") {
      
      df_coef <- data.frame(u0 = attr(VarCorr(boot_fit)[[1]], "stddev")
                            ,resid = sigma(boot_fit)
                            ,intercept = lme4::fixef(boot_fit)[1])
      
    }
    
    # Sample bootstrapped residuals to add back to predictions
    data_new[[id_var]] <- as.character(data_new[[id_var]])
    data_new_filter <- data_new %>%
      filter(get(id_var) %in% i[[id_var]]) 

    eps_sample <- sample(residuals(boot_fit), size = nrow(data_new_filter), replace = TRUE)
    data_new_filter$y_hat <- predict(boot_fit, data_new_filter)
    data_new_filter$y_star <- data_new_filter$y_hat + eps_sample
    
    return(list(coefficients = df_coef
                ,predictions = data_new_filter))           
    
    # Remove objects from memory
    rm(boot_fit
       ,df_coef
       ,eps_sample
       ,data_new_filter)
    
  }
  
  ## Return objects ----
  return(list(coefficients = dplyr::bind_rows(boot_results[,1])
              ,predictions = dplyr::bind_rows(boot_results[,2])
        )
  )
  
  ## Stop cluster
  stopCluster(cluster)
  
}
