#' @name bootstrapFitParallel
#' @author Ed Wilkes
#' 
#' @description fits new rlmer models to bootstrap resampled data in parallel. NB: currently
#' is only coded to work with rlmer model objects.
#' 
#' @param model lmer or rlmer object containing fitted model
#' @param data Data frame containing original data to be sampled
#' @param id_var Vector of ID variables corresponding to nesting levels
#' @param x_var String denoting the fixed effect (if required)
#' @param response String denoting the name of the response variable being modelled
#' @param seed Random seed number
#' @param max_boot Numeric denoting the maximum number of iterations to perform
#' @param type String denoting model type "random" or "mixed"
#' @param data_new Data from which to draw response predictions
#' @param n_cores Numeric denoting number of cores (virtual) to use as workers
#'
#' @return List containing the following:
#'           (i) Data frame of resampled effect estimates
#'           (ii) List of predicted responses for each bootstrapped model
#'
bootstrapFitForEach <- function(model
                                ,data
                                ,id_var
                                ,x_var = NA
                                ,response
                                ,seed
                                ,max_boot
                                ,type
                                ,data_new
                                ,n_cores) {
  
  # Packages
  require(doParallel)
  require(dplyr)
  require(foreach)
  require(iterators)
  require(lme4)
  require(robustlmm)
  
  # Mixed models
  if (type == "mixed") {
    
    form <- as.formula(paste0(response, "~1 +", x_var, " + (1|", id_var, ")"))
    
  } else if (type == "random") {
    
    form <- as.formula(paste0(response, "~1 + (1|", id_var, ")"))
    
  }
  
  # Generate bootstrap samples in list
  set.seed(seed)
  list_samples <- list()
  for (i in 1:max_boot) {
    list_samples[[i]] <- nestedBootSample(model = model
                                          ,data = data
                                          ,id_var = id_var
                                          ,response = response)
  }
  
  # Set up clusters
  cores_max <- detectCores()
  if (n_cores > cores_max) {
    stop("Number of cores detected is lower than n_cores!")
  }
  cluster <- makeCluster(n_cores)
  registerDoParallel(cluster)
  
  # Foreach parallel loop over list_samples
  boot_results <- foreach(i = iter(list_samples)
                          ,.combine = rbind
                          ,.packages = "lme4") %dopar% {
    
    # Sample data
    df_sample <- i
    
    # Model and return coefficients
    boot_fit <- robustlmm::rlmerRcpp(form, data = df_sample)
    df_new_filter <- dplyr::filter(data_new, get(id_var) %in% df_sample[[id_var]]) 
    df_pred <- data.frame(y_hat = predict(boot_fit, df_new_filter)
                          ,x = df_new_filter[[x_var]]
                          ,id_var = df_new_filter[[id_var]])
    
    if (type == "mixed") {
      
      df_coef <- data.frame(u0 = attr(lme4::VarCorr(boot_fit)[[1]], "stddev")
                            ,resid = robustlmm::getME(boot_fit, "sigma")
                            ,intercept = lme4::fixef(boot_fit)[1]
                            ,slope = lme4::fixef(boot_fit)[2])
      
    } else if (type == "random") {
      
      df_coef <- data.frame(u0 = attr(VarCorr(boot_fit)[[1]], "stddev")
                            ,resid = robustlmm::getME(boot_fit, "sigma")
                            ,intercept = lme4::fixef(boot_fit)[1])
      
    }
    
    return(list(coefficients = df_coef
                ,predictions = df_pred)
    )           
    
    rm(boot_fit, df_coef, df_pred, df_sample, df_new_filter)
    
  }
  
  return(list(coefficients = dplyr::bind_rows(boot_results[,1])
              ,predictions = dplyr::bind_rows(boot_results[,2])
        )
  )
  
  stopCluster(cluster)
  
}
