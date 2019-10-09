#' @name nestedBootSample
#' @author Ed Wilkes
#' 
#' @description performs nested bootstrap resampling of a given data frame
#' 
#' @param model lmer or rlmer object containing fitted model
#' @param data Data frame containing original data to be sampled
#' @param id_var Vector of ID variables corresponding to nesting levels
#' @param response String denoting the name of the response variable being modelled
#'
#' @return Data frame containing nested bootstrap sampled data
#'
nestedBootSample <- function(model
                             ,data
                             ,id_var
                             ,response) {
  
  # This code has been adapted from the excellent response from /u/Ben Bolker on this
  # Stackoverflow thread: https://stats.stackexchange.com/questions/231074/confidence-intervals-on-predictions-for-a-non-linear-mixed-model-nlme
  
  # Packages
  require(dplyr)
  
  pp <- predict(model)
  rr <- residuals(model)
  dd <- data.frame(data, y_hat = pp, residual = rr)
  
  ## Sample top-level groups "id_var" with replacement
  u0 <- levels(dd[[id_var]])
  bsamp_1 <- sample(u0, size = length(u0), replace = TRUE)
  
  ## Sample next level (residuals), *within* top-level groups, with replacement
  bsamp_2 <- lapply(bsamp_1,
                    function(x) {
                      ddb <- dd[dd[[id_var]] == x,]
                      ddb$y_star <- ddb$y_hat +
                        sample(ddb$residual, size = nrow(ddb), replace = TRUE)
                      return(ddb)
                    })
  
  ## Flatten bootstrap results into data frame
  res <- bind_rows(bsamp_2)
  return(res)
  
}