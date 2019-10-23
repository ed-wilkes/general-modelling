#' @name clusterBootSample
#' @author Ed Wilkes
#' 
#' @description performs non-parametric bootstrap resampling of a given data frame by
#' cluster resampling only
#' 
#' @param data Data frame containing original data to be sampled
#' @param id_var String denoting the variable to cluster data by
#'
#' @return Data frame containing cluster bootstrap resampled data
#'
clusterBootSample <- function(data, id_var) {
  
  # Packages
  require(dplyr)
  
  # Sample clusters (i.e., id_var)
  u0 <- unique(as.character(data[[id_var]]))
  u0_sample <- sample(u0, size = length(u0), replace = TRUE)
  
  list_u0_sample <- lapply(u0_sample, function(x) {
    df_sample <- as.data.frame(filter(data, get(id_var) == x))
    return(df_sample)
  })
  
  data_boot <- as.data.frame(bind_rows(list_u0_sample))
  return(data_boot)
  
}