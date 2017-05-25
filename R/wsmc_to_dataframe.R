#'@export
wsmc_to_dataframe <- function(results, ...){
  th <- results$thetas_history
  if (!is.null(results$target$parameter_names)){
    parameter_names <- results$target$parameter_names
  } else {
    parameter_names <- paste0("X", 1:ncol(results$thetas_history[[1]]))
  }
  nsteps <- length(th)
  df <- data.frame()
  for (i in 1:(nsteps)){
    df_ <- data.frame(cbind(th[[i]], rep(i, nrow(th[[i]])), rep(as.numeric(results$compute_times)[i], nrow(th[[i]]))))
    names(df_) <- c(parameter_names, "step", "time")
    df <- rbind(df, df_)
  }
  names(df) <- c(parameter_names, "step", "time")
  return(df)
}


# wsmc_to_dataframe <- function(results, parameter_names){
#   th <- results$thetas_history
#   nsteps <- length(th)
#   df <- data.frame()
#   for (i in 1:(nsteps)){
#     df_ <- data.frame(cbind(th[[i]], rep(i, nrow(th[[i]]))))
#     names(df_) <- c(parameter_names, "step")
#     df <- rbind(df, df_)
#   }
#   names(df) <- c(parameter_names, "step")
#   return(df)
# }
