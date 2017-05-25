#'@export
plot_threshold <- function(wsmcresults){
  gthreshold <- qplot(x = 1:(length(wsmcresults$threshold_history)), y = wsmcresults$threshold_history, geom = "line")
  gthreshold <- gthreshold + xlab("step") + ylab("threshold")
  return(gthreshold)
}

#'@export
plot_threshold_time <- function(wsmcresults){
  gthreshold <- qplot(x = wsmcresults$compute_times, y = wsmcresults$threshold_history, geom = "line")
  gthreshold <- gthreshold + xlab("time (s)") + ylab("threshold")
  return(gthreshold)
}

#'@export
plot_ncomputed <- function(wsmcresults){
  ncomputed <- qplot(x = 1:(length(wsmcresults$ncomputed)), y = wsmcresults$ncomputed, geom = "line")
  ncomputed <- ncomputed + xlab("step") + ylab("# distances computed")
  return(ncomputed)
}

#'@export
plot_bivariate <- function(wsmcresults, i1, i2, from_step = 0){
  wsmc.df <- wsmc_to_dataframe(wsmcresults)
  nsteps <- max(wsmc.df$step)
  parameter_names <- names(wsmc.df)[c(i1, i2)]
  names(wsmc.df)[c(i1,i2)] <- c("temp_name_1", "temp_name_2")
  g <- ggplot(wsmc.df %>% filter(step > from_step), aes(x = temp_name_1, y = temp_name_2, colour = step, group = step))
  g <- g + geom_point(alpha = 0.5)
  g <- g + scale_colour_gradient2(midpoint = from_step+floor((nsteps-from_step)/2)) + theme(legend.position = "none")
  g <- g + xlab(parameter_names[i1]) + ylab(parameter_names[i2])
  return(g)
}

#'@export
plot_bivariate_polygon <- function(wsmcresults, i1, i2, from_step = 0){
  wsmc.df <- wsmc_to_dataframe(wsmcresults)
  nsteps <- max(wsmc.df$step)
  parameter_names <- names(wsmc.df)[c(i1, i2)]
  names(wsmc.df)[c(i1,i2)] <- c("temp_name_1", "temp_name_2")
  g <- ggplot(wsmc.df %>% filter(step > from_step), aes(x = temp_name_1, y = temp_name_2, colour = step, group = step))
  g <- g + stat_density_2d(aes(fill = step), geom = "polygon")
  g <- g + scale_colour_gradient2(midpoint = from_step+floor((nsteps-from_step)/2)) + theme(legend.position = "none")
  g <- g + scale_fill_gradient2(midpoint = from_step+floor((nsteps-from_step)/2))
  g <- g + xlab(parameter_names[i1]) + ylab(parameter_names[i2])
  return(g)
}

#'@export
plot_marginal <- function(wsmcresults, i, from_step = 0){
  wsmc.df <- wsmc_to_dataframe(wsmcresults)
  parameter_name <- names(wsmc.df)[i]
  names(wsmc.df)[i] <- c("temp_name_i")
  g <- ggplot(wsmc.df %>% filter(step > from_step), aes(x = temp_name_i, colour = step, group = step)) + geom_density(aes(y = ..density..))
  g <- g + theme(legend.position = "none") + xlab(parameter_name)
  g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
  return(g)
}

#'@export
plot_marginal_time <- function(wsmcresults, i, from_step = 0){
  wsmc.df <- wsmc_to_dataframe(wsmcresults)
  parameter_name <- names(wsmc.df)[i]
  names(wsmc.df)[i] <- c("temp_name_i")
  g <- ggplot(wsmc.df %>% filter(step > from_step), aes(x = temp_name_i, colour = time, group = factor(time))) + geom_density(aes(y = ..density..))
  g <- g + xlab(parameter_name) # +  theme(legend.position = "none")
  g <- g + scale_color_gradient(name = "time (s)", low = rgb(1,0.5,0.5), high = "darkblue")
  return(g)
}

