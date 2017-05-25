#'@export
gandkinversecdf <- function(uniforms, theta){
  return(gandkinversecdf_(uniforms, theta))
}

#'@export
gandkinversecdf_givennormals <- function(normals, theta){
  return(gandkinversecdf_givennormals_(normals, theta))
}

#'@export
gandkcdf <- function(y, theta, maxsteps = 1000, tolerance = 1e-10, lower = 1e-20, upper = 1-1e-20){
  return(gandkcdf_(y, theta, maxsteps, tolerance, lower, upper))
}
