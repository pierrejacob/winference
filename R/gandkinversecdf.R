#'@export
gandkinversecdf <- function(uniforms, theta){
  return(gandkinversecdf_(uniforms, theta))
}

#'@export
gandkinversecdf_givennormals <- function(normals, theta){
  return(gandkinversecdf_givennormals_(normals, theta))
}
