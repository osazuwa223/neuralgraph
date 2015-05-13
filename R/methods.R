#' Get deviance
#' 
#' Deviance is defined as residual sum squares on all the output variables and their predictions
#' 
#' @param g, a fitted signal graph model
#' @export
get_deviance <- function(g){
  sum((unlist(V(g)[is.observed]$output.signal) - unlist(V(g)[is.observed]$observed) )^2)
}