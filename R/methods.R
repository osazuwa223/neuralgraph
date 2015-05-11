#' Get deviance
#' 
#' Deviance is defined as residual sum squares on all the output variables and their predictions
#' 
#' @param g, a fitted signal graph model
#' @export
get_deviance <- function(g){
  sum((unlist(V(g)[type == "output"]$output.signal) - unlist(V(g)[type == "output"]$observed) )^2)
}