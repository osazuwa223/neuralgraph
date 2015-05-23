#' Calculate the signal propagation in the signalgraph model
#' 
#' For a given vertex, the values for the vertex attribute 'input.signal' are calculated as the linear 
#' combination of outputs with edge weights as weights. 'output.signal' is calculated by applying the activation
#' function to the input signal.
#' 
#' @param g, the graph model
#' @param v.index, the index of a vertex
#' @return an updated graph model
calculateVals <- function(g, v.index){ 
  weight_vector <- matrix(E(g)[to(v)]$weight, ncol=1)
  linear_combination <- lucy::parents(g, v) %>%
    ensure_that(length(.) > 0) %>%
    {do.call("cbind", V(g)[.]$output.signal)} %>%
    `%*%`(weight_vector) %>% 
    as.numeric
  V(g)[v]$input.signal <- list(linear.combination) # Add it to input signal attribute
  if(!is.null(g$activation.prime)){
    V(g)[v]$f.prime.input <- list(g$activation.prime(linear.combination)) # Apply f prime and add to the f prime input signal attribute
  }
  output <- g$activation(linear.combination) # apply activation function and add it to activation function attribute
  V(g)[v]$output.signal <- list(output)
  V(g)[v]$output.signal %>% unlist %>% 
    ensure_that(validateVector(.)
  g
}

#' Update the 'signal' at each vertex in a signalgraph object
#' 
#' Updates the 'output.signal' attribute of each vertex in the signalgraph.  Relies on the 'updateVertices' 
#' propagation function in the lucy package, using calculateVals as the callback.
#' @param g a signalgraph object.
#' @return a signal graph object with updated values for the signal attributes.
updateSignals <- function(g){
  updateVertices(getDeterminers = igraph::iparents, callback = calculateVals)
}

#' Update the weight at each edge in a signalgraph object
#' 
#' Updates the 'weight' attribute of each edge in the signalgraph.  Relies on the 'updateedges' 
#' propagation function in the lucy package, using calculateVals as the callback.
#' @param g a signalgraph object.
#' @return a signalgraph object with updated weights.
updateWeights <- function(g){
  updateEdges(g, getDeterminers = getDependentEdges, callback = fitWeightsForEdgeTarget)
}



#' Determine the edges in g whose weights impact the optimization of the weight of edge e 
getDependentEdges <- function(g, e){
  v_trg <- get.edgelist(g)[e, 2] %>%
    {V(g)[.]} 
  v_observed <- V(g)[is.observed]
  dependent.edges <- NULL
  if(!(v_trg == output_v)){
    if(v_observed %in% ichildren(g, v_trg)) {
      dependent_edges <- E(g)[v_trg %->% v_observed]
    }else{
      dependent_edges <- getConnectingEdges(g, v_trg, v_observed)
    }
  }
  dependent.edges
}



#' Fit an initialized network
#' 
#' @param g an igraph object initialized to an unfit signalgraph object
#' @param epsilon when means square area falls below epsilon, stop
#' @param max.iter maximum number of iterations
#' @param verbose if TRUE print messages generated during optimization
fitInitializedNetwork <- function(g, epsilon = 1e-3, max.iter = 100){
  mse <-  getMSE(g)
  for(i %in% 1:max.iter){
    print(i)
    g <- resetUpdateAttributes(g) %>%
      updateEdges %>% 
      updateSignals
    message("First 3 Weights: ", paste(round(E(g)$weight[1:3], 3), collapse =", "), "\n")
    mse <- getMSE(g) %T>% # Get the new MSE
      {message("Mean Squared Error: ", round(., 3), "\n")} %T>% # Message out the error
      {if((mse - .) < epsilon) return(g)} # Stop if improvement is super small.
  }
  g
}

#' Fit a Signalgraph model
#' 
#' Initializes a graph into a signalgraph then fits the model.
#' 
#' @param g igraph object. The vertices must be named. 
#' @param data a data frame. All of the names in the data from must match a vertex name.
#' @param graph_attr list of graph attributes.  Graph attributes include:  
#' \itemize{
#'  \item{penalty}{penalized least squares error penalty parameter value}
#'  \item{activation}{the activation function (this actually is an R function), defaults to logistic.}
#'  \item{activation.prime}{The derivative fo the activation function, used in gradient calculation. Defaults to NULL}
#'  \item{min.max.constraints}{2 element numeric containing the acceptable range for each rate.}
#'  }
#' @param epsilon when means square area falls below epsilon, stop
#' @param max.iter maximum number of iterations
#' @param verbose if TRUE print messages generated during optimization
#' @return A fitted signalgraph object.
#' @export
fitNetwork <- function(g, data, graph_attr, epsilon = 1e-3, max.iter = 100, verbose=FALSE){
  g %>% 
    initializeGraph(graph_attr, n = nrow(data)) %>%
    fitInitializedNetwork(epsilon, max.iter, verbose)
}