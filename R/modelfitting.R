#' Calculate the signal propagation in the signalgraph model
#' 
#' For a given vertex, the values for the vertex attribute 'input.signal' are calculated as the linear 
#' combination of outputs with edge weights as weights. 'output.signal' is calculated by applying the activation
#' function to the input signal.
#' 
#' @param g, the graph model
#' @param v.index, the index of a vertex
#' @return an updated graph model
calculateVals <- function(g, v){ 
  weight_vector <- matrix(E(g)[to(v)]$weight, ncol=1)
  linear_combination <- iparents(g, v) %>%
    ensure_that(length(.) > 0) %>%
    {do.call("cbind", V(g)[.]$output.signal)} %>%
    `%*%`(weight_vector) %>% 
    as.numeric %>%
    ensure_that(checkVector(.))
  V(g)[v]$input.signal <- list(linear_combination) # Add it to input signal attribute
  if(!is.null(g$activation.prime)){
    V(g)[v]$f.prime.input <- list(g$activation.prime(linear_combination)) # Apply f prime and add to the f prime input signal attribute
  }
  output <- g$activation(linear_combination) # apply activation function and add it to activation function attribute
  V(g)[v]$output.signal <- list(output)
  V(g)[v]$output.signal %>% unlist %>% 
  ensure_that(checkVector(.))
  g
}

#' Update the 'signal' at each vertex in a signalgraph object
#' 
#' Updates the 'output.signal' attribute of each vertex in the signalgraph.  Relies on the 'updateVertices' 
#' propagation function in the lucy package, using calculateVals as the callback.
#' @param g a signalgraph object.
#' @return a signal graph object with updated values for the signal attributes.
updateSignals <- function(g){
  updateVertices(g, getDeterminers = lucy::iparents, callback = calculateVals)
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

#' Find edges that affect an edge's optimization
#' 
#' When optimizing weights by back-propagation, before optimizing the objective function for a 
#' given weight, it is necessary to first optimize for all weights that are no paths downstream
#' from the given weight and upstream of an observed variable.  This function identifies all such 
#' downstream edges relative an input edge.
#' 
#' @param g a signalgraph object
#' @param e an edge index in the signal graph
#' @return set of edge indices 
getDependentEdges <- function(g, e){
  # Find the target vertex of e
  trg_vertex <- V(g)[get.edgelist(g)[e, 2]]
  # Find observed nodes that are downstream of the target vertex 
  downstream_observed <- intersect(V(g)[is.observed], getDownstreamNodes(g, trg_vertex)) 
  #Find edges on paths between the target vertex and the downstream observed nodes.
  dependent_edges <- NULL
  for(v in downstream_observed){
    # Use tryCatch here because I want this function to return NULL instead of erroring out 
    # if a node is not downstream.
    connecting_edges <- tryCatch(getConnectingEdges(g, trg_vertex, v),
                                 error = function(e) NULL)
    # tag 'em and bag 'em
    dependent_edges <- c(dependent_edges, connecting_edges)
  }
  dependent_edges %>% ensure_that(!(e %in% .)) # sanity check
}

#' Fit an initialized network
#' 
#' @param g an igraph object initialized to an unfit signalgraph object
#' @param epsilon when means square area falls below epsilon, stop
#' @param max.iter maximum number of iterations
fitInitializedNetwork <- function(g, epsilon = 1e-4, max.iter = 3){
  mse_last <-  getMSE(g)
  for(i in 1:max.iter){
    g <- resetUpdateAttributes(g) %>%
      updateWeights %>% 
      updateSignals
    `if`(ecount(g) > 3,
         message("First 3 Weights: ", paste(round(E(g)$weight[1:3], 3), collapse =", "), "\n"),
         message("Weights: ", paste(round(E(g)$weight, 3), collapse =", "), "\n")
    )
    mse <- getMSE(g) %T>% # Get the new MSE
      {message("Mean Squared Error: ", round(., 6), "\n")} 
    if(abs(mse - mse_last) < epsilon){
        message("Algorithm stabilized in ", i, " iterations.")
        return(g) # Stop if improvement is super small.
      } 
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
#'  \item{L1_pen}{penalized least squares error L1 penalty parameter value}
#'  \item{L2_pen}{penalized least squares error L2 penalty parameter value}
#'  \item{activation}{the activation function (this actually is an R function), defaults to logistic.}
#'  \item{activation.prime}{The derivative fo the activation function, used in gradient calculation. Defaults to NULL}
#'  \item{min.max.constraints}{2 element numeric containing the acceptable range for each rate.}
#'  }
#' @param epsilon when means square area falls below epsilon, stop
#' @param max.iter maximum number of iterations
#' @return A fitted signalgraph object.
#' @export
fitNetwork <- function(g, data, fixed = NULL, graph_attr = NULL, epsilon = 1e-3, max.iter = 3){
  g %>% 
    initializeGraph(data, fixed, graph_attr) %>%
    fitInitializedNetwork(epsilon, max.iter)
}