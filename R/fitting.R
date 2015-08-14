#' Get deviance
#' 
#' Deviance is defined as residual sum squares on all the output variables and their predictions
#' 
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
  #old_output <- unlist(V(g)[v]$output.signal)
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

#' Order edges for ideal weighting
#' Recursively backtracks through structure and orders edges by proximity to the
#' target node, starting with the leaves.
#' @return a numeric array of edge indices
order_edges <- function(g, targets = get_leaves(g)){
  ordering <- as.numeric(E(g)[to(targets)]) 
  sources <- unique(get_edge_vertex(g, ordering, "from"))
  if(length(setdiff(sources, get_roots(g))) > 0){
    ordering <- c(ordering, order_edges(g, sources))
  }
  unique(ordering)
}

#' Update the weight at each edge in a signalgraph object
#' 
#' Updates the 'weight' attribute of each edge in the signalgraph.  Relies on the 'updateedges' 
#' propagation function in the lucy package, using calculateVals as the callback.
#' @param g a signalgraph object.
#' @return a signalgraph object with updated weights.
update_weights <- function(g, verbose = FALSE){
  ordering <- order_edges(g)
  for(e in ordering){
    g <- edge_updater(g, e, get_dependent_edges, fit_weights_for_edge_target, verbose = verbose) 
  }
  if(!all(V(g)$updated)){
    warning("The following were not updated: ", 
            paste(paste(get_edge_vertex(E(g)[!updated]), collapse = "<-"), collapse = ", ")
    )
  }
  g
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
get_dependent_edges <- function(g, e){
  # Find the target vertex of e
  trg_vertex <- V(g)[get.edgelist(g)[e, 2]]
  # Find observed nodes that are downstream of the target vertex 
  downstream_observed <- intersect(V(g)[is.observed], get_downstream_nodes(g, trg_vertex)) 
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
#' @param verbose if TRUE prints info on propagation
#' @export
fit_initialized_sg <- function(g, epsilon = 1e-4, max.iter = 3, verbose = FALSE){
  mse_last <-  getMSE(g)
  for(i in 1:max.iter){
    g <- resetUpdateAttributes(g) %>%
      update_weights(verbose = verbose) %>% 
      update_signals
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
#' Initializes a graph into a signalgraph then fits the model with penalized least squares.
#' The default L2 penalty parameter is .01  
#' 
#' @param g igraph object. The vertices must be named. 
#' @param data a data frame. All of the names in the data from must match a vertex name.
#' @param fixed names of fixed variables in the vertices
#' @param graph_attr list of graph attributes.  Graph attributes include:  
#' \itemize{
#'  \item{L1_pen:}{ penalized least squares error L1 penalty parameter value}
#'  \item{L2_pen:}{ penalized least squares error L2 penalty parameter value}
#'  \item{activation:}{ the activation function (this actually is an R function), defaults to logistic.}
#'  \item{activation.prime:}{ The derivative fo the activation function, used in gradient calculation. Defaults to NULL}
#'  \item{min.max.constraints:}{ 2 element numeric containing the acceptable range for each rate.}
#'  }
#' @param epsilon when means square area falls below epsilon, stop
#' @param max.iter maximum number of iterations
#' @return A fitted signalgraph object.
#' @export
fitNetwork <- function(g, data, fixed = NULL, graph_attr = list(L2_pen = .01), 
                       epsilon = 1e-3, max.iter = 3, verbose = FALSE){
  g %>% 
    initializeGraph(data, fixed, graph_attr) %>%
    fit_initialized_sg(epsilon, max.iter, verbose = verbose)
}




