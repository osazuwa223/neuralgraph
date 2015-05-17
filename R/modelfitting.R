
#' Calculate the Values of Vertices in Neural Network Model
#' 
#' For a given vertex, the values for the vertex attribute 'input.signal' are calculated as the linear 
#' combination of outputs with edge weights as weights. 'output.signal' is calculated by applying the activation
#' function to the input signal.
#' 
#' @param g, the graph model
#' @param v.index, the index of a vertex
#' @return an updated graph model
calculateVals <- function(g, v.index){
  v <- V(g)[v.index]
  v.parents <- iparents(g, v) # Grab the parents of the index
  if(length(v.parents) == 0) stop("Attempting apply activation function
                                  to a node without parents.")
  parent.val.mat <- do.call("cbind", V(g)[v.parents]$output.signal) # Create a matrix with the parent node output values in the colums
  weights <- matrix(E(g)[to(v)]$weight, ncol=1) # Create a verticle weight vector
  linear.combination <- as.numeric(parent.val.mat %*% weights) # Calculate the linear combination and convert to numeric
  V(g)[v]$input.signal <- list(linear.combination) # Add it to input signal attribute
  if(!is.null(g$activation.prime)){
    V(g)[v]$f.prime.input <- list(g$activation.prime(linear.combination)) # Apply f prime and add to the f prime input signal attribute
  }
  output <- g$activation(linear.combination) # apply activation function and add it to activation function attribute
  V(g)[v]$output.signal <- list(output)
  V(g)[v]$output.signal %>% unlist %>% {!isValidV(.)} %>% `if`(stop("When calculating values of vertex ", v.index, " there were NA, infinite, or otherwise invalid values."))
  g
}

#' Determine the edges in g whose weights impact the optimization of the weight of edge e 
getDependentEdges <- function(g, e){
  e <- E(g)[ e]
  v_trg_name <- get.edgelist(g)[e, 2]
  v.trg <- V(g)[v_trg_name] %>% as.numeric
  v.observed <- V(g)[is.observed] %>% as.numeric
  dependent.edges <- NULL
  if(!(v.trg == output.v)){
    if(v.observed %in% ichildren(g, v.trg)) {
      dependent.edges <- E(g)[v.trg %->% v.observed]
    }else{
      dependent.edges <- getConnectingEdges(g, v.trg, v.observed)
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
fitInitializedNetwork <- function(g, epsilon = 1e-3, max.iter = 100, verbose=F){
  e <-  getLoss(g)  / g$n 
  i <- 0
  test <- TRUE
  while(i < max.iter){
    print(i)
    g <- resetUpdateAttributes(g)
    if(verbose){
      g <- updateEdges(g, getDeterminers = getDependentEdges, callback = fitWeightsForEdgeTarget)
    }else{
      g <- suppressMessages(updateEdges(g, getDeterminers = getDependentEdges, callback = fitWeightsForEdgeTarget))
    }
    g <- updateVertices(g, getDeterminers = iparents, callback = calculateVals)
    e.new <- 2 *  getLoss(g) / g$n  
    message("Error: ", round(e.new, 3), "\n")
    message("Weights: ", paste(round(E(g)$weight[1:3], 3), collapse =", "), "\n")
    if((e - e.new) < epsilon) return(g)
    e <- e.new
    i <- i + 1
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