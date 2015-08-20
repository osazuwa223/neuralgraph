##################################################################################################
# High level functions for fitting a signal graph
##################################################################################################

#' Fit a Signalgraph model
#' 
#' fitNetwork initializes a graph into a signalgraph then fits the model with penalized least squares.
#' fit_initialized_sg fits an network already initialized by initializeGraph.  fit_initialized_sg
#' is called by fitNetwork.
#' 
#' @param g igraph object. The vertices must be named. 
#' @param data a data frame. All of the names in the data from must match a vertex name.
#' @param fixed names of fixed variables in the vertices
#' @param graph_attr list of graph attributes.  Graph attributes include:  
#' \itemize{
#'  \item{L1_pen:}{ penalized least squares error L1 penalty parameter value}
#'  \item{L2_pen:}{ penalized least squares error L2 penalty parameter value. The default L2 penalty parameter is .04.}
#'  \item{activation:}{ the activation function (this actually is an R function), defaults to logistic.}
#'  \item{activation.prime:}{ The derivative fo the activation function, used in gradient calculation. Defaults to NULL}
#'  \item{min.max.constraints:}{ 2 element numeric containing the acceptable range for each rate.}
#'  }
#' @param min.iter minimum number of iterations
#' @param max.iter maximum number of iterations
#' @param epsilon after the minimum number of iterations, when change in means squared error 
#' between iterations falls below epsilon, cease optimizing and return current estimates
#' @param verbose if set to true, prints messages of details of optimization.
#' @return A fitted signalgraph object.
#' @export
fitNetwork <- function(g, data, fixed = NULL, graph_attr = list(L2_pen = .04), 
                       min.iter = 2,  max.iter = 5, epsilon = 1e-4, verbose = FALSE){
  g %>% 
    initializeGraph(data, fixed, graph_attr) %>%
    fit_initialized_sg(min.iter = min.iter, max.iter = max.iter, epsilon = epsilon, 
                       verbose = verbose)
}
#' @rdname fitNetwork
#' @export
fit_initialized_sg <- function(g, min.iter = 2, max.iter = 5, epsilon = 1e-4, verbose = FALSE){
  mse <- getMSE(g)
  for(i in 1:max.iter){
    g_last <- g
    g <- resetUpdateAttributes(g_last) %>%
      update_weights(verbose = verbose) %>% 
      ensure(!identical(., g_last), err_desc = "failed to produced updated graph") %>%
      resetUpdateAttributes %>%
      update_signals
    mse_last <- mse
    mse <- getMSE(g) %T>% # Get the new MSE
      {message("Mean Squared Error: ", round(., 6), "\n")} 
    if(i >= min.iter && abs(mse - mse_last) < epsilon){
      message("Change in MSE: ", mse - mse_last)
      message("Algorithm stabilized in ", i, " iterations.")
      return(g) # Stop if improvement is super small.
    }
    gc()
    }
  g
}

##################################################################################################
# Calculation of new fitted values (output.signal) given weights in the graph 
##################################################################################################

#' Propagation of fitted values through the graph.
#' 
#' 'update_signals' uses graph propagation to update the fitted values (the output.signal attribute) at each vertex
#'  in a signalgraph object.  It relies on the graph propagation function 'update_vertices' in the lucy
#'  package, using calculate_signal as the callback.
#'  The 'calculate_signal' callback, when applied to a given vertex, calculates the linear combination
#'  of fitted values from the parent nodes, and then applies the activation function to calculate 
#'  the 'output.signal' attribute.  The linear combination is calculated in the function 'getLinearCombination'.
#'  
#' @param g a signalgraph object.
#' @param v the vertex to be updated given the parents
#' @param wts the weights of the incoming edges given the parents
#' @param model.mat a matrix constructing of the "output.signal" attributes of the parents
#' @return a signal graph object with updated values for the signal attributes.
update_signals <- function(g, verbose = FALSE){
  update_vertices(g, get_determiners = lucy::iparents, callback = calculate_signal, verbose = verbose)
}
#' @rdname update_signals
calculate_signal <- function(g, v){ 
  old_output <- unlist(V(g)[v]$output.signal)
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
#' @rdname update_signals
getLinearCombination <- function(wts, model.mat) as.numeric(model.mat %*% wts)


##################################################################################################
# Traversing the network for fitting weights
##################################################################################################

#' Traversing the network for fit weight sets.
#' 
#' To reduce dimensionality of the problem of fitting all the weights in the graph, for each vertex
#' the weights on the incoming edges of that vertex are fit seperately.  Priority goes to 'fitted' vertices 
#' (vertices that are used to calculate the mean squared error loss function), then vertices that 
#' are upstream by 1 degree, then 2 degrees, etc.  'update_weights' traverses the vertices and 
#' fits them individually. It uses 'vertex_updater" from the lucy package to apply the 'fit_weights_for_node'
#' callback on a vertex.  'vertex_updater' will first check all the downstream vertices for those that
#' are not yet updated (via the 'get_downstream_vertices' callback) and recursively fit those if they 
#' are not yet updated.  However, 'order_vertices' is first applied to minimize the need for this
#' recursion.  
#' @param g a signalgraph object.
#' @param verbose if true prints messages about fitting performance.
#' @return a signalgraph object with updated weights.
update_weights <- function(g, verbose = FALSE){
  ordering <- order_vertices(g)
  for(v in ordering){
    g <- vertex_updater(g, v, get_downstream_vertices, fit_weights_for_node, verbose = verbose) 
  }
  if(!all(V(g)$updated)){
    warning("The following vertices had edges that were not updated: ", 
            paste(V(g)[!updated]$name, collapse = ", ")
    )
  }
  g
}
#' @rdname update_weights
order_vertices <- function(g){
  # Start with roots, then all optimization targets that are not roots
  all_optimized <- intersect(V(g)[is.random], V(g)[is.observed])
  roots <- as.numeric(V(g)[is.leaf])
  starting_nodes <- c(roots, setdiff(all_optimized, roots))
  # Then iterate up
  order_remainder <- function(targets){
    target_parents <- lapply(targets, function(target) iparents(g, target)) %>%
      unlist %>% unique %>% setdiff(starting_nodes)
    ordering <- target_parents
    if(length(setdiff(target_parents, get_roots(g))) > 0){
      ordering <- c(ordering, order_remainder(target_parents)) 
    }
    unique(ordering)
  }
  ordering <- c(starting_nodes, order_remainder(starting_nodes))
  # Remmove roots because they have no incoming edges
  setdiff(ordering, V(g)[is.root])
}
#' @rdname update_weights
get_downstream_vertices <- function(g, v){
  # Find observed vertext that are downstream of the target vertex
  random_and_observed <- intersect(V(g)[is.observed], V(g)[is.random])
  intersect(V(g)[random_and_observed], get_downstream_nodes(g, v)) 
}


##################################################################################################
# Optimization of weights
##################################################################################################

#' Fit incoming edge weights for intervention vertex
#' 
#' Constructs a loss function specific to the vertex, and uses this to optimize the weights on the
#' vertices incoming edges using the current weights as the initial values in the optimization.
#' The incoming edge weights are replaced with the optimized values, and the updated graph object.
#' is returned.
#' 
#' @param g a signal graph object
#' @param v the vertex whose incoming edge weights are to be optimized
#' @seealso \code{\link{fit_weights_for_node}}, \code{\link{getObjective}}, \code{\link{getOptimizationFunction}}, \code{\link{getOptimizationFunctionNG}}, \code{\link{getPrediction}}
fit_weights_for_node <- function(g, v){
  message("Optimizing incoming edges to vertex ", V(g)[v]$name)
  lossFunction <- getObjective(g, v) 
  weights.initial <- E(g)[to(v)]$weight
  optimizer <- getOptimizationFunctionNG(g, lossFunction)
  weights.updated <- optimizer(weights.initial)
  E(g)[to(v)]$weight <- weights.updated
  # Add to the trace
  for(j in 1:length(E(g)[to(v)])){
    path <- E(g)[to(v)][j]$path[[1]]
    updated_path <- c(path, weights.updated[j]) 
    E(g)[to(v)][j]$path <- list(updated_path)
    len <- length(updated_path)
    edge_name <- paste(V(g)[get_edge_vertex(g, E(g)[to(v)][j])]$name, collapse = "->")
    message("Last ", min(len, 5), " estimates for edge ", edge_name, ": ", 
            paste(round(updated_path, 3), collapse = ", "))
  }
  V(g)[v]$updated <- TRUE 
  g
}

#' Closure for generating a penalized least squares loss function for use in optimization
#' 
#' This function is a closure that creates a loss function that takes a weight vector as an argument.
#' It produces a function that is passed to closure that generates the optimization function -- this optimization
#' function takes initial weight value and produces an optimized weight value.
#' The closure takes a graph vertex as an argument.  The algorithm identifies the incoming edges to the vertex,
#' then pulls the weights from those edges, constructing a weight vector.  This weight vector is the argument
#' for the loss function, and is used as the starting values in the optimization. The function calculates squared
#' error loss by first taking this weight vector input and predicting the outcome variable, then calculating loss
#' from the predicted and observed values.  
#' @param initial_graph an initialized signal graph, the prefix 'initial' in 'initial_graph' is in contrast to 
#' candidate graphs against which the initial_graph will be compared in the optimization.
#' @param v the vertex whose weights will be optimized.
#' @seealso \code{\link{fit_weights_for_node}}, \code{\link{getObjective}}, \code{\link{getOptimizationFunction}}, \code{\link{getOptimizationFunctionNG}}, \code{\link{getPrediction}}
getObjective <- function(initial_graph, v){
  if(is.null(initial_graph$L1_pen) || is.null(initial_graph$L2_pen)) stop("Penalty has not been specified.")
  lossFunction <- function(wts){
    candidate_graph <- getPrediction(initial_graph, v, wts)
    random_and_observed <- as.numeric(intersect(V(initial_graph)[is.random], V(initial_graph)[is.observed]))
    prediction <- unlist(V(candidate_graph)[random_and_observed]$output.signal)
    observed <- unlist(V(initial_graph)[random_and_observed]$observed)
    sum((observed - prediction) ^ 2) + 
      initial_graph$L1_pen * sum(abs(E(candidate_graph)$weight)) + 
      initial_graph$L2_pen * sum(E(candidate_graph)$weight ^ 2)
  }
  lossFunction
} 

#' Closure for optimization with BFGS and specified gradient function
#' Called by fit_weights_for_node.  Creates a closure based on the 'optim' function.  Uses
#' BFGS or unbounded BFGS depending on if bounds are provided.
#' @seealso \code{\link{fit_weights_for_node}}, \code{\link{getObjective}}, \code{\link{getOptimizationFunction}}, \code{\link{getOptimizationFunctionNG}}, \code{\link{getPrediction}}
getOptimizationFunction <- function(g, lossFunction, getGradient){
  if(!is.null(g$min.max.constraints)){
    low <- g$min.max.constraints[1] 
    high <- g$min.max.constraints[2]
    names(low) <- names(high) <- NULL
    optimFunction <- function(wts){
      optim(wts, fn = lossFunction, gr = getGradient, method="L-BFGS-B",
            lower=rep(low, length(wts)), upper=rep(high, length(wts)))$par
    }
  }else{
    optimFunction <- function(wts){
      optim(wts, fn = lossFunction, gr = getGradient, method="BFGS")$par
    }
  }
  optimFunction
}

#' Closure for optimization with BFGS, no gradient function specified
#' 
#' @seealso \code{\link{fit_weights_for_node}}, \code{\link{getObjective}}, \code{\link{getOptimizationFunction}}, \code{\link{getOptimizationFunctionNG}}, \code{\link{getPrediction}}
getOptimizationFunctionNG <- function(g, lossFunction){
  if(!is.null(g$min.max.constraints)){
    low <- g$min.max.constraints[1] 
    high <- g$min.max.constraints[2]
    names(low) <- names(high) <- NULL
    optimFunction <- function(wts){
      optim(wts, fn = lossFunction, method="L-BFGS-B",
            lower=rep(low, length(wts)), upper=rep(high, length(wts)))$par
    }
  }else{
    optimFunction <- function(wts){
      optim(wts, fn = lossFunction, method="BFGS")$par
    }
  }
  optimFunction
}

#' Get a new graph used for prediction
#' 
#' A given vertex has a given set of incoming edges, each of these edges has weights.
#' Upon supplying the vertex and weights on its incoming edges, this function updates the graph and
#' returns the updated graph object.
#' 
#' @param g, fitted graph
#' @param v, the vertex whose incoming edges will be updated
#' @param weights, a numeric vector containing the weights for those incoming edges.  
#' @return a new updated graph that can be used to generate a prediction.
#' @seealso \code{\link{fit_weights_for_node}}, \code{\link{getObjective}}, \code{\link{getOptimizationFunction}}, \code{\link{getOptimizationFunctionNG}}, \code{\link{getPrediction}}
getPrediction <- function(g, v, new_weights){
  #message("Prediction function call: candidate weights being propagated forward.")
  observed_and_random <- intersect(V(g)[is.random], V(g)[is.observed])
  prediction_graph <- resetUpdateAttributes(g) %>%
    ensure_that(length(E(.)[to(v)]) > 0, 
                err_desc = "Attempted to update incoming edge weights for parentless node.") %>%
    ensure_that(length(E(.)[to(v)]) == length(new_weights),
                err_desc = "# of weights doesn't match # of incoming edges.") 
  if(!all(E(prediction_graph)[to(v)]$weight == new_weights)){ # TRUE when optim test X0
    E(prediction_graph)[to(v)]$weight <- new_weights
    prediction_graph <- update_signals(prediction_graph) %>%
      ensure_that({ #Make sure all the output signals are valid vectors
        lapply(V(.)[is.observed]$output.signal, checkVector) %>% unlist %>% all
      }, err_desc = "One observed vertex has invalid values in output.signal")# %>%
    #       ensure_that({
    #         all(examine_signal_graph(prediction_graph)$vertices$output.signal[[1]][, observed_and_random] ==
    #               examine_signal_graph(g)$vertices$output.signal[[1]][, observed_and_random])
    #       }, err_desc = "signals not getting propagated.")
  } 
  prediction_graph
}

##################################################################################################
# Implementation of gradient descent for backpropagation.
##################################################################################################

#' Gradient descent on weights
#' 
#' @param initial weight values
#' @param loss function that takes a weight vector as an input and returns a single value
#' @param gradient function that takes weight vector as an input are returns a gradient vector of the same length
#' @param step step size
#' @param maxit maximum number of iterations
#' @param epsilon stop when change in loss output is less than epsilon 10 iterations in a row.
#' @seealso \code{\link{gradientDescent}}, \code{\link{getGradientFunction}}, \code{\link{doChainRule}}
gradientDescent <- function(wts_init, grad, loss, maxit = 100, epsilon = .01){
  i <- 0
  little_movement <- 0
  wts <- wts_init
  step <- .001
  ls_maker <- function(wts, gr) Vectorize(function(step) loss(wts - gr * step)) # a closure for doing line search
  while(i < maxit){
    print(i)
    paste("wts=", paste(round(wts, 3), collapse=" ")) %>% print
    gr <- grad(wts) # get gradient
    paste("gr=", paste(round(gr, 3), collapse=" ")) %>% print
    ls <- ls_maker(wts, gr) # get line search function
    step_new <- tryCatch(optimise(ls, c(0, 100))$minimum, # do line search
                         error = function(e) "failed")
    print(step_new)
    wts_new <- wts - step_new * grad(wts) # get a new weight
    e <- abs(loss(wts) - loss(wts_new)) # Check against the old weight
    if(e < epsilon) little_movement <- little_movement + 1 else little_movement <- 0
    if(little_movement > 10){
      message("little movement")
    } 
    wts <- wts_new
    step <- step_new
    i <- i + 1
  }
  wts
}

#' Gradient generating closure
#' 
#' In the network, any given node v has a set of incoming edges.  Those incoming edges have weights.
#' The fitting procedure fits the weights in the network by seperately these weights for each node.
#' More specifically, a seperate optimization is conducted for each node v, where in the weights on the 
#' incoming edges to v are optimized on the loss function.  
#' 
#' This closure returns a function that computes the gradient of all the weights 
#' on the incoming edges of a node v.  The output function is intended to be passed to an optimization
#' functional.  
#' @seealso \code{\link{gradientDescent}}, \code{\link{getGradientFunction}}, \code{\link{doChainRule}}
getGradientFunction <- function(g, v){
  v.incoming.edges <- E(g)[to(v)]
  edge.names <- paste(v.incoming.edges)
  incoming.edge.count <- length(v.incoming.edges)
  v.observed <- V(g)[is.observed]
  gradientFunction <- function(wts){
    #Calculates the gradient for a set of weights using the doChainRule function
    prediction <- getPrediction(g, v, wts)
    observed <- unlist(v.observed$observed)
    loss.function.derivative <- -2 * (observed - prediction) 
    chain.rule.output <- matrix(NA, nrow = g$n, ncol = incoming.edge.count, dimnames = list(NULL, edge.names))
    for(e.index in v.incoming.edges){
      e <- E(g)[e.index]
      chain.rule.output[, paste(e)] <- doChainRule(g, v.observed, e)
    }
    gradient.output <- colSums(
      apply(chain.rule.output, 2, function(chain.rule.result){
        loss.function.derivative * chain.rule.result
      })
    )
    as.numeric(gradient.output)
  }
  gradientFunction
}

#' Back-propagation of derivative calculation
#' 
#' The getGradientFunction closure creates a gradient function specific to a given vertex.
#' It takes as an argument the weights corresponding to the incoming parents of that vertex.
#' When this weight vector is varied, the predictions on the observed vertices change, and thus
#' the cost function output changes.  The gradient function captures that rate of change.  It relies 
#' on back-propagation, calculating the rates of change for every downstream node that is affected
#' by changing the weight vector.  This function does that back-propagation calculation for each
#' element of the weight vector
#' 
#' @param g a model
#' @param v vertex index whose incoming weights are the argument for the gradient function
#' @param e edge index of the individual element of the weight vector where the derivative is being calculated
#' @return a vector corresponding to the derivative
#' @seealso \code{\link{gradientDescent}}, \code{\link{getGradientFunction}}, \code{\link{doChainRule}}
doChainRule <- function(g, v, e){
  if(!(length(v) == 1 && length(e) == 1)) stop("doChainRule() works on only one vertex and only one edge at a time.") 
  e.src <- getEdgeVertex(g, e, "from")
  if(v == e.src) stop("The chainrule has gone back too far, v: ", v, " e: ", e)
  e.trg <- getEdgeVertex(g, e, "to") #Grab the target vertex
  if(!(e.trg %in%  v || isBDownstreamOfA(g, a = e.trg, b = v))){
    stop("You've attempted to find the gradient of a node's output signal
         w.r.t an edge weight that that does not affect that output signal. v: ", v, " e: ", e)  
  }
  f.prime.input <- unlist(V(g)[v]$f.prime.input)
  #Next check that the edge is an incoming edge to v
  if(e.trg == v){
    output <- unlist(V(g)[e.src]$output.signal) * f.prime.input
  }else{
    connected.nodes <- V(g)[lucy::getConnectingNodes(g, e.trg, v)]
    varying.parents <- V(g)[intersect(lucy::iparents(g, v), connected.nodes)]
    parent.names <- paste(varying.parents)
    v.parents.chain.rule.result <- matrix(NA, nrow = g$n, ncol = length(varying.parents), 
                                          dimnames = list(NULL, parent.names))
    for(v.parent.index in varying.parents){
      v.parent <- V(g)[v.parent.index]
      v.parents.chain.rule.result[, paste(v.parent)] <- doChainRule(g, v.parent, e)
    }
    output <- rowSums(
      apply(v.parents.chain.rule.result, 2, function(parent.result){
        parent.result * f.prime.input
      })
    )
  }
  #message("Completed chainrule calculation for node ", v, " w.r.t edge ", e)
  if(!is.numeric(output)){
    stop("v: ", 
         V(g)[v], 
         ", e: ", 
         E(g)[e], " output: ", 
         paste(round(head(output), 2), collapse=" "), 
         ", fpi: ", 
         paste(round(head(f.prime.input), 2), collapse=" ")
    )
  }
  output
}