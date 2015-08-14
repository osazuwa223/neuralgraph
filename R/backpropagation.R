
#' A Simple Matrix Multiplication to Calculate Linear Inputs
getLinearCombination <- function(wts, model.mat) as.numeric(model.mat %*% wts)

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

#' Gradient descent on weights
#' 
#' @param initial weight values
#' @param loss function that takes a weight vector as an input and returns a single value
#' @param gradient function that takes weight vector as an input are returns a gradient vector of the same length
#' @param step step size
#' @param maxit maximum number of iterations
#' @param epsilon stop when change in loss output is less than epsilon 10 iterations in a row.
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

#' Update the 'signal' at each vertex in a signalgraph object
#' 
#' Updates the 'output.signal' attribute of each vertex in the signalgraph.  Relies on the 'updateVertices' 
#' propagation function in the lucy package, using calculateVals as the callback.
#' @param g a signalgraph object.
#' @return a signal graph object with updated values for the signal attributes.
update_signals <- function(g){
  update_vertices(g, get_determiners = lucy::iparents, callback = calculateVals)
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

#' MSE calculation over random variables in the graph
#' 
#' Calculates the mean squared error for all the random variables in a signalgraph 
#' object, then reports the mean of means.
#' 
#' @param g a signal graph object
#' @export 
getMSE <- function(g){
  response_variables <- intersect(V(g)[is.random], V(g)[is.observed])
  k <- length(V(g)[response_variables])
  observed <- unlist(V(g)[response_variables]$observed)
  prediction <- unlist(V(g)[response_variables]$output.signal)
  out <- sum((observed - prediction) ^ 2) / g$n / k
  out
}

#' Closure for generating a penalized least squares loss function for use in optimization
#' 
#' This function is a closure that creates a loss function that takes a weight vector as an argument.
#' The closure takes a graph vertex as an argument.  The algorithm identifives the incoming edges to the vertex,
#' then pulls the weights from those edges, constructing a weight vector.  This weight vector is the argument
#' for the loss function, and is used as the starting values in the optimization. The function calculates squared
#' error loss by first taking this weight vector input and predicting the outcome variable, then calculating loss
#' from the predicted and observed values.  
#' @param initial_graph an initialized signal graph, the prefix 'initial' in 'initial_graph' is in contrast to 
#' candidate graphs against which the initial_graph will be compared in the optimization.
#' @param v the vertex whose weights will be optimized.
getObjective <- function(initial_graph, v){
  if(is.null(initial_graph$L1_pen) || is.null(initial_graph$L2_pen)) stop("Penalty has not been specified.")
  lossFunction <- function(wts){
    candidate_graph <- getPrediction(initial_graph, v, wts)
    random_and_observed <- as.numeric(intersect(V(g)[is.random], V(g)[is.observed]))
    prediction <- unlist(V(candidate_graph)[random_and_observed]$output.signal)
    observed <- unlist(V(initial_graph)[random_and_observed]$observed)
    sum((observed - prediction) ^ 2) + 
      initial_graph$L1_pen * sum(abs(E(candidate_graph)$weight)) + 
      initial_graph$L2_pen * sum(E(candidate_graph)$weight ^ 2)
  }
  lossFunction
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
#' functional.  I use a closure because the optim functional only accepts a weights are
# 1) Get the parents and parent matrices
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

#' Closure for optimization with BFGS and specified gradient function
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

#' CLosure for optimization with BFGS, no gradient function specified
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

fitWeightsForNode <- function(g, v){
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
  g
}

#' Fit edge weight
#' 
#' Adjust the weight for a given edge by minimizing a loss function.  This function takes a signal edge
#' as its argument.  It identifies the target vertex of the edge, then calls fitWeightsForNode, 
#' which optimizes the loss function over all incoming edges to that target vertex, rather than just the 
#' edge argument alone.  Once the weight has been optimized, the updated attribute is set to true.
#' 
#' @param g a signalgraph object
#' @param e an index for the edge in the signalgraph object
#' @return an updated signalgraph object
#' @seealso fitWeightsForNode
fit_weights_for_edge_target <- function(g, e){
  e <- as.numeric(e)
  old_weight <- E(g)[e]$weight
  edge.target <- get_edge_vertex(g, e, "to")
#   message("Fitting for edge ", 
#           paste(V(g)[get_edge_vertex(g, E(g)[e])]$name, 
#                 collapse = "->"))
  g <- fitWeightsForNode(g, edge.target)
#   new_weight <- E(g)[e]$weight
#   message(E(g)[e]$name, ': Old weight = ', round(old_weight, 4), ', New Weight = ', round(new_weight, 4))
  E(g)[e]$updated <- TRUE
  g
}

#' The logistic function
#' 
#' @export
logistic <- function(z) 1 / (1 + exp(-z))
#' @rdname logistic
#' @export 
logistic_prime <- function(z) logistic(z) * (1 - logistic(z))
