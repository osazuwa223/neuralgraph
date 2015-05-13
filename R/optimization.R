#' The logistic function
#' 
#' @export
logistic <- function(z) 1 / (1 + exp(-z))
#' @rdname logistic
#' @export 
logistic_prime <- function(z) logistic(z) * (1 - logistic(z))

#' A Simple Matrix Multiplication to Calculate Linear Inputs
getLinearCombination <- function(wts, model.mat) as.numeric(model.mat %*% wts)

#' Back-propagation of derivative calculation
#' 
#' The getGradientFunction closure creates a gradient function specific to a given vertex.
#' It takes as an argument the weights corresponding to the incoming parents of that vertex.
#' When this weight vector is varied, the predictions on the graphs's outputs change, and thus
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
    stop("You've attempted to find the gradient of a node's output
         w.r.t an edge weight that that does not affect that output. v: ", v, " e: ", e)  
  }
  f.prime.input <- unlist(V(g)[v]$f.prime.input)
  #Next check that the edge is an incoming edge to v
  if(e.trg == v){
    output <- unlist(V(g)[e.src]$output.signal) * f.prime.input
  }else{
    connected.nodes <- V(g)[getConnectingNodes(g, e.trg, v)]
    varying.parents <- V(g)[intersect(iparents(g, v), connected.nodes)]
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
  prediction.graph <- resetUpdateAttributes(g)
  if(length(E(prediction.graph)[to(v)]) == 0) stop("Attempted to update incoming edge weights for parentless node.")
  if(length(E(prediction.graph)[to(v)]) != length(new_weights)) stop("# of weights doesn't match # of incoming edges.")
  E(prediction.graph)[to(v)]$weight <- new_weights
  prediction.graph <- updateVertices(prediction.graph, getDeterminers = iparents, callback = calculateVals)
  prediction <- unlist(V(prediction.graph)[type == "output"]$output.signal)
  if(!isValid(prediction)) stop("An error occured in predicting vertex ", v)
  prediction.graph
}

getLoss <- function(g){
  observed <- unlist(V(g)[type=="output"]$observed)
  prediction <- unlist(V(g)[type=="output"]$output.signal)
  sum( (observed - prediction) ^ 2)
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
  if(is.null(initial_graph$penalty)) stop("Penalty has not been specified.")
  lossFunction <- function(wts){
    candidate_graph <- getPrediction(initial_graph, v, wts)
    prediction <- unlist(V(candidate_graph)[type=="output"]$output.signal)
    observed <- unlist(V(initial_graph)[type=="output"]$observed)
    sum((observed - prediction) ^ 2) + initial_graph$penalty * sum(E(candidate_graph)$weight ^ 2)
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
  output.node <- V(g)[type=="output"]
  gradientFunction <- function(wts){
    #Calculates the gradient for a set of weights using the doChainRule function
    prediction <- getPrediction(g, v, wts)
    observed <- unlist(output.node$observed)
    loss.function.derivative <- -2 * (observed - prediction) 
    chain.rule.output <- matrix(NA, nrow = g$n, ncol = incoming.edge.count, dimnames = list(NULL, edge.names))
    for(e.index in v.incoming.edges){
      e <- E(g)[e.index]
      chain.rule.output[, paste(e)] <- doChainRule(g, output.node, e)
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
    low <- g$min.max.constraints["min"] 
    high <- g$min.max.constraints["max"]
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
    low <- g$min.max.constraints["min"] 
    high <- g$min.max.constraints["max"]
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
  v <- V(g)[v]
  lossFunction <- getObjective(g, v)
  getGradient <- getGradientFunction(g, v)
  weights.initial <- E(g)[to(v)]$weight
  optimizer <- getOptimizationFunctionNG(g, lossFunction)
  weights.updated <- optimizer(weights.initial)
  E(g)[to(v)]$weight <- weights.updated
  g
}

fitWeightsForEdgeTarget <- function(g, e){
  old_weight <- E(g)[e]$weight
  edge.target <- get.edgelist(g)[e, 2]
  message("Fitting for edge ", E(g)[e]$name)
  g <- fitWeightsForNode(g, edge.target)
  new_weight <- E(g)[e]$weight
  message(E(g)[e]$name, ': Old weight = ', old_weight, ', New Weight = ', new_weight)
  E(g)[e]$updated <- TRUE
  g
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
    if(step_new == "failed") browser()
    wts_new <- wts - step_new * grad(wts) # get a new weight
    e <- abs(loss(wts) - loss(wts_new)) # Check against the old weight
    if(e < epsilon) little_movement <- little_movement + 1 else little_movement <- 0
    if(little_movement > 10){
      message("little movement")
      browser() # Return if we already have little movement
    } 
    wts <- wts_new
    step <- step_new
    i <- i + 1
  }
  wts
}