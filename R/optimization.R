#' The logistic function
#' 
#' @export
logistic <- function(z) 1 / (1 + exp(-z))
#' @rdname logistic
#' @export 
logistic_prime <- function(z) logistic(z) * (1 - logistic(z))

#' A Simple Matrix Multiplication to Calculate Linear Inputs
getLinearCombination <- function(weights, model.mat) as.numeric(model.mat %*% weights)

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
  e.src <- getEdgeVertex(g, e, "from")
  if(v == e.src) stop("The chainrule has gone back too far, v: ", v, " e: ", e)
  e.trg <- getEdgeVertex(g, e, "to") #Grab the target vertex
  if(!(e.trg %in%  v || isBDownstreamOfA(g, a = e.trg, b = v))){
    stop("You've attempted to find the gradient of a node's output
         w.r.t an edge weight that that does not affect that output. v: ", v, " e: ", e)  
  }
  f.prime.input <- unlist(V(g)[v]$f.prime.input)
  #Next check that the edge is not an incoming edge to v
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

#' Get a new predictions for output vertices after reweighting edges
#' 
#' A given vertex has a given set of incoming edges, each of these edges has weights.
#' Upon supplying the vertex and weights on its incoming edges, this function updates the prediction
#' on the output edges.
#' 
#' @param g, fitted graph
#' @param v, the vertex whose incoming edges will be updated
#' @param weights, a numeric vector containing the weights for those incoming edges.  
#' @export   
getPrediction <- function(g, v, new_weights){
  #message("Prediction function call: candidate weights being propagated forward.")
  prediction.graph <- resetUpdateAttributes(g)
  if(length(E(prediction.graph)[to(v)]) == 0) stop("Attempted to update incoming edge weights for parentless node.")
  E(prediction.graph)[to(v)]$weight <- new_weights
  prediction.graph <- updateVertices(prediction.graph, getDeterminers = iparents, callback = calculateVals)
  prediction <- unlist(V(prediction.graph)[type == "output"]$output.signal)
  if(!isValid(prediction)) stop("An error occured in predicting vertex ", v)
  prediction
}

getLoss <- function(g){
  observed <- unlist(V(g)[type=="output"]$observed)
  prediction <- unlist(V(g)[type=="output"]$output.signal)
  .5 * sum( (observed - prediction) ^ 2)
}

#' Create a loss function dependent on a given vertex
#' 
#' This is a closure that creates a loss function that varies given a weight vector.
#' The weights are weight attributes of incoming edges to the vertex given in the argument.
#' Thus for a given vertex, you can inspect how varying the weights that determine the vertex affect loss. 
getLossFunction <- function(g, v){
  #Creating a temporary graph where new weights for v are added, the values are propagated forward.
  #And a new prediction is generated
  lossFunction <- function(weights){
    prediction <- getPrediction(g, v, weights)
    observed <- unlist(V(g)[type=="output"]$observed)
    .5 * sum( (observed - prediction) ^ 2)
  }
  lossFunction
} 

getGradientFunction <- function(g, v){
  # This closure returns a function that computes the gradient of all the weights 
  # on the incoming edges of a node v.  Unlike in multi-layer perceptrons where the
  # weights at each level are optimized together, here the weights for each incoming 
  # edge for each node are optimized together
  # 1) Get the parents and parent matrices
  v.incoming.edges <- E(g)[to(v)]
  edge.names <- paste(v.incoming.edges)
  incoming.edge.count <- length(v.incoming.edges)
  output.node <- V(g)[type=="output"]
  gradientFunction <- function(weights){
    #Calculates the gradient for a set of weights using the doChainRule function
    prediction <- getPrediction(g, v, weights)
    observed <- unlist(output.node$observed)
    loss.function.derivative <- -1 * (observed - prediction) # derivative of .5 * squared loss
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
    gradient.output    
  }
  gradientFunction
}

getOptimizationFunction <- function(g, lossFunction, getGradient){
  if(!is.null(g$min.max.constraints)){
    lower <- rep(g$min.max.constraints["min"], length(weights))
    upper <- rep(g$min.max.constraints["max"], length(weights))
    names(lower) <- names(upper) <- NULL
    optimFunction <- function(weights){
      optim(weights, fn = lossFunction, gr = getGradient, method="L-BFGS-B",
            lower=lower, upper=upper)$par
    }
  }else{
    optimFunction <- function(weights){
      optim(weights, fn = lossFunction, gr = getGradient, method="BFGS")$par
    }
  }
  optimFunction
}

fitWeightsForNode <- function(g, v){
  v <- V(g)[v]
  lossFunction <- getLossFunction(g, v)
  getGradient <- getGradientFunction(g, v)
  weights.initial <- E(g)[to(v)]$weight
  optimizer <- getOptimizationFunction(g, lossFunction, getGradient)
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
  if(new_weight == old_weight) message(E(g)[e]$name, ': Old weight = ', old_weight, ', New Weight = ', new_weight)
  E(g)[e]$updated <- TRUE
  g
}