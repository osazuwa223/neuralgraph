#' Checks the validity of a numeric vertex attribute
#' 
#' @examples
#' g <- get_gate()
#' isValid(unlist(V(g)[type == "output"]$observed))
isValid <- function(item){
  valid <- TRUE
  if(any(is.na(item))) valid <- FALSE
  if(any(is.infinite(item))) valid <- FALSE
  if(length(item) == 0) valid <- FALSE
  valid
}

#' Calculate the Values of Vertices in Neural Network Model
#' 
#' For a given vertex, input and the output signals are calculated.
#' 
#' @param g, the graph model
#' @param v.index, the index of a vertex
#' @return an updated graph model
calculateVals <- function(g, v.index){
  v <- V(g)[v.index]
  v.parents <- iparents(g, v)
  if(length(v.parents) == 0) stop("Attempting apply activation function
                                  to a node without parents.")
  parent.val.mat <- do.call("cbind", V(g)[v.parents]$output.signal)
  weights <- matrix(E(g)[to(v)]$weight, ncol=1)
  linear.combination <- as.numeric(parent.val.mat %*% weights)
  V(g)[v]$input.signal <- list(linear.combination)
  V(g)[v]$f.prime.input <- list(g$activation.prime(linear.combination))
  output <- g$activation(linear.combination)
  V(g)[v]$output.signal <- list(output)
  V(g)[v]$output.signal %>% unlist %>% {!isValid(.)} %>% `if`(stop("Error in output for vertex ", v.index))
  g
}

#' Reset Model Attributes
#' 
#' Resets the attributes of a model, after some of the attributes have been updated.
#' @param g a model
#' @return A model with updated attributes reset to FALSE
resetUpdateAttributes <- function(g){
  V(g)$updated <- FALSE
  E(g)$updated <- FALSE
  V(g)[type %in% c("input", "intercept")]$updated <- TRUE
  g
}

#' Add Nodes Corresponding to Biases
#' 
#' If intercepts (biases) already exist in the graph, no intercepts are added.
addInterceptNodes <- function(g){
  # If intercepts already exist, do nothing
  if("intercept" %in% V(g)$type) return(g)
  non.input.nodes <- V(g)[inDegree(g, V(g)) != 0]
  g.new <- g
  for(v in non.input.nodes){
    v <- as.numeric(v)
    intercept.name <- paste("int", v, sep=".")
    g.new <- g.new + intercept.name #Add the intercept to the graph
    g.new <- initializeVertexVectors(g.new, intercept.name) #Having added the intercept, give it the correct properties
    V(g.new)[intercept.name]$type <- "intercept" #Label the intercept as 'intercept'
    V(g.new)[intercept.name]$output.signal <- list(rep(1, g$n)) #Give the value of 1
    g.new <- g.new + igraph::edge(intercept.name, V(g)[v]$name)
  }
  g.new
}

#' Initialize numerically-valued vertex attributes
#' 
#' For a given vertex, vertex attibutes that take a numeric vector as a value are 
#' initialized with a placeholder. Specifically, the vertex attributes are;
#' \itemize{
#'  \item{input.signal}{The vector of linear combination of values from the parent nodes}
#'  \item{output.signal}{The output of the activation function applied to input.signal}
#'  \item{f.prime.input}{The derivative of the activation function applied to the input.signal}
#'  \item{observed}{Observed values in the data.  Not present for hidden and input nodes.}
#'  }
#' All of these are initialized with NA values
#' 
#' @param g graph model
#' @param v.index vertex index
#' @return and igraph object where the above attributes are initialized to the value \code{list(rep(NA, g$n))}
initializeVertexVectors <- function(g, v.index){
  na.placeholder <- list(rep(NA, g$n)) 
  V(g)[v.index]$input.signal <- na.placeholder
  V(g)[v.index]$f.prime.input <- na.placeholder
  V(g)[v.index]$output.signal <- na.placeholder
  V(g)[v.index]$observed <- na.placeholder
  g
}

#' Simulating starting weights
initializeWeights <- function(g){  
  if(any(is.infinite(g$min.max.constraints)) || is.null(g$min.max.contraints)) {
    E(g)$weight <- rnorm(ecount(g), sd = 3)
  } else {
    E(g)$weight <- runif(ecount(g), min = g$min.max.constraints[1], max = g$min.max.constraints[2])
  }
  g
}


#' Primes a Graph Object for Fitting a Neural Network Model
#' @param g igraph object with vertices corresponding to input nodes, output nodes, and 
#' hidden nodes.
#' @param input.table data frame of the input values
#' @param output.table data frame of the output values
#' @param activation function, the desired activation function
#' @param activation.prime, the Derivative of the desired activation function
#' @param min.max.constraints (optional) numeric containing the limiting range of the estimates
#' 
#' @return A graph with all the attributes needed to fit the neural network model.
#' 
#' @export
initializeGraph <- function(g, input.table, output.table, activation=logistic, 
                            activation.prime=logistic.prime, min.max.constraints=NULL){
  if(length(
    intersect(list.graph.attributes(g), 
                      c("activation", "activation.prime", 
                              "min.max.constraints", "n"))
            ) > 1 ){
    stop("This graph structure seems to have already been updated.")
  }
  g$activation <- activation
  g$activation.prime <- activation.prime
  if(!is.null(min.max.constraints)) names(min.max.constraints) <- c("min", "max")
  g$min.max.constraints <- min.max.constraints 
  g$n <- nrow(output.table)
  if(!is.null(V(g)$type)) stop("Graph vertices already have a type attribute.")
  V(g)$type <- "middle"
  V(g)[names(input.table)]$type <- "input"
  V(g)[names(output.table)]$type <- "output"
  for(v in V(g)){
    g <- initializeVertexVectors(g, v)
  }
  g <- addInterceptNodes(g)
  g <- resetUpdateAttributes(g)
  for(input in names(input.table)){
    V(g)[input]$output.signal <- list(input.table[, input])
  }
  V(g)$observed <- NA
  for(output in names(output.table)){
    V(g)[output]$observed <- list(output.table[, output])
  }
  #Reinitialize names
  g <- nameEdges(g)
  #initialize weights
  g <- initializeWeights(g)
  ##V(g)$intercept <- runif(vcount(g))
  ##V(g)[type == "input"]$intercept <- NA
  g <- updateVertices(g, getDeterminers = iparents, callback = calculateVals)
  g
}

fitInitializedNetwork <- function(g, epsilon, max.iter, verbose=F){
  e <- 2 * getLoss(g)  / g$n  #Multiply * 2 because of .5 coefficient in loss function. Devide by n to get mean error loss
  i <- 1
  test <- TRUE
  while(test){
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
    test <- (e - e.new) > epsilon || i < max.iter
    i <- i + 1
    e <- e.new
  }
  g
}

#' Fit a Neural Network on a Graph Structure
#' 
#' @export
fitNetwork <- function(g, input.table, output.table, epsilon=.05, max.iter=3, 
                       activation=logistic, activation.prime=logistic.prime, 
                       min.max.constraints = NULL, verbose=F){
  g <- initializeGraph(g, input.table, output.table, 
                       activation = activation, 
                       activation.prime = activation.prime,
                       min.max.constraints = min.max.constraints)
  g <- fitInitializedNetwork(g, epsilon, max.iter, verbose)
}

#' Add new data to fitted model in order to do prediction
newDataUpdate <- function(g, input.table){
  for(input in names(input.table)){
    V(g)[input]$output.signal <- list(input.table[, input])
  }
  g <- updateVertices(g, getDeterminers = iparents, callback = calculateVals)
}

#' Get a data frame of output signals   
getDF <- function(g){
  output <- do.call(data.frame, V(g)[type != "intercept"]$output.signal)
  names(output) <- paste(V(g)[type != "intercept"])
  output
}

getDependentEdges <- function(g, e){
  #Determine the edges in g whose weights impact the 
  #optimization of the weight of edge e 
  e <- E(g)[ e]
  v_trg_name <- get.edgelist(g)[e, 2]
  v.trg <- V(g)[v_trg_name] %>% as.numeric
  output.v <- V(g)[type == "output"] %>% as.numeric
  dependent.edges <- NULL
  if(!(v.trg == output.v)){
    if(output.v %in% ichildren(g, v.trg)) {
      dependent.edges <- E(g)[v.trg %->% output.v]
    }else{
      dependent.edges <- getConnectingEdges(g, v.trg, output.v)
    }
  }
  dependent.edges
}