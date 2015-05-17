#' Checks the validity of a numeric vertex attribute
#' 
#' @examples
#' g <- get_gate()
#' isValidV(unlist(V(g)[is.observed]$observed))
isValidV <- function(item){
  valid <- TRUE
  if(any(is.na(item))) valid <- FALSE
  if(any(is.infinite(item))) valid <- FALSE
  if(length(item) == 0) valid <- FALSE
  valid
}

#' Reset Model Attributes
#' 
#' Resets the 'updated' vertex attribute for each vertex. All vertices, except for roots (which includes biases)
#' have this attribute set to FALSE.  The biases have this attribute set to TRUE.
#' @param g a model
#' @return A model with 'updated' attribute for each vertex is reset.
resetUpdateAttributes <- function(g){
  v_updated <- V(g)[V(g)$is.root]$updated
  if(!is.null(v_updated)){
    if(any(!(v_updated))) stop("Roots should not have FALSE value for updated attribute.")
  }
  V(g)$updated <- FALSE
  E(g)$updated <- FALSE
  V(g)[V(g)$is.root]$updated <- TRUE
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

#' Initialize graph attributes of a graph object
#' 
#' Signalgraph stores key parameters used in fitting and other operations on the signal graph as igraph 
#' attributes.  igraph has 3 kinds of attributes, edge attributes, vertex attributes, and graph attributes.
#' This function sets graph attributes.  Key examples include:
#'  \itemize{
#'  \item{penalty}{penalized least squares error penalty parameter value}
#'  \item{activation}{the activation function (this actually is an R function)}
#'  \item{activation.prime}{The derivative fo the activation function, used in gradient calculation.}
#'  \item{min.max.constraints}{2 element numeric containing the acceptable range for each rate.}
#'  }
#' 
#' @param g an igraph object not yet initialized as a signal graph
#' @param graph_attr a list of objects to be used as graph attributes
#' @param n number of rows in the data
initializeGraphAttributes <- function(g, graph_attr, n){
  basic_attributes <- c("activation", "min.max.constraints", "n", "penalty")
  if(any(basic_attributes %in% list.graph.attributes(g))) stop("Input graph has graph attributes reserved for signal graph.")
  for(attrib in names(graph_attr)) g <- set.graph.attribute(g, attrib, graph_attr[[attr]])
  graph_attributes <- list.graph.attributes(g)
  if(!("penalty" %in% graph_attributes)) g$penalty <- 0 # Have a basic penalty.
  if(!("activation" %in% graph_attributes)){
    g$activation <- logistic
    warning("No activation function specified. Defaulting logistic activation.")
  }
  g$n <- n
  g
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
#'  
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

#
addDataToVertices <- function(g, data){
  if(!(all(names(data) %in% V(g)$name))) stop("Data contains variable names that don't match vertices in the graph.")  
  for(item in names(data)){
    data_list <- list(data[, item])
    V(g)[item]$observed <- data_list
    V(g)[item]$is.observed <- TRUE
    if(V(g)[item]$is.root) V(g)[item]$output.signal <- data_list
  }
}

#' Add Nodes Corresponding to Biases
#' 
#' If intercepts (biases) already exist in the graph, no intercepts are added.
addInterceptNodes <- function(g){
  # If intercepts already exist, do nothing
  if(any(V(g)$is.bias)) return(g)
  non.root.nodes <- V(g)[inDegree(g, V(g)) != 0]
  g.new <- g
  for(v in non.root.nodes){
    v <- as.numeric(v)
    intercept.name <- paste("int", v, sep=".")
    g.new <- g.new + intercept.name #Add the intercept to the graph
    g.new <- initializeVertexVectors(g.new, intercept.name) #Having added the intercept, give it the correct properties
    V(g.new)[intercept.name]$is.bias <- TRUE #Label the intercept as 'intercept'
    V(g.new)[intercept.name]$output.signal <- list(rep(1, g$n)) #Give the value of 1
    g.new <- g.new + igraph::edge(intercept.name, V(g)[v]$name)
  }
  g.new
}

#' Initialize boolean vertex attibutes
#' 
#' There are 4 core boolean vertex attributes:
#' \itemize{
#'  \item{is.bias}{TRUE if the vertex is a bias or intercept}
#'  \item{is.observed}{TRUE if the vertex is observed in the data}
#'  \item{is.root}{TRUE if the vertex is a root in the graph}
#'  \item{is.leaf}{TRUE if the vertex is a leaf in the graph}
#'  }
#'  
#'  This function initializes all these vertex attributes for later updating.
initializeVertexBooleans <- function(g){
  V(g)$is.bias <- FALSE
  V(g)$is.observed <- FALSE
  V(g)$is.root <- FALSE
  V(g)[get_roots(g)]$is.root <- TRUE
  V(g)$is.leaf <- FALSE 
  V(g)[get_roots(g)]$is.leaf <- TRUE
}

applyVertexVectors <- function(g) {
  for(v in V(g)) g <- initializeVertexVectors(g, v) 
  g
}

#' Initialize the vertex attributes and values
#' 
#' @seealso initializeVertexBooleans, addInterceptNodes, applyVertexVectors, addDataToVertices, resetUpdateAttributes
initializeVertices <- function(g){
  g %>%
    initializeVertexBooleans %>%
    addInterceptNodes %>%
    applyVertexVectors %>%
    addDataToVertices %>%
    resetUpdateAttributes
}

#' Initalize edge attributes
#' 
#' Name and initialize the edges.
#' 
initializeEdges <- function(g){
  g %>% 
    nameEdges %>% 
    initializeWeights 
}

#' Priming an igraph bbject for fitting a signalgraph model
#' 
#' Any vertex that does not map to a name in the  
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
#' @return A graph with all the attributes needed to fit the neural network model.
initializeGraph <- function(g, data, graph_attr){
  initializeGraphAttributes(g, graph_attr) %>%
    initializeVertices %>%
    addData(data) %>%
    initializeEdges %>%
    updateVertices(getDeterminers = iparents, callback = calculateVals) # update vertex values given weights
}