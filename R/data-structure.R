

#' Check validity of arguments
#' 
#' Check validity of both inputs to constuction of signalgraph object
#' @param g an igraph object, vertices should be named.
#' @param data a dataframe, each variable name should match a vertex name in the graph
#' @return the input graph, ready for passing to a subsequent function
checkArgs <- function(g, data){
  if(!is.directed(g)) stop("Signal graph requires a directed graph")
  if(vcount(g) < 2 || ecount(g) == 0) stop("There has to be at least 2 vertices and 1 directed edge.")
  if(is.null(V(g)$name)) stop("Vertices must be named")
  if(!all(names(data) %in% V(g)$name)) stop("Data contains variables that are not named in the graph.")
  leaves <- V(g)[get_leaves(g)]$name
  roots <- V(g)[get_roots(g)]$name
  if(length(intersect(roots, leaves)) != 0) stop("Detected at least one vertex that is both a root and a leaf. Perhaps there is an unconnected vertex?")
  if(length(setdiff(c(roots, leaves), names(data))) > 0) stop("Graph roots and leaves must be observed in the data.")  
  basic_attributes <- c("activation", "min.max.constraints", "n", "penalty")
  if(any(basic_attributes %in% list.graph.attributes(g))) stop("Input graph has graph attributes reserved for signal graph.")
  g
}

#' Checks the validity of a vertex vector attribute
#' 
#' @examples
#' g <- get_gate()
#' checkVector(unlist(V(g)[is.observed]$observed))
checkVector <- function(item){
  valid <- TRUE
  if(any(is.na(item))) valid <- FALSE
  if(any(is.infinite(item))) valid <- FALSE
  if(length(item) == 0) valid <- FALSE
  valid
}

#' Reset Model Attributes
#' 
#' Resets the 'updated' attribute for each vertex and edge. All vertices, except for roots (which includes biases)
#' have this attribute set to FALSE.  The biases have this attribute set to TRUE -- this attribute should always
#' be true for biases.
#' @param g a igraph object
#' @return A model with 'updated' attribute for each vertex and edge is reset.
resetUpdateAttributes <- function(g){
  root_update <- V(g)[is.root]$updated
  if(!is.null(root_update)){
    if(any(!(root_update))) stop("Roots should not have FALSE value for updated attribute.")
  }
  V(g)$updated <- FALSE
  E(g)$updated <- FALSE
  V(g)[V(g)$is.root]$updated <- TRUE
  g
}

#' Simulating starting weights
#' 
#' Starting weights used when fitting signalgraph objects are simulated from a uniform distribution
#' where the range is given by the 'min.max.constraints' attribute.
#' 
#' @param g igraph object
#' @return a graph with updated weights
initializeWeights <- function(g){  
  if(any(is.infinite(g$min.max.constraints)) || is.null(g$min.max.constraints)) {
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
  for(attrib in names(graph_attr)) g <- set.graph.attribute(g, attrib, graph_attr[[attrib]])
  graph_attributes <- list.graph.attributes(g)
  if(!("penalty" %in% graph_attributes)) g$penalty <- 0 # Have a basic penalty.
  if(!("activation" %in% graph_attributes)){
    g$activation <- logistic
    message("No activation function specified. Defaulting logistic activation.")
  }
  g$n <- n
  g
}

#' Initialize numerically-valued vertex attributes
#' 
#' For each vertex, there are several sets of vectors representing data and calculations
#' made on the data.  These are stored as vertex attributes in the form of numeric objects
#' wrapped in a list so the length is one.  The vector attributes are;
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
initializeVectorAttributesForVertex <- function(g, v.index){
  na.placeholder <- list(rep(NA, g$n)) 
  V(g)[v.index]$input.signal <- na.placeholder
  V(g)[v.index]$f.prime.input <- na.placeholder
  V(g)[v.index]$output.signal <- na.placeholder
  V(g)[v.index]$observed <- na.placeholder
  V(g)
  g
}

#' Add data to vertices
#' 
#' Adds data to the 'observed' attributes.  Any vertex that does not map to a name in the data will be treated 
#' as a hidden node, meaning \code{V(g)[v]$is.observed} will be \code{FALSE}.
#' If a variable in the data is not present amongst the graph vertices, an error is thrown.
#' 
#' @param g igraph object. The vertices must be named. 
#' @param data a data frame. All of the names in the data from must match a vertex name.
#' @return an igraph object where the 
addDataToVertices <- function(g, data){
  for(item in names(data)){
    data_list <- list(data[, item])
    V(g)[item]$observed <- data_list
    V(g)[item]$is.observed <- TRUE 
    if(V(g)[item]$is.root) V(g)[item]$output.signal <- data_list
  }
  # Label the hidden nodes
  hidden_from_data <- V(g)[!is.bias] %>% # pull the non-bias nodes
    {setdiff(.$name, names(data))} # find those that are not in the data  
  V(g)[hidden_from_data]$is.hidden <- TRUE
  # Give 1 value ot the interceopts
  V(g)[is.bias]$output.signal <- list(rep(1, g$n))
  g
}

#' Add Nodes Corresponding to Biases
#' 
#' If biases already exist in the graph, no biases are added.
#' This function uses the igraph method on the `+` generic. In this code it is called
#' explicitly with igraph::`+.igraph` for safety reasons. 
#' 
#' @param g igraph object
#' @return an igraph object
addBiases <- function(g){
  # If biases already exist, do nothing
  if(any(V(g)$is.bias)) return(g)
  non.root.nodes <- V(g)[inDegree(g, V(g)) != 0]
  g.new <- g
  for(v in non.root.nodes){
    bias_name <- paste("bias", v, sep="_") # create a name for the bias, eg. "bias_2"
    g.new <- (g.new + bias_name) %>% # Add a vertex to the graph. The `+.igraph` method adds vertices with `+` primitive.  
      initializeVectorAttributesForVertex(bias_name) #Having added the bias, give it the correct properties
    V(g.new)[bias_name]$is.bias <- TRUE #Label the bias as 'bias'
    V(g.new)[bias_name]$output.signal <- list(rep(1, g$n)) #Give the value of 1
    g.new <- g.new + igraph::edge(bias_name, V(g)[v]$name) #same as 'g.new + igraph::edge(bias_name, V(g)[v]$name)' but safer
  }
  g.new
}

#' Initialize boolean vertex attibutes
#' 
#' There are 4 core boolean vertex attributes:
#' \itemize{
#'  \item{is.bias}{TRUE if the vertex is a bias or bias}
#'  \item{is.observed}{TRUE if the vertex is observed in the data}
#'  \item{is.hidden}{TRUE if the vertex is not an bias and not observed in data}
#'  \item{is.root}{TRUE if the vertex is a root in the graph}
#'  \item{is.leaf}{TRUE if the vertex is a leaf in the graph}
#' }
#'  
#' This function initializes all these vertex attributes for later updating.
#' @return g
initializeVertexBooleans <- function(g){
  V(g)$is.bias <- FALSE
  V(g)[grepl('bias', V(g)$name)]$is.bias <- TRUE # grep the biases by name and label them.
  V(g)$is.observed <- FALSE
  V(g)$is.hidden <- FALSE
  V(g)$is.root <- FALSE
  V(g)[get_roots(g)]$is.root <- TRUE
  V(g)$is.leaf <- FALSE 
  V(g)[get_leaves(g)]$is.leaf <- TRUE
  g
}

#' Add vector attributes to graph vertices
#' 
#' signalgraph uses several vertex attributes that store a list of values.  These are
#' called vector attributes, because the set of values is used as a vector in calculations.
#' @param g igraph object
initializeVertexVectors <- function(g) {
  for(v in V(g)) g <- initializeVectorAttributesForVertex(g, v) 
  g
}

#' Initialize the vertex attributes and values
#' 
#' Initializes the boolean attributes of the vectors, 
#'  \enumerate{
#'   \item Initialize the boolean vertex attributes
#'   \item Add bias nodes nodes
#'   \item Add vector attributes to vertices
#'   \item Populate vector attributes with data
#'   \item Reset the 'update' vertex and edge attributes
#' }
#' @param g an igraph object
initializeVertices <- function(g, data){
  g %>%
    addBiases %>%
    initializeVertexBooleans %>%
    initializeVertexVectors %>%
    addDataToVertices(data) %>%
    resetUpdateAttributes
}

#' Initalize edge attributes
#' 
#' Used in contruction of a signalgraph.  Names and initializes the edges.
#' @param g igraph object. The vertices must be named. 
#' @return g igraph object
initializeEdges <- function(g){
  g %>% 
    nameEdges %>% 
    initializeWeights 
}

#' Create an signal graph object that is unfitted
#' 
#' igraph objects have three kinds attributes; graph attributes, edge attributes, and vertex attributes.
#' This function builds a signalgraph object from an igraph object using these attributes.  First the 
#' graph attributes are added, then vertex attributes.  The model takes a data frame as an input.  
#' The name of each variable in the data must match a vertex name in the graph.  The values for a given variable
#' are added as a vertex attribute to that vertex.  Next, edge weights are added as edge attributes.  Finally, 
#' the weights are updated.
#' 
#' @param g igraph object. The vertices must be named. 
#' @param data a data frame. All of the names in the data from must match a vertex name.
#' @param graph_attr list of graph attributes.  Graph attributes include:  
#' \itemize{
#'  \item{penalty}{penalized least squares error penalty parameter value}
#'  \item{activation}{the activation function (this actually is an R function), defaults to logistic.}
#'  \item{activation.prime}{The derivative of the activation function, used in gradient calculation. Defaults to NULL}
#'  \item{min.max.constraints}{2 element numeric containing the acceptable range for each rate.}
#'  }
#' @return A graph with all the attributes needed to fit the neural network model.
initializeGraph <- function(g, data, graph_attr = NULL){
  g %>%
    checkArgs(data) %>% # Check the arguments
    initializeGraphAttributes(graph_attr, nrow(data)) %>% # Add the graph attributes
    initializeVertices(data) %>% # Add biases and vertex attributes
    initializeEdges %>% # Add edge weights
    updateSignals # update vertex values given weights
}