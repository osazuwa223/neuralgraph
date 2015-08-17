#' Load Markov Blankets
#' 
#' A Markov neighborhood for a given feature is the set of feature that includes that feature itself, and
#' all the nodes with which it shares a causal relationship.  If causality is given by a DAG, then it is
#' equivilent to the node and its Markov blanket.  This function adds this latter Markov blanket + identity 
#' based Markov neighborhood.
#' @param g a signalgraph
#' @export
loadCN <- function(g){
  g <- name_vertices(g)
  va <- list.vertex.attributes(g)
  if("is.root" %in% va){
    g_no_root <- induced.subgraph(g, V(g)[!is.root])
    for(v in V(g_no_root)){
      v_name <- V(g_no_root)[v]$name
      cn_name <- V(g_no_root)[c(imb(g_no_root, v), v)]$name
      V(g)[v_name]$causal_nbr <- list(V(g)[cn_name]) 
    }
  } else {
    for(v in V(g)){
      V(g)[v]$causal_nbr <- list(V(g)[c(imb(g), v)])  
    }
  }
  g
}

#' Check validity of arguments
#' 
#' Check validity of both inputs to constuction of signalgraph object
#' @param g an igraph object, vertices should be named.
#' @param data a dataframe, each variable name should match a vertex name in the graph
#' @return the input graph, ready for passing to a subsequent function
checkArgs <- function(g, data, fixed){
  if(length(setdiff(fixed, names(data)) > 0)) stop("Specified fixed variables not found in the data")
  if(!is.directed(g)) stop("Signal graph requires a directed graph")
  if(vcount(g) < 2 || ecount(g) == 0) stop("There has to be at least 2 vertices and 1 directed edge.")
  if(is.null(V(g)$name)) stop("Vertices must be named")
  if(!all(names(data) %in% V(g)$name)) stop("Data contains variables that are not named in the graph.")
  leaves <- V(g)[get_leaves(g)]$name
  roots <- V(g)[get_roots(g)]$name
  if("is.bias" %in% list.vertex.attributes(g)) roots <- setdiff(roots, V(g)[is.bias]$name)
  if(length(intersect(roots, leaves)) != 0) stop("Detected at least one vertex that is both a root and a leaf. Perhaps there is an unconnected vertex?")
  if(length(setdiff(leaves, names(data))) > 0) stop("Graph leaves must be observed in the data.")  
  basic_attributes <- c("activation", "min.max.constraints", "n", "L1_pen", "L2_pen")
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
    if(!is.null(g$prop_sparse)){
      stopifnot(g$prop_sparse > 0 || g$prop_sparse < 1)
      prop_non_zero <- 1 - g$prop_sparse
      E(g)$weight <-  rep(0, ecount(g))
      for(i in 1:ecount(g)){
        E(g)[i]$path <- list(NULL)
      }
      non_zero_count <- ceiling(prop_non_zero * ecount(g))
      non_zero_edge_index <- E(g)[sample(E(g), non_zero_count)]
      E(g)[non_zero_edge_index]$weight <- rnorm(non_zero_count, sd = 3)
    } else{
      E(g)$weight <- rnorm(ecount(g), sd = 3)
    }
  } else {
    e_min <- g$min.max.constraints[1]
    e_max <- g$min.max.constraints[2]
    if(!is.null(g$prop_sparse)){ # If the prop_sparse attribute is present, then prop_sparse many edges will be initialized at 0
      stopifnot(g$prop_sparse > 0 || g$prop_sparse < 1)
      prop_non_zero <- 1 - g$prop_sparse
      E(g)$weight <-  rep(0, ecount(g))
      non_zero_count <- ceiling(prop_non_zero * ecount(g))
      non_zero_edge_index <- E(g)[sample(E(g), non_zero_count)]
      E(g)[non_zero_edge_index]$weight <- runif(non_zero_count, min = g$min.max.constraints[1], max = g$min.max.constraints[2])
    } else {
      E(g)$weight <- rnorm(ecount(g), sd = 3)
      E(g)$weight[E(g)$weight < g$min.max.constraints[1]] <- g$min.max.constraints[1]
      E(g)$weight[E(g)$weight > g$min.max.constraints[2]] <- g$min.max.constraints[2]
    }  
  }
  for(i in 1:ecount(g)){
    E(g)[i]$path <- list(E(g)[i]$weight)
  }
  g
}

#' Initialize graph attributes of a graph object
#' 
#' Signalgraph stores key parameters used in fitting and other operations on the signal graph as igraph 
#' attributes.  igraph has 3 kinds of attributes, edge attributes, vertex attributes, and graph attributes.
#' This function sets graph attributes.  Key examples include:
#'  \itemize{
#'  \item{L1_pen}{penalized least squares error L1 penalty parameter value}
#'  \item{L2_pen}{penalized least squares error L2 penalty parameter value}
#'  \item{activation}{the activation function (this actually is an R function)}
#'  \item{activation.prime}{The derivative fo the activation function, used in gradient calculation.}
#'  \item{min.max.constraints}{2 element numeric containing the acceptable range for each rate.}
#'  }
#' @param g an igraph object not yet initialized as a signal graph
#' @param graph_attr a list of objects to be used as graph attributes
#' @param n number of rows in the data
initializeGraphAttributes <- function(g, graph_attr, n){
  basic_attributes <- c("activation", "min.max.constraints", "n", "L1_pen", "L2_pen")
  for(attrib in names(graph_attr)) g <- set.graph.attribute(g, attrib, graph_attr[[attrib]])
  graph_attributes <- list.graph.attributes(g)
  if(!("L1_pen" %in% graph_attributes)) g$L1_pen <- 0 # parameter of 0 unless specified.
  if(!("L2_pen" %in% graph_attributes)) g$L2_pen <- 0 # parameter of 0 unless specified.
  if(!("activation" %in% graph_attributes)) g$activation <- logistic
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
addDataToVertices <- function(g, data, fixed){
  for(item in names(data)){
    data_list <- list(data[, item])
    V(g)[item]$observed <- data_list
    V(g)[item]$is.observed <- TRUE 
    if(V(g)[item]$is.root) V(g)[item]$output.signal <- data_list
  }
  for(item in fixed){
    V(g)[fixed]$is.fixed <- TRUE
  }
  nonrandom <- c(V(g)[is.bias], V(g)[is.fixed]) #b iases and fixed variables are non-random
  random <- setdiff(V(g), nonrandom) # everything else is random
  V(g)[random]$is.random <- TRUE # Every thing that is not fixed and is not random is biased.
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
#' @param fixed names of fixed variables in the vertices
#' @return an igraph object
addBiases <- function(g, fixed){
  # If biases already exist, do nothing
  if(any(V(g)$is.bias)) return(g)
  #non_root_nodes <- V(g)[inDegree(g, V(g)) != 0]
  non_root_nodes <- setdiff(V(g)$name, fixed)
  g.new <- g
  for(v in non_root_nodes){
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
#' @return signal graph object
initializeVertexBooleans <- function(g){
  V(g)$is.bias <- FALSE
  V(g)[grepl('bias', V(g)$name)]$is.bias <- TRUE # grep the biases by name and label them.
  V(g)$is.observed <- FALSE
  V(g)$is.random <- FALSE
  V(g)$is.fixed <- FALSE
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
#' @param fixed character array, vertex names of fixed variables 
initializeVertices <- function(g, data, fixed){
  g %>%
    addBiases(fixed) %>%
    initializeVertexBooleans %>%
    initializeVertexVectors %>%
    addDataToVertices(data, fixed) %>%
    resetUpdateAttributes
}

#' Initalize edge attributes
#' 
#' Used in contruction of a signalgraph.  Names and initializes the edges.
#' @param g igraph object. The vertices must be named. 
#' @return g igraph object
initializeEdges <- function(g){
  g %>% 
    name_edges %>% 
    initializeWeights 
}

#' Create an signal graph object that is unfitted
#' 
#' igraph objects have three kinds attributes; graph attributes, edge attributes, and vertex attributes.
#' This function builds a signalgraph object from an igraph object using these attributes.  First the 
#' graph attributes are added, then vertex attributes.  The model takes a data frame as an input. Fixed
#' variables have to be named in the fixed argument, or else variables will be considered random. 
#' The name of each variable in the data must match a vertex name in the graph.  The values for a given variable
#' are added as a vertex attribute to that vertex.  Next, edge weights are added as edge attributes.  Finally, 
#' the weights are updated.
#' 
#' @param g igraph object. The vertices must be named. 
#' @param data a data frame. All of the names in the data from must match a vertex name.
#' @param fixed character array, vertex names of fixed variables in data. Defaults to NULL meaning all variables
#' in the data are treated as random.
#' @param graph_attr list of graph attributes.  Graph attributes include:  
#' \itemize{
#'  \item{L1_pen}{penalized least squares error L1 penalty parameter value}
#'  \item{L2_pen}{penalized least squares error L2 penalty parameter value}
#'  \item{activation}{the activation function (this actually is an R function), defaults to logistic.}
#'  \item{activation.prime}{The derivative of the activation function, used in gradient calculation. Defaults to NULL}
#'  \item{min.max.constraints}{2 element numeric containing the acceptable range for each rate.}
#'  }
#' @return A graph with all the attributes needed to fit the neural network model.
#' @export
initializeGraph <- function(g, data, fixed = NULL, graph_attr = NULL){
  g %>%
    checkArgs(data, fixed) %>% # Check the arguments
    initializeGraphAttributes(graph_attr, nrow(data)) %>% # Add the graph attributes
    initializeVertices(data, fixed) %>% # Add biases and vertex attributes
    initializeEdges %>% # Add edge weights
    update_signals # update vertex values given weights
}