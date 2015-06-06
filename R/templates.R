#' Test Case: Random graph and matching data
#' 
#' Generates an igraph object and a corresponding dataset.  Used to 
#' randomly generate test cases
#' 
#' @param m the number of desired nodes
#' @param k the number of desired observed nodes
#' @param n the number of desired rows in the data
#' @return A list of two elements, graph and data.
rand_case <- function(m, k, n = m * k + m){
  stopifnot(m > k)
  g <- lucy::generateMultiConnectedDAG(m) %>% nameVertices
  roots <- V(g)[get_roots(g)]$name 
  leaves <- V(g)[get_leaves(g)]$name
  num_roots_leaves <- length(c(roots, leaves))
  if(k < length(c(roots, leaves))){
    warning("k too low, setting k such that only roots and leaves are observed.")
    k <- num_roots_leaves
  }
  data <- c(lapply(roots, function(root) runif(n)), 
            lapply(leaves, function(leaf) runif(n))) %>%
    as.data.frame %>%
    `names<-`(c(roots, leaves))
  if(k > num_roots_leaves){
    observed_middle_nodes <- setdiff(V(g)$name, c(leaves, roots)) %>% 
      sample(k - num_roots_leaves)
    data <- lapply(observed_middle_nodes, FUN = function(item) runif(n)) %>%
      as.data.frame %>%
      `names<-`(observed_middle_nodes) %>%
      cbind(data)
  } 
  list(g = g, data = data)
}

#' Generate a random unfit signalgraph object
#' 
#' Uses the generateMultiConnectedDAG in the lucy \url{https://github.com/osazuwa223/lucy} package.  
#' Since the vertices that must have data are the roots and the leaves, (everything in between can be hidden),
#' then a data frame is simulated for only those nodes.
#' 
#' @param m the number of desired nodes
#' @param k the number of desired observed nodes
#' @param n the number of desired rows in the data
#' @return an initialized, but unfit (unoptimized) signal graph model
#' @export
random_unfit_sg <- function(m, k, n = m * k + m){
  g <- rand_case(m, k, n) %>% 
    {initializeGraph(.$g, .$data, graph_attr = list(activation = logistic))}
}

#' Generate a random unfit signalgraph object
#' 
#' Uses the generateMultiConnectedDAG in the lucy \url{https://github.com/osazuwa223/lucy} package.  
#' Since the vertices that must have data are the roots and the leaves, (everything in between can be hidden),
#' then a data frame is simulated for only those nodes.
#' 
#' @param m the number of desired nodes
#' @param k the number of desired observed nodes
#' @param n the number of desired rows in the data
#' @return a signalgraph model
#' @export
random_sg <- function(m, k, n){
  g <- rand_case(m, k, n) %>% 
      {fitNetwork(.$g, .$data, graph_attr = list(activation = logistic))}
}


#' Create an signal graph model of a multi-layer perceptron logic gate (unfit case).
#' 
#' Creates an signal graph model of a multi-layer perceptron logic gate.  The signal graph weights are not fit.
#' This is used primarily for experiments.
#' 
#' @param outputs character vector describing the desired gates.  Defaults to "all", in which case  
#' (AND, OR, NAND, NOR, XOR, and XNOR) are the outputs.  Otherwise any subset of these will work.
#' @param layers a numeric vector where each element corresponds to a layer, and the value of an element
#' is the number of hidden nodes in the layer.  Default is NULL.  If value is null, no hidden layers will be
#' present (just an input and output layer).
#' @return a signal graph model object with unfit weights
#' 
#' @export
get_gate <- function(outputs = "all", layers=NULL){
  gates <- c("AND", "OR", "NAND", "NOR", "XOR", "XNOR")  
  if(!all(outputs %in% c("all", gates))) stop("Invalid set of gates.")
  logic_gates <- expand.grid(list(I1 = c(0, 1), I2 = c(0, 1))) %>% #Create the logic gate table
    transform(AND = (I1 * I2 == 1) * 1, 
              OR = (I1 + I2 > 0) * 1) %>%
    transform(NAND = (!AND) * 1,
              NOR = (!OR) * 1,
              XOR = (I1 + I2 == 1) * 1, 
              XNOR = (I1 == I2) * 1)
  g <- igraphr::mlp_graph(inputs = c("I1", "I2"),
                    outputs = gates, layers = layers) %>%
    initializeGraph(logic_gates[, c("I1", "I2", gates)])
  if(!identical(outputs, "all")){
    output_nodes <- c(V(g)[c("I1", "I2")], V(g)[name %in% outputs])
    exclusion_nodes <- V(g)[is.observed]  %>% setdiff(output_nodes)
    nuisance_biases <- exclusion_nodes %>% # Find the biases for the outputs that will be excluded
      lapply(lucy::iparents, g=g) %>%
      unlist %>%
      intersect(V(g)[is.bias])
    g <- igraph::induced.subgraph(g, setdiff(V(g), c(exclusion_nodes, nuisance_biases)))
  }
  nameEdges(g) 
}



