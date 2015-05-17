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
    initializeGraph(logic_gates[, c("I1", "I2")], logic_gates[, gates])
  if(!identical(outputs, "all")){
    output_nodes <- V(g)[name %in% outputs]
    exclusion_nodes <- V(g)[is.observed]  %>% setdiff(output_nodes)
    nuissance_intercepts <- exclusion_nodes %>% # Find the intercepts for the outputs that will be excluded
      lapply(iparents, g=g) %>%
      unlist %>%
      intersect(V(g)[type== "intercept"])
    g <- igraph::induced.subgraph(g, setdiff(V(g), c(exclusion_nodes, nuissance_intercepts)))
  }
  g
}

#' 
#' 
random_sg <- function(k, n){
  generateMultiConnectedDAG(k) %>% nameVertices
  roots <- V(g)[get_roots(g)]$name 
  leaves <- V(g)[get_leaves(g)]$name
  sim_data <- c(lapply(roots, function(root) runif(n)), 
                lapply(leaves, function(leaf) rep(NA, n))) %>%
    `names<-`(c(roots, leaves)) %>%
    data.frame
  # Initialize the graph, os the output.signals are all there.
  g <- initializeGraph(g, data = sim_data)
}
