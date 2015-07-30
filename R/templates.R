#' Test Case: Random graph and matching data
#' 
#' Generates an igraph object and a corresponding dataset.  A number of unobserved vertices will 
#' be chosen at random.  Used to randomly generate test cases.
#' 
#' @param m the number of desired nodes
#' @param n the number of desired rows in the data
#' @return A list of two elements, graph and data.
#' @export
rand_case <- function(m, n = m * m + m, method = "ordered"){
  g <- lucy::sim_DAG(m, method = method) %>% name_vertices
  roots <- V(g)[get_roots(g)]$name 
  leaves <- V(g)[get_leaves(g)]$name
  num_roots_leaves <- length(c(roots, leaves))
  k <- sample(num_roots_leaves:m, 1) # Number of observed nodes
  .data <- c(lapply(roots, function(root) runif(n)), 
            lapply(leaves, function(leaf) runif(n))) %>%
    as.data.frame %>%
    `names<-`(c(roots, leaves))
  if(k > num_roots_leaves){
    observed_middle_nodes <- setdiff(V(g)$name, c(leaves, roots)) %>% 
      sample(k - num_roots_leaves)
    .data <- lapply(observed_middle_nodes, FUN = function(item) runif(n)) %>%
      as.data.frame %>%
      `names<-`(observed_middle_nodes) %>%
      cbind(.data)
  } 
  list(g = g, data = .data)
}

#' Generate a random unfit signalgraph object
#' 
#' Uses the sim_DAG in the lucy \url{https://github.com/osazuwa223/lucy} package.  
#' Since the vertices that must have data are the roots and the leaves, (everything in between can be hidden),
#' then a data frame is simulated for only those nodes.
#' 
#' @param m the number of desired nodes
#' @param n the number of desired rows in the data
#' @param num_fixed the number of vertices treated as fixed variables.
#' @param no_fixed boolean if TRUE then all vertices are treated as random. Defaults to FALSE
#' @param no_fixed_prob probability of generating a case with no fixed variable nodes
#' @return an initialized, but unfit (unoptimized) signal graph model
#' @export
random_unfit_sg <- function(m, n = m * m + m, no_fixed_prob = .3, ...){
  g <- rand_case(m, n)
  no_fixed <- sample(c(TRUE, FALSE), 1, prob = c(no_fixed_prob, 1 - no_fixed_prob))
  if(no_fixed){
    return(initializeGraph(g$g, g$data, ...))
  }
  roots <- V(g$g)[get_roots(g$g)]$name
  initializeGraph(g$g, g$data, fixed = roots, ...)
}

#' Generate a random fitted signalgraph object
#' 
#' Uses the sim_DAG in the lucy \url{https://github.com/osazuwa223/lucy} package.  
#' The vertices that must have data are the roots and the leaves, (everything in between can be hidden),
#' so data is simulated for only those nodes.
#' 
#' @param m the number of desired nodes
#' @param n the number of desired rows in the data
#' @param ... additional arguments, including graph attributes
#' @param no_fixed boolean if TRUE then all vertices are treated as random. Defaults to FALSE
#' @return a signalgraph object
#' @export
random_sg <- function(m, n, max.iter = 1, no_fixed = FALSE,...){
  g <- rand_case(m, n)
  if(no_fixed){
    return(fitNetwork(g$g, g$data, ...))
  }
  roots <- V(g$g)[get_roots(g$g)]$name
  fitNetwork(g$g, g$data, fixed = roots, ...)
}

#' Generate a random signalgraph object for data simulation
#' 
#' Produces a signalgraph objectwhere the values of the observed data and the fitted data are the same.  
#' This is designed to produce a gold standard, which can simulate data, where upon a new model can be 
#' fit on the simulated data, and the parameters of the standard and learned parameters of the new model
#' can be compared.
#' 
#' @param m the number of desired nodes
#' @param n the number of desired rows in the data
#' @param error_sd desired standard deviation for Gaussian error of logit(x)
#' @param ... arguments past to fitNetwork, including graph attributes
#' @return a signalgraph object
#' @export
sim_system <- function(m, n, error_sd = .2, ...){
  g <- rand_case(m, n) %>%
    {initializeGraph(.$g, .$data, fixed = get_roots(.$g), ...)} %>%
    loadCN
  fitted_vals <- get_fitted(g)
  logit <- function(x) log(x / (1+x))
  for(node in names(fitted_vals)){ 
    if(V(g)[node]$is.observed){
      observed_val <- fitted_vals[, node] %>%
        logit %>%
        {1 / (1 + exp(-1 * (. + rnorm(n, sd = error_sd))))}
      V(g)[node]$observed <- list(observed_val)
    } 
  }
  g
}

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
#
calculateVals <- function(g, v){ 
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

checkVector <- function(item){
  valid <- TRUE
  if(any(is.na(item))) valid <- FALSE
  if(any(is.infinite(item))) valid <- FALSE
  if(length(item) == 0) valid <- FALSE
  valid
}

update_signals <- function(g){
  update_vertices(g, getDeterminers = lucy::iparents, callback = calculateVals)
}

#' Simulate data from a signal graph
#' 
#' Simulates new data for all the fixed variables and propagates those changes
#' through the rest of the system.
#' @param g a signal graph object
#' @param n desired number of rows in output dataframe
#' @param add_error if true, adds noise to each column in output dataframe using the jitter function.
#' @return a dataframe of data
#' @export
sim_system_data <- function(g, n, add_error = FALSE){
  prediction_graph <- resetUpdateAttributes(g)
  fixed_v <- V(prediction_graph)[is.fixed]
  # First sim new values for the output signal
  prediction_graph$n <- n
  for(v in fixed_v){
    V(prediction_graph)[v]$output.signal <- list(runif(n))
  }
  prediction_graph <- update_signals(prediction_graph) 
  # Pull the values into a data 
  .data <- get_fitted(prediction_graph)
  if(add_error) {
    .data <- lapply(.data, jitter) %>%
      as.data.frame %>%
      `names<-`(names(.data))
  }
  .data
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
    initializeGraph(logic_gates[, c("I1", "I2", gates)], fixed = c("I1", "I2"))
  if(!identical(outputs, "all")){
    output_nodes <- c(V(g)[c("I1", "I2")], V(g)[name %in% outputs])
    exclusion_nodes <- V(g)[is.observed]  %>% setdiff(output_nodes)
    nuisance_biases <- exclusion_nodes %>% # Find the biases for the outputs that will be excluded
      lapply(lucy::iparents, g=g) %>%
      unlist %>%
      intersect(V(g)[is.bias])
    g <- igraph::induced.subgraph(g, setdiff(V(g), c(exclusion_nodes, nuisance_biases)))
  }
  name_edges(g) 
}