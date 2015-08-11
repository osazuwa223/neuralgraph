#' Simulate a scale-free signal graph structure from an input structure
#' 
#' Takes an input graph structure, fits a power law model, and simulates a new
#' graph based on the fit.  This function uses the Barabasi-Albert model to 
#' generate graphs implemented in lucy::power_law_sim, however an alteration is made.  
#' The edges of the input graph are reversed before estimating the power law fit, and the 
#' simulated graph edges are reversed again before it is returned from the function.  
#' This way, the graph is characterized not by preferential attachment of incoming edges, but
#' by preferential attachment of outgoing edges.  Further, vertices with 0 degrees are avoided 
#' by having high degree vertices 'donate' edges to degreeless vertices.
#' @param g an igraph object.  
#' @param p number of desired vertices in the output graph.
#' @return A new graph simulated from a power law based on g.
#' @export 
power_signal_graph <- function(g, p){
  g_out <- g %>%
    reverse_edges %>%
    power_sim_no_singleton(p) %>%
    simplify %>% 
    reverse_edges
  V(g_out)$name <- V(g_out) %>% as.numeric %>% as.character
  g_out
}

#' Singleton-free power law sim
#' The objective is to simulate graphs with the same out-degree distribution as an input graph using
#' the Barabasi-Albert algorithm.  However, if there are nodes in the input graph that have no 
#' out-degree, then there is a non-zero probability of simulating a it is possible the simulation
#' will produce a graph with singletons -- nodes with no edges.  To algorithm performs the 
#' simulation normally, the removes the singletons from the graph. Using the resulting graph as a 
#' starting graph, it continues simulating new nodes and edges but now with 0 probability of a node
#' having no outgoing edges, until the desired amount of  nodes is reached.  The output graph will
#' generally not have the same out-degree distribution, but will be close.  The lower the probability 
#' of having 0 outgoing edges, the closer they will be. 
#' @param g the input igraph object. 
#' @param p the desired amount of vertices
#' @return a simulated graph based on g.
power_sim_no_singleton <- function(g, p){
  gamma <- fit_barabasi_power(g)
  out_degree_dist <- igraph::degree.distribution(g, mode = "out")
  sim_graph <- igraph::barabasi.game(p, power = gamma, out.dist = out_degree_dist) 
  singletons <- V(sim_graph)[igraph::degree(sim_graph) == 0]
  if(length(singletons) > 0){
    start_graph <- sim_graph - singletons
    out_degree_dist[1] <- 0
    sim_graph <- igraph::barabasi.game(p, power = gamma, out.dist = out_degree_dist,
                                       start.graph = start_graph)
  }
  sim_graph %>%
    ensure_that(is.dag(.)) %>%
    ensure_that(is.simple(.))
}

#' Test Case: Random graph and matching data
#' 
#' Generates an igraph object and a corresponding dataset.  A number of unobserved vertices will 
#' be chosen at random.  Used to randomly generate test cases.  If supplied an input graph structure,
#' the new structure is simulated using a power law.  Otherwise, the function uses the sim_DAG in 
#' the lucy \url{https://github.com/robertness/lucy} package.
#' 
#' @param p the number of desired nodes
#' @param n the number of desired rows in the data
#' @param input_g an igraph object.  If supplied then the graph is simulated from
#' a power law fit on this input graph. 
#' @param method the fitting method for simulating a directed acyclic graph.  Ignored
#' if input_g is supplied. 
#' @return A list of two elements, graph and data.
#' @seealso power_signal_graph
#' @seealso sim_DAG
#' @export
rand_case <- function(p, n = p * p + p, input_g = NULL, method = "ordered"){
  if(is.null(input_g)){
    g <- lucy::sim_DAG(p, method = method)
  } else {
    g <- power_signal_graph(input_g, p)
  }
  # Generation of data matching the graph
  roots <- V(g)[get_roots(g)]$name 
  leaves <- V(g)[get_leaves(g)]$name
  num_roots_leaves <- length(c(roots, leaves))
  k <- sample(num_roots_leaves:p, 1) # Number of observed nodes
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
  fixed = V(g)[get_roots(g)]$name
  list(g = g, data = .data)
}

#' Generate a random unfit signalgraph object
#'   
#' Since the vertices that must have data are the roots and the leaves, (everything in between can be hidden),
#' then a data frame is simulated for only those nodes.
#' 
#' @param p the number of desired nodes
#' @param n the number of desired rows in the data
#' @param input_g an igraph object.  If supplied then the graph is simulated from
#' a power law fit on this input graph. 
#' @param method the fitting method for simulating a directed acyclic graph.  Ignored
#' if input_g is supplied. 
#' @param no_fixed_prob probability of generating a case with no fixed variable nodes
#' @return an initialized, but unfit (unoptimized) signal graph model
#' @seealso power_signal_graph
#' @seealso sim_DAG
#' @seealso rand_case
#' @export
random_unfit_sg <- function(p, n = p * p + p, input_g = NULL, method = "ordered", no_fixed_prob = .3, ...){
  if(is.null(input_g)){
    case <- rand_case(p, n, method = "method")
  } else {
    case <- rand_case(p, n, input_g = input_g)
  }
  no_fixed <- sample(c(TRUE, FALSE), 1, prob = c(no_fixed_prob, 1 - no_fixed_prob))
  if(no_fixed){
    return(initializeGraph(case$g, case$data, ...))
  }
  roots <- V(case$g)[get_roots(case$g)]$name
  initializeGraph(case$g, case$data, fixed = roots, ...)
}

#' Generate a random fitted signalgraph object
#' 
#' Uses the sim_DAG in the lucy \url{https://github.com/robertness/lucy} package.  
#' The vertices that must have data are the roots and the leaves, (everything in between can be hidden),
#' so data is simulated for only those nodes.
#' 
#' @param p the number of desired nodes
#' @param n the number of desired rows in the data
#' @param ... additional arguments, including graph attributes
#' @param input_g an igraph object.  If supplied then the graph is simulated from
#' a power law fit on this input graph. 
#' @param method the fitting method for simulating a directed acyclic graph.  Ignored
#' if input_g is supplied. 
#' @param no_fixed boolean if TRUE then all vertices are treated as random. Defaults to FALSE
#' @return a signalgraph object
#' @export
random_sg <- function(p, n, max.iter = 1, input_g = NULL, method = "ordered", no_fixed = FALSE,...){
  if(is.null(input_g)){
    case <- rand_case(p, n, method = "method")
  } else {
    case <- rand_case(p, n, input_g = input_g)
  }
  if(no_fixed){
    return(fitNetwork(case$g, case$data, ...))
  }
  roots <- V(case$g)[get_roots(case$g)]$name
  fitNetwork(case$g, case$data, fixed = roots, ...)
}

#' Generate a random signalgraph object for data simulation
#' 
#' Produces a signalgraph objectwhere the values of the observed data and the fitted data are the same.  
#' This is designed to produce a gold standard, which can simulate data, where upon a new model can be 
#' fit on the simulated data, and the parameters of the standard and learned parameters of the new model
#' can be compared.
#' 
#' @param p the number of desired nodes
#' @param n the number of desired rows in the data
#' @param input_g an igraph object.  If supplied then the graph is simulated from
#' a power law fit on this input graph. 
#' @param method the fitting method for simulating a directed acyclic graph.  Ignored
#' if input_g is supplied. 
#' @param ... arguments past to fitNetwork, including graph attributes
#' @return a signalgraph object
#' @seealso power_signal_graph
#' @seealso sim_DAG
#' @seealso rand_case
#' @export
sim_system <- function(p, n, input_g = NULL, method = "ordered", ...){
  if(is.null(input_g)){ 
    case <- rand_case(p, n, method = method)
  } else { 
    case <- rand_case(p, n, input_g = input_g)
  }
  g <- case %>% 
    {initializeGraph(.$g, .$data, fixed = get_roots(.$g), ...)} %>%
    loadCN
  fitted_vals <- get_fitted(g)
  for(node in names(fitted_vals)){ 
    if(V(g)[node]$is.observed){
      observed_val <- fitted_vals[, node]
      V(g)[node]$observed <- list(observed_val)
    } 
  }
  g
}

#' Generate data from a simulated signal graph
#' Pull data frame from a simulated signal graph and add error to simulate realistic data.
#' @param g the output of sim_system
#' @param factor numeric quantifying the amount of error to add to each random element in the 
#' data frame (see factor argument in base function 'jitter').
#' @return a data frame that can be used as data from the simulated system.
#' @export
sim_data_from_system <- function(g, factor){
  observed_random <- intersect(V(g)[is.observed], V(g)[is.random])
  recover_design(g) %>%
    add_error(V(g)[observed_random]$name, factor = 1000)
}

#' 
#' Regenerate the data values in a signal graph
#' 
#' Simulates new data for all the fixed variables and propagates those changes
#' through the rest of the system.
#' @param g a signal graph object
#' @param n desired number of rows in output dataframe
#' @return a dataframe of data
#' @export
regenerate_graph_data <- function(g, n){
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

#' Add Gaussian error
#' 
#' In the creation of a data frame from a simulated system, 
#' this adds Gaussian error to the variables of a table and rescales them 
#' to between 0 and 1.
#' @param .data a data frame
#' @param vars names or indices of the data.frame to add the error. 
#' @param noise the desired amount of standard error in the noise to add to the variables.
#' Note that after noise is added the var is rescaled to 0 and 1.
add_gauss_error <- function(.data, vars, noise){
  for(v in vars){
    .data[, v] <- .data[, v] + rnorm(nrow(.data), sd = noise) %>%
      {(. - min(.)/(max(.) - min(.)))}
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

#' Add error to a specific columns in a table
#' @param df a data frame 
#' @param vars the names or indices of columns 
#' @param factor numeric for factor argument in jitter function
add_error <- function(df, vars, factor = 1){
  for(v in vars){
    df_v <- df[, v]
    df_v <- jitter(df_v, factor)
    df[, v] <- (df_v - min(df_v)) / (max(df_v) - min(df_v)) 
  } 
  df
}