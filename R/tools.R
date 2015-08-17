#' The logistic function
#' The logistic function and its derivative (used in gradient calculations).
#' @export
logistic <- function(z) 1 / (1 + exp(-z))
#' @rdname logistic
#' @export 
logistic_prime <- function(z) logistic(z) * (1 - logistic(z))

#' MSE calculation
#' 
#' vertexMSE calculates the mean squared error for a given fitted variables in a signalgraph 
#' object. vertexMSEs calculates MSE for each fitted vertices and returns an array of values.
#' getMSE calculates the mean squared error for all the fitted variables in a signalgraph 
#' object, then reports the mean of means.
#' 
#' @param g a signal graph object
#' @param v a vertex id in the signal graph
#' @export 
getMSE <- function(g){
  response_variables <- intersect(V(g)[is.random], V(g)[is.observed])
  k <- length(V(g)[response_variables])
  observed <- unlist(V(g)[response_variables]$observed)
  prediction <- unlist(V(g)[response_variables]$output.signal)
  sum((observed - prediction) ^ 2) / g$n / k
}
#' @rdname getMSE
#' @export
vertexMSE <- function(g, v){
  if(!V(g)[v]$is.random) stop("vertex is not a random variable.")
  (sum(unlist(V(g)[v]$observed) - unlist(V(g)[v]$output.signal))^2) / g$n
}
#' @rdname getMSE
#' @export
vertexMSEs <- function(g){
  observed_and_random <- intersect(V(g)[is.observed], V(g)[is.random])
  lapply(V(g)[observed_and_random], function(v) vertexMSE(g, v)) %>%
    unlist %>%
    structure(names = V(g)[observed_and_random]$name)
}

#' Rescale the column of a data frame to between 0 and 1
rescale_df <- function(df){
  observed.list <- lapply(df, function(col){
    x.max <- max(col)
    x.min <- min(col)
    new.vals <- (col - x.min)/(x.max - x.min)
    list(new.vals = new.vals, x.min = x.min, x.max = x.max)
  })
  new.df <- do.call("cbind", lapply(observed.list, function(item) item$new.vals))
  new.df <- as.data.frame(new.df)
  min.max.list <- lapply(observed.list, function(item){
    c(item$x.min, item$x.max)
  })
  list(df = new.df, min.max = min.max.list)
}

#' Get a data frame of output signals   
#' @param g A fitted graph
#' @return a data frame containing fitted values
#' @seealso recover_design
#' @export
get_fitted <- function(g){
  V(g)[!is.bias]$output.signal %>%
    data.frame %>% 
    `names<-`(V(g)[!is.bias]$name)
}

#' Pull observed data from a signalgraph object
#' 
#' Pull the observed data (AKA examples) for the observed vertices in a signal graph model into a data frame.
#' 
#' @param g A fitted graph
#' @return a data frame, essentially the design matrix.
#' @seealso get_fitted
#' @export
recover_design <- function(g){
  V(g)[is.observed]$observed %>%
    data.frame %>% 
    `names<-`(V(g)[is.observed]$name) # Put into a data frme
}

#' A helper function used in summarizing signal graph structures
#' @seealso examine_signal_graph
format_vertex_list <- function(output){
  if(is.list(output)){
    if(length(output) > 1){
      output.list.item <- as.data.frame(do.call("cbind", lapply(output, head)))
      if(ncol(output.list.item) == vcount(g)){
        names(output.list.item) <- V(g)
      }
      output <- list(output.list.item)
    } else {
      output <- list(head(unlist(output)))
    }
  }
  output
}

#' Summarize the values of an signal graph model
#' 
#' A summary function that prints out the values of various graph, edge, and vertex attributes in formatted lists.
#'  
#' @param g a signal graph object
#' @export   
examine_signal_graph <- function(g){
  igraphr::examineGraph(g, formatVertexAttr=format_vertex_list)
}

#' Perform long tests
#'  
#' Select whether or not to perform more time consuming tests by setting opt to TRUE or FALSE
long_test <- function(option){
  if(!option){
    eval(quote(testthat::skip("Skipping time consuming test.")), parent.frame())
  } 
}

#' Visualize a signal graph object
#' 
#' Fixed nodes are displayed as orange.  Of the two kinds of random nodes, 
#' observed nodes are displayed in green, hidden nodes are displayed in blue.
#' For random observed nodes, the width of node borders corresponds to error.
#' 
#' @export
sg_viz <- function(g, main = NULL, sub = NULL, show_biases = FALSE){
  if(!show_biases)  g <- igraph::induced.subgraph(g, V(g)[!is.bias])
  fill <- structure(rep("white", vcount(g)), names = V(g)$name)
  fill[V(g)$is.observed] <- "light green"
  fill[V(g)$is.root] <- "dark orange"
  fill[V(g)$is.hidden] <- "blue"
  if("is.bias" %in% list.vertex.attributes(g)) fill[V(g)$is.bias] <- "grey"
  col <- structure(rep("black", vcount(g)), names = V(g)$name)
  observed_and_random <- intersect(V(g)[is.observed], V(g)[is.random])
  col[observed_and_random] <- "dark red"
  lwd <- structure(rep(1, vcount(g)), names = V(g)$name)
  mses <- vertexMSEs(g) 
  lwd_vals <- 8 * (mses / max(mses))^3 + 1
  lwd[observed_and_random] <- lwd_vals
  node_list <- list(fill = fill, col = col, lwd = lwd) 
  g_out <- g %>% name_vertices %>% # Give vertices names if they do not have nay 
    igraph.to.graphNEL(.) %>% # convert to a graphNEL
    {Rgraphviz::layoutGraph(.)} %>% # lay the graph out
    {graph::`nodeRenderInfo<-`(., node_list)} # add the node annotation
  # if(!is.null(main)) graph::graph.par(list(graph = list(main = main))) # Add a title if one is given
  graph::graph.par(list(graph = list(main = main, sub = sub))) # Add a title if one is given
  Rgraphviz::renderGraph(g_out) # Render the graph
}

#' Convert logistic activation parameters to parameters in u / (1 + u) function
logistic_to_positive <- function(g, v){
  weight_vector <- matrix(E(g)[to(v)]$weight, ncol=1)
  parents <- iparents(g, v)
  parent_mat <- parents %>%
    ensure_that(length(.) > 0) %>%
    {do.call("cbind", V(g)[.]$output.signal)}
  box_mat <- t(parent_mat) %*% parent_mat
  inv_box <- solve(box_mat)
  linear_combination <-  parent_mat %>%
    `%*%`(weight_vector) %>% 
    as.numeric %>%
    ensure_that(checkVector(.))
  new_weights <- inv_box %*% t(parent_mat) %*% exp(linear_combination) %>%
    as.numeric %>%
    structure(names = V(g)[parents]$name) %>%
    ensure_that(. > 0)
  new_weights
}

#' Pull igraph structure from a signal graph
#' @param g a signalgraph object
#' @export
get_structure <- function(g){
  g %>%
    {induced.subgraph(., V(g)[!is.bias])} %>%
    get.edgelist %>% # Pull out the edges
    graph.edgelist # Add them back in
}




