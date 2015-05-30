

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
#' @export
sg_viz <- function(g, main = NULL, show_biases = TRUE){
  if(!show_biases)  g <- igraph::induced.subgraph(g, V(g)[!is.bias])
  col <- structure(rep("white", vcount(g)), names = V(g)$name)
  col[V(g)$is.observed] <- "light green"
  col[V(g)$is.hidden] <- "blue"
  if("is.bias" %in% list.vertex.attributes(g)) col[V(g)$is.bias] <- "grey"
  node_list <- list(fill = col) 
  g_out <- g %>% nameVertices %>% # Give vertices names if they do not have nay 
    igraph.to.graphNEL(.) %>% # convert to a graphNEL
    {Rgraphviz::layoutGraph(.)} %>% # lay the graph out
    {graph::`nodeRenderInfo<-`(., node_list)} # add the node annotation
  if(!is.null(main)) graph::graph.par(list(graph = list(main = main))) # Add a title if one is given
  Rgraphviz::renderGraph(g_out) # Render the graph
  graph::graph.par(list(graph = list(main = ""))) # Reset graph parameters
}

