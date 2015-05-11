#' Rescale the column of a data frame to between 0 and 1
rescale_df <- function(df){
  output.list <- lapply(df, function(col){
    x.max <- max(col)
    x.min <- min(col)
    new.vals <- (col - x.min)/(x.max - x.min)
    list(new.vals = new.vals, x.min = x.min, x.max = x.max)
  })
  new.df <- do.call("cbind", lapply(output.list, function(item) item$new.vals))
  new.df <- as.data.frame(new.df)
  min.max.list <- lapply(output.list, function(item){
    c(item$x.min, item$x.max)
  })
  list(df = new.df, min.max = min.max.list)
}

#' Pull fitted values from a signal graph into a data frame
#' 
#' @param g A fitted graph
#' @return a data frame containing fitted values
#' @seealso recover_design
#' @export
get_fitted <- function(g){
  vertices_of_interest <- V(g)[!(type %in% "intercept")]
  vertices_of_interest$output.signal %>% #Pull everything but the intercepts
{do.call("rbind", .)} %>% t %>% as.data.frame %>% `names<-`(vertices_of_interest$name) # Put into a data frme
}


#' Pull observed data for the inputs and outputs of a signal graph model into a data frame
#' 
#' Pull the observed data (AKA examples) for the input and output vertices of a signal graph model into a data frame.
#' 
#' @param g A fitted graph
#' @return a data frame, essentially the design matrix.
#' @seealso get_fitted
#' @export
recover_design <- function(g){
  c(V(g)[type == "input"]$output.signal, 
    V(g)[type == "output"]$observed) %>%
    data.frame %>% 
    `names<-`(c(V(g)[type == "input"]$name, V(g)[type == "output"]$name)) # Put into a data frme
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
