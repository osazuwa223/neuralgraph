#' K order node selection
#' 
#' Ranks vertices in the signal graph by a given ranking method.  Chooses the K
#' vertices with the top K highest ranking.  Returns the union of the causal 
#' neighborhoods of those K vertices.
#' @export
causal_prioritization <- function(g, k){
  mse_set_list <- vertexMSEs(g) %>%
    sort(decreasing = TRUE) %>%
    names %>%
    {structure(V(g)[.]$causal_nbr, names = .)} %>%
    lapply(function(set) {set$name})
  if(k > length(mse_set_list)) warning("k is larger than the number of random observed vertices.")
  betweenness_set_list <- g %>% 
    induced.subgraph(V(.)[!is.bias]) %>%
    betweenness(weights = rep(1, ecount(.))) %>%
    sort(decreasing = TRUE) %>% 
    .[. > 0] %>%
    names %>%
    {structure(V(g)[.]$causal_nbr, names = .)} %>%
    lapply(function(set) {set$name})
  match_sets(mse_set_list, betweenness_set_list, k) %>%
    `names<-`(c("mse_set", "btw_set"))
}

#' Intersection of top k sets
#' 
#' \emph{k_order_union} takes a list of sets of items where the order of the list implies rank.
#' It returns a single set corresponding to the union
#' of the first k sets, meaning the sets with the 1st to kth highest rankings. 
#' 
#' \emph{match_sets} takes two lists of sets, and first finds the union set of the top k sets of the first list.
#' Then it finds a union set of the second list such that the two union sets are matched in length.
#' First, it iteratively increases the union set of the second list until the second union set is the same size
#' as the first. If on the jth iteration the size of the second union set is smaller
#' than the first, then the algorithm iterates and the union of the first j + 1 sets is found.
#' If on the jth iteration the second union set is larger than the first, then the second union set 
#' is truncated to be the same size as the first set in a way that emphasizes the difference between the two union sets.  
#' If the second union set is smaller and there is no j + 1 set in the second list, then the first
#' union set is truncated to the same size as the second, again in a way that emphasizes their differences.
#' @param set_list a list of arrays, and it is assumed the ordering of the list
#' corresponds to ranking (eg first item on the list has the highest rank.)
#' @param k the order of the union, i.e. 'find the union of the first k items in the list'
#' @param set_1 first list
#' @param set_2 second list
#' @return an array of list items
k_order_union <- function(set_list, k){
  set_list[1:k] %>%
    unlist %>%
    unique
}
#' @rdname k_order_union
match_sets <- function(list_1, list_2, k){
  union_set_1 <- k_order_union(list_1, k)
  m <- length(list_2)
  for(j in 1:m){
    union_set_2 <- k_order_union(list_2, j)
    if(length(union_set_2) >= length(union_set_1)){
      set_intersection <- intersect(union_set_2, union_set_1)
      set_diff <- setdiff(union_set_2, union_set_1)
      union_set_2 <- c(set_diff, set_intersection)[1:length(union_set_1)]
      break
    } else if(j == m){
      set_intersection <- intersect(union_set_1, union_set_2)
      set_diff <- setdiff(union_set_1, union_set_2)
      union_set_1 <- c(set_diff, set_intersection)[1:length(union_set_2)]
      break
    }
  }
  list(set_1 = union_set_1, set_2 = union_set_2)
}

