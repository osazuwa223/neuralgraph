#' K order node selection
#' 
#' Ranks vertices in the signal graph by a given ranking method.  Chooses the K
#' vertices with the top K highest ranking.  Returns the union of the causal 
#' neighborhoods of those K vertices.
#' @export
causal_prioritization <- function(g, k, method = "wtd_betweenness"){
  if(method == "wtd_betweenness"){
    wts <- E(g)$weight
    if(any(wts < 0)) wts <- abs(wts) 
    ranking <- betweenness(g, weights = wts)[V(g)[!is.bias]] %>%
      sort(decreasing = TRUE)
    top_k <- names(ranking[1:k])
    final_set <- V(g)[top_k]$causal_nei %>% unlist %>% unique %>% sort
  } else if(method == "betweenness"){
    ranking <- betweenness(g, weights = rep(1, ecount(g)))[V(g)[!is.bias]] %>%
      sort(decreasing = TRUE)
    top_k <- names(ranking[1:k])
    final_set <- V(g)[top_k]$causal_nei %>% unlist %>% unique %>% sort
  } 
  final_set
}