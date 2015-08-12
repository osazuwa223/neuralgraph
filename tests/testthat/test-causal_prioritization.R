option <- FALSE
context("Causal Prioritization")
test_that("loadCN does not include roots.",{
  g <- random_unfit_sg(10, 100)
  g <- loadCN(g)
  lapply(V(g), function(v){
    if(!V(g)[v]$is.root){
      example_cn <- V(g)[unlist(V(g)[v]$causal_nbr)]
      expect_true(all(!example_cn$is.root))
    }
  })
})

test_that("causal_selection gets the desired nodes", {
  skip("skipping causal selection")
  g <- random_unfit_sg(10, 100) 
  g <- loadCN(g)
  k <- 1
  betweenness(g, weights = abs(E(g)$weight))[V(g)[!is.bias]] %>%
    sort(decreasing = TRUE) %>% # rank
    .[1:k] %>% # pull the top k
    names %>% # pull the names of the top k
    {V(g)[.]$causal_nbr} %>% # Pull out the causal neighborhood
    unlist %>% # pool them together
    unique %>% # find union
    sort %>% # sort the name
    expect_identical(sort(causal_prioritization(g, 1, method = "wtd_betweenness"))) # compare
  betweenness(g, weights = rep(1, ecount(g)))[V(g)[!is.bias]] %>%
    sort(decreasing = TRUE) %>%
    .[1:k] %>%
    names %>%
    {V(g)[.]$causal_nbr} %>%
    unlist %>%
    unique %>%
    sort %>%
    expect_identical(sort(causal_prioritization(g, 1, method = "betweenness"))) # compare  
})

test_that("k_order_intersection works on trivial examples ", {
  list1 <- list(c(1, 4, 5, 6), c(4, 5, 6, 7), c(7, 8))
  expect_identical(k_order_union(list1, 2), c(1, 4, 5, 6, 7))
  expect_identical(k_order_union(list1, 3), c(1, 4, 5, 6, 7, 8))
})

test_that("match_set works as expect on a trival examples", {
  list1 <- list(c(1, 4, 5, 6), c(4, 5, 6, 7), c(7, 8, 11, 34, 25, 22))
  list2 <- list(c(5, 6), c(4, 6, 7), c(7, 8), c(8, 9))
  match1 <- match_sets(list1, list2, 1)
  match2 <- match_sets(list1, list2, 2)
  match3 <- match_sets(list1, list2, 3)
  expect_identical(match1, list(set_1 = k_order_union(list1, 1), set_2 = c(7, 5, 6, 4)))
  expect_identical(match2, list(set_1 = k_order_union(list1, 2), set_2 = c(8, 5, 6, 4, 7)))
  expect_identical(match_sets(list1, list2, 3), list(set_1 = c(1, 11, 34, 25, 22, 4), 
                                                     set_2 = k_order_union(list2, 4)))
})
