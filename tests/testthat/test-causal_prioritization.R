context("causal prioritization")
devtools::load_all("../../R/optimization.R")
devtools::load_all("../../R/tools.R")
option <- FALSE
#devtools::load_all("R/optimization.R")
#devtools::load_all("R/templates.R")
test_that("causal_selection gets the desired nodes", {
  g <- random_unfit_sg(10, 100) 
  g <- loadMB(g)
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
  
  