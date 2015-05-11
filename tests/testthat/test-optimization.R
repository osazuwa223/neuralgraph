devtools::load_all("../../R/optimization.R")
#devtools::load_all("R/optimization.R")
library(plyr)
library(dplyr)
context("optimization")

# Testing using Titanic3 data
data(titanic3)
titan <- filter(titanic3, !is.na(age), !is.na(survived), !is.na(fare)) %>% #Not worrying about NA vals
  mutate(survived = as.numeric(survived)) %>%
  select(age, survived, fare) %>%
  rescale_df %$% #Note, everything is rescaled to between 0 and 1
  df

g <- mlp_graph("age", "survived") %>%
  initializeGraph(input.table = select(titan, age), 
                  output.table = select(titan, survived)) %>%
  {induced.subgraph(., V(.)[c("age", "survived")])} %>% # Having removed the intercept, I need to reupdate 
  resetUpdateAttributes %>%
  updateVertices(getDeterminers = iparents, callback = calculateVals)


###########################################################################################
# Test that functions used in optimization are working correctly.
test_that("prediction from logistic function works as expected.", {
  logistic(.5 * V(g)["age"]$output.signal[[1]]) %>% # calculation of prediction w/ weight of 5 via logistic function
  identical(getPrediction(g, V(g)["survived"], .5)) %>%# compared to algorithms generation of prediction
  expect_true
})

test_that("doChainRule in layer-free univarite produces logistic_prime(input * weight) * input", {
  weight <- E(g)[to("survived")]$weight
  age_out <- V(g)["age"]$output.signal %>% unlist
  expected <- logistic_prime(weight * age_out) * age_out 
  chain_rule_output <- doChainRule(g, V(g)["survived"], E(g)[to("survived")])
  expect_equal(expected, chain_rule_output)
})

test_that("getGradient in layer-free univarite produces -(Y - f(input)) * logistic_prime(input * weight) * input", {
  weight <- E(g)[to("survived")]$weight
  age_out <- V(g)["age"]$output.signal %>% unlist
  Y <- V(g)["survived"]$observed %>% unlist
  Y_hat <- V(g)["survived"]$output.signal %>% unlist
  expected <- sum(-(Y - Y_hat) * logistic_prime(weight * age_out) * age_out)
  gradient_output <- getGradientFunction(g, V(g)["survived"])(weight) %>% as.numeric
  expect_equal(expected, gradient_output)
})

g <- mlp_graph(c("age", "fare"), "survived") %>%
  initializeGraph(input.table = select(titan, age, fare), 
                  output.table = select(titan, survived))

test_that("prediction from logistic function works as expected.", {
  linear_combination <- V(g)[c("age", "fare", "int.3")]$output.signal %>% 
    {do.call("cbind", .)} %>%
    `%*%`(matrix(c(-.5, -.5, -.5), ncol = 1)) %>% 
    as.numeric
  logistic(linear_combination) %>% # calculation of prediction w/ weight of 5 via logistic function
    identical(getPrediction(g, V(g)["survived"], c(-.5, -.5, -.5))) %>%# compared to algorithms generation of prediction
    expect_true
})

test_that("hand calculation of doChainRule in layer-free multivarite produces same results as package functions", {
  weights <- E(g)[to("survived")]$weight
  age_out <- V(g)["age"]$output.signal %>% unlist
  fare_out <- V(g)["fare"]$output.signal %>% unlist
  linear_combo <- V(g)[c("age", "fare", "int.3")]$output.signal %>% 
    {do.call("cbind", .)} %>%
    `%*%`(matrix(weights, ncol = 1))
  expected <- cbind(age_weight = logistic_prime(linear_combo) * age_out, 
                  fare_weight = logistic_prime(linear_combo) * fare_out,
                  int.3_weight = logistic_prime(linear_combo))
  chain_rule_output <- sapply(E(g)[to("survived")], doChainRule, g = g, v = V(g)["survived"])
  expect_equal(colSums(expected), colSums(chain_rule_output))
  #######
  Y <- V(g)["survived"]$observed %>% unlist
  Y_hat <- V(g)["survived"]$output.signal %>% unlist
  expected <- apply(chain_rule_output, 2, function(item) -(Y - Y_hat) * item) %>% colSums
  gradient_output <- getGradientFunction(g, V(g)["survived"])(weights) %>% as.numeric
  expect_equal(expected, gradient_output)
})

# Now complex graphs
g <- mlp_graph(c("age", "fare"), "survived", layers = c(2, 2)) %>%
  initializeGraph(input.table = select(titan, age, fare), 
                  output.table = select(titan, survived))

test_that("doChainRule errors in the case when the vertex value does not depend on the weight", {
  expect_error(doChainRule(g, v = V(g)["H21"], e = E(g)["H11" %->% "H22"]),
               "\n  You've attempted to find the gradient of a node's output\n         w.r.t an edge weight that that does not affect that output. v: 5 e: 7\n")
               
})