devtools::load_all("../../R/optimization.R")
devtools::load_all("../../R/tools.R")
#devtools::load_all("R/optimization.R")
library(plyr)
library(dplyr)
context("optimization")

# Testing using Titanic3 data
data(titanic3)
titan <- dplyr::filter(titanic3, !is.na(age), !is.na(survived), !is.na(fare)) %>% #Not worrying about NA vals
  {dplyr::mutate(., survived = as.numeric(survived))} %>%
  {dplyr::select(., age, survived, fare)} %>%
  rescale_df %$% #Note, everything is rescaled to between 0 and 1
  df

###########################################################################################
# Univariate no intercept, no hidden layer case
g <- mlp_graph("age", "survived") %>%
  initializeGraph(input.table = select(titan, age), 
                  output.table = select(titan, survived)) %>%
  {induced.subgraph(., V(.)[c("age", "survived")])} %>% # Having removed the intercept, I need to reupdate 
  resetUpdateAttributes %>%
  updateVertices(getDeterminers = iparents, callback = calculateVals)


test_that("prediction from logistic function works as expected.", {
  logistic(.5 * V(g)["age"]$output.signal[[1]]) %>% # calculation of prediction w/ weight of 5 via logistic function
  identical(getPrediction(g, V(g)["survived"], .5)) %>%# compared to algorithms generation of prediction
  expect_true
})

###########################################################################################
# Multivariate with hidden layers
g <- mlp_graph(c("age", "fare"), "survived", layers = c(3, 2)) %>%
  initializeGraph(input.table = select(titan, age, fare), 
                  output.table = select(titan, survived))

test_that("prediction from logistic function works as expected.", {
  linear_combination <- V(g)[c("age", "fare", "int.3")]$output.signal %>% 
    {do.call("cbind", .)} %>%
    `%*%`(matrix(c(-.5, -.5, -.5), ncol = 1)) %>% 
    as.numeric
  g <- resetUpdateAttributes(g)
  E(g)[to("H11")]$weight <- c(-.5, -.5, -.5)
  g <- updateVertices(g, getDeterminers = iparents, callback = calculateVals)
  logistic(linear_combination) %>% # calculation of prediction w/ weight of 5 via logistic function
    identical(unlist(V(g)["H11"]$output.signal)) %>%# compared to algorithms generation of prediction
    expect_true
})

