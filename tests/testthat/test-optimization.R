devtools::load_all("../../R/optimization.R")
devtools::load_all("../../R/tools.R")
#devtools::load_all("R/optimization.R")
library(plyr)
library(dplyr)
context("optimization")

# Testing using Titanic3 data
data(titanic3)
titan <- filter(titanic3, !is.na(age), !is.na(survived), !is.na(fare)) %>% #Not worrying about NA vals
  mutate(survived = as.numeric(survived)) %>%
  select(., age, survived, fare) %>%
  rescale_df %$% #Note, everything is rescaled to between 0 and 1
  df

###########################################################################################
# Univariate no bias, no hidden layer case
g <- mlp_graph("age", "survived") %>%
  initializeGraph(select(titan, age, survived), fixed = "age") %>%
  {induced.subgraph(., V(.)[c("age", "survived")])} %>% # Having removed the bias, I need to reupdate 
  resetUpdateAttributes %>%
  updateSignals

test_that("prediction from logistic function works as expected.", {
  prediction <- unlist(V(getPrediction(g, V(g)["survived"], .5))["survived"]$output.signal)
  logistic(.5 * V(g)["age"]$output.signal[[1]]) %>% # calculation of prediction w/ weight of 5 via logistic function
    expect_equal(prediction)# compared to algorithms generation of prediction  
})

###########################################################################################
# Multivariate with hidden layers
g <- mlp_graph(c("age", "fare"), "survived", layers = c(3, 2)) %>%
  initializeGraph(select(titan, age, fare, survived), fixed = c("age", "fare"))

test_that("prediction from logistic function works as expected.", {
  linear_combination <- V(g)[c("age", "fare", "bias_H11")]$output.signal %>% 
    {do.call("cbind", .)} %>%
    `%*%`(matrix(c(-.5, -.5, -.5), ncol = 1)) %>% 
    as.numeric
  g <- resetUpdateAttributes(g)
  E(g)[to("H11")]$weight <- c(-.5, -.5, -.5)
  g <- updateSignals(g)
  logistic(linear_combination) %>% # calculation of prediction w/ weight of 5 via logistic function
    expect_equal(unlist(V(g)["H11"]$output.signal)) # compared to algorithms generation of prediction
})

###########################################################################################
# Testing penalized loss
test_that("penalized loss is greater than unpenalized loss for the same graph",{
  g_pen <- g
  g_pen$L1_pen <- .05
  g_pen$L2_pen <- .05
  g$L1_pen<- 0
  g$L2_pen<- 0
  v <- V(g)["survived"]
  pen_loss <- getObjective(g_pen, v)
  loss <- getObjective(g, v)
  expect_more_than(pen_loss(E(g)[to(v)]$weight), loss(E(g)[to(v)]$weight))
})

