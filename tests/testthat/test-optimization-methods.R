# Testing and evaluation of performance of various optimization techniques
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

g <- mlp_graph("age", "survived") %>%
  initializeGraph(input.table = select(titan, age), 
                  output.table = select(titan, survived)) %>%
  {induced.subgraph(., V(.)[c("age", "survived")])} %>% # Having removed the intercept, I need to reupdate 
  resetUpdateAttributes %>%
  updateVertices(getDeterminers = iparents, callback = calculateVals)


###########################################################################################
# Tests of optimization of logistic activation :
#   * logistic regression of 'survived' vs 'age' in the Titanic3 data
#   * no intercept -- one dimensional weight
#   * no hidden layers -- straight forward logistic regression on squared error loss
#   * loss function is convex in the weight
#   * age is rescaled to between 0 and 1


test_that("in single variate no intercept case, pure loss minimization gets the same results as nls", {
  expected <- nls(survived ~ logistic(age * w), data = titan, 
                  start = list(w = E(g)[to("survived")]$weight)) %>%
    coef %>% as.numeric
  # Test with just optimize
  output <- getLossFunction(g, "survived") %>% # Generate the loss function
    optimize(c(-5, 5)) %$% # Optimize it. True minimum based on data is about -.95
    minimum
  expect_equal(expected, output, tolerance = .02)
})

###########################################################################################
# Next expanding to multivariate.  Again, the goal is to compare to nls().

g <- mlp_graph(c("age", "fare"), "survived") %>%
  initializeGraph(input.table = select(titan, age, fare), output.table = select(titan, survived))

test_that("in multivariate case, pure loss minimization gets the same results as nls", {
  expected <- nls(survived ~ logistic(w0 + w1 * age + w2 * fare), data = titan, 
                  start = as.list(structure(E(g)[to("survived")]$weight, names = c("w1", "w2", "w0")))) %>%
    coef %>% # Get the coefficients 
    as.numeric
  output <- getLossFunction(g, "survived") %>% # Generate the loss function 
    {optim(E(g)$weight, .)} %$%
    par 
  expect_equal(expected, output, tolerance = .1)
})

###########################################################################################
# Next expanding to the case of multiple hidden layers
g <- mlp_graph(c("age", "fare"), "survived", c(5, 5)) %>%
  initializeGraph(input.table = select(titan, age, fare), 
                  output.table = select(titan, survived))

test_that("a limited BFGS optimization reduces loss for the output node.", {
  get_loss <- getLossFunction(g, "survived")
  initial_weights <- E(g)[to("survived")]$weight
  updated_weights <- optim(initial_weights, get_loss, control = list(maxit = 30))$par #Just allowing for 10 iterations
  expect_true(get_loss(initial_weights) > get_loss(updated_weights)) 
})

test_that("a limited BFGS optimization reduces loss for an arbitrary hidden node.", {
  get_loss <- getLossFunction(g, "H11")
  initial_weights <- E(g)[to("H11")]$weight
  updated_weights <- optim(initial_weights, get_loss, control = list(maxit = 30))$par #Just allowing for 10 iterations
  expect_true(get_loss(initial_weights) > get_loss(updated_weights)) 
})
