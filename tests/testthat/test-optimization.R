devtools::load_all("R/optimization.R")
library(plyr)
library(dplyr)
context("optimization")

# Here I wish to check optimization with logistic activation and gradient works in the case
# of the loss function being convex in a one-dimensional weight variable.  
# So I build a 0-hidden layer model of logistic regression of 'survived' vs 'age' in the 
# Titanic3 data.  No intercept is given so the weight is one-dimensional.  'age' is 
# rescaled to between 0 and 1. 
data(titanic3)
titan <- filter(titanic3, !is.na(age), !is.na(survived)) %>% #Not worrying
  mutate(survived = as.numeric(survived)) %>%
  select(age, survived, fare) %>%
  rescale_df %$%
  df

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

test_that("gradient is positive if greater than minimum, negative if less than minimum",{
  min_val <- getLossFunction(g, "survived") %>% # Generate the loss function
    optimize(c(-5, 5)) %$% # Optimize it
    minimum # Pull out the minimum
  get_gradient <- getGradientFunction(g, "survived") # Generate the gradient
  starting_vals <- runif(100, -1, 1) # Generate a bunch of starting vals
  lapply(starting_vals, function(starting_val){ # Test if the sign if the gradient is correct
    if(starting_val < min_val){
      get_gradient(starting_val) %>% expect_less_than(0)
    } else {
      get_gradient(starting_val) %>% expect_more_than(0)
    }
  })
})

test_that("gradient is near 0 when loss is minimized",{
  val_at_min <- getLossFunction(g, "survived") %>% # Generate the loss function
    optimize(c(-5, 5)) %$% # Optimize it
    minimum # Pull out the weight value at the minimum
  get_gradient <- getGradientFunction(g, "survived") # Generate the gradient
  expect_equal(get_gradient(val_at_min ), 0, tolerance = .1)
})

test_that("we get the same results as glm", {
  expected <- glm(survived ~ 0 + age, family = binomial(link = "logit"), data = titan) %>%
    coef %>% as.numeric
  # Test with just optimize
  output <- getLossFunction(g, "survived") %>% # Generate the loss function
    optimize(c(-5, 5)) %$% # Optimize it. True minimum based on data is about -.95
    minimum
  expect_equal(expected, output, tolerance = .02)
})

test_that("numeric integral of gradient looks like loss", {
  get_loss <- getLossFunction(g, "survived")
  # Numeric results are the most correct around the minimum
  get_gradient <- getGradientFunction(g, "survived")
  # The numeric integral of the gradiet should approximate loss (at least close to the minimum ~ -.95)
  integral <- integrate(Vectorize(get_gradient), -1, -.9)$value
  difference <- get_loss(-.9) - get_loss(-1)
  expect_equal(integral, difference, tolerance = .01)
})

test_that("gradient should be giving expected values as numeric derivative",{
  get_loss <- Vectorize(getLossFunction(g, "survived"))
  get_grad <- Vectorize(getGradientFunction(g, "survived"))
  val <- -0.9468121
  expected_deriv <- get_grad(val) %>% as.numeric
  numeric_deriv <- numericDeriv(quote(get_loss(val)), "vals") %>% attr("gradient") %>% diag 
  expect_equal(numeric_deriv, expected_deriv, tolerance = .1)  
})

# Next expanding the analysis to the multivariate case.  The goal is to assure the logistic gradient
# and loss functions still behave well in the multivariate-convex palce.
g <- mlp_graph("age", "survived") %>%
  initializeGraph(input.table = titan[, "age", drop = F], 
                  output.table = titan[, "survived", drop = F]) 

# Same tests as before except adding an intercept
# This should yield a two dimensional gradient, so comparison to optim instead of optimize will have to be used.
# In this case however, the weights are still identifiable

test_that("prediction from logistic function works as expected.", {
  linear_combination <- .5 * V(g)["age"]$output.signal[[1]] + .5 * 1 # .5 * 1 is for the intercept
  logistic(linear_combination) %>% # calculation of prediction w/ weight of 5 via logistic function
    identical(getPrediction(g, V(g)["survived"], c(.5, .5))) %>%# compared to algorithms generation of prediction
    expect_true
})

test_that("doChainRule in layer-free univarite produces logistic_prime(input * weight) * input", {
  weights <- E(g)[to("survived")]$weight
  age_out <- V(g)[c("age")]$output.signal %>% unlist
  linear_combo <- V(g)[c("age", "int.2")]$output.signal %>% 
    {do.call("cbind", .)} %>%
    `%*%`(matrix(weights, ncol = 1))
  expected <- cbind(age_weight = logistic_prime(linear_combo) * age_out, 
                    int.2_weight = logistic_prime(linear_combo))
  chain_rule_output <- sapply(E(g)[to("survived")], doChainRule, g = g, v = V(g)["survived"])
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

test_that("gradient is near 0 when loss is minimized",{
  val_at_min <- getLossFunction(g, "survived") %>% # Generate the loss function
    optimize(c(-5, 5)) %$% # Optimize it
    minimum # Pull out the weight value at the minimum
  get_gradient <- getGradientFunction(g, "survived") # Generate the gradient
  expect_equal(get_gradient(val_at_min ), 0, tolerance = .1)
})


test_that("we get the same results as glm", {
  expected <- glm(survived ~ age, family = binomial(link = "logit"), data = titan) %>%
    coef %>% as.numeric %>% sort
  output <- getLossFunction(g, "survived") %>% # Generate the loss function
    {optim(E(g)$weight, ., getGradientFunction(g, "survived"))} %$%
    par %>% sort
  expect_equal(expected, output, tolerance = .1)
})





context("complex graphs")

test_that("doChainRule produces 0 in the case when the vertex value does not depend on the weight", {})

test_that("a limited optimization reduces loss.", {
  g <- get_gate("AND", layers = c(1, 2))
  get_loss <- getLossFunction(g, "AND")
  initial_weights <- E(g)[to("AND")]$weight
  updated_weights <- optim(initial_weights, get_loss, control = list(maxit = 40))$par #Just allowing for 10 iterations
  expect_true(get_loss(initial_weights) > get_loss(updated_weights)) 
})

test_that("use of gradient should at least speed things up.", {
  g <- get_gate("AND", layers = c(2, 1))
  get_loss <- getLossFunction(g, "AND")
  get_gradient <- getGradientFunction(g, "AND")
  initial_weights <- E(g)[to("AND")]$weight
  time_no_grad <- system.time(optim(initial_weights, get_loss, method = "BFGS", 
                                    control = list(maxit = 10)))[["elapsed"]] #Just allowing for 40 iterations
  time_w_grad <- system.time(optim(initial_weights, get_loss, gr = get_gradient, method = "BFGS", 
                                   control = list(maxit = 10)))[["elapsed"]] #Just allowing for 40 iterations
  expect_true(time_w_grad < time_no_grad)    
})


test_that("a gradient descent step reduces loss in a MLP case.", {
  g <- get_gate(outputs = "AND", layers = c(2, 3))
  v <- V(g)["AND"]
  loss <- getLossFunction(g, v)
  gradient <- getGradientFunction(g, v)
  weights_initial <- E(g)[to(v)]$weight
  weights_optimized <- optim(weights_initial, fn = loss, gr = gradient, method="BFGS")$par
  expect_true(loss(weights_optimized) < loss(weights_initial))
})





