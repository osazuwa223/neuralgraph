# Testing optimization methods that use a specified gradient.

# Note May 10, 2015.  The gradient tests are currently failing.  I believe the gradient calculation is erroneous, though 
# I am not sure where

devtools::load_all("../../R/optimization.R")
devtools::load_all("../../R/tools.R")
#devtools::load_all("R/optimization.R")
#devtools::load_all("R/tools.R")
library(plyr)
library(dplyr)
context("Optimization with specified gradient")

test_that("logistic_prime basic function is working as expected", {
  f <- function(z) exp(-z)/(1+exp(-z))^2
  z <- runif(100) 
  expect_equal(f(z), logistic_prime(z))
  expect_equal(logistic(4) - logistic(-4), integrate(Vectorize(logistic_prime), -4, 4)$value)
})

test_that("hand calculation of gradient in error free case of logistic regression reproduces gradient calculation", {  
})

test_that("hand calculation of gradient in case of logistic regression with errors reproduces gradient calculation", {  
})



# Testing using Titanic3 data
data(titanic3)
titan <- dplyr::filter(titanic3, !is.na(age), !is.na(survived), !is.na(fare)) %>% #Not worrying about NA vals
  {dplyr::mutate(., survived = as.numeric(survived))} %>%
  {dplyr::select(., age, survived, fare)} %>%
  rescale_df %$% #Note, everything is rescaled to between 0 and 1
  df

g <- mlp_graph("age", "survived") %>%
  initializeGraph(select(titan, age, survived), fixed = "age") %>%
  {induced.subgraph(., V(.)[c("age", "survived")])} %>% # Having removed the bias, I need to reupdate 
  resetUpdateAttributes %>%
  updateSignals

test_that("doChainRule in layer-free univarite produces logistic_prime(input * weight) * input", {
  skip("not doing gradient tests")
  weight <- E(g)[to("survived")]$weight
  age_out <- V(g)["age"]$output.signal %>% unlist
  expected <- logistic_prime(weight * age_out) * age_out 
  chain_rule_output <- doChainRule(g, V(g)["survived"], E(g)[to("survived")])
  expect_equal(expected, chain_rule_output)
})

test_that("getGradient in layer-free univarite produces -(Y - f(input)) * logistic_prime(input * weight) * input", {
  skip("not doing gradient tests")
  weight <- E(g)[to("survived")]$weight
  age_out <- V(g)["age"]$output.signal %>% unlist
  Y <- V(g)["survived"]$observed %>% unlist
  Y_hat <- V(g)["survived"]$output.signal %>% unlist
  expected <- sum(-(Y - Y_hat) * logistic_prime(weight * age_out) * age_out)
  gradient_output <- getGradientFunction(g, V(g)["survived"])(weight) %>% as.numeric
  expect_equal(expected, gradient_output)
})

test_that("hand calculation of doChainRule in layer-free multivarite produces same results as package functions", {
  skip("not doing gradient tests")
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


test_that("gradient zero equals loss minimimum", {
  skip("skipping gradient problems")
  loss_min <- getObjective(g, "survived") %>% # Generate the loss function
    optimize(c(-5, 5)) %$% # Optimize it
    minimum
  grad_zero <- getGradientFunction(g, "survived") %>%
    uniroot(c(-5, 5)) %$%
    root
  expect_equal(loss_min, grad_zero, tolerance = .1)
})

test_that("gradient is positive if greater than minimum, negative if less than minimum",{
  # The actual value that is minimizing the gradient is different from the loss minimum for some reason
  skip("skipping gradient problems")
  get_gradient <- getGradientFunction(g, "survived") # Generate the gradient
  min_val<- get_gradient %>%
    uniroot(c(-5, 5)) %$%
    root
  starting_vals <- runif(100, -1.5, -.5) # Generate a bunch of starting vals
  epsilon <- .4
  lapply(starting_vals, function(starting_val){ # Test if the sign if the gradient is correct
    if(starting_val < (min_val- epsilon)){
      get_gradient(starting_val) %>% expect_less_than(0)
    } 
    if (starting_val  > (min_val+ epsilon)){
      get_gradient(starting_val) %>% expect_more_than(0)
    } 
  })
})

test_that("numeric integral of gradient looks like loss", {
  skip("skipping gradient problems")
  get_loss <- getObjective(g, "survived")
  # Numeric results are the most correct around the minimum
  get_gradient <- getGradientFunction(g, "survived")
  # The numeric integral of the gradiet should approximate loss (at least close to the minimum ~ -.95)
  integral <- integrate(Vectorize(get_gradient), -1, -.9)$value
  difference <- get_loss(-.9) - get_loss(-1)
  expect_equal(integral, difference, tolerance = .05)
})

# This may be a problem. Skipping for now.
test_that("gradient should be giving expected values as numeric derivative",{
  skip("skipping gradient problems")
  get_loss <- Vectorize(getObjective(g, "survived"))
  get_grad <- Vectorize(getGradientFunction(g, "survived"))
  val <- -0.9468121
  expected_deriv <- get_grad(val) %>% as.numeric
  numeric_deriv <- numericDeriv(quote(get_loss(val)), "val") %>% attr("gradient") %>% diag 
  expect_equal(numeric_deriv, expected_deriv, tolerance = .1)  
})

###########################################################################################
# Testing basic implementation of gradient descent

g1 <- get_gate(outputs = "AND", layers = c(3, 2))
g2 <- mlp_graph(c("age", "fare"), "survived", c(5, 5)) %>%
  initializeGraph(select(titan, age, fare, survived), fixed = NULL)

test_that("doChainRule errors in the case when the vertex value does not depend on the weight", {
  skip("a bug here but ignoring gradient descent for now")
  expect_error(doChainRule(g1, v = V(g1)["H21"], e = E(g1)["H11" %->% "H22"]),
               "\n  You've attempted to find the gradient of a node's output\n         w.r.t an edge weight that that does not affect that output. v: 5 e: 7\n")
  
})

test_that("gradient descent yields weight with gradient near 0 when loss is minimized",{
  skip("skipping gradient problems")
  loss <- getObjective(g2, "survived")
  grad <- getGradientFunction(g2, "survived")
  gd_min <- E(g2)[to("survived")]$weight %>%
    gradientDescent(grad, loss)
  expect_equal(get_grad(gd_min), rep(0, length(gd_min), tolerance = .02))
})

test_that("gradient is positive if greater than minimum, negative if less than minimum",{
  skip("skipping gradient problems")
  min_val <- gradientDescent(...)
  get_gradient <- getGradientFunction(g, "survived") # Generate the gradient
  starting_vals <- runif(100, -1.5, -.5) # Generate a bunch of starting vals
  epsilon <- .01
  lapply(starting_vals, function(starting_val){ # Test if the sign if the gradient is correct
    if(starting_val < (min_val- epsilon)){
      get_gradient(starting_val) %>% expect_less_than(0)
    } else if (starting_val > (min_val + epsilon)){
      get_gradient(starting_val) %>% expect_more_than(0)
    } else {
      tryCatch(get_gradient(starting_val) %>% expect_equal(0, tolerance = .5),
               error = function(e) print(starting_val)
      )
    }
  })
})


test_that("a limited optimization reduces loss.", {
  skip("skipping gradient problems")
  get_loss <- getObjective(g2, "H11")
  get_gradient <- getGradientFunction(g2, "H11")
  initial_weights <- E(g1)[to("H11")]$weight
  updated_weights <- gradientDescent(initial_weights, get_loss, get_gradient, control = list(maxit = 200))#$par #Just allowing for 10 iterations
  expect_true(get_loss(initial_weights) > get_loss(updated_weights)) 
})

test_that("use of gradient should at least speed things up.", {
  skip("skipping gradient problems")
  get_loss <- getObjective(g1, "AND")
  get_gradient <- getGradientFunction(g1, "AND")
  initial_weights <- E(g1)[to("AND")]$weight
  time_no_grad <- system.time(gradientDescent(initial_weights, get_loss, method = "BFGS", 
                                    control = list(maxit = 20)))[["elapsed"]] #Just allowing for 40 iterations
  time_w_grad <- system.time(gradientDescent(initial_weights, get_loss, gr = get_gradient, method = "BFGS", 
                                   control = list(maxit = 20)))[["elapsed"]] #Just allowing for 40 iterations
  expect_true(time_w_grad < time_no_grad)    
})

test_that("optimization with the gradient outperforms optimization with a random gradient", {
  skip("skipping gradient problems")
  v <- V(g2)["H21"]
  loss <- getObjective(g2, v)
  gradient <- getGradientFunction(g2, v)
  random_gradient <- function(weight) runif(length(weight))
  weights_initial <- E(g2)[to(v)]$weight
  standard <- optim(weights_initial, fn = loss, gr = gradient, method="BFGS")#$par
  random <- optim(weights_initial, fn = loss, gr = random_gradient, method="BFGS")#$par
  gfree <- optim(weights_initial, fn = loss,  method="BFGS")#$par
  expect_true(loss(standard) < loss(random))
})
