devtools::load_all("R/optimization.R")

context("optimization: testing logistic functions and prediction")

test_that("logistic_prime basic function is working as expected", {
  f <- function(z) exp(-z)/(1+exp(-z))^2
  z <- runif(100) 
  expect_equal(f(z), logistic_prime(z))
  expect_equal(logistic(4) - logistic(-4), integrate(Vectorize(logistic_prime), -4, 4)$value)
})

test_that("after the weights of a given vertex has changed, the prediction should change",{
  g <- get_gate("AND")
  original <- V(g)[type=="output"]$output.signal %>% unlist
  updated <- getPrediction(g, V(g)["AND"], runif(3))
  expect_true(!identical(original, updated))
})

context("optimization - univariate case on titanic model")
# Some tests on the titanic 3 data.  A 0-hidden layer model is fit so the cost function is 
# convex in the weights.  The model is univariate and without intercepts (biases) so that 
# the gradient has one dimension.  This is to check basic expectations of optimization behavior.
data(titanic3)
titan <- titanic3[!is.na(titanic3$age), ]
titan <- titan[!is.na(titan$survived), ] 
titan$survived <- as.numeric(titan$survived)
titan <- titan[, c("age", "survived")]
titan <- rescale_df(titan)[[1]]
# Creating a simple logistic regression on non-rescaled age, with no intercept.  
# So this should yield simple one dimensional gradient.
g <- mlp_graph("age", "survived") %>%
  initializeGraph(input.table = titan[, "age", drop = F], 
                  output.table = titan[, "survived", drop = F]) %>%

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

test_that("we get the same results as glm", {
  expected <- glm(survived ~ 0 + age, family = binomial(link = "logit"), data = titan) %>%
    coef %>% as.numeric
  # Test with just optimize
  output <- getLossFunction(g, "survived") %>% # Generate the loss function
    optimize(c(-5, 5)) %$% # Optimize it
    minimum
  expect_equal(expected, output, tolerance = .02)
})


test_that("loss is close to numeric integral of the gradient.", {
  get_loss <- getLossFunction(g, "survived")
  get_gradient <- getGradientFunction(g, "survived") 
  integral <- integrate(Vectorize(get_gradient), .1, .2)$value
  difference <- get_loss(.2) - get_loss(.1)
  expect_equal(integral, difference, tolerance = .01)
  # There is too much of a discrepency here.  
})

test_that("gradient should be giving expected values as numeric derivative",{
  get_loss <- Vectorize(getLossFunction(g, "survived"))
  get_grad <- Vectorize(getGradientFunction(g, "survived"))
  starting_vals <- runif(5, -1, 1)
  expected_deriv <- get_grad(starting_vals) %>% as.numeric
  numeric_deriv <- numericDeriv(quote(get_loss(starting_vals)),"starting_vals") %>% attr("gradient") %>% diag 
  expect_equal(numeric_deriv, expected_deriv, tolerance = .1)
  # There is too much of a discrepency here.  
})

context("optimization: univariate case with intercept")

test_that("doChainRule produces 0 in the case when the vertex value does not depend on the weight", {})

context("complex graphs")

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
  time_no_grad <- system.time(optim(initial_weights, get_loss,method = "BFGS", 
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

test_that("an update to an edge should result in a new edge weight (at least in early iterations of fitting", {
})




