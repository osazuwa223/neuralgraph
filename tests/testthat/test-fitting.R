library(plyr, quietly = TRUE)
library(dplyr, quietly = TRUE)
option=TRUE
data(mapk_g)
data(titanic3)
titan <- filter(titanic3, !is.na(age), !is.na(survived), !is.na(fare)) %>% #Not worrying about NA vals
  mutate(survived = as.numeric(survived)) %>%
  select(., age, survived, fare) %>%
  rescale_df %$% #Note, everything is rescaled to between 0 and 1
  df
################################################################################
context("Basic sanity checks in model fitting")
################################################################################
test_that("despite not working a single node graph input, regression on a constant works if the constant
          is given explicity.", {
            long_test(option)
            g_const <- graph.empty(2) %>% name_vertices %>% `+`(edge(c("2", "1"))) 
            data2 <- data.frame(runif(10), rep(1, 10)) %>% `names<-`(c("1", "2"))
            fitNetwork(g_const, data2, max.iter = 1)
          })

test_that("in an initialized neural net style model where the predicted and observed output are 
          exactly the same, the model results of fitting the model should be a graph with 0 error 
          and unchanged weights.", {
            long_test(option)
            set.seed(30)
            g <- random_unfit_sg(10, 3)
            #Set the observed values exactly to the predicted values
            V(g)[is.observed]$observed <- V(g)[is.observed]$output.signal
            # After fitting I expect no changes
            g2 <- fit_initialized_sg(g, .05, max.iter = 1)
            # The difference between predicted and observed should still be 0
            expect_equal(getMSE(g2), 0, .001)
            # The weights should be unchanged
            expect_equal(E(g2)$weight, E(g)$weight, tolerance= .01)
          })

################################################################################
context("Single predictor, no middle layers or intercept, on Titanic data and AND gate")
################################################################################
g <- mlp_graph("age", "survived") %>%
  initializeGraph(select(titan, age, survived), fixed = "age") %>%
  {induced.subgraph(., V(.)[c("age", "survived")])} %>% # Having removed the bias, I need to reupdate 
  resetUpdateAttributes %>%
  update_signals

test_that("prediction from logistic function works as expected.", {
  prediction <- unlist(V(getPrediction(g, V(g)["survived"], .5))["survived"]$output.signal)
  logistic(.5 * V(g)["age"]$output.signal[[1]]) %>% # calculation of prediction w/ weight of 5 via logistic function
    expect_equal(prediction)# compared to algorithms generation of prediction  
})

test_that("pure loss minimization gets the same results as nls", {
  expected <- nls(survived ~ logistic(age * w), data = titan, 
                  start = list(w = E(g)[to("survived")]$weight)) %>%
    coef %>% as.numeric
  # Test with just optimize
  output <- getObjective(g, "survived") %>% # Generate the loss function
    optimize(c(-5, 5)) %$% # Optimize it. True minimum based on data is about -.95
    minimum
  expect_equal(expected, output, tolerance = .02)
})


################################################################################
context("Multiple predictors, no layers or intercept, Titanic data ")
################################################################################
g <- mlp_graph(c("age", "fare"), "survived") %>%
  initializeGraph(select(titan, age, fare, survived), fixed = c("age", "fare"))

test_that("in multivariate case, pure loss minimization with optim gets the same results as nls", {
  set.seed(2702) # Sometimes you get a singular matrix, setting the seed fixes that.
  expected <- nls(survived ~ logistic(w0 + w1 * age + w2 * fare), data = titan, 
                  start = as.list(structure(E(g)[to("survived")]$weight, names = c("w1", "w2", "w0")))) %>%
    coef %>% # Get the coefficients 
    as.numeric
  output <- getObjective(g, "survived") %>% # Generate the loss function 
    {optim(E(g)$weight, .)} %$%
    par 
  expect_equal(expected, output, tolerance = .1)
})

################################################################################
context("Simple multi-layer models on gates")
################################################################################
test_that("reasonable prediction on one-layer gate model", {
  long_test(option)
  set.seed(31)
  g <- mlp_graph(c("I1", "I2"), "AND", 2)
  system <- list(I1 = c(1, 0), I2 = c(1, 0)) %>%
    expand.grid %>%
    `names<-`(c("I1", "I2")) %>%
  {dplyr::mutate(., AND = I1 * I2)}
  #igraphviz(g)
  fit <- fitNetwork(g, system, fixed = c("I1", "I2"), max.iter = 2, 
                  graph_attr = list(L2_pen = 0, min.max.constraints = c(-8,8)))
  expect_equal(unlist(V(fit)["AND"]$output.signal), 
               unlist(V(fit)["AND"]$observed), tolerance = .1)
})

test_that("reasonable prediction on multi-layer gate model", {
  g <- mlp_graph(c("I1", "I2"), c("AND", "OR", "NOR"), c(3, 2, 4))
  system <- list(I1 = c(1, 0), I2 = c(1, 0)) %>%
    expand.grid %>%
    `names<-`(c("I1", "I2")) %>%
    {dplyr::mutate(., AND = I1 * I2)} %>%
    mutate(AND = I1 * I2, OR = (I1 + I2 > 0) * 1, NOR = (I1 + I2 == 0) * 1)
  fit <- fitNetwork(g, system, fixed = c("I1", "I2"), max.iter = 2, 
                    graph_attr = list(L2_pen = 0, min.max.constraints = c(-8,8)))
  expect_equal(unlist(V(fit)["AND"]$output.signal), 
               unlist(V(fit)["AND"]$observed), tolerance = .1)
})

################################################################################
context("Single variate response with multiple-layers on Titanic data")
################################################################################
test_that("multi-layer model has less loss than nls given it has more parameters.", {
  long_test(option)
  set.seed(33)
  g <- mlp_graph(c("age", "fare"), "survived", c(2, 1)) %>%
  {initializeGraph(., data = dplyr::select(titan, age, fare, survived), fixed = c("age", "fare"))}
  fit <- fit_initialized_sg(g,  epsilon = .01, max.iter = 2)
  nls_fit <-  nls(survived ~ logistic(w0 + w1 * age + w2 * fare), data = titan, 
                start = as.list(structure(E(g)[to("H11")]$weight, names = c("w1", "w2", "w0"))))
  expect_less_than(vertexMSE(fit, V(g)["survived"]), deviance(nls_fit))
})

################################################################################
context("Evaluating denser MLP on the Titanic data")
################################################################################
g <- mlp_graph(c("age", "fare"), "survived", layers = c(3, 2)) %>%
  initializeGraph(select(titan, age, fare, survived), fixed = c("age", "fare"))
test_that("prediction from logistic function works as expected.", {
  linear_combination <- V(g)[c("age", "fare", "bias_H11")]$output.signal %>% 
    {do.call("cbind", .)} %>%
    `%*%`(matrix(c(-.5, -.5, -.5), ncol = 1)) %>% 
    as.numeric
  g <- resetUpdateAttributes(g)
  E(g)[to("H11")]$weight <- c(-.5, -.5, -.5)
  g <- update_signals(g)
  logistic(linear_combination) %>% # calculation of prediction w/ weight of 5 via logistic function
    expect_equal(unlist(V(g)["H11"]$output.signal)) # compared to algorithms generation of prediction
})

g <- mlp_graph(c("age", "fare"), "survived", c(5, 5)) %>%
  initializeGraph(select(titan, age, fare, survived), fixed = c("age", "fare"))
test_that("a limited BFGS optimization reduces loss for the output node.", {
  get_loss <- getObjective(g, "survived")
  initial_weights <- E(g)[to("survived")]$weight
  updated_weights <- optim(initial_weights, get_loss, control = list(maxit = 30))$par #Just allowing for 10 iterations
  expect_true(get_loss(initial_weights) > get_loss(updated_weights)) 
})

test_that("a limited BFGS optimization reduces loss for an arbitrary hidden node.", {
  get_loss <- getObjective(g, "H11")
  initial_weights <- E(g)[to("H11")]$weight
  updated_weights <- optim(initial_weights, get_loss, control = list(maxit = 30))$par #Just allowing for 10 iterations
  expect_true(get_loss(initial_weights) > get_loss(updated_weights)) 
})

################################################################################
context("Evaluating fits on less structured DAGs")
################################################################################
test_that("A simulated DAG approach doesn't fail", {
  long_test(option)
  set.seed(100)
  g <- random_unfit_sg(25, 100, input_g = mapk_g)
  g$L2_pen <- 0.0389
  #sg_viz(g)
  fit <- fit_initialized_sg(g, max.iter = 1, verbose = TRUE)  
  #sg_viz(fit)
  expect_true(getMSE(g) > getMSE(fit))
  fit2 <- fit_initialized_sg(fit, max.iter = 1, verbose = TRUE)
  sg_viz(fit2)
  expect_true(getMSE(fit) > getMSE(fit2))
})

test_that("If a downstream vertex is not optimized, it gets optimized first", {
  long_test(option)
  set.seed(100)
  g <- random_unfit_sg(10, 100, method = "ic-dag")
  #sg_viz(g)
  #Here optimizing 7 should not occur until 10, 9, and 5 are fit, but not 8
  g$L2_pen <- 0.0389
  # start of fit_initialized_sg
  g <- resetUpdateAttributes(g)
  # start of update_weights, where vertices are normalis
  expect_message(
    vertex_updater2(g, 7, get_downstream_vertices, fit_weights_for_node, verbose = TRUE),
    "propagating to vertices 5, 9, 10"
  ) # not sure how to exclude 8
  }
})

################################################################################
context("Introduction of penalty terms")
################################################################################
g_structure <- mlp_graph(c("age", "fare"), "survived", c(4, 3)) %>%
  initializeGraph(select(titan, age, fare, survived), fixed = c("age", "fare"))
test_that("penalized loss is greater than unpenalized loss for the same graph",{
  g <- g_structure
  g_pen <- g
  g_pen$L1_pen <- .04
  g_pen$L2_pen <- .04
  g$L1_pen<- 0
  g$L2_pen<- 0
  v <- V(g)["survived"]
  pen_loss <- getObjective(g_pen, v)
  loss <- getObjective(g, v)
  expect_more_than(pen_loss(E(g_pen)[to(v)]$weight), loss(E(g)[to(v)]$weight))
})

test_that("L1 norm with super high penalty should bring weights from non-bias terms to 0", {
  long_test(option)
  set.seed(30)
  g_pen <- g_structure
  g_pen$L1_pen <- 50
  g_pen <- fit_initialized_sg(g_pen, epsilon = .01, max.iter = 1, verbose =TRUE)
  expect_equal(sum(round(E(g_pen)[from(V(g_pen)[!is.bias])]$weight, 4)), 0)  
})


################################################################################
context("Testing gradient-related calculations")
################################################################################
test_that("logistic_prime basic function is working as expected", {
  skip("skipping gradient related problems")
  f <- function(z) exp(-z)/(1+exp(-z))^2
  z <- runif(100) 
  expect_equal(f(z), logistic_prime(z))
  expect_equal(logistic(4) - logistic(-4), integrate(Vectorize(logistic_prime), -4, 4)$value)
})

g <- mlp_graph("age", "survived") %>%
  initializeGraph(select(titan, age, survived), fixed = "age") %>%
  {induced.subgraph(., V(.)[c("age", "survived")])} %>% # Having removed the bias, I need to reupdate 
  resetUpdateAttributes %>%
  update_signals

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
  initializeGraph(select(titan, age, fare, survived), fixed = c("age", "fare"))

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
