library(plyr, quietly = TRUE)
library(dplyr, quietly = TRUE)
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
data(titanic3)
titan <- filter(titanic3, !is.na(age), !is.na(survived), !is.na(fare)) %>% #Not worrying about NA vals
  mutate(survived = as.numeric(survived)) %>%
  select(., age, survived, fare) %>%
  rescale_df %$% #Note, everything is rescaled to between 0 and 1
  df
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

test_that("in multivariate case, pure loss minimization gets the same results as nls", {
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
  fit <- fitNetwork(g, system, fixed = c("I1", "I2"), max.iter = 2, 
                  graph_attr = list(L2_pen = 0, min.max.constraints = c(-20,20)))
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
context("Evaluating denser networks on the Titanic data")
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

test_that("for a dense network, edges should update at each step", {
  expect_true(FALSE)
})


################################################################################
context("Introduction of penalty terms")
################################################################################
g_structure <- mlp_graph(c("age", "fare"), "survived", c(4, 3))
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

test_that("L1 norm with super high penalty should bring weights from non-bias terms to 0", {
  long_test(option)
  set.seed(30)
  g_pen <- {fitNetwork(g_structure, dplyr::select(titan, age, fare, survived),
                       fixed = c("age", "fare"),
                       graph_attr = list(L1_pen = 100), epsilon = .01, max.iter = 1)}
  expect_equal(sum(round(E(g_pen)[from(V(g_pen)[!is.bias])]$weight, 4)), 0)  
})

test_that("L2 norm has less sum squares of fitted weight than unpenalized.", {
  long_test(option)
  set.seed(30)
  g_structure <- mlp_graph(c("age", "fare"), "survived", c(4, 3))
  g_no_pen <- {fitNetwork(g_structure, dplyr::select(titan, age, fare, survived), 
                          fixed = c("age", "fare"), epsilon = .01, max.iter = 1)}
  g_pen <- {fitNetwork(g_structure, dplyr::select(titan, age, fare, survived), 
                       fixed = c("age", "fare"),
                       graph_attr = list(L2_pen = 3), 
                       epsilon = .01, max.iter = 1)}
  expect_less_than(sum(E(g_pen)$weight^2), sum(E(g_no_pen)$weight^2))  
})

################################################################################
context("Comparing one dimensional gradient function to calculated gradient")
################################################################################
test_that("logistic_prime basic function is working as expected", {
  f <- function(z) exp(-z)/(1+exp(-z))^2
  z <- runif(100) 
  expect_equal(f(z), logistic_prime(z))
  expect_equal(logistic(4) - logistic(-4), integrate(Vectorize(logistic_prime), -4, 4)$value)
})