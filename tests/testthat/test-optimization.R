context("optimization")

test_that("after the weights of a given vertex has changed, the prediction should change",{
  g <- get_gate("AND")
  original <- V(g)[type=="output"]$output.signal %>% unlist
  updated <- getPrediction(g, V(g)["AND"], runif(3))
  !identical(original, updated) %>% expect_true
})


test_that("a limited optimization reduces loss.", {
  g <- get_gate("AND", layers = c(1, 2))
  get_loss <- getLossFunction(g, "AND")
  initial_weights <- E(g)[to("AND")]$weight
  updated_weights <- optim(initial_weights, get_loss, control = list(maxit = 40))$par #Just allowing for 10 iterations
  get_loss(initial_weights) > get_loss(updated_weights) %>%
    expect_true
})

test_that("a limited optimizaton with gradient reduces loss.", {
  g <- get_gate("AND", layers = c(2, 1))
  get_loss <- getLossFunction(g, "AND")
  get_gradient <- getGradientFunction(g, "AND")
  initial_weights <- E(g)[to("AND")]$weight
  updated_weights_w_grad <- optim(initial_weights, get_loss, gr = get_gradient, method = "BFGS", 
                                  control = list(maxit = 40))$par #Just allowing for 10 iterations
  get_loss(updated_weights_no_grad) > get_loss(updated_weights_w_grad) %>%
    expect_true
})

test_test("gradient should be giving expected values",{
  g <- get_gate("AND")
  initial_weights <- E(g)[to("AND")]$weight
  get_loss <- getLossFunction(g, "AND")
  get_gradient <- getGradientFunction(g, "AND") 
  #Integrate in the first dimension
  first_gradient <- function(w) get_gradient(c(w, initial_weights[2:3]))[1] %>% as.numeric
  integrate(Vectorize(first_gradient), -10, initial_weights[1])
  get_loss(initial_weights) - get_loss(c(-10, initial_weights[2:3]))
})

test_that("a gradient descent step reduces loss in a MLP case.", {
  `names<-`(c("I1", "I2")) %>%
    mutate(AND = I1 * I2)
  g <- mlpgraph(c("I1", "I2"), c(3, 2), c("AND")) %>% #Use a 2 layer MLP
    initializeGraph(input.table = system[, c("I1", "I2")], 
                    output.table = system[, "AND", drop = F])
  v <- V(g)["AND"]
  #loss <- getLossFunction(g, v)
  # getPrediction is the problem
  loss <- function(weights){
    output_vertex <- V(g)[type=="output"]
    #prediction <- getPrediction(g, v, weights)
    observed <- unlist(output_vertex$observed)
    .5 * sum( (observed - prediction) ^ 2)
  }
  gradient <- getGradientFunction(g, v)
  weights_initial <- E(g)[to(v)]$weight
  weights_optimized <- optim(weights_initial, fn = loss, gr = gradient, method="BFGS")$par
  # ls  
  expect_true(loss(weights_optimized) < loss(weights_initial))
})

test_that("an update to an edge should result in a new edge weight (at least in early iterations of fitting", {
})




mtcars2 <- rescale_df(mtcars)$df
Y <- mtcars2[, 1]
p <- 3
X <- mtcars2[, 2:(p+1)] %>% as.matrix
w <- runif(p)
logistic
loss <- function(w){
  linear_combo <- X %*% matrix(w, ncol = 1) %>% as.numeric
  output <- logistic(linear_combo)
  sum(.5 * (Y - (nonlinear_output))^2)
}
grad <- function(w){
  linear_combo <- X %*% matrix(w, ncol = 1) %>% as.numeric
  output <- logistic(linear_combo)
  derivative <- logistic.prime
  colSums((Y - (linear_combo)) * (-X)) %>%
    `names<-`(paste("w", 1:length(w), sep=""))
}

#Perhaps genarlize to any gradient?
get_first_gradient <- function(fixed){
  Vectorize(function(w) grad(c(w, fixed))[1])
}

test_that("assuming logistic activation, gradient of the loss function numerically itegrates to the loss function.", {
  fixed <- structure(runif(p - 1), names = paste("w", 2:(p), sep="")) #fixed dimensions
  rng <- c(lower = -2, upper = 2) 
  first_gradient <- get_first_gradient(fixed)
  integral_output <- integrate(first_gradient, rng[1], rng[2])$value
  regular_output <- loss(c(rng[2], fixed)) - loss(c(rng[1], fixed))
  expect_equal(integral_output, regular_output , tolerance = .002)
})

