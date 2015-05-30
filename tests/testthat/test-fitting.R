devtools::load_all("../../R/tools.R")
#devtools::load_all("R/tools.R")
option <- FALSE
context("Models gets reasonable fits")

test_that("despite not working a single node graph input, regression on a constant works if the constant
          is given explicity.", {
            g_const <- graph.empty(2) %>% nameVertices %>% `+`(edge(c("2", "1"))) 
            data2 <- data.frame(runif(10), rep(1, 10)) %>% `names<-`(c("1", "2"))
            fitNetwork(g_const, data2)
          })

test_that("in an initialized neural net style model where the predicted and observed output are 
          exactly the same, the model results of fitting the model should be a graph with 0 error 
          and unchanged weights.", {
  long_test(option)
  set.seed(30)
  g <- random_unfit_sg(8)
  #Set the observed values exactly to the predicted values
  V(g)[is.observed]$observed <- V(g)[is.observed]$output.signal
  output.df[, paste(outputs)] <- unlist(V(g)[is.observed]$observed)
  # After fitting I expect no changes
  g2 <- fitInitializedNetwork(g, .05, 3)
  # The difference between predicted and observed should still be 0
  expect_equal(getLoss(g2), 0, .001)
  # The weights should be unchanged
  expect_equal(E(g2)$weight, E(g)$weight, tolerance= .01)
})

test_that("I can get descent results on a linear model with hidden nodes", {})

### On titanic
data(titanic3)
titan <- dplyr::filter(titanic3, !is.na(age), !is.na(survived), !is.na(fare)) %>% #Not worrying about NA vals
  {dplyr::mutate(., survived = as.numeric(survived))} %>%
  {dplyr::select(., age, survived, fare)} %>%
  rescale_df %$% #Note, everything is rescaled to between 0 and 1
  df

test_that("multi-layer model has less loss than nls given it has more parameters.", {
  long_test(option)
  set.seed(30)
  g <- mlp_graph(c("age", "fare"), "survived", c(2, 1)) %>%
    {initializeGraph(., data = dplyr::select(titan, age, fare, survived))}
  fit <- fitInitializedNetwork(g,  epsilon = .01, verbose = T)
  nls_fit <-  nls(survived ~ logistic(w0 + w1 * age + w2 * fare), data = titan, 
                  start = as.list(structure(E(g)[to("H11")]$weight, names = c("w1", "w2", "w0"))))
  expect_less_than(get_deviance(fit), deviance(nls_fit))
})

test_that("model should perform a reasonable MLP prediction on a toy problem with single output.", {
  # In the future, expand to multivariate case
  #g <- mlpgraph(c("I1", "I2"), c(3, 2, 4), c("AND", "OR", "NOR"))
  long_test(option)
  set.seed(31)
  g <- mlp_graph(c("I1", "I2"), "AND", 2)
  system <- list(I1 = c(1, 0), I2 = c(1, 0)) %>%
    expand.grid %>%
    `names<-`(c("I1", "I2")) %>%
    {dplyr::mutate(., AND = I1 * I2)}
    #mutate(AND = I1 * I2, OR = (I1 + I2 > 0) * 1, NOR = (I1 + I2 == 0) * 1)
  fit <- fitNetwork(g, input.table = dplyr::select(system, I1, I2), 
                    output.table = dplyr::select(system, AND), verbose = T)
  expect_equal(unlist(V(fit)["AND"]$output.signal), unlist(V(fit)["AND"]$observed), tolerance = .1)
})

test_that("penalized least squares has less sum squares of fitted weight than unpenalized.", {
  long_test(option)
  set.seed(30)
  g_structure <- mlp_graph(c("age", "fare"), "survived", c(4, 3))
  g_no_pen <- {fitNetwork(g_structure, input.table = dplyr::select(titan, age, fare), 
     output.table = dplyr::select(titan, survived), epsilon = .01, verbose = T)}
  g_pen <- {fitNetwork(g_structure, input.table = dplyr::select(titan, age, fare), 
     output.table = dplyr::select(titan, survived), penalty = .05, epsilon = .01, verbose = T)}
  expect_less_than(sum(E(g_pen)$weight^2), sum(E(g_no_pen)$weight^2))  
})


