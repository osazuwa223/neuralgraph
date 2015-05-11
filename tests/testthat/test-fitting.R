devtools::load_all("../../R/tools.R")
#devtools::load_all("R/tools.R")
context("Models gets reasonable fits")

opt <-  T # Option to skip tests where the network is fit
test_fit <- function(){
  if(opt == T) skip("Skipping a fit")
}

test_that("in an initialized model where the predicted and observed output are exactly the same, 
the model results of fitting the model should be a graph with 0 error and unchanged weights.", {
  #test_fit()
  g <- generateMultiConnectedDAG(8)
  inputs <- V(g)[igraph::degree(g, mode="in") == 0]
  outputs <- V(g)[igraph::degree(g, mode="out") == 0]
  input.val.list <- lapply(inputs, function(input) runif(1000))
  input.df <- data.frame(input.val.list)
  names(input.df) <- inputs
  output.list <- lapply(outputs, function(output) rep(NA, 1000))
  output.df <- data.frame(output.list)
  names(output.df) <- outputs
  g <- initializeGraph(g, input.table = input.df, 
                       output.table = output.df)
  #Set the observed values exactly to the predicted values
  V(g)[type == "output"]$observed <- V(g)[type == "output"]$output.signal
  output.df[, paste(outputs)] <- unlist(V(g)[type == "output"]$observed)
  # After fitting I expect no changes
  g2 <- fitInitializedNetwork(g, .05, 3)
  # The difference between predicted and observed should still be 0
  sum((unlist(V(g)[type == "output"]$observed) - 
         unlist(V(g)[type == "output"]$output.signal))^2) %>%
    expect_equal(0)
  # The weights should be unchanged
  expect_equal(E(g2)$weight, E(g)$weight)
})

### On titanic
data(titanic3)
titan <- filter(titanic3, !is.na(age), !is.na(survived), !is.na(fare)) %>% #Not worrying about NA vals
  mutate(survived = as.numeric(survived)) %>%
  select(age, survived, fare) %>%
  rescale_df %$% #Note, everything is rescaled to between 0 and 1
  df

test_that("multi-layer model has less loss than nls given it has more parameters.", {
  g <- mlp_graph(c("age", "fare"), "survived", c(2, 1)) %>%
    initializeGraph(input.table = select(titan, age, fare), 
                    output.table = select(titan, survived))
  fit <- fitInitializedNetwork(g,  epsilon = .01, verbose = T)
  nls_fit <-  nls(survived ~ logistic(w0 + w1 * age + w2 * fare), data = titan, 
                  start = as.list(structure(E(g)[to("H11")]$weight, names = c("w1", "w2", "w0"))))
  expect_less_than(get_deviance(fit), deviance(nls_fit))
})

test_that("model should perform a reasonable MLP prediction on a toy problem with single output.", {
  # In the future, expand to multivariate case
  #g <- mlpgraph(c("I1", "I2"), c(3, 2, 4), c("AND", "OR", "NOR"))
  g <- mlp_graph(c("I1", "I2"), "AND", 2)
  system <- list(I1 = c(1, 0), I2 = c(1, 0)) %>%
    expand.grid %>%
    `names<-`(c("I1", "I2")) %>%
    mutate(AND = I1 * I2)
    #mutate(AND = I1 * I2, OR = (I1 + I2 > 0) * 1, NOR = (I1 + I2 == 0) * 1)
  fit <- fitNetwork(g, input.table = select(system, I1, I2), 
                    output.table = select(system, AND), verbose = T)
  expect_equal(unlist(V(fit)["AND"]$output.signal), unlist(V(fit)["AND"]$observed), tolerance = .1)
})
