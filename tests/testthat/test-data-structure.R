devtools::load_all("../../R/optimization.R")
devtools::load_all("../../R/tools.R")
option <- TRUE
#devtools::load_all("R/optimization.R")
context("Signal graph data structure")

test_that("initializeGraph returns a graph structure ready for fitting.", {
  set.seed(21)
  g <- generateMultiConnectedDAG(5)
  inputs <- c("2", "3")
  outputs <- c("4")
  input.table <- as.data.frame(matrix(runif(1000 * 2), ncol = 2,
                                      dimnames = list(NULL, inputs)))
  names(input.table) <- inputs
  output.table <- input.table %>%
    as.matrix %>% 
    {. %*% matrix(c(1.7, .8), ncol = 1)} %>%
    {(function(x) x / (1 + x))(.)} %>%
    data.frame %>%
    `names<-`("4")
  g <- initializeGraph(g, input.table, output.table, 
                       activation = logistic,
                       activation.prime = logistic_prime,
                       min.max.constraints = c(min = -Inf, max = Inf))
  # Fails if you initialize twice
  expect_error(
    initializeGraph(g, input.table, output.table,
                    activation = logistic, 
                    activation.prime = logistic_prime,
                    min.max.constraints = c(min = -Inf, max = Inf)),
    "This graph structure seems to have already been updated."
  )
  # Need weight, name, and updated edge attributes.
  list.edge.attributes(g) %>%
    identical(c("weight", "name", "updated")) %>% 
    expect_true 
  # updated attribute should all be false
  expect_true(!E(g)$updated %>% all)
  # Need the following vertex attributes 
  list.vertex.attributes(g) %>%
    identical(c("name", "type", "input.signal", "f.prime.input",
                "output.signal", "observed", "updated")) %>%
    expect_true
  # The input.signal, f.prime.input, output.signal, and observed values should
  # should be numerics.
  lapply(
    list(
      V(g)[1]$input.signal[[1]],
      V(g)[1]$f.prime.input[[1]],
      V(g)[1]$output.signal[[1]],
      V(g)[type == "output"]$observed[[1]]
    ), 
    function(item){
      (length(item) > 0) %>%
        expect_true
    }
  )
  # Graph attributes should be activation, activation.prime, min.max.constraints
  list.graph.attributes(g) %>%
    identical(c("activation", "activation.prime", "min.max.constraints", "n")) %>%
    expect_true
})

test_that("fitNetwork returns a graph structure", {
  long_test(option)
  set.seed(21)
  g <- generateMultiConnectedDAG(5)
  inputs <- c("2", "3")
  outputs <- c("4")
  input.table <- as.data.frame(matrix(runif(1000 * 2), ncol = 2,
                                      dimnames = list(NULL, inputs)))
  names(input.table) <- inputs
  output.table <- input.table %>%
    as.matrix %>% 
    {. %*% matrix(c(1.7, .8), ncol = 1)} %>%
    {(function(x) x / (1 + x))(.)} %>%
    data.frame %>%
    `names<-`("4")
  g <- fitNetwork(g, input.table, output.table, 
                  activation = logistic,
                  activation.prime = logistic_prime,
                  min.max.constraints = c(min = -Inf, max = Inf),
                  verbose=T)
  expect_true(class(g) == "igraph")
})


test_that("after the weights of a given vertex has changed, the prediction should change",{
  g <- get_gate("AND")
  original <- V(g)[type=="output"]$output.signal %>% unlist
  updated <- getPrediction(g, V(g)["AND"], runif(3))
  expect_true(!identical(original, updated))
})

test_that("prediction should never produce NA or other invalid values, it should rather error out", {
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
  g <- updateVertices(g, getDeterminers = iparents, callback = calculateVals)
  #Every thing that is not an input or an intercept should work.
  V(g)[!(type %in% c("input", "intercept"))]$output.signal %>% lapply(isValid) %>% lapply(expect_true)
})

test_that("after a pass at fitting, all edges have been traversed.", {
  long_test(option)
  g1 <- get_gate("AND", c(3, 2))
  g2 <- updateEdges(g1, getDeterminers = getDependentEdges, callback = fitWeightsForEdgeTarget)
  # We know an edge is traversed when fitWeightsForEdgeTarget sets 'updated' attribute to 'TRUE'.
  c(!E(g1)$updated,  # All edges should initially be unupdated
    E(g2)$updated) %>% # All edges should finally be updated
    all %>%
    expect_true
})

test_that("test that if the updated status of intercepts/biases and input nodes ever change, an error is thrown", {
  long_test(option)
  data(titanic3)
  titan <- dplyr::filter(titanic3, !is.na(age), !is.na(survived), !is.na(fare)) %>% #Not worrying about NA vals
    {dplyr::mutate(., survived = as.numeric(survived))} %>%
    {dplyr::select(., age, survived, fare)} %>%
    rescale_df %$% #Note, everything is rescaled to between 0 and 1
    df
  g <- mlp_graph(c("age", "fare"), "survived", c(2, 1)) %>%
    {initializeGraph(., input.table = dplyr::select(titan, age, fare), 
                    output.table = dplyr::select(titan, survived))}
  V(g)[type == "input"]$updated <- FALSE
  expect_error(fitInitializedNetwork(g,  epsilon = .01, verbose = T), "Inputs or biases had FALSE for updated.")
})
