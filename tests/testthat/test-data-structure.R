devtools::load_all("../../R/optimization.R")
devtools::load_all("../../R/tools.R")
option <- FALSE
#devtools::load_all("R/optimization.R")
context("Signal graph data structure")

test_that("an error is thrown if the input graph does not have named vertices", {})
test_that("an error is thrown if the input data is not a named data frame or list", {})


test_that("signal graph objects have vertex attributes 'is.observed', 'is.root', 'is.leaf', and 'is.bias'", {
  sg_list <- list(gate = get_gate(layers = c(2, 3)), random = random_sg(4, 100))
  attribs <- c("is.observed", "is.root", "is.leaf", "is.bias")
  lapply(sg_list, function(g){
    expect_true(all(attribs %in% list.vertex.attributes(g)))
  })
})

test_that("if there are variables in the data that are not named nodes in the graph, an error is thrown",{
  
  
})
#is.observed, - 
#is.root, -input signal is null
#is.leaf leaves must be observed
#is.bias,

test_that("the biases are FALSE for is.observed, TRUE for is.bias, and TRUE for is.root.", {
  
})

test_that("", {
  
})

test_that("initializeGraph returns a graph structure ready for fitting.", {
  set.seed(21)  
  g <- random_sg(4, 5) # This generates a graph, and calls 'initializeGraph in the final step.
  # Fails if you initialize twice
  g_data <- recover_design(g)
  expect_error(
    initializeGraph(g, g_data),
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
      V(g)[is.observed]$observed[[1]]
    ), 
    function(item){
      (length(item) > 0) %>%
        expect_true
    }
  )
  # Graph attributes should be activation, activation.prime, min.max.constraints
  list.graph.attributes(g) %>%
    identical(c("penalty", "activation", "activation.prime", "min.max.constraints", "n")) %>%
    expect_true
})

test_that("fitNetwork returns a graph structure", {
  long_test(option)
  set.seed(21)
  g <- random_sg(3, 2)
  observed <- recover_design(g)
  g <- fitNetwork(g, observed, verbose=T)
  expect_true(class(g) == "igraph")
})


test_that("after the weights of a given vertex has changed, the prediction should change",{
  g <- get_gate("AND")
  original <- V(g)[is.observed]$output.signal %>% unlist
  updated <- getPrediction(g, V(g)["AND"], runif(3))
  expect_true(!identical(original, updated))
})

test_that("prediction should never produce NA or other invalid values, it should rather error out", {
  g <- random_sg(3, 2)
  observed <- recover_design(g)
  g <- updateVertices(g, getDeterminers = iparents, callback = calculateVals)
  #Every thing that is not an input or an intercept should work.
  V(g)[!is.bias]$output.signal %>% lapply(isValidV) %>% lapply(expect_true)
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

test_that("test that if the updated status of intercepts/roots and input nodes ever change, an error is thrown", {
  long_test(option)
  data(titanic3)
  titan <- dplyr::filter(titanic3, !is.na(age), !is.na(survived), !is.na(fare)) %>% #Not worrying about NA vals
    {dplyr::mutate(., survived = as.numeric(survived))} %>%
    {dplyr::select(., age, survived, fare)} %>%
    rescale_df %$% #Note, everything is rescaled to between 0 and 1
    df
  g <- mlp_graph(c("age", "fare"), "survived", c(2, 1)) %>%
    {initializeGraph(., data = dplyr::select(titan, age, fare, survived))}
  V(g)[is.root]$updated <- FALSE
  expect_error(fitInitializedNetwork(g,  epsilon = .01, verbose = T), "Inputs or biases had FALSE for updated.")
})
