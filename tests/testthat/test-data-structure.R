library(dplyr)
devtools::load_all("../../R/optimization.R")
devtools::load_all("../../R/tools.R")
option <- TRUE
#devtools::load_all("R/optimization.R")
#devtools::load_all("R/templates.R")
context("Signal graph data structure")

test_that("initializeGraph generates a signalgraph with all the desirable attributes", {
  case <- rand_case(8, 3)
  g <- case$g
  test_data <- case$data
  g <- initializeGraph(g, test_data)
  expect_true("igraph" %in% class(g))
  c("name", "input.signal", "f.prime.input", "output.signal", "observed", "is.bias", "is.observed", "is.hidden", "is.root", "is.leaf", "updated") %in% 
    list.vertex.attributes(g) %>% 
    all %>%
    expect_true 
  c("name", "weight", "updated") %in% 
    list.edge.attributes(g) %>%
    all %>%
    expect_true
})

none <- function(x) !all(x)
test_that("observed nodes have correct attributes", {
  g <- random_unfit_sg(8, 3)
  observed <- V(g)[is.observed]
  observed$input.signal %>% lapply(is.na) %>% unlist %>% any %>% expect_true # assuming there should be at least one non-root
  observed$is.observed %>% all %>% expect_true # by definition
  observed$is.hidden %>% none %>% expect_true # by definition
  observed$is.bias %>% none %>% expect_true # biases are not counted among the 'observed'
  observed$is.root %>% any %>% expect_true # there must be at least one observed node that is a root
  observed$is.leaf %>% any %>% expect_true # there must be at least one observed node that is a leaf
  observed$updated %>% all %>% expect_true # intialize graph should have made 
})
test_that("hidden nodes have correct attributes", {
  g <- random_unfit_sg(8, 3)
  hidden <- V(g)[is.hidden]
  hidden$input.signal %>% lapply(is.na) %>% unlist %>% none %>% expect_true # input.signal attribute should always have an input signal
  hidden$is.observed %>% none %>% expect_true # Not observed by definition
  hidden$is.hidden %>% all %>% expect_true # Hidden by definition
  hidden$is.bias %>% none %>% expect_true # Not biases by definition
  hidden$is.root %>% none %>% expect_true # Hidden nodes cannot be roots
  hidden$is.leaf %>% none %>% expect_true # Hidden nodes cannot be leaves 
  hidden$updated %>% all %>% expect_true # always have updated set to TRUE
})
test_that("biases have correct attributes", {
  g <- random_unfit_sg(8, 3)
  biases <- V(g)[is.bias]
  biases$input.signal %>% lapply(is.na) %>% unlist %>% all %>% expect_true # input.signal attribute should be NA
  biases$is.observed %>% none %>% expect_true # Not observed
  biases$is.hidden %>% none %>% expect_true # Not hidden
  biases$is.bias %>% all %>% expect_true # obviously
  biases$is.root %>% all %>% expect_true # is a root
  biases$is.leaf %>% none %>% expect_true # is not a leaf
  biases$updated %>% all %>% expect_true # always have updated set to TRUE
})
test_that("root nodes have correct attributes", {
  g <- random_unfit_sg(8, 3)
  roots <- V(g)[is.root]
  roots$input.signal %>% lapply(is.na) %>% unlist %>% all %>% expect_true # roots can't have input.signal 
  roots$is.observed %>% any %>% expect_true # Roots must have at least on non-bias, which must be observed
  roots$is.hidden %>% none %>% expect_true # Roots cannot be hidden
  roots$is.bias %>% any %>% expect_true # There must at least be one root that is a bias
  roots$is.root %>% all %>% expect_true # obviously
  roots$is.leaf %>% none %>% expect_true # Leaves cannot be roots
  roots$updated %>% all %>% expect_true # always have updated set to TRUE
})
test_that("leaf nodes have correct attributes", {
  g <- random_unfit_sg(8, 3)
  leaves <- V(g)[is.leaf]
  leaves$input.signal %>% lapply(is.na) %>% unlist %>% none %>% expect_true # Leaves must have input.signal attribute 
  leaves$is.observed %>% all %>% expect_true # Leaves must be observed
  leaves$is.hidden %>% none %>% expect_true # Leaves cannot be hidden
  leaves$is.bias %>% none %>% expect_true # Leaves cannot be biases
  leaves$is.root %>% none %>% expect_true # Leaves cannot be roots
  leaves$is.leaf %>% all %>% expect_true # obviously
  leaves$updated %>% all %>% expect_true # always have updated set to TRUE
})
test_that("an error is thrown if there are not at least 2 layers (roots and leaves must be mutually exclusive),
          part of this means it will not work on a single node graph input.", {
            g1 <- graph.empty(1) %>% nameVertices
            data1 <- data.frame(runif(10)) %>% `names<-`("1")
            g2 <- graph.empty(2) %>% nameVertices
            data2 <- data.frame(runif(10), runif(10)) %>% `names<-`(c("1", "2"))
            error = "There has to be at least 2 vertices and 1 directed edge."
            expect_error(initializeGraph(g1, data1), error)
            expect_error(initializeGraph(g2, data2), error)
})
test_that("recover_design gets the data back", {
  case <- rand_case(8, 3)
  g <- case$g
  test_data <- case$data
  g %>%
    initializeGraph(test_data) %>%
    recover_design %>%
    `[`(names(test_data)) %>% # sort names
    expect_identical(test_data)
})
test_that("for initializeGraph an error is thrown if the input graph does not have named vertices", {
  case <- rand_case(8, 3)
  g <- case$g
  test_data <- case$data  
  remove.vertex.attribute(g, "name") %>% 
    {expect_error(initializeGraph(., test_data), "Vertices must be named.")}
})
test_that("an error is thrown if the data contains variables not in the graph", {
  case <- rand_case(8, 3)
  g <- case$g
  test_data <- case$data
  test_data %>% 
    transform(new_var = runif(nrow(test_data))) %>%
    {expect_error(initializeGraph(g, .), "\n  Data contains variables that are not named in the graph.\n")}
})
test_that("if a vertex is not in a variable in the data, it becomes a hidden variable", {
  case <- rand_case(8, 3)
  g <- case$g
  test_data <- case$data
  hidden <- setdiff(V(g)$name, names(test_data))
  g %>%
    initializeGraph(test_data) %>%
    {V(.)[hidden]$is.hidden} %>%
    all %>% 
    expect_true
})
test_that("if a root or a leaf is not observed in the data, an error is thrown", {
  case <- rand_case(8, 3)
  g <- case$g
  test_data <- case$data %>%
    {`[`(.,setdiff(names(.), V(g)[get_roots(g)[1]]))}
  expect_error(initializeGraph(g, test_data), "Graph roots and leaves must be observed in the data.")
})
test_that("initializeEdges provides a graph with 'weight' edge attribute and if the attribute
          already exists, they should have been changed.", {
            case <- rand_case(8, 3)
            g <- case$g
            w1 <- E(g)$weight
            g <- initializeEdges(g)
            w2 <- E(g)$weight
            expect_true(!identical(w1, w2))
          })
test_that("resetUpdateAttributes changes updated structure of all nodes EXCEPT root nodes", { 
  case <- rand_case(8, 3)
  g <- case$g
  test_data <- case$data
  g %>%
    initializeGraph(test_data) %>%
    resetUpdateAttributes %T>%
    {V(.)[is.root]$updated %>%
      all %>%
      expect_true} %T>%
    {V(.)[!is.root]$updated %>%
      all %>%
      expect_false}
})
test_that("initializeGraph returns a graph structure ready for fitting.", {
  long_test(option)
  g <- random_sg(5, 4, 10) # This generates a graph, and calls 'initializeGraph in the final step.
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
    identical(c("L1_pen", "L2_pen", "activation", "activation.prime", "min.max.constraints", "n")) %>%
    expect_true
})

test_that("fitNetwork returns a graph structure", {
  long_test(option)
  g <- random_sg(3, 2, 10)
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

test_that("test that if the updated status of biases/roots and input nodes ever change, an error is thrown", {
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
