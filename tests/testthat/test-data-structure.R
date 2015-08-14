library(dplyr)
option <- FALSE
source("test_dir/test_helpers.R")

################################################################################
context("Data structure")
################################################################################
test_that("loadCN produces Markov blankets without bias vertices", {
  case <- rand_case(8, 3)
  g <- initializeGraph(case$g, case$data) 
  g_test <- induced.subgraph(g, V(g)[!is.bias]) %>% loadCN
  g <- loadCN(g)
  for(v in V(g_test)$name){
    g_cn <- as.numeric(unlist(V(g)[v]$causal_nbr))
    g_test_cn <- as.numeric(c(imb(g_test, V(g_test)[v]), V(g_test)[v]))
    expect_equal(g_cn , g_test_cn)
  }
})
test_that("initializeGraph generates a signalgraph with all the desirable attributes", {
  case <- rand_case(8, 3)
  g <- case$g
  test_data <- case$data
  g <- initializeGraph(g, test_data)
  expect_true("igraph" %in% class(g))
  c("name", "input.signal", "f.prime.input", "output.signal", "observed", "is.bias", "is.observed", 
    "is.random", "is.fixed",  "is.hidden", "is.root", "is.leaf", "updated") %in% 
    list.vertex.attributes(g) %>% 
    all %>%
    expect_true 
  c("name", "weight", "updated") %in% 
    list.edge.attributes(g) %>%
    all %>%
    expect_true
})

test_that("observed nodes have correct attributes", {
  g <- random_unfit_sg(8, 10)
  observed <- V(g)[is.observed]
  observed %>% receives_input(g) %>% expect_one_or_more # some must be at least 
  #one random var which must recieve input
  observed$is.observed %>% expect_all # by definition
  observed$is.hidden %>% expect_none # by definition
  observed$is.bias %>% expect_none # biases are not counted among the 'observed'
  observed$is.root %>% unconstrained  # depends how many are fixed or random
  observed$is.leaf %>% expect_one_or_more # there must be at least one observed node that is a leaf
  observed$is.random %>% expect_one_or_more # At least one observed variable must be random.
  observed$is.fixed %>% expect_not_all # There can be zero fixed or some but not all.
  #observed$is.fixed # Having no fixed variables is possible 
})
test_that("hidden nodes have correct attributes", {
  g <- random_unfit_sg(8, 3)
  hidden <- V(g)[is.hidden]
  hidden %>% receives_input(g) %>% expect_all # input.signal attribute should always have an input signal
  hidden$is.observed %>% expect_none # Not observed by definition
  hidden$is.hidden %>% expect_all # Hidden by definition
  hidden$is.bias %>% expect_none # Not biases by definition
  hidden$is.root %>% expect_none # Hidden nodes cannot be roots
  hidden$is.leaf %>% expect_none # Hidden nodes cannot be leaves 
  hidden$is.random %>% expect_all
  hidden$is.fixed %>% expect_none
})
test_that("biases have correct attributes", {
  g <- random_unfit_sg(8, 3)
  biases <- V(g)[is.bias]
  biases %>% receives_input(g) %>% expect_none # input.signal attribute should be NA
  biases$is.observed %>% expect_none # Not observed
  biases$is.hidden %>% expect_none # Not hidden
  biases$is.bias %>% expect_all # obviously
  biases$is.root %>% expect_all # is a root
  biases$is.leaf %>% expect_none # is not a leaf
  biases$is.random %>% expect_none
  biases$is.fixed %>% expect_none
})
test_that("root nodes have correct attributes", {
  g <- random_unfit_sg(8, 3)
  roots <- V(g)[is.root]
  roots %>% receives_input(g) %>% expect_none # roots can't have input.signal 
  roots$is.observed %>% expect_not_all # because biases are roots and are not observed
  roots$is.hidden %>% expect_none # Roots cannot be hidden
  roots$is.bias %>% any %>% expect_true # There must at least be one root that is a bias
  roots$is.root %>% expect_all # obviously
  roots$is.leaf %>% expect_none # Leaves cannot be roots
  roots$is.random %>% expect_none # Roots are either fixed or bias and cannot be random
  roots$is.fixed %>% expect_not_all # There must be at least one bias.
})
test_that("leaf nodes have correct attributes", {
  g <- random_unfit_sg(8, 3)
  leaves <- V(g)[is.leaf]
  leaves %>% receives_input(g) %>% expect_all # Leaves must have input.signal attribute 
  leaves$is.observed %>% expect_all # Leaves must be observed
  leaves$is.hidden %>% expect_none # Leaves cannot be hidden
  leaves$is.bias %>% expect_none # Leaves cannot be biases
  leaves$is.root %>% expect_none # Leaves cannot be roots
  leaves$is.leaf %>% expect_all # obviously
  leaves$is.random %>% expect_all # all the leaves must be random
  leaves$is.fixed %>% expect_none # none of the leaves can be fixed
})
test_that("random nodes have correct attributes", {
  g <- random_unfit_sg(8, 3)
  randoms <- V(g)[is.random]
  randoms %>% receives_input(g) %>% expect_all # random values must have input.signal 
  #hidden check omitted, some, or all can be hidden
  randoms$is.observed %>% expect_one_or_more
  randoms$is.hidden %>% expect_not_all # random nodes shouldn't all be hidden
  randoms$is.bias %>% expect_none # random nodes can't be biases
  randoms$is.root %>% expect_none # random nodes can't be roots
  randoms$is.leaf %>% expect_one_or_more # random nodes must be leaves
  randoms$is.random %>% expect_all # by definition
  randoms$is.fixed 
})
test_that("fixed nodes have correct attributes", {
  g <- random_unfit_sg(8, 3)
  fixed <- V(g)[is.fixed]
  fixed %>% receives_input(g) %>% expect_none # Fixed variables are roots and thus
  # cannot recieve signal
  fixed$is.observed %>% expect_all #  observed by definition
  fixed$is.hidden %>% expect_none # fixed variable nodes cannot be hidden
  fixed$is.bias %>% expect_none # distinguishing fixed variables from bias constant but technically
  # a bias can be considered a fixed variable
  fixed$is.root %>% expect_all # fixed variable nodes must be roots
  fixed$is.leaf %>% expect_none # fixed variable nodes cannot be leaves
  fixed$is.random %>% expect_none
  fixed$is.fixed %>% expect_all
})

test_that("MLP has expected boolean attributes.", {
  data(titanic3)
  titan <- dplyr::filter(titanic3, !is.na(age), !is.na(survived), !is.na(fare)) %>% #Not worrying about NA vals
    {dplyr::mutate(., survived = as.numeric(survived))} %>%
    {dplyr::select(., age, survived, fare)} %>%
    rescale_df %$% #Note, everything is rescaled to between 0 and 1
    df
  g <- mlp_graph(c("age", "fare"), "survived", c(2, 1)) %>%
    {initializeGraph(., data = dplyr::select(titan, age, fare, survived), fixed = c("age", "fare"))} 
  V(g)[is.observed]$name %>% expect_equal(c("age", "fare", "survived"))
  V(g)[is.hidden]$name %>% expect_equal(c("H11", "H12", "H21"))
  V(g)[is.bias]$name %>% expect_equal(c("bias_H11", "bias_H12", "bias_H21", "bias_survived"))
  V(g)[is.root]$name %>% expect_equal(c("age", "fare", "bias_H11", "bias_H12", "bias_H21", "bias_survived"))
  V(g)[is.leaf]$name %>% expect_equal("survived")
  V(g)[is.random]$name %>% expect_equal(c("H11", "H12", "H21", "survived"))
  V(g)[is.fixed]$name %>% expect_equal(c("age", "fare"))
})
test_that("fixed variables become roots, random variables become not-roots, and get biases ", {
  expect_equal(names(formals(initializeGraph)), c("g", "data", "fixed", "graph_attr"))
  expect_equal(names(formals(fitNetwork)), c("g", "data", "fixed", "graph_attr", "epsilon", "max.iter"))
  data(mtcars)
  my_data <- mtcars %>% 
    rescale_df %>%
    .[[1]]
  fixed <- c("vs", "am", "gear")
  g <- lucy::sim_DAG(8)
  random <- setdiff(names(my_data), fixed)
  V(g)$name <- random
  g <- g + igraph::vertices(fixed)
  g <- g + igraph::edges(c("vs", "mpg", "vs", "hp", "am", "carb", "gear", "drat", "gear", "qsec"))
  g <- initializeGraph(g, my_data, fixed = fixed)
  expect_all(!V(g)[random]$is.root)
  expect_all(V(g)[random] == V(g)[is.random])
  expect_true(all(V(g)[fixed]$is.root))
})
test_that("an error is thrown if there are not at least 2 layers (roots and leaves must be mutually exclusive),
          part of this means it will not work on a single node graph input.", {
            g1 <- graph.empty(1) %>% name_vertices
            data1 <- data.frame(runif(10)) %>% `names<-`("1")
            g2 <- graph.empty(2) %>% name_vertices
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
    {expect_error(initializeGraph(., test_data, fixed = NULL), "Vertices must be named.")}
})
test_that("an error is thrown if the data contains variables not in the graph", {
  case <- rand_case(8, 3)
  g <- case$g
  test_data <- case$data
  test_data %>% 
    transform(new_var = runif(nrow(test_data))) %>%
    {expect_error(initializeGraph(g, .), "\n  Data contains variables that are not named in the graph.\n")}
})
test_that("an error is thrown if the 'fixed' argument contains variables not in the data", {
  case <- rand_case(8, 3)
  g <- case$g
  test_data <- case$data
  fixed <- c(names(test_data)[1:3], "foo")
  expect_error(initializeGraph(g, test_data, fixed = fixed), 
               "\n  Specified fixed variables not found in the data\n")
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

test_that("initializeGraph makes the output signal of fixed variables equivilent to their observed values",{
  case <- rand_case(15, 3)
  fixed <- get_roots(case$g)
  g <- initializeGraph(case$g, case$data,fixed = fixed) 
  lapply(fixed, function(v){
    expect_identical(V(g)[v]$observed, V(g)[v]$output.signal)
  })
})

test_that("if a leaf is not observed in the data, an error is thrown", {
  case <- rand_case(8, 3)
  g <- case$g
  test_data <- case$data %>%
    {`[`(.,setdiff(names(.), V(g)[get_leaves(g)[1]]))}
  expect_error(initializeGraph(g, test_data), "Graph leaves must be observed in the data.")
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
    initializeGraph(test_data, fixed = NULL) %>%
    resetUpdateAttributes %T>%
    {V(.)[is.root]$updated %>%
      all %>%
      expect_true} %T>%
    {V(.)[!is.root]$updated %>%
      all %>%
      expect_false}
})

test_that("initializeGraph returns a graph structure ready for fitting.", {
  g <- random_unfit_sg(10, 2) # This generates a graph, and calls 'initializeGraph in the final step.
  # Fails if you initialize twice
  g_data <- recover_design(g)
  expect_error(
    initializeGraph(g, g_data),
    "Input graph has graph attributes reserved for signal graph."
  )
  # Need weight, name, and updated edge attributes.
  list.edge.attributes(g) %>%
    identical(c("weight", "name", "updated")) %>% 
    expect_true 
  # updated attribute should all be false
  expect_true(!E(g)$updated %>% all)
  # Need the following vertex attributes 
  setdiff(list.vertex.attributes(g), c("name", "input.signal", "f.prime.input", "output.signal",
                                       "observed", "is.bias", "is.observed", "is.hidden", 
                                       "is.root", "is.leaf", "is.random", "is.fixed", "updated")) %>%
  {length(.) == 0} %>%
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
  setdiff(list.graph.attributes(g), c("activation", "L1_pen", "L2_pen", "activation", 
                              "activation.prime", "min.max.constraints", "n")) %>%
    {length(.) == 0} %>%
    expect_true
})

test_that("adjusting 'min.max.constraint' argument gives broader initial weights.", {
  case <- rand_case(5)
  .data <- case$data
  case$g %>% 
    initializeGraph(.data, graph_attr = list(min.max.constraints = c(-100, 100))) %>%
    {E(.)$weight}    
})

test_that("fitNetwork returns a graph structure", {
  long_test(option)
  g <- random_unfit_sg(3, 2)
  observed <- recover_design(g)
  g_structure <- get_structure(g)
  g_fit <- fitNetwork(g_structure, observed, max.iter = 1)
  expect_true(class(g_fit) == "igraph")
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
  g2 <- update_edges(g1, get_determiners = get_dependent_edges, 
                     callback = fit_weights_for_edge_target)
  # We know an edge is traversed when fit_weights_for_edge_target sets 'updated' attribute to 'TRUE'.
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
    {initializeGraph(., data = dplyr::select(titan, age, fare, survived), fixed = c("age", "fare"))}
  V(g)[is.root]$updated <- FALSE
  expect_error(fit_initialized_sg(g,  epsilon = .01, max.iter = 1), "Roots should not have FALSE value for updated attribute.")
})

test_that("model reproduces logistic activation connection between nodes.", { 
  g <- data.frame(from = c("A", "B", "C", "D"),
                  to = c("C", "C", "E", "E")) %>%
    graph.data.frame 
  w.0_c = .2; w.a_c = .3; w.b_c = .5; w.0_e = -.2; w.c_e = -4; w.d_4 = 5
  .data <- data.frame(A = runif(10),
                      B = runif(10),
                      D = runif(10)) %>%
    mutate(C = logistic(w.0_c + w.a_c * A + w.b_c * B ),
           E = logistic(w.0_e + w.c_e * C + w.d_4 * D))
  fit <- fitNetwork(g, .data, graph_attr = c(L1_pen = 0, L2_pen = 0))
  expect_equal(getMSE(fit), 0)  
})

test_that("Randomly selected weights are within min-max constraints", {
  case <- rand_case(20)
  .data <- case$data
  case$g %>%
    initializeGraph(.data, graph_attr = list(min.max.constraints = c(-8, 8))) %>%
    {E(.)$weight} %>%
    summary %T>%
    {expect_less_than(.["Max."], 8)} %T>%
    {expect_more_than(.["Min."], -8)}
})