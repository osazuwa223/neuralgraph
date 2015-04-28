context("Neural network implementation")

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
                  activation.prime = logistic.prime,
                  min.max.constraints = c(min = -Inf, max = Inf))
 # Fails if you initialize twice
 expect_error(
   initializeGraph(g, input.table, output.table,
                       activation = logistic, 
                       activation.prime = logistic.prime,
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
                  activation.prime = logistic.prime,
                  min.max.constraints = c(min = -Inf, max = Inf),
                  verbose=T)
  expect_true(class(g) == "igraph")
})

test_that("in an initialized model where the predicted and observed output 
are exactly the same, the model results of fitting the model should 
be a graph with 0 error and unchanged weights.", {
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
  E(g2)$weight %>% identical(E(g)$weight) %>% expect_true  
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
  system <- list(c(1, 0), c(1, 0)) %>%
    expand.grid %>%
    `names<-`(c("I1", "I2")) %>%
    mutate(AND = I1 * I2)
  g1 <- mlpgraph(c("I1", "I2"), c(3, 2), c("AND")) %>% #Use a 2 layer MLP
    initializeGraph(input.table = system[, c("I1", "I2")], 
                    output.table = system[, "AND", drop = F])
  g2 <- updateEdges(g1, getDeterminers = getDependentEdges, callback = fitWeightsForEdgeTarget)
  # We know an edge is traversed when fitWeightsForEdgeTarget sets 'updated' attribute to 'TRUE'.
  c(!E(g1)$updated,  # All edges should initially be unupdated
    E(g2)$updated) %>% # All edges should finally be updated
    all %>%
    expect_true
})

#Check out functionals deriv, 
test_that("a gradient descent step reduces loss with numeric derivatives.", {
  system <- list(c(1, 0), c(1, 0)) %>%
    expand.grid %>%
    `names<-`(c("I1", "I2")) %>%
    mutate(AND = I1 * I2)
  g <- matrix(c("I1", "AND",
                "I2", "AND"), byrow = T, ncol = 2) %>%
    graph.data.frame %>%
    initializeGraph(input.table = system[, c("I1", "I2")], 
                    output.table = system[, "AND", drop = F]) 
})

test_that("a gradient descent step reduces loss in a MLP case.", {
  system <- list(c(1, 0), c(1, 0)) %>%
    expand.grid %>%
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
  stop()
})

test_that("model should perform a reasonable MLP prediction on a toy problem with single output.", {
  require(dplyr)
  #g <- mlpgraph(c("I1", "I2"), c(3, 2, 4), c("AND", "OR", "NOR"))
  g <- mlpgraph(c("I1", "I2"), c(3, 2), c("AND"))
  #igraphviz(g)
  system <- list(c(1, 0), c(1, 0)) %>%
    expand.grid %>%
    `names<-`(c("I1", "I2")) %>%
    mutate(AND = I1 * I2)
    #mutate(AND = I1 * I2, OR = (I1 + I2 > 0) * 1, NOR = (I1 + I2 == 0) * 1)
  fit <- fitNetwork(g, input.table = system[, c("I1", "I2")], 
                    output.table = system[, "AND", drop = F], 
                    epsilon = 2, verbose = T)
              matrix(ncol = 2, byrow = T) %>%
              graph.edgelist
            list(c(1, 0), c(1, 0)) %>%
              expand.grid %>%
              `names<-`(c("A", "B")) %>%
              mutate(AND = A * B, OR = (A + B > 0) * 1, NOR = (A + B == 0) * 1)
            
            inputs <- V(g)[igraph::degree(g, mode="in") == 0]
            outputs <- V(g)[igraph::degree(g, mode="out") == 0]
            input.df <- lapply(c("A", "B"), function(input) c(1,0)) %>%
              expand.grid %>%
              `names<-`(V(g)[inputs]$name) %>% 
              mutate(AND = ifelse())
            
            
            input.df <- data.frame(input.val.list)
            names(input.df) <- inputs
            output.list <- lapply(outputs, function(output) rep(NA, 1000))
            output.df <- data.frame(output.list)
            names(output.df) <- outputs
            g <- initializeGraph(g, input.table = input.df, 
                                 output.table = output.df)
            })

test_that("for a simple muli-node input, single-node output graph, calculateVals should reproduce
simple arithmetic.", {
  stop()
})

test_that("no unexpected input to resetUpdateAttributes, specifically one where all nodes or
edges are updated == TRUE or all are FALSE", {
  stop()
})

test_that("adding one bias node per non-input node should have same parameter count as a comparable
multi-layer perceptron.", {
  stop()
})

test_that("values and updated status for input nodes and bias nodes should never change", {
  stop()
})

test_that("works when no derivitive of activation function is provided.", {
  stop()
})

test_that("fitNetwork returns a MSPR on the infert dataset that is close to 
          that of a network of the same shape fit with the neuralnet package", {
  stop()
})

test_that("fitting of the model works like it would in lm or glm", {
  stop()
})

test_that("fitted and predict works like they would work in lm or glm (see getDF and newDataUpdate)", {
  stop()
})

test_that("fetching of residuals work like they would in lm or glm", {
  stop()
})