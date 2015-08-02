# devtools::load_all("../../R/tools.R")
# devtools::load_all("../../R/templates.R")
#devtools::load_all("R/tools.R")
#devtools::load_all("R/templates.R")
context("templates for commmonly used signal graph models")
library(dplyr, quietly = TRUE)

test_that("get_gate with all outputs works as expected", {
  gates <- c("AND", "OR", "NAND", "NOR", "XOR", "XNOR") 
  g <- get_gate(layers = c(2, 3))
  V(g)[is.observed]$name %>%
    setdiff(c(gates, "I1", "I2")) %>%
    {length(.) == 0 } %>%
    expect_true  
})

test_that("get_gate with partial outputs works as expected", {
  gates <- c("AND", "NAND", "XNOR") 
  g <- get_gate(outputs = gates, layers = c(2, 3))
  V(g)[is.observed]$name %>%
    setdiff(c(gates, "I1", "I2")) %>%
    {length(.) == 0 } %>%
    expect_true  
})

test_that("get_gate produces the expected outputs on 2 input system", {  
  # THis is what the output data frame should look like
  logic_gates <- expand.grid(list(I1 = c(0, 1), I2 = c(0, 1))) %>% 
    {dplyr::mutate(., AND = (I1 * I2 == 1) * 1, #dplyr has to be prefixed here for the code to work for some reason.
           OR = (I1 + I2 > 0) * 1 ,
           NAND = (!AND) * 1,
           NOR = (!OR) * 1,
           XOR = (I1 + I2 == 1) * 1, 
           XNOR = (I1 == I2) * 1)}
  # So recreate this data frame from the function output and compare
  get_gate(layers = c(3, 3)) %>% 
    recover_design %>% #The outputs in the design should be the same as the logic_gates table
    identical(logic_gates) %>%
    expect_true
})

test_that("get_gate replicates a hand made version", {
  system <- expand.grid(list(I1 = c(0, 1), I2 = c(0, 1))) %>% 
    mutate(AND = (I1 * I2 == 1) * 1)
  g1 <- get_gate("AND", c(3, 2))
  g2 <- mlp_graph(c("I1", "I2"), "AND", c(3, 2)) %>% #Use a 2 layer MLP
    initializeGraph(select(system, I1, I2, AND), fixed = c("I1", "I2"))
  expect_equal(V(g1)$name, V(g2)$name)
})

test_that("sim_system produces a model that looks like a fitted model", {
  g <- sim_system(5, 100)
  g_new <- fit_initialized_sg(g)
  expect_identical(get_fitted(g), get_fitted(g_new))
})

test_that("sim_system_data produces data gaussian error with gaussian error", {
  gauss_error <- function(x, err){
    x <- x + rnorm(length(x), sd = err)
    (x - min(x)) / max(x)
  }
  expect_true(FALSE)
})

test_that("sim_system enables control of edge density", {
  1:10 %>%# number of graphs to generate
    lapply(sim_system(10, edge_density = .3)) %>%
    lapply(edge_density) %>%
    unlist %>%
    mean %>%
    expect_equal(.3)
})

test_that("sim_system enables control of the proportion of edge weights that are 0", {
  edge_sparcity <- function(g) sum(E(g)weight == 0) / ecount(g)
  1:10 %>%# number of graphs to generate
    lapply(sim_system(10, edge_sparcity = .3)) %>%
    lapply(edge_sparcity) %>%
    unlist %>%
    mean %>%
    expect_equal(.3)
})


test_that("sim_system_data produces data that follows the graph mapping", {
  set.seed(14)
  g <- sim_system(10, 100)
  .data <- sim_system_data(g)
  wts <- E(g)[to("3")]
  expect_identical(.data[, "3"] <- .data[, c("1", "8", "9", "10")])
})

source_gist("https://gist.github.com/robertness/ecb79071301924f43d09")
mapk_g <- g
test_that("power_signal_graph works on the mapk network", {
    #Code for grabbing MapK phosphorylations from KEGG
    g <- power_signal_graph(mapk_g, 25)
    V(g)[igraph::degree(g) == 0] %>%
      length %>%
      expect_equal(0)
})

expect_not_too_ungaussian <- function(x){
  x <- (x - mean(x)) / sd(x)
  test <- ks.test(x, pnorm)
  expect_true(test$p.value > .2)
}

test_that("sim_system_from_graph produces a graph that can sim data", {
  g <- sim_system_from_graph(mapk_g, 25, 100, error_sd = .2)
  expect_identical(vcount(g), 25)
  observed_and_random <- intersect(V(g)[is.observed], V(g)[is.random])
  lapply(V(g)[observed_and_random], function(v){
    err <- unlist(V(g)[v]$output.signal) - unlist(V(g)[v]$observed) 
    expect_true(sum(err) != 0) # There should be some error 
    expect_not_too_ungaussian(err)
  })
})
