
context("templates for commmonly used signal graph models")
library(dplyr, quietly = TRUE)
source("test_dir/test_helpers.R")


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

context("Simulating a system and data from the system")

# When simulating a system, we want a gold standard that is identical 
# to a fitted model in terms of data structure.  The simulation simulates
# weights, generates fitted 'output.signal' values based on those weights, then 
# enters values for observed values from those weight.  

data(mapk_g)
test_that("sim_system produces observed values that are the same as the 
          'fitted' (output.signal) values.", {
            list(g1 = sim_system(10, 100), 
                 g2 = sim_system(10, 100, mapk_g)) %>%
              lapply(function(g){
                observed <- recover_design(g) 
                fitted <- get_fitted(g)[, names(observed)] 
                expect_identical(observed, fitted)
              })
          })

test_that("we can add Gaussian error to get realistic data.", {
  fixed <- V(g)[is.fixed]
  observed_and_random <- intersect(V(g)[is.observed], V(g)[is.random])
  unrandom_data <- sim_system(10, 100) %>%
    recover_design 
  random_data <- unrandom_data %>%
    add_gauss_error(V(g)[observed_and_random])
  for(v in names(unrandom_data)){
    if(v %in% fixed$name){
      expect_equal(rep(0, nrow(unrandom_data)), abs(unrandom_data[,v] - random_data[, v]))
    }else{
      expect_true(its_not_too_ungaussian(abs(unrandom_data[,v] - random_data[, v])))
    }
  }
  
  
  x <- seq(from = 0, to = 1, by = .0001) %>% # Typically we scale to {0, 1}. Suppose {0, 1} were the original scale
    {. + rnorm(length(.), .2)} %>% # Then this is adding Gaussian error on the original scale
    {(. - min(.))/(max(.) - min(.))}  # An this is the scaling step.
  # The simulations add error in this way.
  # The error can still be seen as Gaussian if the {0, 1} bounds are not viewed as proper
  # bounds but just the consequence of a linear transformation that enables model fitting.
  expect_equal(.2, add_gauss_error(x, .2), tolerance = 0.02)
})

test_that("we can simulate a dataset with error", {
  # Data simulation pulls the fitted values, adds error, and returns the observed values. 
  list(g1 = sim_system(10, 1000), 
       g2 = sim_system(10, 1000, mapk_g)) %>%
    lapply(function(g){
      fixed <- V(g)[is.fixed]
      observed_random <- intersect(V(g)[is.observed], V(g)[is.random])
      sim_data <- sim_data_from_graph(g, gauss, gauss_sd = .2)
      expect_true(all(names(sim_data) %in% V(g)[is.observed]$name))
      test_data <- g %>%
        recover_design %$% 
        observed_random %>%        
        add_gauss_error(.2)
      expect_equal(mean(sim_data[, observed_random][, 1]), mean(test_data[, 1]))
    })
})

test_that("sim_system enables control of the proportion of edge weights that are 0", {
  edge_sparcity <- function(g) sum(E(g)$weight == 0) / ecount(g)
  1:100 %>%# number of graphs to generate
    lapply(sim_system(10, edge_sparcity = .3)) %>%
    lapply(edge_sparcity) %>%
    unlist %>%
    mean %>%
    expect_equal(.3, tolerance = .05)
})