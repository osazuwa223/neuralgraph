devtools::load_all("../../R/tools.R")
devtools::load_all("../../R/templates.R")
#devtools::load_all("R/tools.R")
#devtools::load_all("R/templates.R")
context("templates for commmonly used signal graph models")
library(dplyr)




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

test_that("sim_system produces a model that looks like a fitted model", {})
test_that("the fitted values and the observed values are the same", {})

