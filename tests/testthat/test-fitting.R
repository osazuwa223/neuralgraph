# devtools::load_all("../../R/tools.R")
#devtools::load_all("R/tools.R")
option <- FALSE
library(dplyr, quietly = TRUE)
test_that("model reproduces deterministic connection between nodes.", { 
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


context("Models gets reasonable fits")

test_that("despite not working a single node graph input, regression on a constant works if the constant
          is given explicity.", {
            long_test(option)
            g_const <- graph.empty(2) %>% name_vertices %>% `+`(edge(c("2", "1"))) 
            data2 <- data.frame(runif(10), rep(1, 10)) %>% `names<-`(c("1", "2"))
            fitNetwork(g_const, data2, max.iter = 1)
          })

test_that("in an initialized neural net style model where the predicted and observed output are 
          exactly the same, the model results of fitting the model should be a graph with 0 error 
          and unchanged weights.", {
  long_test(option)
  set.seed(30)
  g <- random_unfit_sg(10, 3)
  #Set the observed values exactly to the predicted values
  V(g)[is.observed]$observed <- V(g)[is.observed]$output.signal
  # After fitting I expect no changes
  g2 <- fit_initialized_sg(g, .05, max.iter = 1)
  # The difference between predicted and observed should still be 0
  expect_equal(getMSE(g2), 0, .001)
  # The weights should be unchanged
  expect_equal(E(g2)$weight, E(g)$weight, tolerance= .01)
})

test_that("I can get descent results on a linear model with hidden nodes", {
  case <- rand_case(20)
  .data <- case$data
  case$g %>%
    initializeGraph(.data, graph_attr = list(min.max.constraints = c(-100, 100))) %>%
    {E(.)$weight} %>%
    summary %T>%
    {expect_less_than(.["Max."], 100)} %T>%
    {expect_more_than(.["Min."], -100)}
})

### On titanic
data(titanic3)
titan <- dplyr::filter(titanic3, !is.na(age), !is.na(survived), !is.na(fare)) %>% #Not worrying about NA vals
  {dplyr::mutate(., survived = as.numeric(survived))} %>%
  {dplyr::select(., age, survived, fare)} %>%
  rescale_df %$% #Note, everything is rescaled to between 0 and 1
  df

test_that("multi-layer model has less loss than nls given it has more parameters.", {
  long_test(option)
  set.seed(33)
  g <- mlp_graph(c("age", "fare"), "survived", c(2, 1)) %>%
    {initializeGraph(., data = dplyr::select(titan, age, fare, survived), fixed = c("age", "fare"))}
  fit <- fit_initialized_sg(g,  epsilon = .01, max.iter = 2)
  nls_fit <-  nls(survived ~ logistic(w0 + w1 * age + w2 * fare), data = titan, 
                  start = as.list(structure(E(g)[to("H11")]$weight, names = c("w1", "w2", "w0"))))
  expect_less_than(vertexMSE(fit, V(g)["survived"]), deviance(nls_fit))
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
  fit <- fitNetwork(g, system, fixed = c("I1", "I2"), max.iter = 2, 
                    graph_attr = list(L2_pen = 0, min.max.constraints = c(-20,20)))
  expect_equal(unlist(V(fit)["AND"]$output.signal), unlist(V(fit)["AND"]$observed), tolerance = .1)
})

test_that("L1 norm with super high penalty should bring weights from non-bias terms to 0", {
  long_test(option)
  set.seed(30)
  g_structure <- mlp_graph(c("age", "fare"), "survived", c(4, 3))
  g_pen <- {fitNetwork(g_structure, dplyr::select(titan, age, fare, survived),
                       fixed = c("age", "fare"),
                       graph_attr = list(L1_pen = 100), epsilon = .01, max.iter = 1)}
  expect_equal(sum(round(E(g_pen)[from(V(g_pen)[!is.bias])]$weight, 4)), 0)  
})


test_that("L2 norm has less sum squares of fitted weight than unpenalized.", {
  long_test(option)
  set.seed(30)
  g_structure <- mlp_graph(c("age", "fare"), "survived", c(4, 3))
  g_no_pen <- {fitNetwork(g_structure, dplyr::select(titan, age, fare, survived), 
                          fixed = c("age", "fare"), epsilon = .01, max.iter = 1)}
  g_pen <- {fitNetwork(g_structure, dplyr::select(titan, age, fare, survived), 
                       fixed = c("age", "fare"),
                       graph_attr = list(L2_pen = 3), 
                       epsilon = .01, max.iter = 1)}
  expect_less_than(sum(E(g_pen)$weight^2), sum(E(g_no_pen)$weight^2))  
})




