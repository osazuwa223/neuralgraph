# context("when using NLS for fitting")
# test_that("a collider with no hidden nodes maps to a child to parent formula", {
#   set.seed(10)
#   g <- random_unfit_sg(3, 10) %>%
#     load_formulas
#   expected_formula <- formula("v2 ~ g$activation(w.1_2 * v1 + w.3_2 * v3 + wb_2 *bias_2)")
#   expect_identical(V(g)[2]$formula, expected_formula)
#   expect_true(is.na(V(g)[c("1", "2", "bias_2")]$formula)
# })
# test_that("nls can work with the simple collider formula", {
#   set.seed(10)
#   g <- random_unfit_sg(3, 10) %>%
#     load_formulas
#   expected_formula <- formula("v2 ~ g$activation(w.1_2 * v1 + w.3_2 * v3 + w.b_2 * bias_2)")
#   .data <- recover_design(g) %>%
#     `names<-`(paste0("v", names(.))) %>%
#     mutate(bias_2 = 1)
#   starting_weights <- E(g)[to(2)]$weight %>%
#     as.list %>%
#     structure(names = c("w.1_2", "w.3_2", "w.b_2"))
#   try_nls <- function(.formula, .data, start){
#     tryCatch(nls(.formula, .data, start = .start), 
#              error = function(e) {
#                tryCatch(nls(.formula, .data, start = start, algorithm = "port"),
#                         error = function(e) "try simulating from the starting_weights.")
#              })
#   }
#   fit <- try_nls(expected_formula, .data, start = starting_weights) 
#   expect_true(is.numeric(coef(fit)))
# })
# test_that("adding in one hidden node to the v structure gives the correct formula", {
#   set.seed(17)
#   g <- random_unfit_sg(4, 10) %>%
#     load_formulas 
#   expected_formula <- formula("v4 ~ g$activation(w.2_4 * v2 + w.b_4 * bias_4 + g$activation(w.3_1 * v3 + w.b_1 * bias_1 ))")
#   expect_identical(V(g)[2]$formula, expected_formula)
#   expect_true(is.na(V(g)[c("2", "1", "bias_4", "bias_1")]$formula)
#   sg_viz(g, show_biases = T)
# })
# test_that("nls can work with the recursion implied by the hidden node.", {
#   set.seed(17)
#   g <- sim_system(4, 10, .2) %>%
#     load_formulas
#   expected_formula <- formula("v4 ~ g$activation(w.2_4 * v2 + w.b_4 * bias_4 + g$activation(w.3_1 * v3 + w.b_1 * bias_1 ))")
#   .data <- recover_design(g) %>%
#     `names<-`(paste0("v", names(.))) %>%
#     mutate(bias_4 = 1, bias_1 = 1)
#   starting_weights <- c(E(g)[to("4")]$weight, E(g)[to("1")]$weight) %>%
#     as.list %>%
#     structure(names = c("w.1_4", "w.2_4", "w.b_4", "w.3_1", "w.b_1")) 
#   starting_weights <- lapply(starting_weights, jitter)
#   try_nls <- function(.formula, .data, start){
#     tryCatch(nls(.formula, .data, start = start), 
#              error = function(e) nls(.formula, .data, start = start, algorithm = "port"))
#   }
#   fit <- try_nls(expected_formula, .data, start = starting_weights) 
#   expect_true(is.numeric(coef(fit)))
# })
# 
#   
# test_that("get the correct formula for a single output MLP", {})
# 
# 
# 
# #generation of the design dfs with biases
# #dynamic construction of formula
# #function for creating starting_weight list