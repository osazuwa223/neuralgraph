###################################
# test_that("we can recover the dynamics of a simple linear system.", {
#   w0 <- 2; w1 <- 2.5; v0 <- -1.3; v1 <- 3.2;
#   f <- function(u) u / (1 + u)
#   simple_system <- data.frame(X1 = runif(10)) %>%
#     {dplyr::mutate(., 
#                   X2 = f(w0 + w1 * X1),
#                   X3 = f(v0 + v1 * X2))}
#   g <- data.frame(from = c("X1", "X2"), to = c("X2", "X3")) %>%
#     graph.data.frame 
#   g <- initializeGraph(g, simple_system, simple_system, penalty = .01, activation = f, activation.prime = NULL)
#   
# })