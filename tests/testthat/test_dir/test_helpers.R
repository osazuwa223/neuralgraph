expect_none <- function(x) expect_true(sum(x) == 0)
expect_one_or_more <- function(x) expect_true(any(x))  
expect_all <- function(x) expect_true(all(x))
expect_not_all <- function(x) expect_true(!all(x))
unconstrained <- function(x) expect_true(TRUE)
receives_input <- function(v_set, g) {
  V(g)[v_set]$input.signal %>%
    lapply(function(array) sum(is.na(array)) == 0) %>% 
    unlist
}
is_not_too_ungaussian <- function(x){
  x <- (x - mean(x)) / sd(x)
  test <- ks.test(x, pnorm)
  expect_true(test$p.value > .2)
}
