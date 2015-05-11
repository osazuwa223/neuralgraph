library(dplyr)
library(magrittr)
transformations <- data.frame(
  A = runif(5),
  B = runif(5),
  Z = rnorm(5)) %>%
  mutate(
    A_squared = A ^ 2,
    A_root = A ^ .5,
    A_log = log(A),
    A_inverse = A ^ -1,
    AB_linear = .2 * A + .5 * B,
    AB_lin_err = AB_linear + Z,
    AB_rescaled = (AB_linear) /(1 + (AB_linear)), 
    AB_scale1_nl_error =  (AB_lin_err) /(1 + (AB_lin_err)),
    AB_logistic = exp(AB_linear) /(1 + exp(AB_linear)),
    AB_logistic_err =  exp(AB_lin_err) / (1 + exp(AB_lin_err))
  )
devtools::use_data(transformations)

