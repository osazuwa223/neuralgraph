historic_fits <- list(
  fitA1 = list(
    fit = fit1,
    notes = "First fit based on only 500 cells, and a graph structure (that I am changing as of today,
    Aug 29 2015, to remove edges that add phosphates but inhibit).  No L1 penalty, L2 is .01, constraints
    were -50 to 50. mmse of 0.0309.  Errors generally high, bulk on PIP2 with 1.0624 mse."
    ),
  fitA2 = list(
    fit = fit11,
    notes = "Another iteration on fitA1. mmse of .0301.  Bulk of error on Raf with .0187 mse."
    ),
  fitA3 = list(
    fit = fit111,
    notes = "Another 3 iterations on fitA2.  mmse of .028.  Bulk of error on Mek with .0036 mse.")
)
use_data(historic_fits, overwrite = RUE)
