#' Create a loss function dependent on a given vertex and penalizes the graph based on
#' penalizd least squares.
#' 
#' This is a closure that creates a loss function that varies given a weight vector.
#' The weights are weight attributes of incoming edges to the vertex given in the argument.
#' Thus for a given vertex, you can inspect how varying the weights that determine the vertex affect loss. 