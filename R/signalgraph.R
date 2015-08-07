#' Titanic Data
#' 
#' Records of the survival of Titanic Passengers and such information as sex, age, fare each person paid, 
#' number of parents/children aboard, number of siblings or spouses aboard, passenger class and other fields.
#' 
#' Variable descriptions and how the data set was constructed: \url{http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/titanic3info.txt}
#' 
#' @source \url{http://biostat.mc.vanderbilt.edu/DataSets} 
#' 
#' @docType data
#' @keywords datasets
#' @name titanic3
#' @usage data(titanic3)
#' @format A data frame with 14 variables and 1309 rows
NULL

#' Transformations
#' 
#' A artificial dataset for experimental purposes.  Contains two uniform random variables and 
#' various transformations on those variables.  Some of the transformations have error added.
#' 
#'  \itemize{
#'  \item A, simulated with runif(5)
#'  \item B, simulated with runif(5)
#'  \item A_squared = A^2
#'  \item A_root = A ^ .5,
#'  \item A_log = log(A),
#'  \item A_inverse = A ^ -1,
#'  \item AB_linear = .2 * A + .5 * B,
#'  \item AB_lin_err = AB_linear + Z,
#'  \item AB_rescaled = (AB_linear) /(1 + (AB_linear)), 
#'  \item AB_scale1_nl_error =  (AB_lin_err) /(1 + (AB_lin_err)),
#'  \item AB_logistic = exp(AB_linear) /(1 + exp(AB_linear)),
#'  \item AB_logistic_err = exp(AB_lin_err) /(1 + exp(AB_lin_err))
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name transformations
#' @usage data(transformations)
#' @format A data frame with 12 variables and 5 rows
NULL

#' KEGG MAPK Signaling Pathway Phosphorylation Network
#' 
#' A network of all the phosphorylation edges pulled from the MAPK Signaling Pathway on KEGG.
#' Edges from the human map ("hsa") are pulled.
#' \url{http://www.kegg.jp/kegg-bin/show_pathway?map=hsa04010&show_description=show}
#' 
#' The script for generating the graph is available here:
#' \url{https://gist.github.com/robertness/ecb79071301924f43d09}
#' 
#' @docType data
#' @keywords datasets
#' @name mapk_g
#' @usage data(mapk_g)
#' @format An igraph object with 60 vertices and 76 edges.  Vertex attributes include "label"
#' ("MAP4K4", "PLA2G4B", "STMN1"...) and "kid" for KEGG ids ("hsa:9448", "hsa:100137049", "hsa:3925",...) 
NULL
