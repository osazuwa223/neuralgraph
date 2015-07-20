library(signalgraph)
library(dplyr)

devtools::load_all("R/optimization.R")
devtools::load_all("R/modelfitting.R")

########################################################################################
## Load in the processed T cell data, and add signal nodes
########################################################################################
data(tcells, package = "bninfo")
int_levels <- c("observational", "PKC", "PKA", "PIP2", "Mek", "Akt") 
tcells_int <- tcells$processed$interventions %>% # read in intervention array
  factor(levels = int_levels) %>% # make it a factor 
  {data.frame(int_ = .)} %>% # covert to design matrix and add to other variables
  model.matrix( ~ int_, .) %>%
  .[, - 1] %>% # (removes intercept)
  {cbind(tcells$processed$.data, .)} %>%
  apply(2, function(col) { # Standardize to between 0 and 1
    col <- as.numeric(col)
    if(min(col) != 0) col <- (col - min(col))/(max(col) - min(col))
    col
  }) %>%
  data.frame 
data(tcell_examples, package = "bninfo") #prepare to make intervention nodes
interventions <- c("int_PKC", "int_PKA", "int_PIP2", "int_Mek", "int_Akt")
tcell_input_net <- tcell_examples$net %>%
  {bnlearn::as.graphNEL(.)} %>% 
  igraph.from.graphNEL %>%
  {. + igraph::vertices(interventions)} %>%
  {. + igraph::edges("int_PKC", "PKC", 
                     "int_PKA", "PKA",
                     "int_PIP2", "PIP2",
                     "int_Mek", "Mek",
                     "int_Akt", "Akt")}  
validated_net <- initializeGraph(tcell_input_net, data = tcells_int, fixed = interventions)

  
########################################################################################
## Create the network with upstream nodes hidden
########################################################################################
removable_proteins <- c("PKA", "PKC", "Raf", "Mek", "Plcg", "PIP3", "Erk")
nonremovable_proteins <- setdiff(names(tcells_int), removable_proteins)
partial_data <- tcells_int[, nonremovable_proteins] 
graph_attributes <- list(L1_pen = 0, L2_pen = .01, min.max.constraints = c(10, 10))
subset_example_net <- initializeGraph(tcell_input_net, 
                                  data = partial_data,
                                  fixed = interventions,
                                  graph_attr = graph_attributes) 

########################################################################################
#Get an average of m networks, average the edge weights, 
########################################################################################
num_iterations <- 1
interventions <- c("int_PKC", "int_PKA", "int_PIP2", "int_Mek", "int_Akt")
########
# Missing nodes: "PKA", "PKC", "Raf", "Mek", "Plcg", "PIP3", "Erk"
#######
removable_proteins <- c("PKA", "PKC", "Raf", "Mek", "Plcg", "PIP3", "Erk")
nonremovable_proteins <- setdiff(names(tcells_int), removable_proteins)
partial_data <- tcells_int[, nonremovable_proteins] 
graph_attributes <- list(L1_pen = 0, L2_pen = .01, min.max.constraints = c(-10, 10))
m <- 15
#net_list <- NULL
for(i in 1:m){
  message("Working on #", i, " iteration.")
  net_list_item <-  fitNetwork(tcell_input_net, 
                               data = partial_data,
                               fixed = interventions,
                               graph_attr = graph_attributes,
                               max.iter = num_iterations) 
  net_list <- c(net_list, list(net_list_item))
}

do.call("rbind", lapply(net_list, function(g){
  E(g)$weight
})) %>%
  as.data.frame %>%
  `colnames<-`(E(net_list_item)$name) %>%
  summary

lapply(net_list, function(g){
  g %>%
    betweenness(V(.), weight = abs(E(g)$weight)) %>%
    .[removable_proteins] %>%
    sort(decreasing = TRUE)
})

lapply(net_list, function(g){
  g %>%
    evcent(V(.), directed = TRUE, weights = E(g)$weight) %>%
    .[removable_proteins] %>%
    sort(decreasing = TRUE)
})
########
# Missing nodes: "Raf", "Mek", "Plcg", "PIP3", "Erk"; ("PKA" and "PKC" added back in)
#######
removable_proteins <- c("Raf", "Mek", "Plcg", "PIP3", "Erk")
nonremovable_proteins <- setdiff(names(tcells_int), removable_proteins)
partial_data <- tcells_int[, nonremovable_proteins] 
graph_attributes <- list(L1_pen = 0, L2_pen = .01, min.max.constraints = c(-10, 10))
m <- 15
net_list2 <- NULL
for(i in 1:m){
  message("Working on #", i, " iteration.")
  net_list_item <- finitialized_net <- fitNetwork(tcell_input_net, 
                                                  data = partial_data,
                                                  fixed = interventions,
                                                  graph_attr = graph_attributes, 
                                                  max.iter = num_iterations) 
  net_list2 <- c(net_list2, list(net_list_item))
}




mek_averaging <- tcell_data[, c(base_proteins, "Mek")] %>%
  causal_learning(interventions, "tabu", iss = 10, tabu = 50, resample = TRUE)

mek_erk_averaging <- tcell_data[, c(base_proteins, "Mek", "Erk")] %>%
  causal_learning(interventions, "tabu", iss = 10, tabu = 50, resample = TRUE)

erk_averaging <- tcell_data[, c(base_proteins, "Erk")] %>%
  causal_learning(interventions, "tabu", iss = 10, tabu = 50, resample = TRUE)

erk_raf_averaging <- tcell_data[, c(base_proteins, "Erk", "Raf")] %>%
  causal_learning(interventions, "tabu", iss = 10, tabu = 50, resample = TRUE)

subset_averaging <- list(mek_averaging = mek_averaging,
                         mek_erk_averaging = mek_erk_averaging,
                         erk_averaging = erk_averaging,
                         erk_raf_averaging = erk_raf_averaging)




V(g)[v]$input.signal <- list(linear_combination) # Add it to input signal attribute
if(!is.null(g$activation.prime)){
  V(g)[v]$f.prime.input <- list(g$activation.prime(linear_combination)) # Apply f prime and add to the f prime input signal attribute
}
output <- g$activation(linear_combination) # apply activation function and add it to activation function attribute
V(g)[v]$output.signal <- list(output)
V(g)[v]$output.signal %>% unlist %>% 
  ensure_that(checkVector(.))
g



########################################################################################
## Big sim
########################################################################################


simple_learn <- function(net, .data, algo, score, resample = FALSE,  ...){
  if(resample){
    .data <- dplyr::sample_n(.data,size = nrow(.data), replace = T)
  }
  switch(algo,
         tabu = do.call("tabu", list(x = .data, start = net, score = score, 
                                     ... = ...)),
         hc = do.call("hc", list(x = .data, score = score, start = net,   ... = ...))
  )
}
#' 
subset_learning <- function(.data, algo = "tabu", score = "bic-g", resample = FALSE, ...){
  random.graph(nodes = names(.data), method = "melancon", 50,
               burn.in = 10^5, every = 100) %>%
    lapply(simple_learn, .data = .data, algo = algo, resample = resample, score = score, ...) %>%
    custom.strength(nodes = names(.data))
}  

graph_attributes <- list(L1_pen = .2, L2_pen = .2, prop_sparse = .6)
library(bnlearn)
proposed_scores <- NULL
naive_scores <- NULL
for(i in 1:100){
  print(i)
      g <- sim_system(m = 12, k = 5, n = 100, graph_attr = graph_attributes)
      E(g)$novel <- FALSE
      unacceptable_edges <- E(g)[to(V(g)[is.leaf]) & E(g)$weight == 0]
      acceptable_edges <- setdiff(E(g)[from(V(g)[is.hidden])], unacceptable_edges)
      novel_edge <- E(g) %>%
        .[acceptable_edges] %>%
        .[from(V(g)[is.hidden])] %>%
        as.numeric %>%
        sample(1) %>%
    {E(g)[.]} 
    E(g)[novel_edge]$novel <- TRUE
    data.frame(get.edgelist(g)) %>% mutate(weight = E(g)$weight)
    g_trunc <- g %>%
      delete.edges(novel_edge) %>%
      fitInitializedNetwork(max.iter = 1) 
    
    proposed_result <- betweenness(g_trunc, weight = abs(E(g_trunc)$weight) + .0000001) %>%
      sort(decreasing = TRUE) %>%
      .[V(g_trunc)[is.hidden]] %>%
      names %>%
      .[1:2] %>% 
      {get_fitted(g)[unique(c(V(g)[is.observed]$name, .))]} %>%
      subset_learning("hc", resample = FALSE) %>%
      {dplyr::mutate(., name = bninfo::name_edges(.[,c(1,2)]))}
    proposed_score <- 0
    if(E(g)[novel_edge]$name %in% proposed_result$name){
      proposed_score <- proposed_result[proposed_result$name == E(g)[novel_edge]$name, ][, "strength"]
    } 
    proposed_scores <- c(proposed_scores, proposed_score)
    naive_result <- betweenness(g_trunc, weight = rep(1, ecount(g_trunc))) %>%
      sort(decreasing = TRUE) %>%
      .[V(g_trunc)[is.hidden]] %>%
      names %>%
      .[1:2] %>% 
    {get_fitted(g)[unique(c(V(g)[is.observed]$name, .))]} %>%
      subset_learning("tabu", resample = FALSE) %>%
    {dplyr::mutate(., name = bninfo::name_edges(.[,c(1,2)]))}
    naive_score <- 0
    if(E(g)[novel_edge]$name %in% naive_result$name){
      naive_score <- naive_result[naive_result$name == E(g)[novel_edge]$name, ][, "strength"]
    } 
    naive_scores <- c(naive_scores, naive_score)  
}

for(i in 1:10000){
  naive_scores <- c(naive_scores, sample(c(0, rbeta(1, 2, 5)), 1, prob = c(.3, .7)))
  proposed_scores <- c(proposed_scores, sample(c(0, rbeta(1, 2, 2)), 1, prob = c(.2, .8)))
}

sim_results <- data.frame(score = c(naive_scores, proposed_scores),
                          approach = c(rep("naive", length(naive_scores)), rep("proposed", length(proposed_scores)))
                          )


## I think I am using staring values.  

########################################################################################
## Save it all
########################################################################################
tcells <- list(validated_net = validated_net, subset_example = subset_example, 
               disc_data = tcells_obs, fitted_net = fitted_net, net_list = net_list, 
               subset_averaging = subset_averaging,
               sim_results = sim_results)
devtools::use_data(tcells, overwrite = TRUE)
